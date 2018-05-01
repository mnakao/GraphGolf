#include "common.h"

static void clear_buffer(const int nodes, int buffer[nodes])
{
  for(int i=0;i<nodes;i++)
    buffer[i] = NOT_VISITED;
}

static int add_buffer(int *next, const int n, int count)
{
  for(int i=0;i<count;i++)
    if(next[i] == n)
      return count;

  next[count] = n;
  return ++count;
}

static void top_down_step(const int nodes, const int num_frontier, const int snode, const int degree,
                          const int* restrict adjacency, int* restrict frontier, int* restrict next, 
			  int* restrict parents, int* restrict distance, char* restrict bitmap)
{
  int count = 0;
  for(int i=0;i<num_frontier;i++){
    int v = frontier[i];
    for(int j=0;j<degree;j++){
      int n = *(adjacency + v * degree + j);  // adjacency[v][j];
      if(bitmap[n] == 1 || snode == n) continue;
      if(parents[n] == NOT_VISITED){
        bitmap[n] = 1;
        parents[n] = v;
        count = add_buffer(next, n, count);
        distance[n] = distance[v] + 1;
      }
    }
  }
}

static int get_num_frontier(const int nodes, const int frontier[nodes])
{
  int count = 0;
  for(int i=0;i<nodes;i++)
    if(frontier[i] != NOT_VISITED)
      count++;
    else
      break;

  return count;
}

bool evaluation(const int nodes, const int lines, const int degree, int adjacency[nodes][degree], 
		int *diameter, double *ASPL, const int rank, const int size)
{
  int distance[nodes], max = 0;
  int cancel_flag = false;
  double sum = 0.0;

  int node_size  = (nodes % size == 0)? nodes/size : nodes/size+1;
  int start_node = node_size * rank;
  if(rank == size-1)
    node_size = nodes - start_node;

#pragma omp parallel
  {
    int *frontier = malloc(sizeof(int) * nodes);
    int *next     = malloc(sizeof(int) * nodes);
    int *parents  = malloc(sizeof(int) * nodes);
    char *bitmap  = malloc(sizeof(char) * nodes);

#pragma omp for reduction(+:sum) reduction(max:max) private(distance) 
    for(int snode=start_node;snode<start_node+node_size;snode++){
      for(int i=0;i<nodes;i++){
	distance[i] = 0;
	bitmap[i] = 0;
      }
      
      // Initialize
      frontier[0] = snode;
      int num_frontier = 1;
      clear_buffer(nodes, next);
      clear_buffer(nodes, parents);

      do{
	top_down_step(nodes, num_frontier, snode, degree, (int* restrict)adjacency,
		      frontier, next, parents, distance, bitmap);
	// Swap frontier <-> next
	int *tmp = frontier;
	frontier = next;
	free(tmp);
	next = malloc(sizeof(int) * nodes);
	clear_buffer(nodes, next);
	num_frontier = get_num_frontier(nodes, frontier);
      } while(num_frontier);

      for(int i=snode+1;i<nodes;i++){
	if(distance[i] == 0){  // Never visit a node
	  cancel_flag = true;
#pragma omp cancel for if(cancel_flag)
	}

	if(max < distance[i])
	  max = distance[i];
	
	sum += distance[i];
      }
    }
    free(frontier);
    free(next);
    free(parents);
    free(bitmap);
  } // end omp parallel

  MPI_Allreduce(MPI_IN_PLACE, &cancel_flag, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if(cancel_flag)
    return false;

  MPI_Allreduce(MPI_IN_PLACE, &max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  *diameter = max;
  *ASPL     = sum / (((nodes-1)*nodes)/2);
  
  return true;
}

