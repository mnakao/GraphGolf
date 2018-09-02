#include "common.h"

static void clear_buffer(const int n, int buffer[n])
{
  if(n<THRESHOLD){
    for(int i=0;i<n;i++)
      buffer[i] = NOT_VISITED;
  }
  else{  // When using if-clause, performance becomes slow
#pragma omp parallel for
    for(int i=0;i<n;i++)
      buffer[i] = NOT_VISITED;
  }
}

static int add_buffer(int *next, const int n, int count)
{
  for(int i=0;i<count;i++)
    if(next[i] == n)
      return count;

  next[count] = n;
  return ++count;
}

static int top_down_step(const int nodes, const int num_frontier, const int snode, const int degree,
			  const int* restrict adjacency, int* restrict frontier, int* restrict next, 
			  int* restrict parents, int* restrict distance, char* restrict bitmap)
{
  int count = 0;
#pragma omp parallel
  {
    int  local_count    = 0;
    int* local_frontier = (int*)malloc(sizeof(int) * nodes);
#pragma omp for
     for(int i=0;i<num_frontier;i++){
      int v = frontier[i];
      for(int j=0;j<degree;j++){
	int n = *(adjacency + v * degree + j);  // adjacency[v][j];
	if(bitmap[n] == 1 || snode == n) continue;
	if(parents[n] == NOT_VISITED){
	  bitmap[n]   = 1;
	  parents[n]  = v;
	  local_count = add_buffer(local_frontier, n, local_count);
	  distance[n] = distance[v] + 1;
	}
      }
     }
#pragma omp critical
    {
      memcpy(&next[count], local_frontier, local_count*sizeof(int));
      count += local_count;
    }
    free(local_frontier);
  }
  return count;
}

bool evaluation(const int nodes, int based_nodes, const int groups, const int lines, const int degree,
		int adjacency[nodes][degree], int *diameter, double *ASPL, const int added_centers)
{
  timer_start(TIMER_BFS);
  bool reeached = true;
  int distance[nodes], max = 0;
  double sum = 0.0;
  
  int first_task  = based_nodes;
  int second_task = added_centers - 1;
  int whole_task  = (added_centers)? first_task+second_task : first_task;
  int per_task    = (whole_task%size==0)? whole_task/size : whole_task/size+1;
  int start_task  = per_task*rank;
  if(whole_task <= start_task)
    per_task = 0;
  else if(whole_task <= start_task + per_task)
    per_task = whole_task - start_task;
  
  int *frontier = malloc(sizeof(int) * nodes);
  int *next     = malloc(sizeof(int) * nodes);
  int *parents  = malloc(sizeof(int) * nodes);
  char *bitmap  = malloc(sizeof(char) * nodes);

  //// First Search ////
  int tmp_per_task = (first_task < start_task)? 0 : per_task;
  tmp_per_task = MIN(tmp_per_task, first_task-start_task);
  if(tmp_per_task < 0)
    tmp_per_task = 0;
  
  for(int snode=start_task;snode<start_task+tmp_per_task;snode++){
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
      num_frontier = top_down_step(nodes, num_frontier, snode, degree, (int* restrict)adjacency,
				   frontier, next, parents, distance, bitmap);
  
      // Swap frontier <-> next
      int *tmp = frontier;
      frontier = next;
      free(tmp);
      next = malloc(sizeof(int) * nodes);
      clear_buffer(nodes, next);
    } while(num_frontier);
    
    for(int i=snode+1;i<groups*based_nodes;i++){
      if(distance[i] == 0)  // Never visit a node
	reeached = false;
      
      if(max < distance[i])
	max = distance[i];
      
      sum += distance[i] * (groups - i/based_nodes);
      // sum += distance[i];
    } 
    
    for(int i=groups*based_nodes;i<nodes;i++){
      if(distance[i] == 0)  // Never visit a node
	reeached = false;
      
      if(max < distance[i])
	max = distance[i];
      
      sum += distance[i] * groups;
    }
  }

  //// Second Search (only if add_centers exists) ////
  if(start_task > first_task)
    start_task = based_nodes * groups + (start_task-first_task);
  else
    start_task = based_nodes * groups;
  
  tmp_per_task = per_task - tmp_per_task;
  for(int snode=start_task;snode<start_task+tmp_per_task;snode++){
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
      num_frontier = top_down_step(nodes, num_frontier, snode, degree, (int* restrict)adjacency,
				   frontier, next, parents, distance, bitmap);
	  
      // Swap frontier <-> next
      int *tmp = frontier;
      frontier = next;
      free(tmp);
      next = malloc(sizeof(int) * nodes);
      clear_buffer(nodes, next);
    } while(num_frontier);
    
    for(int i=snode+1;i<nodes;i++){
      if(distance[i] == 0)  // Never visit a node
	reeached = false;
      
      if(max < distance[i])
	max = distance[i];
      
      sum += distance[i];
    }
  }
  
  free(frontier);
  free(next);
  free(parents);
  free(bitmap);

  MPI_Allreduce(MPI_IN_PLACE, &reeached, 1, MPI_BYTE, MPI_BAND, MPI_COMM_WORLD);
  if(! reeached){
    timer_stop(TIMER_BFS);
    return false;
  }

  MPI_Allreduce(MPI_IN_PLACE, &max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  *diameter = max;
  *ASPL     = (sum) / (((nodes-1)*nodes)/2);
  //  printf("%d\n", (int)sum);
  timer_stop(TIMER_BFS);
  return true;
}

