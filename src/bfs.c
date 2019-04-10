#include "common.h"

#ifdef _OPENMP
static int top_down_step(const int level, const int nodes, const int num_frontier,
                         const int degree, const int* restrict adjacency, int* restrict frontier,
                         int* restrict next, unsigned char* restrict bitmap)
{
  int count = 0;
  int local_frontier[nodes];
#pragma omp parallel private(local_frontier)
  {
    int local_count = 0;
#pragma omp for nowait
     for(int i=0;i<num_frontier;i++){
       int v = frontier[i];
       for(int j=0;j<degree;j++){
         int n = *(adjacency + v * degree + j);  // adjacency[v][j];
         if(bitmap[n] == NOT_VISITED){
           bitmap[n] = level;
           local_frontier[local_count++] = n;
         }
       }
     }  // end for i
#pragma omp critical
     {
       memcpy(&next[count], local_frontier, local_count*sizeof(int));
       count += local_count;
     }
  }
  return count;
}
#else
static int top_down_step(const int level, const int nodes, const int num_frontier,
                         const int degree, const int* restrict adjacency, int* restrict frontier,
                         int* restrict next, unsigned char* restrict bitmap)
{
  int count = 0;
  for(int i=0;i<num_frontier;i++){
    int v = frontier[i];
    for(int j=0;j<degree;j++){
      int n = *(adjacency + v * degree + j);  // int n = adjacency[v][j];
      if(bitmap[n] == NOT_VISITED){
        bitmap[n] = level;
        next[count++] = n;
      }
    }
  }

  return count;
}
#endif

bool evaluation(const int nodes, int based_nodes, const int groups, const int lines, const int degree,
		int adjacency[nodes][degree], int *diameter, double *ASPL, const int added_centers)
{
  timer_start(TIMER_BFS);

  unsigned char *bitmap = malloc(sizeof(unsigned char) * nodes);
  int *frontier = malloc(sizeof(int));
  int *next     = malloc(sizeof(int) * nodes);
  bool reached  = true;
  double sum    = 0.0;
  *diameter     = 0;

  for(int s=rank;s<based_nodes;s+=size){
    int num_frontier = 1, level = 0;
    for(int i=0;i<nodes;i++)
      bitmap[i] = NOT_VISITED;

    frontier[0] = s;
    bitmap[s]   = level;

    while(1){
      num_frontier = top_down_step(level++, nodes, num_frontier, degree,
				   (int *)adjacency, frontier, next, bitmap);
      if(num_frontier == 0) break;
  
      int *tmp = frontier;
      frontier = next;
      free(tmp);
      next = malloc(sizeof(int) * nodes);
    }

    *diameter = MAX(*diameter, level-1);
       
    for(int i=s+1;i<nodes;i++){
      if(bitmap[i] == NOT_VISITED)
	reached = false;

      if(i < groups*based_nodes)
	sum += (bitmap[i] + 1) * (groups - i/based_nodes);
      else
	sum += (bitmap[i] + 1) * groups; // for add_centers
    }
  }

  // for add_centers
  for(int s=based_nodes*groups+rank;s<nodes;s+=size){
    int num_frontier = 1, level = 0;
    for(int i=0;i<nodes;i++)
      bitmap[i] = NOT_VISITED;
    
    frontier[0] = s;
    bitmap[s]   = level;
    
    while(1){
      num_frontier = top_down_step(level++, nodes, num_frontier, degree,
				   (int *)adjacency, frontier, next, bitmap);
      if(num_frontier == 0) break;
	  
      int *tmp = frontier;
      frontier = next;
      free(tmp);
      next = malloc(sizeof(int) * nodes);
    }
    
    *diameter = MAX(*diameter, level-1);
	
    for(int i=s+1;i<nodes;i++){
      if(bitmap[i] == NOT_VISITED)
	reached = false;
      
      sum += bitmap[i] + 1;
    }
  }

  free(bitmap);
  free(frontier);
  free(next);

  MPI_Allreduce(MPI_IN_PLACE, &reached, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
  if(!reached){
    timer_stop(TIMER_BFS);
    return false;
  }

  MPI_Allreduce(MPI_IN_PLACE, diameter, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &sum,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  *ASPL = sum / ((((double)nodes-1)*nodes)/2);

  timer_stop(TIMER_BFS);
  return true;
}
