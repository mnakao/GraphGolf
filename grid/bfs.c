#include "common.h"

#ifdef _OPENMP
static int top_down_step(const int level, const int nodes, const int num_frontier,
                         const int degree, const int* restrict adjacency, int* restrict frontier,
                         int* restrict next, int* restrict distance, char* restrict bitmap)
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
           bitmap[n]   = VISITED;
	   distance[n] = level;
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
                         int* restrict next, int* restrict distance, char* restrict bitmap)
{
  int count = 0;
  for(int i=0;i<num_frontier;i++){
    int v = frontier[i];
    for(int j=0;j<degree;j++){
      int n = *(adjacency + v * degree + j);  // int n = adjacency[v][j];
      if(bitmap[n] == NOT_VISITED){
	bitmap[n]   = VISITED;
	distance[n] = level;
        next[count++] = n;
      }
    }
  }

  return count;
}
#endif

bool evaluation(const int nodes, const int lines, const int degree,
		int adjacency[nodes][degree], int *diameter, double *ASPL)
{
  timer_start(TIMER_BFS);

  char *bitmap  = malloc(sizeof(char) * nodes);
  int *frontier = malloc(sizeof(int)  * nodes);
  int *distance = malloc(sizeof(int)  * nodes);
  int *next     = malloc(sizeof(int)  * nodes);
  bool reached  = true;
  double sum    = 0.0;
  *diameter     = 0;

  for(int s=rank;s<nodes;s+=size){
    int num_frontier = 1, level = 0;
    for(int i=0;i<nodes;i++)
      bitmap[i] = NOT_VISITED;

    frontier[0] = s;
    distance[s] = level;
    bitmap[s]   = VISITED;

    while(1){
      num_frontier = top_down_step(level++, nodes, num_frontier, degree,
				   (int *)adjacency, frontier, next, distance, bitmap);
      if(num_frontier == 0) break;
  
      int *tmp = frontier;
      frontier = next;
      next     = tmp;
    }

    *diameter = MAX(*diameter, level-1);
       
    for(int i=s+1;i<nodes;i++){
      if(bitmap[i] == NOT_VISITED)
	reached = false;

      sum += (distance[i] + 1);
    }
  }

  free(bitmap);
  free(frontier);
  free(distance);
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
