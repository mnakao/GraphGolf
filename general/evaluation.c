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

static bool bfs(const int nodes, int based_nodes, const int groups, const int lines, const int degree,
		int adjacency[nodes][degree], int *diameter, double *ASPL, const int added_centers)
{
  char *bitmap  = malloc(sizeof(char) * nodes);
  int *frontier = malloc(sizeof(int));
  int *distance = malloc(sizeof(int) * nodes);
  int *next     = malloc(sizeof(int) * nodes);
  bool reached  = true;
  double sum    = 0.0;
  *diameter     = 0;

  for(int s=rank;s<based_nodes;s+=size){
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
      free(tmp);
      next = malloc(sizeof(int) * nodes);
    }

    *diameter = MAX(*diameter, level-1);
       
    for(int i=s+1;i<nodes;i++){
      if(bitmap[i] == NOT_VISITED)
	reached = false;

      if(i < groups*based_nodes)
	sum += (distance[i] + 1) * (groups - i/based_nodes);
      else
	sum += (distance[i] + 1) * groups; // for add_centers
    }
  }

  // for add_centers
  for(int s=based_nodes*groups+rank;s<nodes;s+=size){
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
      free(tmp);
      next = malloc(sizeof(int) * nodes);
    }
    
    *diameter = MAX(*diameter, level-1);
	
    for(int i=s+1;i<nodes;i++){
      if(bitmap[i] == NOT_VISITED)
	reached = false;
      
      sum += distance[i] + 1;
    }
  }

  free(bitmap);
  free(frontier);
  free(distance);
  free(next);

  MPI_Allreduce(MPI_IN_PLACE, &reached, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
  if(!reached){
    timer_stop(TIMER_APSP);
    return false;
  }

  MPI_Allreduce(MPI_IN_PLACE, diameter, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &sum,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  *ASPL = sum / ((((double)nodes-1)*nodes)/2);

  return true;
}


static bool matrix_op(const int nodes, const int degree, const int* restrict adjacency,
		      const int groups, int *diameter, double *ASPL)
{
  unsigned int elements = (nodes/groups+(UINT64_BITS-1))/UINT64_BITS;
  unsigned int chunk = (elements+(size-1))/size;
  size_t s = nodes * chunk * sizeof(uint64_t);
  uint64_t* A = malloc(s);  // uint64_t A[nodes][chunk];
  uint64_t* B = malloc(s);  // uint64_t B[nodes][chunk];
  int parsize = (elements+(chunk-1))/chunk;
  double sum = 0.0;

  clear_buffers(A, B, nodes * chunk);
  *diameter = 1;
  for(int t=rank;t<parsize;t+=size){
    uint64_t kk, l;
    clear_buffers(A, B, nodes * chunk);
    for(l=0; l<UINT64_BITS*chunk && UINT64_BITS*t*chunk+l<nodes/groups; l++){
      unsigned int offset = (UINT64_BITS*t*chunk+l)*chunk+l/UINT64_BITS;
      A[offset] = B[offset] = (0x1ULL<<(l%UINT64_BITS));
    }

    for(kk=0;kk<nodes;kk++){
#pragma omp parallel for
      for(int i=0;i<nodes;i++)
        for(int j=0;j<degree;j++){
          int n = *(adjacency + i * degree + j);  // int n = adjacency[i][j];
          for(int k=0;k<chunk;k++)
            B[i*chunk+k] |= A[n*chunk+k];
        }

      uint64_t num = 0;
#pragma omp parallel for reduction(+:num)
      for(int i=0;i<chunk*nodes;i++)
        num += POPCNT(B[i]);

      if(num == (uint64_t)nodes*l) break;

      // swap A <-> B
      uint64_t* tmp = A;
      A = B;
      B = tmp;

      sum += ((double)nodes * l - num) * groups;
    }
    *diameter = MAX(*diameter, kk+1);
  }
  MPI_Allreduce(MPI_IN_PLACE, diameter, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  sum += (double)nodes * (nodes - 1);

  if(*diameter > nodes){
    //     PRINT_R0("This graph is not connected graph.\n");
    return false;
  }

  *ASPL = sum / (((double)nodes-1)*nodes);
  free(A);
  free(B);

  return true;
}

static bool matrix_op_low_mem(const int nodes, const int degree, const int* restrict adjacency,
                              const int groups, int *diameter, double *ASPL)
{
  unsigned int elements = (nodes/groups+(UINT64_BITS-1))/UINT64_BITS;
  size_t s = nodes * CHUNK * sizeof(uint64_t);
  uint64_t* A = malloc(s);  // uint64_t A[nodes][CHUNK];
  uint64_t* B = malloc(s);  // uint64_t B[nodes][CHUNK];
  int parsize = (elements + CHUNK - 1)/CHUNK;
  double sum = 0.0;

  *diameter = 1;
  for(int t=rank;t<parsize;t+=size){
    unsigned int kk, l;
    clear_buffers(A, B, nodes * CHUNK);
    for(l=0; l<UINT64_BITS*CHUNK && UINT64_BITS*t*CHUNK+l<nodes/groups; l++){
      unsigned int offset = (UINT64_BITS*t*CHUNK+l)*CHUNK+l/UINT64_BITS;
      A[offset] = B[offset] = (0x1ULL<<(l%UINT64_BITS));
    }

    for(kk=0;kk<nodes;kk++){
#pragma omp parallel for
      for(int i=0;i<nodes;i++)
        for(int j=0;j<degree;j++){
          int n = *(adjacency + i * degree + j);  // int n = adjacency[i][j];
          for(int k=0;k<CHUNK;k++)
            B[i*CHUNK+k] |= A[n*CHUNK+k];
        }

      uint64_t num = 0;
#pragma omp parallel for reduction(+:num)
      for(int i=0;i<CHUNK*nodes;i++)
        num += POPCNT(B[i]);

      if(num == (uint64_t)nodes*l) break;

      // swap A <-> B
      uint64_t* tmp = A;
      A = B;
      B = tmp;

      sum += ((double)nodes * l - num) * groups;
    }
    *diameter = MAX(*diameter, kk+1);
  }
  MPI_Allreduce(MPI_IN_PLACE, diameter, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  sum += (double)nodes * (nodes - 1);

  if(*diameter > nodes){
    //     PRINT_R0("This graph is not connected graph.\n");
    return false;
  }

  *ASPL = sum / (((double)nodes-1)*nodes);

  free(A);
  free(B);

  return true;
}

bool evaluation(const int nodes, int based_nodes, const int groups, const int lines, const int degree,
		int adjacency[nodes][degree], int *diameter, double *ASPL, const int added_centers, const int algo)
{
  timer_start(TIMER_APSP);
  
  bool ret;
  if(algo == BFS)
    ret = bfs(nodes, based_nodes, groups, lines, degree, adjacency, diameter, ASPL, added_centers);
  else if(algo == MATRIX_OP)
    ret = matrix_op(nodes, degree, (int *)adjacency, groups, diameter, ASPL);
  else // (algo == MATRIX_OP_LOW_MEM)
    ret = matrix_op_low_mem(nodes, degree, (int *)adjacency, groups, diameter, ASPL);
  
  timer_stop(TIMER_APSP);
  
  return ret;
}
