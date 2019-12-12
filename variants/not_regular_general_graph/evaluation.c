#include "common.h"

#ifdef _OPENMP
static int top_down_step(const int level, const int nodes, const int num_frontier,
                         const int degree, const int* restrict adj, int* restrict frontier,
                         int* restrict next, int* restrict distance,
			 const int *restrict num_degrees, char* restrict bitmap)
{
  int count = 0;
  int local_frontier[nodes];
#pragma omp parallel private(local_frontier)
  {
    int local_count = 0;
#pragma omp for nowait
     for(int i=0;i<num_frontier;i++){
       int v = frontier[i];
       for(int j=0;j<num_degrees[v];j++){
         int n = *(adj + v * degree + j);  // adj[v][j];
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
                         const int degree, const int* restrict adj, int* restrict frontier,
                         int* restrict next, int* restrict distance,
			 const int *restrict num_degrees, char* restrict bitmap)
{
  int count = 0;
  for(int i=0;i<num_frontier;i++){
    int v = frontier[i];
    for(int j=0;j<num_degrees[v];j++){
      int n = *(adj + v * degree + j);  // int n = adj[v][j];
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
		const int* restrict adj, int* restrict diam, double* restrict ASPL,
		const int num_degrees[nodes], const int added_centers)
{
  char *bitmap  = malloc(sizeof(char) * nodes);
  int *frontier = malloc(sizeof(int)  * nodes);
  int *distance = malloc(sizeof(int)  * nodes);
  int *next     = malloc(sizeof(int)  * nodes);
  bool reached  = true;
  double sum    = 0.0;
  *diam         = 0;

  for(int s=rank;s<based_nodes;s+=procs){
    int num_frontier = 1, level = 0;
    for(int i=0;i<nodes;i++)
      bitmap[i] = NOT_VISITED;

    frontier[0] = s;
    distance[s] = level;
    bitmap[s]   = VISITED;

    while(1){
      num_frontier = top_down_step(level++, nodes, num_frontier, degree,
				   adj, frontier, next, distance, num_degrees, bitmap);
      if(num_frontier == 0) break;
  
      int *tmp = frontier;
      frontier = next;
      next     = tmp;
    }

    *diam = MAX(*diam, level-1);
       
    for(int i=s+1;i<nodes;i++){
      if(bitmap[i] == NOT_VISITED)
	reached = false;

      if(i < groups*based_nodes)
	sum += (distance[i] + 1) * (groups - i/based_nodes);
      else
	sum += (distance[i] + 1) * groups; // for added_centers
    }
  }

  if(added_centers){
    int start_rank = based_nodes % procs;
    int start_node = based_nodes*groups+rank-start_rank;
    if(start_node < based_nodes*groups) start_node += procs;
    for(int s=start_node;s<nodes;s+=procs){
      int num_frontier = 1, level = 0;
      for(int i=0;i<nodes;i++)
	bitmap[i] = NOT_VISITED;
      
      frontier[0] = s;
      distance[s] = level;
      bitmap[s]   = VISITED;
      
      while(1){
	num_frontier = top_down_step(level++, nodes, num_frontier, degree, adj,
				     frontier, next, distance, num_degrees, bitmap);
	if(num_frontier == 0) break;
	
	int *tmp = frontier;
	frontier = next;
	next     = tmp;
      }
    
      *diam = MAX(*diam, level-1);
      
      for(int i=s+1;i<nodes;i++){
	if(bitmap[i] == NOT_VISITED)
	  reached = false;
	
	sum += distance[i] + 1;
      }
    }
  }
  
  free(bitmap);
  free(frontier);
  free(distance);
  free(next);

  MPI_Allreduce(MPI_IN_PLACE, &reached, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
  if(reached){
    MPI_Allreduce(MPI_IN_PLACE, diam, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    *ASPL = sum / ((((double)nodes-1)*nodes)/2);
  }
  else{
    *diam = INT_MAX;
    *ASPL = DBL_MAX;
    timer_stop(TIMER_APSP);
  }

  return reached;
}


static bool matrix_op(const int nodes, const int based_nodes, const int degree,
		      const int* restrict adj, const int groups, int* restrict diam,
		      double* restrict ASPL, const int* restrict num_degrees, const int added_centers)
{
  unsigned int elements = (based_nodes+(UINT64_BITS-1))/UINT64_BITS;
  unsigned int chunk    = (elements+(procs-1))/procs;
  size_t s    = nodes*chunk*sizeof(uint64_t);
  uint64_t* A = malloc(s);  // uint64_t A[nodes][chunk];
  uint64_t* B = malloc(s);  // uint64_t B[nodes][chunk];
  int parsize = (elements+(chunk-1))/chunk;
  double sum  = 0.0;

  *diam = 1;
  for(int t=rank;t<parsize;t+=procs){
    uint64_t kk, l;
    clear_buffers(A, B, nodes*chunk);
    for(l=0; l<UINT64_BITS*chunk && UINT64_BITS*t*chunk+l<based_nodes; l++){
      unsigned int offset = (UINT64_BITS*t*chunk+l)*chunk+l/UINT64_BITS;
      A[offset] = B[offset] = (0x1ULL<<(l%UINT64_BITS));
    }

    for(kk=0;kk<nodes;kk++){
#pragma omp parallel for
      for(int i=0;i<nodes;i++)
        for(int j=0;j<num_degrees[i];j++){
          int n = *(adj + i * degree + j);  // int n = adj[i][j];
          for(int k=0;k<chunk;k++)
            B[i*chunk+k] |= A[n*chunk+k];
        }

      uint64_t num1 = 0, num2 = 0;
#pragma omp parallel for reduction(+:num1)
      for(int i=0;i<based_nodes*groups*chunk;i++)
        num1 += POPCNT(B[i]);

#pragma omp parallel for reduction(+:num2)
      for(int i=based_nodes*groups*chunk;i<nodes*chunk;i++)
	num2 += POPCNT(B[i]);

      if(num1+num2 == (uint64_t)nodes*l) break;

      // swap A <-> B
      uint64_t* tmp = A;
      A = B;
      B = tmp;

      sum += ((double)based_nodes*groups * l - num1) * groups;
      sum += ((double)added_centers      * l - num2) * groups * 2;
    }
    *diam = MAX(*diam, kk+1);
  }

  if(added_centers){
    elements = (added_centers+(UINT64_BITS-1))/UINT64_BITS;
    chunk    = (elements+(procs-1))/procs;
    parsize  = (elements+(chunk-1))/chunk;

    int s = based_nodes % procs;
    int new_rank = (rank - s >= 0)? rank-s : rank-s+procs;
    for(int t=new_rank;t<parsize;t+=procs){
      uint64_t kk, l;
      clear_buffers(A, B, nodes*chunk);
      for(l=0; l<UINT64_BITS*chunk && UINT64_BITS*t*chunk+l<added_centers; l++){
	unsigned int offset = (UINT64_BITS*t*chunk+l+(nodes-added_centers))*chunk+l/UINT64_BITS;
	A[offset] = B[offset] = (0x1ULL<<(l%UINT64_BITS));
      }
      
      for(kk=0;kk<nodes;kk++){
#pragma omp parallel for
	for(int i=0;i<nodes;i++)
	  for(int j=0;j<degree;j++){
	    int n = *(adj + i * degree + j);  // int n = adj[i][j];
	    for(int k=0;k<chunk;k++)
	      B[i*chunk+k] |= A[n*chunk+k];
	  }
	
	uint64_t num = 0;
#pragma omp parallel for reduction(+:num)
	for(int i=based_nodes*groups*chunk;i<nodes*chunk;i++)
	  num += POPCNT(B[i]);

	if(num == (uint64_t)added_centers*l) break;

	// swap A <-> B
	uint64_t* tmp = A;
	A = B;
	B = tmp;
	
	sum += ((double)added_centers * l - num);
      }
      *diam = MAX(*diam, kk+1);
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, diam, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  sum += (double)nodes * (nodes - 1);

  free(A);
  free(B);

  if(*diam < nodes){
    *ASPL = sum / (((double)nodes-1)*nodes);
    return true;
  }
  else{
    *diam = INT_MAX;
    *ASPL = DBL_MAX;
    return false;
  }
}

static bool matrix_op_mem_saving(const int nodes, const int based_nodes, const int degree,
				 const int* restrict adj, const int groups, int* restrict diam,
				 double* restrict ASPL, const int* restrict num_degrees, const int added_centers)
{
  unsigned int elements = (based_nodes+(UINT64_BITS-1))/UINT64_BITS;
  size_t s    = nodes*CHUNK*sizeof(uint64_t);
  uint64_t* A = malloc(s);  // uint64_t A[nodes][CHUNK];
  uint64_t* B = malloc(s);  // uint64_t B[nodes][CHUNK];
  int parsize = (elements+(CHUNK-1))/CHUNK;
  double sum  = 0.0;

  *diam = 1;
  for(int t=rank;t<parsize;t+=procs){
    unsigned int kk, l;
    clear_buffers(A, B, nodes*CHUNK);
    for(l=0; l<UINT64_BITS*CHUNK && UINT64_BITS*t*CHUNK+l<based_nodes; l++){
      unsigned int offset = (UINT64_BITS*t*CHUNK+l)*CHUNK+l/UINT64_BITS;
      A[offset] = B[offset] = (0x1ULL<<(l%UINT64_BITS));
    }

    for(kk=0;kk<nodes;kk++){
#pragma omp parallel for
      for(int i=0;i<nodes;i++)
        for(int j=0;j<num_degrees[i];j++){
          int n = *(adj + i * degree + j);  // int n = adj[i][j];
          for(int k=0;k<CHUNK;k++)
            B[i*CHUNK+k] |= A[n*CHUNK+k];
        }

      uint64_t num1 = 0, num2 = 0;
#pragma omp parallel for reduction(+:num1)
      for(int i=0;i<based_nodes*groups*CHUNK;i++)
        num1 += POPCNT(B[i]);

#pragma omp parallel for reduction(+:num2)
      for(int i=based_nodes*groups*CHUNK;i<nodes*CHUNK;i++)
	num2 += POPCNT(B[i]);
      
      if(num1+num2 == (uint64_t)nodes*l) break;

      // swap A <-> B
      uint64_t* tmp = A;
      A = B;
      B = tmp;

      sum += ((double)based_nodes*groups * l - num1) * groups;
      sum += ((double)added_centers      * l - num2) * groups * 2;
    }
    *diam = MAX(*diam, kk+1);
  }

  if(added_centers){
    elements = (added_centers+(UINT64_BITS-1))/UINT64_BITS;
    parsize  = (elements+(CHUNK-1))/CHUNK;

    int s = based_nodes % procs;
    int new_rank = (rank - s >= 0)? rank-s : rank-s+procs;
    for(int t=new_rank;t<parsize;t+=procs){
      unsigned int kk, l;
      clear_buffers(A, B, nodes*CHUNK);
      for(l=0; l<UINT64_BITS*CHUNK && UINT64_BITS*t*CHUNK+l<added_centers; l++){
	unsigned int offset = (UINT64_BITS*t*CHUNK+l+(nodes-added_centers))*CHUNK+l/UINT64_BITS;
	A[offset] = B[offset] = (0x1ULL<<(l%UINT64_BITS));
      }
      
      for(kk=0;kk<nodes;kk++){
#pragma omp parallel for
	for(int i=0;i<nodes;i++)
	  for(int j=0;j<degree;j++){
	    int n = *(adj + i * degree + j);  // int n = adj[i][j];
	    for(int k=0;k<CHUNK;k++)
	      B[i*CHUNK+k] |= A[n*CHUNK+k];
	  }

	uint64_t num = 0;
#pragma omp parallel for reduction(+:num)
	for(int i=based_nodes*groups*CHUNK;i<nodes*CHUNK;i++)
	  num += POPCNT(B[i]);

	if(num == (uint64_t)added_centers*l) break;
	
	// swap A <-> B
	uint64_t* tmp = A;
	A = B;
	B = tmp;
	
	sum += ((double)added_centers * l - num);
      }
      *diam = MAX(*diam, kk+1);
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, diam, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  sum += (double)nodes * (nodes - 1);

  free(A);
  free(B);
  
  if(*diam < nodes){
    *ASPL = sum / (((double)nodes-1)*nodes);
    return true;
  }
  else{
    *diam = INT_MAX;
    *ASPL = DBL_MAX;
    return false;
  }
}

bool evaluation(const int nodes, int based_nodes, const int groups, const int lines, const int degree,
		int* restrict adj, int* restrict diam, double* restrict ASPL, const int added_centers,
		const int* restrict num_degrees, const int algo)
{
  timer_start(TIMER_APSP);
  
  bool ret;
  if(algo == BFS)
    ret = bfs(nodes, based_nodes, groups, lines, degree, adj, diam, ASPL, num_degrees, added_centers);
  else if(algo == MATRIX_OP)
    ret = matrix_op(nodes, based_nodes, degree, adj, groups, diam, ASPL, num_degrees, added_centers);
  else // (algo == MATRIX_OP_MEM_SAVING)
    ret = matrix_op_mem_saving(nodes, based_nodes, degree, adj, groups, diam, ASPL, num_degrees, added_centers);
  
  timer_stop(TIMER_APSP);
  
  return ret;
}
