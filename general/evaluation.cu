#include "common.h"
uint64_t *A_dev, *B_dev;
uint64_t *result, *result_dev;
int *adjacency_dev;

__global__ void matrix_op_init_dev(uint64_t* __restrict__ A, uint64_t* __restrict__ B,
                                   const int based_nodes, const int t, const int chunk)
{
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid<UINT64_BITS*chunk && UINT64_BITS*t*chunk+tid<based_nodes) {
    unsigned int offset = (UINT64_BITS*t*chunk+tid)*chunk+tid/UINT64_BITS;
    A[offset] = B[offset] = (0x1ULL<<(tid%UINT64_BITS));
    tid += blockDim.x * gridDim.x;
  }
}

__global__ void clear_buffers_dev(uint64_t* __restrict__ A, uint64_t* __restrict__ B,
				  const int length)
{
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid<length) {
    A[tid] = B[tid] = 0;
    tid += blockDim.x * gridDim.x;
  }
}

__global__ static void matrix_op_dev(const uint64_t* __restrict__ A, uint64_t* __restrict__ B,
				     const int* __restrict__ adjacency, const int nodes,
				     const int degree, const unsigned int elements)
{
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < nodes*elements) {
    int i = tid / elements;
    int k = tid % elements;
    uint64_t tmp = B[tid];
    for(int j=0;j<degree;j++){
      int n = *(adjacency + i * degree + j);  // int n = adjacency[i][j];
      tmp |= A[n*elements+k];
    }
    B[tid] = tmp;
    tid += blockDim.x * gridDim.x;
  }
}

__global__ static void popcnt_dev(const uint64_t* __restrict__ B, const int nodes, 
				  const unsigned int elements, uint64_t* __restrict__ result)
{
  __shared__ uint64_t cache[THREADS];
  int cacheIndex = threadIdx.x;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  uint64_t num = 0;
  while (tid < elements*nodes) {
    num += POPCNT(B[tid]);
    tid += blockDim.x * gridDim.x;
  }
  cache[cacheIndex] = num;
  __syncthreads();

  int i = blockDim.x/2;
  while (i != 0){
    if (cacheIndex < i)
      cache[cacheIndex] += cache[cacheIndex+i];
    __syncthreads();
    i /= 2;
  }

  if(cacheIndex == 0)
    result[blockIdx.x] = cache[0];
}

extern "C" bool matrix_op(const int nodes, const int degree, const int based_nodes,
			  const int* __restrict__ adjacency, const int groups,
			  int *diameter, double *ASPL)
{
  unsigned int elements = (based_nodes+UINT64_BITS-1)/UINT64_BITS;
  unsigned int chunk = (elements+(procs-1))/procs;
  int parsize = (elements + chunk - 1)/chunk;

  double sum = 0.0;
  *diameter = 1;
  cudaMemcpy(adjacency_dev, adjacency, sizeof(int)*nodes*degree, cudaMemcpyHostToDevice);
  
  for(int t=rank;t<parsize;t+=procs){
    unsigned int kk, l;
    for(l=0; l<UINT64_BITS*chunk && UINT64_BITS*t*chunk+l<based_nodes; l++){}
    
    clear_buffers_dev  <<< BLOCKS, THREADS >>> (A_dev, B_dev, nodes*chunk);
    matrix_op_init_dev <<< BLOCKS, THREADS >>> (A_dev, B_dev, based_nodes, t, chunk);
  
    for(kk=0;kk<nodes;kk++){
      matrix_op_dev <<< BLOCKS, THREADS >>> (A_dev, B_dev, adjacency_dev,
					     nodes, degree, chunk);
      popcnt_dev <<< BLOCKS, THREADS >>> (B_dev, nodes, chunk, result_dev);

      cudaMemcpy(result, result_dev, sizeof(uint64_t)*BLOCKS, cudaMemcpyDeviceToHost);
      uint64_t num = 0;
      for (int i=0;i<BLOCKS;i++)
	num += result[i];

      if(num == (uint64_t)nodes*l) break;

      // swap A <-> B
      uint64_t* tmp = A_dev;
      A_dev = B_dev;
      B_dev = tmp;
    
      sum += ((double)nodes * l - num) * groups;
    }
    *diameter = MAX(*diameter, kk+1);
  }
  MPI_Allreduce(MPI_IN_PLACE, diameter, 1, MPI_INT,    MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &sum,     1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  sum += (double)nodes * (nodes - 1);
  *ASPL = sum / (((double)nodes-1)*nodes);

  if(*diameter > nodes){
    //    ERROR("This graph is not connected graph.\n");
    return false;
  }
  
  return true;
}

extern "C" bool matrix_op_mem_saving(const int nodes, const int degree, const int based_nodes,
				     const int* __restrict__ adjacency,
				     const int groups, int *diameter, double *ASPL)
{
  unsigned int elements = (based_nodes+UINT64_BITS-1)/UINT64_BITS;
  int parsize = (elements+(CHUNK-1))/CHUNK;

  double sum = 0.0;
  *diameter = 1;
  cudaMemcpy(adjacency_dev, adjacency, sizeof(int)*nodes*degree, cudaMemcpyHostToDevice);

  for(int t=rank;t<parsize;t+=procs){
    unsigned int kk, l;
    for(l=0; l<UINT64_BITS*CHUNK && UINT64_BITS*t*CHUNK+l<based_nodes; l++){}
    
    clear_buffers_dev  <<< BLOCKS, THREADS >>> (A_dev, B_dev, nodes*CHUNK);
    matrix_op_init_dev <<< BLOCKS, THREADS >>> (A_dev, B_dev, based_nodes, t, CHUNK);

    for(kk=0;kk<nodes;kk++){
      matrix_op_dev <<< BLOCKS, THREADS >>> (A_dev, B_dev, adjacency_dev, nodes, degree, CHUNK);
      popcnt_dev <<< BLOCKS, THREADS >>> (B_dev, nodes, CHUNK, result_dev);
      
      cudaMemcpy(result, result_dev, sizeof(uint64_t)*BLOCKS, cudaMemcpyDeviceToHost);
      uint64_t num = 0;
      for (int i=0;i<BLOCKS;i++)
	num += result[i];

      if(num == (uint64_t)nodes*l) break;

      // swap A <-> B
      uint64_t* tmp = A_dev;
      A_dev = B_dev;
      B_dev = tmp;

      sum += ((double)nodes * l - num) * groups;
    }
    *diameter = MAX(*diameter, kk+1);
  }
  MPI_Allreduce(MPI_IN_PLACE, diameter, 1, MPI_INT,    MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &sum,     1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  sum += (double)nodes * (nodes - 1);
  *ASPL = sum / (((double)nodes-1)*nodes);

  if(*diameter > nodes){
    //    ERROR("This graph is not connected graph.\n");
    return false;
  }

  return true;
}

extern "C" void init_matrix_dev(const int nodes, const int degree, const int based_nodes,
				const int algo)
{
  cuInit(0);
  int gpus = -1;
  cudaGetDeviceCount(&gpus);
  cudaSetDevice(rank%gpus);
  unsigned int elements = (based_nodes+UINT64_BITS-1)/UINT64_BITS;
  size_t s = (algo == MATRIX_OP)? (elements+procs-1)/procs : CHUNK;
  s *= nodes * sizeof(uint64_t);

  cudaMalloc((void**)&A_dev, s);
  cudaMalloc((void**)&B_dev, s);
  cudaHostAlloc((void**)&result, BLOCKS*sizeof(uint64_t), cudaHostAllocDefault);
  cudaMalloc((void**)&result_dev,      sizeof(uint64_t)*BLOCKS);
  cudaMalloc((void**)&adjacency_dev,   sizeof(int)*nodes*degree);
}

extern "C" void finalize_matrix_dev()
{
  cudaFree(A_dev);
  cudaFree(B_dev);
  cudaFreeHost(result);
  cudaFree(result_dev);
  cudaFree(adjacency_dev);
}

extern "C" bool evaluation(const int nodes, const int based_nodes, const int groups, const int lines, 
			   const int degree, int* __restrict__ adjacency, int* __restrict__ diameter, 
			   double* __restrict__ ASPL, const int added_centers, const int algo)
{
  timer_start(TIMER_APSP);
  bool ret;
  if(algo == MATRIX_OP)
    ret = matrix_op(nodes, degree, based_nodes, adjacency, groups, diameter, ASPL);
  else // algo == MATRIX_OP_MEM_SAVING
    ret = matrix_op_mem_saving(nodes, degree, based_nodes, adjacency, groups, diameter, ASPL);
  
  timer_stop(TIMER_APSP);
  return ret;
}
