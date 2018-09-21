#include "common.h"
// Ref. https://github.com/chaihf/BFS-OpenMP/blob/master/bfs/bfs_up_down_omp.cpp

static void clear_buffer(const int n, int buffer[n], int value)
{
#pragma omp parallel for
  for(int i=0;i<n;i++)
    buffer[i] = value;
}

static int top_down_step_thread(const int nodes, const int num_frontier, const int snode, const int degree,
				const int* restrict adjacency, int* restrict frontier, int* restrict next, 
				int* restrict parents, char* restrict bitmap)
{
  int count = 0;
#pragma omp parallel
  {
    int local_count = 0;
    int *local_frontier = malloc(nodes * sizeof(int));  // (num_frontier*degree)/threads * sizeof(int)
#pragma omp for nowait
     for(int i=0;i<num_frontier;i++){
       int v = frontier[i];
       for(int j=0;j<degree;j++){
	 int n = *(adjacency + v * degree + j);  // adjacency[v][j];
	 if(bitmap[n] == 1 || snode == n) continue;
	 if(__sync_bool_compare_and_swap(&parents[n], NOT_VISITED, parents[v]+1)){
	   bitmap[n] = 1;
	   local_frontier[local_count++] = n;
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

static int top_down_step_nothread(const int nodes, const int num_frontier, const int snode, const int degree,
				  const int* restrict adjacency, int* restrict frontier, int* restrict next, 
				  int* restrict parents, char* restrict bitmap)
{
  int count = 0;
  for(int i=0;i<num_frontier;i++){
    int v = frontier[i];
    for(int j=0;j<degree;j++){
      int n = *(adjacency + v * degree + j);  // adjacency[v][j];
      if(bitmap[n] == 1 || snode == n) continue;
      if(parents[n] == NOT_VISITED){
        bitmap[n] = 1;
        parents[n] = parents[v] + 1;
	next[count++] = n;
      }
    }
  }
  return count;
}

bool evaluation(const int nodes, int based_nodes, const int groups, const int lines, const int degree,
		int adjacency[nodes][degree], int *diameter, double *ASPL, const int added_centers)
{
  timer_start(TIMER_BFS);
  bool reeached = true;
  int max_diam = 0;
  double sum = 0.0;

  int (*top_down_step)(const int, const int, const int, const int,
		       const int* restrict, int* restrict, int* restrict,
		       int* restrict, char* restrict);
  if(threads == 1)
    top_down_step = top_down_step_nothread;
  else
    top_down_step = top_down_step_thread;

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
    // Initialize
    frontier[0] = snode;
    int num_frontier = 1;
    int tmp_diam = 0;
    clear_buffer(nodes, parents, NOT_VISITED);
    memset(bitmap, 0, nodes);

    do{
      num_frontier = top_down_step(nodes, num_frontier, snode, degree, (int *)adjacency,
				   frontier, next, parents, bitmap);
  
      // Swap frontier <-> next
      int *tmp = frontier;
      frontier = next;
      free(tmp);
      next = malloc(sizeof(int) * nodes);
      if(num_frontier != 0) tmp_diam++;
    } while(num_frontier);

    if(max_diam < tmp_diam)
      max_diam = tmp_diam;
       
    for(int i=snode+1;i<groups*based_nodes;i++){
      if(parents[i] == NOT_VISITED)  // Never visit a node
	reeached = false;
      
      sum += (parents[i] + 1) * (groups - i/based_nodes);
      // sum += distance[i];
    } 
    
    for(int i=groups*based_nodes;i<nodes;i++){
      if(parents[i] == NOT_VISITED)  // Never visit a node
	reeached = false;
      
      sum += (parents[i] + 1) * groups;
    }
  }

  //// Second Search (only if add_centers exists) ////
  if(start_task > first_task)
    start_task = based_nodes * groups + (start_task-first_task);
  else
    start_task = based_nodes * groups;
  
  tmp_per_task = per_task - tmp_per_task;
  for(int snode=start_task;snode<start_task+tmp_per_task;snode++){
    // Initialize
    frontier[0] = snode;
    int num_frontier = 1;
    int tmp_diam = 0;
    clear_buffer(nodes, parents, NOT_VISITED);
    memset(bitmap, 0, nodes);
    
    do{
      num_frontier = top_down_step(nodes, num_frontier, snode, degree, (int *)adjacency,
				   frontier, next, parents, bitmap);
	  
      // Swap frontier <-> next
      int *tmp = frontier;
      frontier = next;
      free(tmp);
      next = malloc(sizeof(int) * nodes);
      if(num_frontier != 0) tmp_diam++;
    } while(num_frontier);
    
    if(max_diam < tmp_diam)
      max_diam = tmp_diam;
	
    for(int i=snode+1;i<nodes;i++){
      if(parents[i] == NOT_VISITED)  // Never visit a node
	reeached = false;
      
      sum += parents[i] + 1;
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

  int all_max_diam = 0;
  MPI_Reduce(&max_diam, &all_max_diam, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  *diameter = all_max_diam;
  *ASPL     = (sum) / (((nodes-1)*nodes)/2);
  //  printf("%d\n", (int)sum);

  timer_stop(TIMER_BFS);
  return true;
}

