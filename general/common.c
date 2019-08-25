#include "common.h"

void print_adj(const int nodes, const int degree, const int adj[nodes][degree])
{
  for(int i=0;i<nodes;i++){
    printf("%3d : ", i);
    for(int j=0;j<degree;j++)
      printf("%3d", adj[i][j]);
    printf("\n");
  }
}

void print_edge(const int nodes, const int degree, const int edge[nodes*degree/2][2])
{
  for(int i=0;i<nodes*degree/2;i++)
    printf("%d %d\n", edge[i][0], edge[i][1]);
}

void clear_buffer(int *buffer, const int n)
{
#pragma omp parallel for
  for(int i=0;i<n;i++)
    buffer[i] = 0;
}

void clear_buffers(uint64_t* restrict A, uint64_t* restrict B, const int s)
{
#pragma omp parallel for
  for(int i=0;i<s;i++)
    A[i] = B[i] = 0;
}

void edge_copy(int *restrict buf1, const int *restrict buf2, const int n)
{
#pragma omp parallel for
  for(int i=0;i<n;i++)
    buf1[i] = buf2[i];
}

int getRandom(const int max)
{
  return (int)(random()*((double)max)/(1.0+RAND_MAX));
}

static int get_end_edge(const int start_edge, const int groups, int (*edge)[2], 
			int line[groups], const int based_nodes)
{
  for(int i=0;i<groups;i++)
    for(int j=0;j<2;j++)
      if(edge[line[i]][j] % based_nodes != start_edge)
	return edge[line[i]][j] % based_nodes;

  return start_edge;
}

bool has_duplicated_edge(const int e00, const int e01, const int e10, const int e11)
{
  return ((e00 == e10 && e01 == e11) || (e00 == e11 && e01 == e10));
}

void swap(int *a, int *b)
{
  int tmp = *a;
  *a = *b;
  *b = tmp;
}

int order(int nodes, const int a, const int b, const int added_centers)
{
  if(!added_centers && nodes%2 == 0 && (a-b)%(nodes/2) == 0) return MIDDLE;
  if(a >= nodes-added_centers || b >= nodes-added_centers)   return MIDDLE;
  if(added_centers && (nodes-added_centers)%2 == 0 && (a-b)%((nodes-added_centers)/2)==0) return MIDDLE;

  if(added_centers) nodes -= added_centers;
  if(a < nodes/2.0){
    if(a > b) return LEFT;
    return (a+nodes/2.0 > b)? RIGHT : LEFT;
  }
  else{
    if(a < b) return RIGHT;
    return (a-nodes/2.0 > b)? RIGHT : LEFT;
  }
}

bool check_loop(const int lines, int edge[lines][2])
{
  timer_start(TIMER_CHECK);
  bool flag = true;

#pragma omp parallel for
  for(int i=0;i<lines;i++)
    if(edge[i][0] == edge[i][1])
      flag = false;

  timer_stop(TIMER_CHECK);
  return flag;
}

bool check_duplicate_edge(const int lines, int edge[lines][2])
{
  timer_start(TIMER_CHECK);
  bool flag = true;

  for(int i=0;i<2;i++)
    for(int j=2;j<lines;j++)
      if(has_duplicated_edge(edge[i][0], edge[i][1], edge[j][0], edge[j][1]))
        flag = false;
  
  timer_stop(TIMER_CHECK);
  return flag;
}

bool check_duplicate_current_edge(const int lines, const int groups, const int line[groups],
                                  int (*edge)[2], int tmp_edge[groups][2], const int original_groups,
				  const int nodes, const int added_centers)
{
  timer_start(TIMER_CHECK);
  int based_lines = lines/original_groups;
  int opt = (groups == original_groups)? 1 : 2;  // 1g-opt : 2g-opt
  bool flag = true;

  if(original_groups%2 == 1 && opt == 1){
    int tmp = line[0]%based_lines;

#pragma omp parallel for
    for(int i=rank;i<based_lines;i+=size)
      if(i != tmp)
	for(int j=0;j<groups;j++)
	  if(has_duplicated_edge(edge[i][0], edge[i][1], tmp_edge[j][0], tmp_edge[j][1]))
	    flag = false;
  }
  else if(opt == 2){
    int tmp0 = line[0]%based_lines;
    int tmp1 = line[1]%based_lines;

#pragma omp parallel for
    for(int i=rank;i<based_lines;i+=size)
      if(i != tmp0 && i != tmp1)
	for(int j=0;j<groups;j++)
	  if(has_duplicated_edge(edge[i][0], edge[i][1], tmp_edge[j][0], tmp_edge[j][1]))
	    flag = false;
  }
  else{ 
    assert(original_groups%2 == 0 && opt == 1);
    int tmp = line[0]%based_lines;
    if(distance(nodes, tmp_edge[0][0], tmp_edge[0][1], added_centers) != (nodes-added_centers)/2){

#pragma omp parallel for
      for(int i=rank;i<based_lines;i+=size)
	if(i != tmp)
	  for(int j=0;j<groups;j++)
	    if(has_duplicated_edge(edge[i][0], edge[i][1], tmp_edge[j][0], tmp_edge[j][1]))
	      flag = false;
    }
    else{

#pragma omp parallel for
      for(int i=rank;i<lines;i+=size)
        if(i%based_lines != tmp)
          for(int j=0;j<groups;j++)
            if(has_duplicated_edge(edge[i][0], edge[i][1], tmp_edge[j][0], tmp_edge[j][1]))
	      flag = false;
    }
  }
  
  MPI_Allreduce(MPI_IN_PLACE, &flag, 1,  MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
  timer_stop(TIMER_CHECK);
  return flag;
}

bool edge_1g_opt(int (*edge)[2], const int nodes, const int lines, const int degree, const int based_nodes, const int based_lines, 
		 const int groups, const int start_line, const int added_centers, int* restrict adjacency, const int ii)
{
  //  printf("1g-opt\n");
  if(groups == 1) // assert ?
    return true;

  if(edge[start_line][0] >= nodes-added_centers || edge[start_line][1] >= nodes-added_centers)
    return false;

  int line[groups], tmp_edge[groups][2], pattern;
  for(int i=0;i<groups;i++)
    line[i] = start_line % based_lines + i * based_lines;

  int start_edge = edge[line[0]][0] % based_nodes;
  int end_edge   = get_end_edge(start_edge, groups, edge, line, based_nodes);
  if(end_edge == start_edge){
    /* In n = 9, g = 4,
       edge[line[:]][:] = {1, 28}, {10, 1}, {19, 10}, {28, 19};
    */
    return false;
  }

  int diff = edge[line[0]][0] - edge[line[0]][1];
  while(1){
    pattern = (groups%2 == 0)? getRandom(groups+1) : getRandom(groups);
    if(pattern == groups){
      for(int i=0;i<groups/2;i++){
	tmp_edge[i][0] = start_edge + based_nodes * i;
	tmp_edge[i][1] = tmp_edge[i][0] + (nodes-added_centers)/2;
      }
      for(int i=groups/2;i<groups;i++){
	tmp_edge[i][0] = end_edge + based_nodes * (i-groups/2);
	tmp_edge[i][1] = tmp_edge[i][0] + (nodes-added_centers)/2;
      }
    }
    else{
      for(int i=0;i<groups;i++)
	tmp_edge[i][0] = start_edge + based_nodes * i;
      
      tmp_edge[0][1] = end_edge + based_nodes * pattern;
      for(int i=1;i<groups;i++){
	int tmp = tmp_edge[0][1] + based_nodes * i;
	tmp_edge[i][1] = (tmp < nodes-added_centers)? tmp : tmp - (nodes-added_centers);
      }
    }
    if(diff != (tmp_edge[0][0] - tmp_edge[0][1])) break;
  }

  assert(check_loop(groups, tmp_edge));
  assert(check_duplicate_edge(groups, tmp_edge));
  if(!check_duplicate_current_edge(lines, groups, line, edge, tmp_edge, groups, nodes, added_centers))
    return false;

  for(int i=0;i<groups;i++)
    if(order(nodes, tmp_edge[i][0], tmp_edge[i][1], added_centers) == RIGHT)
      swap(&tmp_edge[i][0], &tmp_edge[i][1]);  // RIGHT -> LEFT
  
  //  if(ii==48){
  /*
  printf("patten = %d\n", pattern);
  printf("nodes = %d, groups = %d\n", nodes, groups);
  for(int i=0;i<groups;i++)
    printf("line[%d] -> (%d, %d)\n", i, edge[line[i]][0], edge[line[i]][1]);
  printf("s = %d, e = %d\n", start_edge, end_edge);
  for(int i=0;i<groups;i++)
    printf(" tmp_edge = %d %d\n", tmp_edge[i][0], tmp_edge[i][1]);
  */
  //  }

  // Change a part of adj.
  int hash[nodes];
#pragma omp parallel for
  for(int i=0;i<groups;i++){
    int x0, x1;
    int y0 = edge[line[i]][0];
    int y1 = edge[line[i]][1];

    for(x0=0;x0<degree;x0++)
      if(adjacency[y0*degree+x0] == y1){
	hash[y0] = x0;
	break;
      }

    for(x1=0;x1<degree;x1++)
      if(adjacency[y1*degree+x1] == y0){
	hash[y1] = x1;
	break;
      }

    if(x0 == degree || x1 == degree)
      ERROR(":%d %d\n", x0, x1);
  }

  //  print_adj(nodes, degree, adjacency);
  //  printf("--\n");
  for(int i=0;i<groups;i++){
    int tmp0 = tmp_edge[i][0];
    int tmp1 = tmp_edge[i][1];
    adjacency[tmp0*degree+hash[tmp0]] = tmp1;
    adjacency[tmp1*degree+hash[tmp1]] = tmp0;
    //    printf("adj[%d][%d] = %d\n", tmp0, hash[tmp0], tmp1);
    //    printf("adj[%d][%d] = %d\n", tmp1, hash[tmp1], tmp0);
  }

  // Set vertexs
#pragma omp parallel for
  for(int i=0;i<groups;i++){
    edge[line[i]][0] = tmp_edge[i][0];
    edge[line[i]][1] = tmp_edge[i][1];
  }

  /*
  if(pattern == 0){
  printf("new:\n");
  for(int i=0;i<nodes;i++){
    printf("%2d:", i);
    for(int j=0;j<degree;j++)
      printf("%3d", adjacency[i*degree+j]);
    printf("\n");
  }

  create_adjacency(nodes, lines, degree, (const int (*)[2])edge, (int *)adjacency);

  printf("correct:\n");
  for(int i=0;i<nodes;i++){
    printf("%2d:", i);
    for(int j=0;j<degree;j++)
      printf("%3d", adjacency[i*degree+j]);
    printf("\n");
  }

  EXIT(0);
  }*/

  return true;
}
