#include "common.h"

void printb(const uint64_t v)
{
  uint64_t mask = 0x1ULL << (sizeof(v) * CHAR_BIT - 1);
  int sum = 0;
  do{
    putchar(mask & v ? '1' : '0');
    sum++;
    if(sum%8==0) putchar(',');
  } while (mask >>= 1);
}

void create_adjacency(const int nodes, const int lines, const int degree,
                      const int edge[lines][2], int adjacency[nodes][degree])
{
  int count[nodes];
  for(int i=0;i<nodes;i++)
    count[i] = 0;

  for(int i=0;i<lines;i++){
    int n1 = edge[i][0];
    int n2 = edge[i][1];
    adjacency[n1][count[n1]++] = n2;
    adjacency[n2][count[n2]++] = n1;
  }
}

bool has_duplicated_vertex(const int e00, const int e01, const int e10, const int e11)
{
  return (e00 == e10 || e01 == e11 || e00 == e11 || e01 == e10);
}

bool has_duplicated_edge(const int e00, const int e01, const int e10, const int e11)
{
  return (e00 == e10 && e01 == e11) || (e00 == e11 && e01 == e10);
}

int getRandom(const int max)
{
  return (int)(random()*((double)max)/(1.0+RAND_MAX));
}

int WIDTH(const int v, const int height)
{
  return v/height;
}

int HEIGHT(const int v, const int height)
{
  return v%height;
}

int ROTATE(const int v, const int nodes, const int groups, const int degree)
{
  if(groups != 2 && groups != 4)
    ERROR("Invalid groups\n");

  if((groups == 2 && degree != 180) || 
     (groups == 4 && (degree != 90 && degree != 180 && degree != 270)))
    ERROR("Invalid degree\n");
  
  int new_v = v + (nodes*degree/360);
  return (new_v < nodes)? new_v : new_v - nodes;
}

void output_edge(const int lines, const int edge[lines*2], const int height)
{
  for(int i=0;i<lines;i++)
    printf("%d,%d %d,%d\n", WIDTH(edge[i*2], height), HEIGHT(edge[i*2], height),
	   WIDTH(edge[i*2+1], height), HEIGHT(edge[i*2+1], height));
}

void copy_edge(int *restrict buf1, const int *restrict buf2, const int n)
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int i=0;i<n;i++)
    buf1[i] = buf2[i];
}

void swap(int *a, int *b)
{
  int tmp = *a;
  *a = *b;
  *b = tmp;
}

bool check_loop(const int lines, const int edge[lines][2])
{
  timer_start(TIMER_CHECK);
  bool flag = true;

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int i=0;i<lines;i++)
    if(edge[i][0] == edge[i][1])
      flag = false;

  timer_stop(TIMER_CHECK);

  if(flag == false){
    for(int i=0;i<lines;i++)
      if(edge[i][0] == edge[i][1]){
	printf("%d: %d %d <--\n", i, edge[i][0], edge[i][1]);
      }
      else{
	printf("%d: %d %d\n", i, edge[i][0], edge[i][1]);
      }
  }
  
  return flag;
}

bool check_duplicate_all_edge(const int lines, const int edge[lines][2])
{
  timer_start(TIMER_CHECK);
  bool flag = true;

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int i=0;i<lines;i++)
    for(int j=i+1;j<lines;j++)
      if(has_duplicated_edge(edge[i][0], edge[i][1], edge[j][0], edge[j][1])){
	printf("%d %d %d %d\n", edge[i][0], edge[i][1], edge[j][0], edge[j][1]);
        flag = false;
      }

  timer_stop(TIMER_CHECK);
  return flag;
}

bool check_duplicate_tmp_edge(const int g_opt, const int groups, int tmp_edge[groups*g_opt][2])
{
  timer_start(TIMER_CHECK);
  bool flag = true;

  for(int i=0;i<g_opt;i++){
    int tmp[2] = {tmp_edge[i][0], tmp_edge[i][1]};
    for(int j=g_opt;j<groups*g_opt;j++)
      if(has_duplicated_edge(tmp[0], tmp[1], tmp_edge[j][0], tmp_edge[j][1]))
	flag = false;
  }

  timer_stop(TIMER_CHECK);
  return flag;
}

bool check_duplicate_current_edge(const int lines, const int edge[lines][2], const int tmp_lines,
				  const int tmp_edge[tmp_lines][2], const int tmp_line[2],
				  const int groups, const int g_opt, const bool is_center)
{
  timer_start(TIMER_CHECK);
  int based_lines = lines/groups;
  bool flag = true;

  if(g_opt == D_2G_OPT){
    int tmp_line0 = tmp_line[0]%based_lines;
    int tmp_line1 = tmp_line[1]%based_lines;

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int i=rank;i<based_lines;i+=procs)
      if(i != tmp_line0 && i != tmp_line1)
        for(int j=0;j<tmp_lines;j++)
          if(has_duplicated_edge(edge[i][0], edge[i][1], tmp_edge[j][0], tmp_edge[j][1]))
            flag = false;
  }
  else if(g_opt == D_1G_OPT){
    int tmp_line0 = tmp_line[0]%based_lines;
    if(! is_center){
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(int i=rank;i<based_lines;i+=procs)
	if(i != tmp_line0)
	  for(int j=0;j<tmp_lines;j++)
	    if(has_duplicated_edge(edge[i][0], edge[i][1], tmp_edge[j][0], tmp_edge[j][1]))
	      flag = false;
    }
    else{
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(int i=rank;i<lines;i+=procs)
        if(i%based_lines != tmp_line0)
          for(int j=0;j<tmp_lines;j++)
            if(has_duplicated_edge(edge[i][0], edge[i][1], tmp_edge[j][0], tmp_edge[j][1]))
              flag = false;
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &flag, 1,  MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);   
  timer_stop(TIMER_CHECK);
  return flag;
}

void create_rotate_table(const int height, const int width, const int groups, int *table)
{
  int nodes = height * width;
  int based_nodes = nodes / groups;
  int based_height = height / 2;
  
  if(groups == 1){
    for(int i=0;i<based_nodes;i++)
      table[i] = i;
  }
  else if(groups == 2){
    for(int i=0;i<based_nodes;i++){
      int j = (i/based_height) * height + (i%based_height);
      int w = WIDTH (j, height);
      int h = HEIGHT(j, height);
      table[i] = j;
      table[i+based_nodes] = (width-w-1)*height + (height-h-1);
    }
  }
  else{
    for(int i=0;i<based_nodes;i++){
      int j = (i/based_height) * height + (i%based_height);
      int w = WIDTH (j, height);
      int h = HEIGHT(j, height);
      table[i] = j;
      table[i+based_nodes  ] = h*height + (height-w-1);
      table[i+based_nodes*2] = (height-w-1)*height + (height-h-1);
      table[i+based_nodes*3] = (height-h-1)*height + w; 
    }
  }
}
