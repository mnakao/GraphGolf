#include "common.h"

void edge_copy(int *restrict buf1, const int *restrict buf2, const int n)
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int i=0;i<n;i++)
    buf1[i] = buf2[i];
}

int getRandom(const int max)
{
  return (int)(random()*((double)max)/(1.0+RAND_MAX));
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
      if(has_duplicated_edge(edge[i][0], edge[i][1], edge[j][0], edge[j][1]))
        flag = false;

  timer_stop(TIMER_CHECK);
  return flag;
}

bool check_duplicate_current_edge(const int lines, const int edge[lines][2], const int selected_line[2], int tmp_edge[2][2])
{
  timer_start(TIMER_CHECK);
  bool flag = true;

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int i=0;i<lines;i++)
    if(i != selected_line[0] && i != selected_line[1])
      for(int j=0;j<2;j++)
	if(has_duplicated_edge(edge[i][0], edge[i][1], tmp_edge[j][0], tmp_edge[j][1]))
	  flag = false;

  timer_stop(TIMER_CHECK);
  return flag;
}

bool edge_1g_opt(int (*edge)[2], const int nodes, const int lines, const int degree, const int based_nodes,
		 const int based_lines, const int groups, const int start_line, const int ii)
{
  if(groups == 1) // assert ?
    return true;

  int tmp_line[groups], tmp_edge[groups][2], pattern;
  for(int i=0;i<groups;i++)
    tmp_line[i] = start_line % based_lines + i * based_lines;
#if 0  
  int start_edge = edge[tmp_line[0]][0] % based_nodes;
  //  int end_edge   = get_end_edge(start_edge, groups, edge, tmp_line, based_nodes);
  //  if(end_edge == start_edge){
    /* In n = 9, g = 4,
       edge[tmp_line[:]][:] = {1, 28}, {10, 1}, {19, 10}, {28, 19};
    */
  //    return false;
  //  }

  int diff = edge[tmp_line[0]][0] - edge[tmp_line[0]][1];
  while(1){
    pattern = getRandom(groups+1);
    if(pattern == groups){
      for(int i=0;i<groups/2;i++){
        tmp_edge[i][0] = start_edge + based_nodes * i;
        tmp_edge[i][1] = tmp_edge[i][0] + nodes/2;
      }
      for(int i=groups/2;i<groups;i++){
        tmp_edge[i][0] = end_edge + based_nodes * (i-groups/2);
        tmp_edge[i][1] = tmp_edge[i][0] + nodes/2;
      }
    }
    else{
      for(int i=0;i<groups;i++)
        tmp_edge[i][0] = start_edge + based_nodes * i;

      tmp_edge[0][1] = end_edge + based_nodes * pattern;
      for(int i=1;i<groups;i++){
        int tmp = tmp_edge[0][1] + based_nodes * i;
        tmp_edge[i][1] = (tmp < nodes)? tmp : tmp - nodes;
      }
    }
    if(diff != (tmp_edge[0][0] - tmp_edge[0][1])) break;
  }

  assert(check_loop(groups, tmp_edge));
  assert(check_duplicate_tmp_edge(1, groups, tmp_edge));
  if(!check_duplicate_current_edge(lines, groups, tmp_line, edge, tmp_edge, groups, 1, (pattern == groups)))
    return false;
  
  for(int i=0;i<groups;i++)
    if(order(nodes, tmp_edge[i][0], tmp_edge[i][1]) == RIGHT)
      swap(&tmp_edge[i][0], &tmp_edge[i][1]);  // RIGHT -> LEFT
	  
  // Set vertexs
#pragma omp parallel for
  for(int i=0;i<groups;i++){
    edge[tmp_line[i]][0] = tmp_edge[i][0];
    edge[tmp_line[i]][1] = tmp_edge[i][1];
  }
#endif
  return true;
}
