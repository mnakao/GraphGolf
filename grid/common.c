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
