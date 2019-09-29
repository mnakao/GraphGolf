#include "common.h"

int getRandom(const int max)
{
  return (int)(random()*((double)max)/(1.0+RAND_MAX));
}

bool has_duplicated_edge(const int e00, const int e01, const int e10, const int e11)
{
  return ((e00 == e10 && e01 == e11) || (e00 == e11 && e01 == e10));
}

bool check_symmetric_edge(const int lines, const int edge[lines][2], const int height,
			  const int width, const int based_height, const int groups)
{
  assert(lines%groups == 0);
  int based_lines = lines / groups;
  int w[2], h[2], tmp_edge[2];

  if(groups == 2){
    for(int i=0;i<based_lines;i++){
      for(int j=0;j<2;j++){
        w[j] = edge[i][j]/height;
        h[j] = edge[i][j]%height;
	tmp_edge[j] = (width-w[j]-1)*height + (height-h[j]-1);
      }
      if(!has_duplicated_edge(edge[based_lines+i][0], edge[based_lines+i][1], tmp_edge[0], tmp_edge[1]))
	return false;
    }
  }
  else if(groups == 4){
    for(int i=0;i<based_lines;i++){
      for(int j=0;j<2;j++){
	w[j] = edge[i][j]/height;
	h[j] = edge[i][j]%height;
      }
      
      for(int j=0;j<2;j++) tmp_edge[j] = h[j]*height + (height-w[j]-1);
      if(!has_duplicated_edge(edge[based_lines+i][0], edge[based_lines+i][1], tmp_edge[0], tmp_edge[1]))
        return false;

      for(int j=0;j<2;j++) tmp_edge[j] = (height-w[j]-1)*height + (height-h[j]-1);
      if(!has_duplicated_edge(edge[based_lines*2+i][0], edge[based_lines*2+i][1], tmp_edge[0], tmp_edge[1]))
	return false;

      for(int j=0;j<2;j++) tmp_edge[j] = (height-h[j]-1)*height + w[j];
      if(!has_duplicated_edge(edge[based_lines*3+i][0], edge[based_lines*3+i][1], tmp_edge[0], tmp_edge[1]))
	return false;
    }
  }
  
  return true;
}

void output_edge(const int lines, const int edge[lines][2], const int height)
{
  for(int i=0;i<lines;i++)
    printf("%d,%d %d,%d\n",
	   edge[i][0]/height, edge[i][0]%height, edge[i][1]/height, edge[i][1]%height);
}

void edge_copy(int *restrict buf1, const int *restrict buf2, const int n)
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
  
  int tmp_line[groups], tmp_edge[groups][2];
  for(int i=0;i<groups;i++){
    int t = start_line + based_lines * i;
    tmp_line[i] = (t < lines)? t : t-lines;
  }
  
  int pattern = getRandom(groups);
  if(groups == 2){
    if(pattern == 0){
      tmp_edge[0][0] = edge[tmp_line[0]][0]; tmp_edge[0][1] = edge[tmp_line[1]][0];
      tmp_edge[1][0] = edge[tmp_line[0]][1]; tmp_edge[1][1] = edge[tmp_line[1]][1];
    }
    else{ // pattern == 1
      tmp_edge[0][0] = edge[tmp_line[0]][0]; tmp_edge[0][1] = edge[tmp_line[1]][1];
      tmp_edge[1][0] = edge[tmp_line[0]][1]; tmp_edge[1][1] = edge[tmp_line[1]][0];
    }
  }
  else{ // groups == 4
    if(pattern == 0){
      tmp_edge[0][0] = edge[tmp_line[0]][0]; tmp_edge[0][1] = edge[tmp_line[1]][1];
      tmp_edge[1][0] = edge[tmp_line[1]][0]; tmp_edge[1][1] = edge[tmp_line[2]][1];
      tmp_edge[2][0] = edge[tmp_line[2]][0]; tmp_edge[2][1] = edge[tmp_line[3]][1];
      tmp_edge[3][0] = edge[tmp_line[3]][0]; tmp_edge[3][1] = edge[tmp_line[0]][1];
    }
    else if(pattern == 1){
      tmp_edge[0][0] = edge[tmp_line[0]][0]; tmp_edge[0][1] = edge[tmp_line[2]][1];
      tmp_edge[1][0] = edge[tmp_line[1]][0]; tmp_edge[1][1] = edge[tmp_line[3]][1];
      tmp_edge[2][0] = edge[tmp_line[2]][0]; tmp_edge[2][1] = edge[tmp_line[0]][1];
      tmp_edge[3][0] = edge[tmp_line[3]][0]; tmp_edge[3][1] = edge[tmp_line[1]][1];
    }
    else if(pattern == 2){
      tmp_edge[0][0] = edge[tmp_line[0]][0]; tmp_edge[0][1] = edge[tmp_line[3]][1];
      tmp_edge[1][0] = edge[tmp_line[1]][0]; tmp_edge[1][1] = edge[tmp_line[0]][1];
      tmp_edge[2][0] = edge[tmp_line[2]][0]; tmp_edge[2][1] = edge[tmp_line[1]][1];
      tmp_edge[3][0] = edge[tmp_line[3]][0]; tmp_edge[3][1] = edge[tmp_line[2]][1];
    }
    else{ // pattern == 3
      tmp_edge[0][0] = edge[tmp_line[0]][0]; tmp_edge[0][1] = edge[tmp_line[2]][0];
      tmp_edge[1][0] = edge[tmp_line[1]][0]; tmp_edge[1][1] = edge[tmp_line[3]][0];
      tmp_edge[2][0] = edge[tmp_line[0]][1]; tmp_edge[2][1] = edge[tmp_line[2]][1];
      tmp_edge[3][0] = edge[tmp_line[1]][1]; tmp_edge[3][1] = edge[tmp_line[3]][1];
    }
  }

  for(int i=0;i<groups;i++){
    edge[tmp_line[i]][0] = tmp_edge[i][0];
    edge[tmp_line[i]][1] = tmp_edge[i][1];
  }

  return true;
}
