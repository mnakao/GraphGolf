#include "common.h"

int getRandom(const int max)
{
  return (int)(random()*((double)max)/(1.0+RAND_MAX));
}

bool has_duplicated_edge(const int e00, const int e01, const int e10, const int e11)
{
  return ((e00 == e10 && e01 == e11) || (e00 == e11 && e01 == e10));
}

int WIDTH(const int v, const int height)
{
  return v/height;
}

int HEIGHT(const int v, const int height)
{
  return v%height;
}

int ROTATE(const int v, const int height, const int degree)
{
  if(degree != 90 && degree != 180 && degree != 270)
    ERROR("Invalid degree\n");

  int w = WIDTH (v, height);
  int h = HEIGHT(v, height);
  
  if(degree == 90)       return h*height + (height-w-1);
  else if(degree == 180) return (height-w-1)*height + (height-h-1);
  else                   return (height-h-1)*height + w; // degree == 270
}

bool check_symmetric_edge(const int lines, const int edge[lines][2], const int height,
			  const int width, const int based_height, const int groups)
{
  assert(lines%groups == 0);
  int tmp_edge[2], based_lines = lines / groups;

  if(groups == 2){
    for(int i=0;i<based_lines;i++){
      for(int j=0;j<2;j++)
	tmp_edge[j] = ROTATE(edge[i][j], height, 180);

      if(!has_duplicated_edge(edge[based_lines+i][0], edge[based_lines+i][1], tmp_edge[0], tmp_edge[1]))
	return false;
    }
  }
  else if(groups == 4){
    // 90 degrees
    for(int i=0;i<based_lines;i++){
      for(int j=0;j<2;j++)
	tmp_edge[j] = ROTATE(edge[i][j], height, 90);
      if(!has_duplicated_edge(tmp_edge[0], tmp_edge[1], edge[based_lines+i][0], edge[based_lines+i][1])){
	if(!( WIDTH (edge[based_lines+i][0], height) + WIDTH (edge[based_lines+i][1], height) == (width-1) &&
	      HEIGHT(edge[based_lines+i][0], height) + HEIGHT(edge[based_lines+i][1], height) == (height-1))){
	  printf("A i=%d: %d,%d-%d,%d %d,%d-%d,%d\n", i,
		 WIDTH(edge[based_lines+i][0], height), HEIGHT(edge[based_lines+i][0], height),
		 WIDTH(edge[based_lines+i][1], height), HEIGHT(edge[based_lines+i][1], height),
		 WIDTH(tmp_edge[0], height), HEIGHT(tmp_edge[0], height),
		 WIDTH(tmp_edge[1], height), HEIGHT(tmp_edge[1], height));
	  return false;
	}
      }

      // 180 degrees
      for(int j=0;j<2;j++)
	tmp_edge[j] = ROTATE(edge[i][j], height, 180);
      if(!has_duplicated_edge(tmp_edge[0], tmp_edge[1], edge[based_lines*2+i][0], edge[based_lines*2+i][1])){
	if(!( WIDTH (edge[based_lines*2+i][0], height) + WIDTH (edge[based_lines*2+i][1], height) == (width-1) &&
	      HEIGHT(edge[based_lines*2+i][0], height) + HEIGHT(edge[based_lines*2+i][1], height) == (height-1))){
	  printf("B i=%d: %d,%d-%d,%d %d,%d-%d,%d\n", i, 
		 WIDTH(edge[based_lines*2+i][0], height), HEIGHT(edge[based_lines*2+i][0], height),
		 WIDTH(edge[based_lines*2+i][1], height), HEIGHT(edge[based_lines*2+i][1], height),
		 WIDTH(tmp_edge[0], height), HEIGHT(tmp_edge[0], height),
		 WIDTH(tmp_edge[1], height), HEIGHT(tmp_edge[1], height));
	  return false;
	}
      }

      // 270 degrees
      for(int j=0;j<2;j++)
	tmp_edge[j] = ROTATE(edge[i][j], height, 270);
      if(!has_duplicated_edge(tmp_edge[0], tmp_edge[1], edge[based_lines*3+i][0], edge[based_lines*3+i][1])){
	if(!( WIDTH (edge[based_lines*3+i][0], height) + WIDTH (edge[based_lines*3+i][1], height) == (width-1) &&
	      HEIGHT(edge[based_lines*3+i][0], height) + HEIGHT(edge[based_lines*3+i][1], height) == (height-1))){
	  printf("C i=%d: %d,%d-%d,%d %d,%d-%d,%d\n", i,
		 WIDTH(edge[based_lines*3+i][0], height), HEIGHT(edge[based_lines*3+i][0], height),
		 WIDTH(edge[based_lines*3+i][1], height), HEIGHT(edge[based_lines*3+i][1], height),
		 WIDTH(tmp_edge[0], height), HEIGHT(tmp_edge[0], height),
		 WIDTH(tmp_edge[1], height), HEIGHT(tmp_edge[1], height));
	  return false;
	}
      }
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
		 const int based_lines, const int height, const int groups, const int start_line, const int ii)
{
  //  printf("edge_1g_opt\n");
  if(groups == 1) // assert ?
    return true;
  
  int tmp_line[groups], tmp_edge[groups][2];
  for(int i=0;i<groups;i++){
    int t = start_line + based_lines * i;
    tmp_line[i] = (t < lines)? t : t-lines;
  }

  //  for(int i=0;i<groups;i++)
  //    printf("line[%d] = %d\n", i, tmp_line[i]);

  //  for(int i=0;i<groups;i++)
  //    printf("Before: %d,%d-%d,%d\n",
  //	   WIDTH (edge[tmp_line[i]][0], height),
  //	   HEIGHT(edge[tmp_line[i]][0], height),
  //	   WIDTH (edge[tmp_line[i]][1], height),
  //	   HEIGHT(edge[tmp_line[i]][1], height));

  if(groups == 2){
    if(getRandom(2) == 0){
      tmp_edge[0][0] = edge[tmp_line[0]][0]; tmp_edge[0][1] = edge[tmp_line[1]][0];
      tmp_edge[1][0] = edge[tmp_line[0]][1]; tmp_edge[1][1] = edge[tmp_line[1]][1];
    }
    else{
      tmp_edge[0][0] = edge[tmp_line[0]][0]; tmp_edge[0][1] = edge[tmp_line[1]][1];
      tmp_edge[1][0] = edge[tmp_line[1]][0]; tmp_edge[1][1] = edge[tmp_line[0]][1];
    }
  }
  else{ // groups == 4
    if(edge[tmp_line[0]][0] == edge[tmp_line[groups-1]][1])  // A cycle is composed of four edges.
      return false;

    int s = edge[tmp_line[0]][0];
    int e = edge[tmp_line[2]][1];
    while(1){
      int pattern = getRandom(groups+1);
      if(pattern != groups){
	tmp_edge[0][0] = s;
	tmp_edge[1][0] = ROTATE(s, height, 90);
	tmp_edge[2][0] = ROTATE(s, height, 180);
	tmp_edge[3][0] = ROTATE(s, height, 270);
	
	if(pattern == 0){
	  tmp_edge[0][1] = e;
	  tmp_edge[1][1] = ROTATE(e, height, 90);
	  tmp_edge[2][1] = ROTATE(e, height, 180);
	  tmp_edge[3][1] = ROTATE(e, height, 270);
	}
	else if(pattern == 1){
	  tmp_edge[0][1] = ROTATE(e, height, 270);
	  tmp_edge[1][1] = e;
	  tmp_edge[2][1] = ROTATE(e, height, 90);
	  tmp_edge[3][1] = ROTATE(e, height, 180);
	}
	else if(pattern == 2){
	  tmp_edge[0][1] = ROTATE(e, height, 180);
	  tmp_edge[1][1] = ROTATE(e, height, 270);
	  tmp_edge[2][1] = e;
	  tmp_edge[3][1] = ROTATE(e, height, 90);
	}
	else{ // pattern == 3
	  tmp_edge[0][1] = ROTATE(e, height, 90);
	  tmp_edge[1][1] = ROTATE(e, height, 180);
	  tmp_edge[2][1] = ROTATE(e, height, 270);
	  tmp_edge[3][1] = e;
	}
      }
      else{ // pattern == groups
	tmp_edge[0][0] = s;
	tmp_edge[0][1] = ROTATE(s, height, 180);
	tmp_edge[1][0] = ROTATE(s, height, 90);
	tmp_edge[1][1] = ROTATE(s, height, 270);
	//
	tmp_edge[2][0] = e;
	tmp_edge[2][1] = ROTATE(e, height, 180);
	tmp_edge[3][0] = ROTATE(e, height, 90);
	tmp_edge[3][1] = ROTATE(e, height, 270);
      }
      if(e != tmp_edge[2][1]) break;
    }
  }

  for(int i=0;i<groups;i++){
    edge[tmp_line[i]][0] = tmp_edge[i][0];
    edge[tmp_line[i]][1] = tmp_edge[i][1];
  }

  //  for(int i=0;i<groups;i++)
  //    printf("After : %d,%d-%d,%d\n",
  //	   WIDTH (edge[tmp_line[i]][0], height),
  //	   HEIGHT(edge[tmp_line[i]][0], height),
  //	   WIDTH (edge[tmp_line[i]][1], height),
  //	   HEIGHT(edge[tmp_line[i]][1], height));
    
  assert(check_symmetric_edge(lines, edge, height, 10, 5, groups));
  return true;
}
