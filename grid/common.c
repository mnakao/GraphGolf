#include "common.h"

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

int ROTATE(const int v, const int height, const int width,
	   const int groups, const int degree)
{
  if(groups != 2 && groups != 4)
    ERROR("Invalid groups\n");

  int w = WIDTH (v, height);
  int h = HEIGHT(v, height);
  if(groups == 2){
    if(degree != 180)
       ERROR("Invalid degree\n");
    
    // degree == 180
    return (width-w-1)*height + (height-h-1);
  }
  else{ // groups == 4
    if(degree != 90 && degree != 180 && degree != 270)
      ERROR("Invalid degree\n");
    
    if(degree == 90)       return h*height + (height-w-1);
    else if(degree == 180) return (height-w-1)*height + (height-h-1);
    else                   return (height-h-1)*height + w; // degree == 270
  }
}

bool check_degree(const int nodes, const int lines, int edge[lines][2])
{
  int n[nodes];
  for(int i=0;i<nodes;i++)
    n[i] = 0;
  
  for(int i=0;i<lines;i++){
    n[edge[i][0]]++;
    n[edge[i][1]]++;
  }
  
  int degree = 2 * lines / nodes;
  for(int i=0;i<nodes;i++)
    if(degree != n[i])
      return false;

  return true;
}

/*
// diameter is not implemented
bool check_vector(const int groups, const int lines, const int height, const int edge[lines][2])
{
  if(groups == 1) return true;
  
  int based_lines = lines/groups;
  int vec_w[lines], vec_h[lines];
  for(int i=0;i<lines;i++){
    vec_w[i] = WIDTH (edge[i][0], height) - WIDTH (edge[i][1], height);
    vec_h[i] = HEIGHT(edge[i][0], height) - HEIGHT(edge[i][1], height);
  }

  if(groups == 2){
    for(int i=0;i<based_lines;i++)
      if(!(vec_w[i] == -1 * vec_w[based_lines+i] && vec_h[i] == -1 * vec_h[based_lines+i]))
	return false;
  }
  else{ // groups == 4
    for(int i=0;i<based_lines;i++){
      if(!(vec_w[i] == -1 * vec_h[based_lines+i] && vec_h[i] == vec_w[based_lines+i])){
	printf("%d,%d %d,%d\n", WIDTH (edge[i][0], height),HEIGHT(edge[i][0], height),WIDTH (edge[i][1], height),HEIGHT(edge[i][1], height));
	printf("%d,%d %d,%d\n", WIDTH (edge[based_lines+i][0], height),HEIGHT(edge[based_lines+i][0], height),
	       WIDTH (edge[based_lines+i][1], height),HEIGHT(edge[based_lines+i][1], height));
	printf("A\n");
        return false;
      }
      if(!(vec_w[i] == -1 * vec_w[based_lines*2+i] && vec_h[i] == -1 * vec_h[based_lines*2+i])){
	printf("%d,%d %d,%d\n", WIDTH (edge[i][0], height),HEIGHT(edge[i][0], height),WIDTH (edge[i][1], height),HEIGHT(edge[i][1], height));
	printf("%d,%d %d,%d\n", WIDTH (edge[based_lines*2+i][0], height),HEIGHT(edge[based_lines*2+i][0], height),
               WIDTH (edge[based_lines*2+i][1], height),HEIGHT(edge[based_lines*2+i][1], height));
	printf("B\n");
	return false;
      }
      if(!(vec_w[i] == vec_h[based_lines*3+i] && vec_h[i] == -1 * vec_w[based_lines*3+i])){
	printf("%d %d : %d %d\n", vec_w[i], vec_w[based_lines+i], vec_h[i], vec_h[based_lines+i]);
	printf("C\n");
        return false;
      }
    }
  }
  
  return true;
}
*/

bool check_symmetric_edge(const int lines, const int edge[lines][2], const int height,
			  const int width, const int based_height, const int groups)
{
  assert(lines%groups == 0);
  int tmp_edge[2], based_lines = lines / groups;

  if(groups == 2){
    for(int i=0;i<based_lines;i++){
      for(int j=0;j<2;j++)
	tmp_edge[j] = ROTATE(edge[i][j], height, width, groups, 180);
      
      if(!has_duplicated_edge(edge[based_lines+i][0], edge[based_lines+i][1], tmp_edge[0], tmp_edge[1]))
	if(!( WIDTH (edge[based_lines+i][0], height) + WIDTH (edge[based_lines+i][1], height) == (width-1) &&
              HEIGHT(edge[based_lines+i][0], height) + HEIGHT(edge[based_lines+i][1], height) == (height-1))){
	  printf("i=%d: %d,%d-%d,%d %d,%d-%d,%d\n", i,
                 WIDTH(edge[based_lines+i][0], height), HEIGHT(edge[based_lines+i][0], height),
                 WIDTH(edge[based_lines+i][1], height), HEIGHT(edge[based_lines+i][1], height),
                 WIDTH(tmp_edge[0], height), HEIGHT(tmp_edge[0], height),
                 WIDTH(tmp_edge[1], height), HEIGHT(tmp_edge[1], height));
	  return false;
	}
    }
  }
  else if(groups == 4){
    // 90 degrees
    for(int i=0;i<based_lines;i++){
      for(int j=0;j<2;j++)
	tmp_edge[j] = ROTATE(edge[i][j], height, width, groups, 90);
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
	tmp_edge[j] = ROTATE(edge[i][j], height, width, groups, 180);
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
	tmp_edge[j] = ROTATE(edge[i][j], height, width, groups, 270);
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
    printf("%d,%d %d,%d\n", WIDTH(edge[i][0], height), HEIGHT(edge[i][0], height),
	   WIDTH(edge[i][1], height), HEIGHT(edge[i][1], height));
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

bool edge_1g_opt(int (*edge)[2], const int nodes, const int lines, const int degree, const int based_nodes,
		 const int based_lines, const int height, const int width, const int groups, const int start_line,
		 const long long ii)
{
  assert(groups != 1);
  
  int tmp_line[groups], tmp_edge[groups][2], pattern;
  for(int i=0;i<groups;i++){
    int t = start_line + based_lines * i;
    tmp_line[i] = (t < lines)? t : t - lines;
  }
#ifdef _DEBUG_MSG
  printf("edge_1g_opt: ii = %lld\n", ii);
  for(int i=0;i<groups;i++)
    printf("line[%d] = %d\n", i, tmp_line[i]);

  for(int i=0;i<groups;i++)
    printf("Before: %d,%d-%d,%d\n",
	   WIDTH (edge[tmp_line[i]][0], height),
  	   HEIGHT(edge[tmp_line[i]][0], height),
  	   WIDTH (edge[tmp_line[i]][1], height),
  	   HEIGHT(edge[tmp_line[i]][1], height));
#endif

  if(edge[tmp_line[0]][0] == edge[tmp_line[groups-1]][1] ||
     edge[tmp_line[0]][1] == edge[tmp_line[groups-1]][0])  // A cycle is composed of four edges.
    return false;
  
  int s = edge[tmp_line[0]][0];
  int e = edge[tmp_line[groups/2]][1];
  while(1){
    pattern = getRandom(groups+1);
    if(groups == 2){
      if(pattern != groups){
	tmp_edge[0][0] = s;
	tmp_edge[1][0] = ROTATE(s, height, width, groups, 180);

	if(pattern == 0){
	  tmp_edge[0][1] = e;
	  tmp_edge[1][1] = ROTATE(e, height, width, groups, 180);
	}
	else{ // pattern == 1
	  tmp_edge[0][1] = ROTATE(e, height, width, groups, 180);
	  tmp_edge[1][1] = e;
	}
      }
      else{ // pattern == groups
	tmp_edge[0][0] = s;
	tmp_edge[0][1] = ROTATE(s, height, width, groups, 180);
	tmp_edge[1][0] = ROTATE(e, height, width, groups, 180);
	tmp_edge[1][1] = e;
      }
    }
    else{ // groups == 4
      if(pattern != groups){
	tmp_edge[0][0] = s;
	tmp_edge[1][0] = ROTATE(s, height, width, groups, 90);
	tmp_edge[2][0] = ROTATE(s, height, width, groups, 180);
	tmp_edge[3][0] = ROTATE(s, height, width, groups, 270);
	
	if(pattern == 0){
	  tmp_edge[0][1] = e;
	  tmp_edge[1][1] = ROTATE(e, height, width, groups, 90);
	  tmp_edge[2][1] = ROTATE(e, height, width, groups, 180);
	  tmp_edge[3][1] = ROTATE(e, height, width, groups, 270);
	}
	else if(pattern == 1){
	  tmp_edge[0][1] = ROTATE(e, height, width, groups, 270);
	  tmp_edge[1][1] = e;
	  tmp_edge[2][1] = ROTATE(e, height, width, groups, 90);
	  tmp_edge[3][1] = ROTATE(e, height, width, groups, 180);
	}
	else if(pattern == 2){
	  tmp_edge[0][1] = ROTATE(e, height, width, groups, 180);
	  tmp_edge[1][1] = ROTATE(e, height, width, groups, 270);
	  tmp_edge[2][1] = e;
	  tmp_edge[3][1] = ROTATE(e, height, width, groups, 90);
	}
	else{ // pattern == 3
	  tmp_edge[0][1] = ROTATE(e, height, width, groups, 90);
	  tmp_edge[1][1] = ROTATE(e, height, width, groups, 180);
	  tmp_edge[2][1] = ROTATE(e, height, width, groups, 270);
	  tmp_edge[3][1] = e;
	}
      }
      else{ // pattern == groups
	tmp_edge[0][0] = s;
	tmp_edge[0][1] = ROTATE(s, height, width, groups, 180);
	tmp_edge[1][0] = ROTATE(s, height, width, groups, 90);
	tmp_edge[1][1] = ROTATE(s, height, width, groups, 270);
	//
	tmp_edge[2][0] = e;
	tmp_edge[2][1] = ROTATE(e, height, width, groups, 180);
	tmp_edge[3][0] = ROTATE(e, height, width, groups, 90);
	tmp_edge[3][1] = ROTATE(e, height, width, groups, 270);
      }
    }
    if(e != tmp_edge[groups/2][1]) break;
  }

  /*
  if(groups == 2){
    int vec_w[groups], vec_h[groups];
    for(int i=0;i<groups;i++){
      vec_w[i] = WIDTH (tmp_edge[i][0], height) - WIDTH (tmp_edge[i][1], height);
      vec_h[i] = HEIGHT(tmp_edge[i][0], height) - HEIGHT(tmp_edge[i][1], height);
    }

    for(int i=1;i<groups;i++){
      if(vec_w[0] == vec_w[i] && vec_h[0] == vec_h[i]){
	swap(&tmp_edge[i][0], &tmp_edge[i][1]);
	printf("AAA\n");
      }
    }
  }
  */

  if(!check_duplicate_current_edge(lines, edge, groups, tmp_edge, tmp_line, groups, D_1G_OPT, (pattern==groups)))
     return false;

  for(int i=0;i<groups;i++){
    edge[tmp_line[i]][0] = tmp_edge[i][0];
    edge[tmp_line[i]][1] = tmp_edge[i][1];
  }
#ifdef _DEBUG_MSG
  for(int i=0;i<groups;i++)
    printf("After : %d,%d-%d,%d\n",
  	   WIDTH (edge[tmp_line[i]][0], height),
  	   HEIGHT(edge[tmp_line[i]][0], height),
  	   WIDTH (edge[tmp_line[i]][1], height),
  	   HEIGHT(edge[tmp_line[i]][1], height));
#endif

  return true;
}
