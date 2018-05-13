#include "common.h"

int getRandom(const int max)
{
  return (int)(rand()*((double)max)/(1.0+RAND_MAX));
}

static int get_end_edge(const int start_edge, const int groups, int (*edge)[2], 
			int line[groups], const int based_nodes)
{
  for(int i=0;i<groups;i++)
    for(int j=0;j<2;j++)
      if(edge[line[i]][j] % based_nodes != start_edge){
	return edge[line[i]][j] % based_nodes;
      }

  return -1;
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

int order(const int nodes, const int a, const int b)
{
  if((a-b)%(nodes/2) == 0) return MIDDLE;

  if(a < nodes/2){
    if(a > b) return LEFT;
    return (a+nodes/2 > b)? RIGHT : LEFT;
  }
  else{
    if(a < b) return RIGHT;
    return (0 <= b && b < a - nodes/2)? RIGHT : LEFT;
  }
}


bool check_loop(const int lines, int (*edge)[2])
{
  for(int i=0;i<lines;i++)
    if(edge[i][0] == edge[i][1])
      return false;

  return true;
}

bool check_duplicate_edge(const int lines, int (*edge)[2])
{
  for(int i=0;i<lines;i++)
    for(int j=i+1;j<lines;j++)
      if(has_duplicated_edge(edge[i][0], edge[i][1], edge[j][0], edge[j][1]))
        return false;

  return true;
}

static bool target_line(const int i, const int groups, const int line[groups])
{
  for(int j=0;j<groups;j++)
    if(i == line[j])
      return false;

  return true;
}

bool check_duplicate_current_edge(const int lines, const int groups, const int line[groups],
                                  int (*edge)[2], int tmp_edge[groups][2])
{
  for(int i=0;i<lines;i++)
    if(target_line(i, groups, line))
      for(int j=0;j<groups;j++)
        if(has_duplicated_edge(edge[i][0], edge[i][1], tmp_edge[j][0], tmp_edge[j][1]))
          return false;

  return true;
}

bool edge_exchange_among_groups(const int based_nodes, const int based_lines, int (*edge)[2], 
				const int groups, const int start_line)
{
  if(groups == 1)
    return true;

  int line[groups], tmp_edge[groups][2];
  int nodes = based_nodes * groups;
  int lines = based_lines * groups;

  for(int i=0;i<groups;i++){
    line[i] = start_line + i * based_lines;
    assert(0 <= line[i] && line[i] < lines);
  }

  int start_edge = edge[line[0]][0] % based_nodes;
  int end_edge   = get_end_edge(start_edge, groups, edge, line, based_nodes);
  if(end_edge == -1) return false;
  int diff       = edge[line[0]][0] - edge[line[0]][1];

  while(1){
    int pattern = getRandom(groups+1);

    if(pattern == groups){
      for(int i=0;i<groups/2;i++){
	tmp_edge[i][0] = edge[line[i]][0];
	int tmp = edge[line[i]][0] + (groups/2) * based_nodes;
	tmp_edge[i][1] = (tmp >= nodes)? tmp - nodes : tmp;
      }
      for(int i=groups/2;i<groups;i++){
	tmp_edge[i][0] = edge[line[i]][1];
	int tmp = edge[line[i]][1] - (groups/2) * based_nodes;
	tmp_edge[i][1] = (tmp < 0)? tmp + nodes : tmp;
      }
    }
    else{
      for(int i=0;i<groups;i++)
	tmp_edge[i][0] = start_edge + based_nodes * i;
      
      tmp_edge[0][1] = end_edge + based_nodes * pattern;
      for(int i=1;i<groups;i++){
	int tmp = tmp_edge[0][1] + based_nodes * i;
	tmp_edge[i][1] = (tmp >= nodes)? tmp - nodes : tmp;
      }
    }

    if(diff != tmp_edge[0][0] - tmp_edge[0][1]) break;
  }

  if(!check_loop(groups, tmp_edge))             return false;
  if(!check_duplicate_edge(groups, tmp_edge))   return false;
  if(!check_duplicate_current_edge(lines, groups, line, edge, tmp_edge))
    return false;

  for(int i=0;i<groups;i++)
    if(order(nodes, tmp_edge[i][0], tmp_edge[i][1]) == RIGHT)
      swap(&tmp_edge[i][0], &tmp_edge[i][1]);  // RIGHT -> LEFT
  
  // Set vertexs
  for(int i=0;i<groups;i++){
    edge[line[i]][0] = tmp_edge[i][0];
    edge[line[i]][1] = tmp_edge[i][1];
  }

  return true;
}
