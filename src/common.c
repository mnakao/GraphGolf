#include "common.h"

int getRandom(const int max)
{
  return (int)(rand()*((double)max)/(1.0+RAND_MAX));
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

  int patterns = (groups%2==0)? groups : groups - 1;
  int tmp_edge[groups][2];
  int nodes = based_nodes * groups;
  int lines = based_lines * groups;
  int line[groups];
  line[0] = start_line;
  for(int i=1;i<groups;i++)
    line[i] = line[0] + i * based_lines;

  int pattern = getRandom(patterns);
  if(groups%2 == 1) pattern++;

  if(pattern == 0){
    for(int i=0;i<groups/2;i++){
      tmp_edge[i][0] = edge[line[i]][0];
      tmp_edge[i][1] = edge[line[i]][0] + (groups/2) * based_nodes;
    }
    for(int i=groups/2;i<groups;i++){
      tmp_edge[i][0] = edge[line[i]][1];
      tmp_edge[i][1] = edge[line[i]][1] - (groups/2) * based_nodes;
    }
  }
  else{
    for(int i=0;i<groups;i++){
      tmp_edge[i][0] = edge[line[i]][0];
      tmp_edge[i][1] = edge[line[i]][1] + based_nodes * pattern;
    }
  }

  for(int i=0;i<groups;i++){
    if(tmp_edge[i][1] < 0)      tmp_edge[i][1] += nodes;
    if(tmp_edge[i][1] >= nodes) tmp_edge[i][1] -= nodes;
  }

  if(!check_loop(groups, tmp_edge))             return false;
  if(!check_duplicate_edge(groups, tmp_edge))   return false;
  if(!check_duplicate_current_edge(lines, groups, line, edge, tmp_edge))
    return false;

  // Set vertexs
  for(int i=0;i<groups;i++){
    edge[line[i]][0] = tmp_edge[i][0];
    edge[line[i]][1] = tmp_edge[i][1];
  }

  return true;
}
