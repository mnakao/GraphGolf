#include "common.h"

static double uniform_rand()
{
  return ((double)random()+1.0)/((double)RAND_MAX+2.0);
}

static void print_result_header()
{
  printf("   Times\t    Temp\tCurrent ASPL (GAP)\tBest ASPL (GAP)\t\t");
  printf("Current Dia. (GAP)\tBest Dia. (GAP)\n");
}

static void print_results(const int num, const double temp, const double current_ASPL, 
			  const double best_ASPL, const double low_ASPL,
			  const int current_diam, const int best_diam, const int low_diam)
{
  printf("%8d\t%f\t", num, temp);
  printf("%f ( %f )\t%f ( %f )\t%d ( %d )\t\t\t%d ( %d )\n",
	 current_ASPL, current_ASPL-low_ASPL, best_ASPL, best_ASPL-low_ASPL,
	 current_diam, current_diam-low_diam, best_diam, best_diam-low_diam);
}  

void create_adjacency(const int nodes, const int lines, const int degree, 
		      int edge[lines][2], int adjacency[nodes][degree])
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

#define CENTER_VERTEX -1
int distance(const int nodes, const int a, const int b, const int center_flag, const int center_vertex)
{
  if(center_flag && (a == center_vertex || b == center_vertex)) return CENTER_VERTEX;
  int v = MAX(a, b) - MIN(a, b);

  if(center_flag){
    //    double m = (nodes-1.0)/2.0;
    //    return (v < m)? v : nodes-1-v
    return (v < nodes/2)? v : nodes-1-v;
  }
  else{
    //    return (v < nodes/2)? v : nodes-v;
    int max  = MAX(a, b);
    int min  = MIN(a, b);
    int tmp1 = max - min;
    int tmp2 = min + nodes - max;
    return MIN(tmp1, tmp2);
  }
}

bool check(const int rank, const int nodes, const int based_nodes, const int lines,
	   const int degree, const int groups, int edge[lines][2], 
	   const int center_flag, const int add_degree_to_center, const int ii)
{
  int based_lines = lines/groups;
  int center_vertex = nodes - 1;

  //  printf("--\n");
  //  for(int i=0;i<lines;i++)
  //     printf("%d %d\n", edge[i][0], edge[i][1]);
    /*
    if(order(nodes, edge[i][0], edge[i][1], center_flag) == RIGHT)
      printf("%d %d : RIGHT\n", edge[i][0], edge[i][1]);
    else if(order(nodes, edge[i][0], edge[i][1], center_flag) == LEFT)
      printf("%d %d : LEFT\n", edge[i][0], edge[i][1]);
    else
    printf("%d %d : MIDDLE\n", edge[i][0], edge[i][1]);*/
  
  for(int i=0;i<based_lines;i++){
    for(int j=1;j<groups;j++){
      int k = j * based_lines + i;
      if(distance(nodes, edge[i][0], edge[i][1], center_flag, center_vertex) != distance(nodes, edge[k][0], edge[k][1], center_flag, center_vertex)){
	PRINT_R0("check 1: %d\n", ii);
	PRINT_R0("edge[%d][0] = %d : edge[%d][1] = %d d=%d\n", i, edge[i][0], i, edge[i][1], distance(nodes, edge[i][0], edge[i][1], center_flag, center_vertex));
	PRINT_R0("edge[%d][0] = %d : edge[%d][1] = %d d=%d\n", k, edge[k][0], k, edge[k][1], distance(nodes, edge[k][0], edge[k][1], center_flag, center_vertex));
	return false;
      }
    }
  }

  for(int i=0;i<based_lines;i++){
    for(int j=1;j<groups;j++){
      int k = j * based_lines + i;
      if(order(nodes, edge[i][0], edge[i][1], center_flag) != order(nodes, edge[k][0], edge[k][1], center_flag)){
	PRINT_R0("check 2 : %d\n", ii);
	PRINT_R0("edge[%d][0] = %d : edge[%d][1] = %d %d\n", i, edge[i][0], i, edge[i][1], order(nodes, edge[i][0], edge[i][1], center_flag));
	PRINT_R0("edge[%d][0] = %d : edge[%d][1] = %d %d\n", k, edge[k][0], k, edge[k][1], order(nodes, edge[k][0], edge[k][1], center_flag));
	return false;
      }
    }
  }

  if(!check_duplicate_edge(lines, edge)){
    PRINT_R0("check 3\n");
    return false;
  }
  
  for(int i=0;i<based_lines;i++){
    if(order(nodes, edge[i][0], edge[i][1], center_flag) != MIDDLE)
      for(int j=1;j<groups;j++){
	int k = j * based_lines + i;
	int tmp0 = edge[k][0] - edge[k-based_lines][0];
	int tmp1 = edge[k][1] - edge[k-based_lines][1];
	if(center_flag){
	  tmp0 = (tmp0 < 0)? tmp0+nodes-1 : tmp0;
	  tmp1 = (tmp1 < 0)? tmp1+nodes-1 : tmp1;
	}
	else{
	  tmp0 = (tmp0 < 0)? tmp0 + nodes : tmp0;
	  tmp1 = (tmp1 < 0)? tmp1 + nodes : tmp1;
	}
	if(tmp0 != based_nodes || tmp1 != based_nodes){
	  PRINT_R0("check 4: %d\n", ii);
	  PRINT_R0("The different group relationship\n");
	  PRINT_R0("edge[%d][0]-edge[%d][0] = %d - %d = %d\n", k, k-based_lines, edge[k][0], edge[k-based_lines][0], tmp0);
	  PRINT_R0("edge[%d][1]-edge[%d][1] = %d - %d = %d\n", k, k-based_lines, edge[k][1], edge[k-based_lines][1], tmp1);
	  return false;
	}
      }
  }

  return true;
}

bool has_duplicated_vertex(const int e00, const int e01, const int e10, const int e11)
{
  return (e00 == e10 || e01 == e11 || e00 == e11 || e01 == e10);
}

static void set_lines(const int opt, int line[2], const int lines, const int groups, const int based_nodes, 
		      int edge[lines][2], const int total_distance[based_nodes])
{
  if(opt == 0 || getRandom(2) == 0){
    while(1){
      line[0] = getRandom(lines);
      line[1] = getRandom(lines);
      if(line[0] != line[1]) return;
    }
  }

  int based_lines = lines/groups, tmp_edge[based_lines], num = 0;
  double priority[based_lines], max = -1.0;

#if 0
  for(int i=0;i<based_nodes;i++)
    printf("%d : %d\n", i, total_distance[i]);
  printf("--\n");
  for(int i=0;i<based_lines;i++){
    int v0 = edge[i][0] % based_nodes;
    int v1 = edge[i][1] % based_nodes;
    printf("%2d : %2d %2d : %d %d : %d\n", i, v0, v1,
           total_distance[v0], total_distance[v1], 
	   total_distance[v0] + total_distance[v1]);
  }
#endif

  for(int i=0;i<based_lines;i++){
    int v0 = edge[i][0] % based_nodes;
    int v1 = edge[i][1] % based_nodes;
    priority[i] = total_distance[v0] + total_distance[v1];
  }

  for(int i=0;i<based_lines;i++){
    if(priority[i] > max){
      max = priority[i];
      tmp_edge[0] = i;
      num = 1;
    }
    else if(priority[i] == max)
      tmp_edge[num++] = i;
  }

  if(num >= 2){
    while(1){
      line[0] = tmp_edge[getRandom(num)] + based_lines * getRandom(groups);
      line[1] = tmp_edge[getRandom(num)] + based_lines * getRandom(groups);
      if(line[0] != line[1]) return;
    }
  }

  // num == 1
  max = -1.0;
  for(int i=0;i<based_lines;i++){
    if(i != tmp_edge[0]){
      if(priority[i] > max){
	max = priority[i];
	tmp_edge[1] = i;
	num = 2;
      }
      else if(priority[i] == max)
	tmp_edge[num++] = i;
    }
  }

  while(1){
    line[0] = tmp_edge[0] + based_lines * getRandom(groups);
    line[1] = tmp_edge[getRandom(num-1)+1] + based_lines * getRandom(groups);
    if(line[0] != line[1]) return;
  }
}

static void edge_exchange(const int nodes, const int lines, const int groups, const int based_nodes,
			  int edge[lines][2], const int total_distance[based_nodes], const int opt,
			  const int center_flag, const int ii)
{
  int line[groups*2], tmp_edge[groups*2][2];
  int based_lines = lines / groups;
  int center_vertex = nodes - 1;
  
  while(1){
    while(1){
      set_lines(opt, line, lines, groups, based_nodes, edge, total_distance);
      if((line[0] - line[1]) % based_lines == 0 || 
	 has_duplicated_vertex(edge[line[0]][0], edge[line[0]][1],
			       edge[line[1]][0], edge[line[1]][1])){
	if(groups == 1) continue;
	else{
	  if(edge_1g_opt(edge, nodes, based_nodes, based_lines, groups, line[0], center_flag)) return;
	  else continue;
	}
      }
      else break;
    }

    bool flag0 = (distance(nodes, edge[line[0]][0], edge[line[0]][1], center_flag, center_vertex) == nodes/2);
    bool flag1 = (distance(nodes, edge[line[1]][0], edge[line[1]][1], center_flag, center_vertex) == nodes/2);
    bool diameter_flag = ((flag0 || flag1) && groups%2 == 0);

    if(diameter_flag){
      if(edge_1g_opt(edge, nodes, based_nodes, based_lines, groups, line[0], center_flag)) return;
      else continue;
    }
    else{ // 2g_opt
      for(int i=1;i<groups;i++){
	int tmp0 = line[0] + based_lines * i;
	int tmp1 = line[1] + based_lines * i;
	line[0+2*i] = (tmp0 >= lines)? tmp0 - lines : tmp0;
	line[1+2*i] = (tmp1 >= lines)? tmp1 - lines : tmp1;
      }

      for(int i=0;i<groups*2;i++)
	for(int j=0;j<2;j++)
	  tmp_edge[i][j] = edge[line[i]][j];

      if(getRandom(2) == 0){
	for(int i=0;i<groups;i++)
	  swap(&tmp_edge[i*2][1], &tmp_edge[i*2+1][1]);
      }
      else{
	for(int i=0;i<groups;i++)
	  swap(&tmp_edge[i*2][1], &tmp_edge[i*2+1][0]);
      }

      assert(check_loop(groups*2, tmp_edge));
      if(!check_duplicate_edge(groups*2, tmp_edge)) continue;
      if(!check_duplicate_current_edge(lines, groups*2, line, edge, tmp_edge, groups, nodes, center_flag))
	continue;

      for(int i=0;i<groups*2;i++)
	if(order(nodes, tmp_edge[i][0], tmp_edge[i][1], center_flag) == RIGHT)
	  swap(&tmp_edge[i][0], &tmp_edge[i][1]); // RIGHT -> LEFT

      for(int i=0;i<groups*2;i++){
	edge[line[i]][0] = tmp_edge[i][0];
	edge[line[i]][1] = tmp_edge[i][1];
      }

      break;
    }
  }
}

static bool accept(const double ASPL, const double current_ASPL, const double temp, const int nodes, const int groups,
		   const bool hill_climbing_flag, const bool detect_temp_flag, double *max_diff_energy)
{
#if 0
  static double max = 100000;
  double tmp = fabs(((current_ASPL-ASPL)*nodes*(nodes-1))/2);
  if(max > tmp && tmp != 0){
    max = tmp;
    printf("%f\n", tmp);
    if(tmp == 1) exit(0);
  }
#endif

  if(ASPL <= current_ASPL) return true;
  if(hill_climbing_flag)   return false; // Only accept when ASPL <= current_ASPL.

  double diff = (double)((current_ASPL-ASPL)*nodes*(nodes-1))/groups;

  if(detect_temp_flag)
    *max_diff_energy = MAX(*max_diff_energy, -1.0 * diff);

  return (exp(diff/temp) > uniform_rand())? true : false;
}

long long sa(const int nodes, const int lines, const int degree, const int groups, double temp, 
	     const long long ncalcs, const double cooling_rate,  const int low_diam,  const double low_ASPL, 
	     const bool hill_climbing_flag, const bool detect_temp_flag, double *max_diff_energy,
	     int edge[lines][2], int *diam, double *ASPL, const int rank, const int size, const int opt, const int cooling_cycle,
	     const int center_flag, const int add_degree_to_center, const int based_nodes)
{
  int current_edge[lines][2], best_edge[lines][2], total_distance[based_nodes];
  long long i;
  edge_copy((int *)best_edge, (int *)edge, lines*2);

  // Create adjacency matrix
  int (*adjacency)[degree] = malloc(sizeof(int)*nodes*degree); // int adjacency[nodes][degree];
  create_adjacency(nodes, lines, degree, edge, adjacency);

  evaluation(nodes, based_nodes, groups, lines, degree, adjacency, diam, ASPL, total_distance, rank, size, opt, center_flag);
  double current_ASPL = *ASPL;
  double best_ASPL    = *ASPL;
  int current_diam    = *diam;
  int best_diam       = *diam;
  int print_interval  = (ncalcs/NUM_OF_PROGRESS == 0)? 1 : ncalcs/NUM_OF_PROGRESS;
  if(rank == 0 && !detect_temp_flag)
    print_result_header();

  for(i=0;i<ncalcs;i++){
    if(i % print_interval == 0 && rank == 0 && !detect_temp_flag)
      print_results(i, temp, current_ASPL, best_ASPL, low_ASPL, 
		    current_diam, best_diam, low_diam);
    while(1){
      edge_copy((int *)current_edge, (int *)edge, lines*2);
      edge_exchange(nodes, lines, groups, based_nodes, current_edge, total_distance, opt, center_flag, (int)i);
      assert(check(rank, nodes, based_nodes, lines, degree, groups, current_edge, center_flag, add_degree_to_center, (int)i));
      create_adjacency(nodes, lines, degree, current_edge, adjacency);
      if(evaluation(nodes, based_nodes, groups, lines, degree, adjacency, diam, ASPL, total_distance, rank, size, opt, center_flag)) break;
    }

    if(accept(*ASPL, current_ASPL, temp, nodes, groups, 
	      hill_climbing_flag, detect_temp_flag, max_diff_energy)){
      current_ASPL = *ASPL;
      current_diam = *diam;
      edge_copy((int *)edge, (int *)current_edge, lines*2);
      if(best_ASPL > current_ASPL){
	edge_copy((int *)best_edge, (int *)current_edge, lines*2);
	best_ASPL = current_ASPL;
      }
      best_diam = MIN(best_diam, current_diam);
      if(best_ASPL == low_ASPL){
	if(rank == 0 && !detect_temp_flag){
	  print_results(i, temp, current_ASPL, best_ASPL, low_ASPL, 
			current_diam, best_diam, low_diam);
	  printf("---\nFound optimum solution.\n");
	}
	break;
      }
    }

    if((i+1)%cooling_cycle == 0)
      temp *= cooling_rate;
  }

  *ASPL = best_ASPL;
  *diam = best_diam;
  edge_copy((int *)edge, (int *)best_edge, lines*2);
  free(adjacency);

  return i;
}

#define ESTIMATED_TIMES 5
double estimated_elapse_time(const long long ncals, const int nodes, const int based_nodes, const int lines, const int degree,
			     const int groups, int edge[lines][2], const int rank, const int size, const int opt,
			     const int center_flag, const int add_degree_to_center)
{
  int diam;    // Not use
  double ASPL; // Not use
  int (*adjacency)[degree] = malloc(sizeof(int)*nodes*degree); // int adjacency[nodes][degree];
  int total_distance[based_nodes], tmp_edge[lines][2];

  timer_start(TIMER_ESTIMATED);
  for(int i=0;i<ESTIMATED_TIMES;i++){
    edge_copy((int *)tmp_edge, (int *)edge, lines*2);
    edge_exchange(nodes, lines, groups, based_nodes, tmp_edge, total_distance, opt, center_flag, (int)i);
    assert(check(rank, nodes, based_nodes, lines, degree, groups, tmp_edge, center_flag, add_degree_to_center, (int)i));
    create_adjacency(nodes, lines, degree, tmp_edge, adjacency);
    evaluation(nodes, based_nodes, groups, lines, degree, adjacency, &diam, &ASPL, total_distance, rank, size, opt, center_flag);
  }
  timer_stop(TIMER_ESTIMATED);
  free(adjacency);

  return timer_read(TIMER_ESTIMATED)/ESTIMATED_TIMES;
}

// This function is mainly useful when groupe is 1.
void check_current_edge(const int nodes, const int degree, const int lines, const int groups,
			const int based_nodes, int edge[lines][2], const double low_ASPL,
			const int rank, const int size, const int center_flag)
{
  int diam;    // Not use
  double ASPL;
  int (*adjacency)[degree] = malloc(sizeof(int)*nodes*degree); // int adjacency[nodes][degree];
  int total_distance[based_nodes];

  create_adjacency(nodes, lines, degree, edge, adjacency);
  if(! evaluation(nodes, based_nodes, groups, lines, degree, adjacency, &diam, &ASPL, total_distance, rank, size, 0, center_flag))
    ERROR("The input file has a node which is never reached by another node.\n");

  if(ASPL == low_ASPL)
    END("The input file has already optimum solution.\n");

  free(adjacency);
}
