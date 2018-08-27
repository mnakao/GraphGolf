#include "common.h"

static double uniform_rand()
{
  return ((double)random()+1.0)/((double)RAND_MAX+2.0);
}

static void print_result_header()
{
  PRINT_R0("   Times\t    Temp\tCur. ASPL GAP\t\tBest ASPL GAP\t\t");
  PRINT_R0("Cur. Dia. GAP\t\tBest Dia. GAP\tAccept Rate\n");
}

static void print_results(const int num, const double temp, const double current_ASPL, 
			  const double best_ASPL, const double low_ASPL, const int current_diam,
			  const int best_diam, const int low_diam, const int accepts, const int rejects)
{
  PRINT_R0("%8d\t%f\t", num, temp);
  PRINT_R0("%f ( %f )\t%f ( %f )\t%d ( %d )\t\t\t%d ( %d )\t\t",
	   current_ASPL, current_ASPL-low_ASPL, best_ASPL, best_ASPL-low_ASPL,
	   current_diam, current_diam-low_diam, best_diam, best_diam-low_diam);
  if(num != 0)
    PRINT_R0("%.4f ( %d / %d )\n", (double)accepts/(accepts+rejects), accepts, (accepts+rejects));
  else
    PRINT_R0("-\n");
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

bool check(const int nodes, const int based_nodes, const int lines,
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

static void edge_exchange(const int nodes, const int lines, const int groups, const int based_nodes,
			  int edge[lines][2], const int center_flag, const int ii)
{
  int line[groups*2], tmp_edge[groups*2][2];
  int based_lines = lines / groups;
  int center_vertex = nodes - 1;
  
  while(1){
    while(1){
      while(1){
	line[0] = getRandom(lines);
	line[1] = getRandom(lines);
	if(line[0] != line[1]) break;
      }
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
		   const bool hill_climbing_flag, const bool detect_temp_flag, const int i, double *max_diff_energy,
		   long long *total_accepts, int *accepts, int *rejects)
{
#if 0
  static double max = 100000;
  //  double tmp = fabs(((current_ASPL-ASPL)*nodes*(nodes-1))/2);
  double tmp = fabs(((current_ASPL-ASPL)*nodes*(nodes-1))/groups);
  if(max > tmp && tmp != 0){
    max = tmp;
    printf("%f\n", tmp);
    if(tmp == 1) exit(0);
  }
#endif

  if(ASPL <= current_ASPL){
    *accepts += 1;
    if(i > SKIP_ACCEPTS) *total_accepts +=1;
    return true;
  }
  if(hill_climbing_flag){ // Only accept when ASPL <= current_ASPL.
    *rejects += 1;
    return false;
  }

  double diff = (double)((current_ASPL-ASPL)*nodes*(nodes-1))/groups;

  if(detect_temp_flag)
    *max_diff_energy = MAX(*max_diff_energy, -1.0 * diff);

  if(exp(diff/temp) > uniform_rand()){
    *accepts += 1;
    if(i > SKIP_ACCEPTS) *total_accepts +=1;
    return true;
  }
  else{
    *rejects += 1;
    return false;
  }
}

long long sa(const int nodes, const int lines, const int degree, const int groups, double temp, 
	     const long long ncalcs, const double cooling_rate,  const int low_diam,  const double low_ASPL, 
	     const bool hill_climbing_flag, const bool detect_temp_flag, double *max_diff_energy,
	     int edge[lines][2], int *diam, double *ASPL, const int cooling_cycle, const int center_flag,
	     const int add_degree_to_center, const int based_nodes, long long *total_accepts)
{
  int current_edge[lines][2], best_edge[lines][2];
  long long i;
  int accepts = 0, rejects = 0;
  edge_copy((int *)best_edge, (int *)edge, lines*2);

  // Create adjacency matrix
  int (*adjacency)[degree] = malloc(sizeof(int)*nodes*degree); // int adjacency[nodes][degree];
  create_adjacency(nodes, lines, degree, edge, adjacency);

  evaluation(nodes, based_nodes, groups, lines, degree, adjacency, diam, ASPL, center_flag);
  double current_ASPL = *ASPL;
  double best_ASPL    = *ASPL;
  int current_diam    = *diam;
  int best_diam       = *diam;
  int print_interval  = (ncalcs/NUM_OF_PROGRESS == 0)? 1 : ncalcs/NUM_OF_PROGRESS;
  if(rank == 0 && !detect_temp_flag)
    print_result_header();

  for(i=0;i<ncalcs;i++){
    if(i % print_interval == 0 && !detect_temp_flag){
      print_results(i, temp, current_ASPL, best_ASPL, low_ASPL, 
		    current_diam, best_diam, low_diam, accepts, rejects);
      accepts = 0;
      rejects = 0;
    }
    
    while(1){
      edge_copy((int *)current_edge, (int *)edge, lines*2);
      edge_exchange(nodes, lines, groups, based_nodes, current_edge, center_flag, (int)i);
      assert(check(nodes, based_nodes, lines, degree, groups, current_edge, center_flag, add_degree_to_center, (int)i));
      create_adjacency(nodes, lines, degree, current_edge, adjacency);
      if(evaluation(nodes, based_nodes, groups, lines, degree, adjacency, diam, ASPL, center_flag)) break;
    }

    if(accept(*ASPL, current_ASPL, temp, nodes, groups, hill_climbing_flag,
	      detect_temp_flag, i, max_diff_energy, total_accepts, &accepts, &rejects)){
      current_ASPL = *ASPL;
      current_diam = *diam;
      edge_copy((int *)edge, (int *)current_edge, lines*2);
      if(best_ASPL > current_ASPL){
	edge_copy((int *)best_edge, (int *)current_edge, lines*2);
	best_ASPL = current_ASPL;
      }
      best_diam = MIN(best_diam, current_diam);
      if(best_ASPL == low_ASPL){
	if(!detect_temp_flag){
	  print_results(i, temp, current_ASPL, best_ASPL, low_ASPL, 
			current_diam, best_diam, low_diam, accepts, rejects);
	  PRINT_R0("---\nFound optimum solution.\n");
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
			     const int groups, int edge[lines][2], const int center_flag, const int add_degree_to_center)
{
  int diam;    // Not use
  double ASPL; // Not use
  int (*adjacency)[degree] = malloc(sizeof(int)*nodes*degree); // int adjacency[nodes][degree];
  int tmp_edge[lines][2];

  timer_start(TIMER_ESTIMATED);
  for(int i=0;i<ESTIMATED_TIMES;i++){
    edge_copy((int *)tmp_edge, (int *)edge, lines*2);
    edge_exchange(nodes, lines, groups, based_nodes, tmp_edge, center_flag, (int)i);
    assert(check(nodes, based_nodes, lines, degree, groups, tmp_edge, center_flag, add_degree_to_center, (int)i));
    create_adjacency(nodes, lines, degree, tmp_edge, adjacency);
    evaluation(nodes, based_nodes, groups, lines, degree, adjacency, &diam, &ASPL, center_flag);
  }
  timer_stop(TIMER_ESTIMATED);
  free(adjacency);

  return timer_read(TIMER_ESTIMATED)/ESTIMATED_TIMES;
}

// This function is mainly useful when groupe is 1.
void check_current_edge(const int nodes, const int degree, const int lines, const int groups,
			const int based_nodes, int edge[lines][2], const double low_ASPL, const int center_flag)
{
  int diam;    // Not use
  double ASPL;
  int (*adjacency)[degree] = malloc(sizeof(int)*nodes*degree); // int adjacency[nodes][degree];

  create_adjacency(nodes, lines, degree, edge, adjacency);
  if(! evaluation(nodes, based_nodes, groups, lines, degree, adjacency, &diam, &ASPL, center_flag))
    ERROR("The input file has a node which is never reached by another node.\n");

  if(ASPL == low_ASPL)
    END("The input file has already optimum solution.\n");

  free(adjacency);
}
