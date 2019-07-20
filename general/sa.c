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

static void print_results(const long long num, const double temp, const double current_ASPL, 
			  const double best_ASPL, const double low_ASPL, const int current_diam,
			  const int best_diam, const int low_diam, const long long accepts, const long long rejects)
{
  PRINT_R0("%8lld\t%f\t", num, temp);
  PRINT_R0("%f ( %f )\t%f ( %f )\t%d ( %d )\t\t\t%d ( %d )\t\t",
	   current_ASPL, current_ASPL-low_ASPL, best_ASPL, best_ASPL-low_ASPL,
	   current_diam, current_diam-low_diam, best_diam, best_diam-low_diam);
  if(num != 0)
    PRINT_R0("%.4f ( %lld / %lld )\n", (double)accepts/(accepts+rejects), accepts, (accepts+rejects));
  else
    PRINT_R0("-\n");
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

#define CENTER_VERTEX -1
int distance(int nodes, const int a, const int b, const int added_centers)
{
  if(a >= nodes-added_centers || b >= nodes-added_centers) return CENTER_VERTEX;

  int v = MAX(a, b) - MIN(a, b);
  if(added_centers) nodes -= added_centers;
  return (v < nodes/2.0)? v : nodes-v;
}

bool check(const int nodes, const int based_nodes, const int lines, const int degree, const int groups,
	   int edge[lines][2], const int added_centers, const int adjacency[nodes][degree], const int ii)
{
  bool flag = true;
  int based_lines = lines/groups;

#pragma omp parallel for
  for(int i=0;i<based_lines;i++){
    for(int j=1;j<groups;j++){
      int k = j * based_lines + i;
      if(distance(nodes, edge[i][0], edge[i][1], added_centers) != distance(nodes, edge[k][0], edge[k][1], added_centers)){
	PRINT_R0("check 1: %d\n", ii);
	PRINT_R0("edge[%d][0] = %d : edge[%d][1] = %d d=%d\n", i, edge[i][0], i, edge[i][1], distance(nodes, edge[i][0], edge[i][1], added_centers));
	PRINT_R0("edge[%d][0] = %d : edge[%d][1] = %d d=%d\n", k, edge[k][0], k, edge[k][1], distance(nodes, edge[k][0], edge[k][1], added_centers));
	flag = false;
      }
    }
  }
  
#pragma omp parallel for
  for(int i=0;i<based_lines;i++){
    for(int j=1;j<groups;j++){
      int k = j * based_lines + i;
      if(order(nodes, edge[i][0], edge[i][1], added_centers) != order(nodes, edge[k][0], edge[k][1], added_centers)){
	PRINT_R0("check 2 : %d\n", ii);
	PRINT_R0("edge[%d][0] = %d : edge[%d][1] = %d %d\n", i, edge[i][0], i, edge[i][1], order(nodes, edge[i][0], edge[i][1], added_centers));
	PRINT_R0("edge[%d][0] = %d : edge[%d][1] = %d %d\n", k, edge[k][0], k, edge[k][1], order(nodes, edge[k][0], edge[k][1], added_centers));
	flag = false;
      }
    }
  }

  if(!check_duplicate_edge(lines, edge)){
    PRINT_R0("check 3\n");
    return false;
  }

#pragma omp parallel for
  for(int i=0;i<based_lines;i++){
    if(order(nodes, edge[i][0], edge[i][1], added_centers) != MIDDLE)
      for(int j=1;j<groups;j++){
	int k = j * based_lines + i;
	int tmp0 = edge[k][0] - edge[k-based_lines][0];
	int tmp1 = edge[k][1] - edge[k-based_lines][1];
	if(added_centers){
	  tmp0 = (tmp0 < 0)? tmp0+nodes-added_centers : tmp0;
	  tmp1 = (tmp1 < 0)? tmp1+nodes-added_centers : tmp1;
	}
	else{
	  tmp0 = (tmp0 < 0)? tmp0 + nodes : tmp0;
	  tmp1 = (tmp1 < 0)? tmp1 + nodes : tmp1;
	}
	if(tmp0 != based_nodes || tmp1 != based_nodes){
	  PRINT_R0("check 4: %d\n", ii);
	  PRINT_R0("The different group relationship\n");
	  PRINT_R0("edge[%d][0]-edge[%d][0] = %d - %d = %d != %d\n", k, k-based_lines, edge[k][0], edge[k-based_lines][0], tmp0, based_nodes);
	  PRINT_R0("edge[%d][1]-edge[%d][1] = %d - %d = %d != %d\n", k, k-based_lines, edge[k][1], edge[k-based_lines][1], tmp1, based_nodes);
	  flag = false;
	}
      }
  }

  if(adjacency != NULL){
    int (*tmp_adjacency)[degree] = malloc(sizeof(int)*nodes*degree);
    create_adjacency(nodes, lines, degree, edge, tmp_adjacency);
  
    int sum[2] = {0,0};
    for(int i=0;i<nodes;i++){
      for(int j=0;j<degree;j++){
	sum[0] += adjacency[i][j];
	sum[1] += tmp_adjacency[i][j];
      }
      if(sum[0] != sum[1]){
	PRINT_R0("Eroor 5 %d %d\n", sum[0], sum[1]);
	flag = false;
	break;
      }
    }
    free(tmp_adjacency);
  }
  
  return flag;
}

bool has_duplicated_vertex(const int e00, const int e01, const int e10, const int e11)
{
  return (e00 == e10 || e01 == e11 || e00 == e11 || e01 == e10);
}

static void edge_exchange(const int nodes, const int lines, const int groups, const int degree, const int based_nodes,
			  int edge[lines][2], const int added_centers, const int ii)
{
  int line[groups*2], tmp_edge[groups*2][2];
  int based_lines = lines / groups;

  while(1){
    while(1){
      while(1){
	line[0] = getRandom(lines);
	line[1] = getRandom(lines);
	if(line[0] != line[1]) break;
      }
      if(has_duplicated_vertex(edge[line[0]][0], edge[line[0]][1], edge[line[1]][0], edge[line[1]][1])){
        continue;
      }
      else if((line[0] - line[1]) % based_lines == 0){
	if(edge_1g_opt(edge, nodes, lines, degree, based_nodes, based_lines, groups, line[0], added_centers))
	  return;
	else continue;
      }
      else break;
    }

    bool flag0 = (distance(nodes, edge[line[0]][0], edge[line[0]][1], added_centers) == (nodes-added_centers)/2);
    bool flag1 = (distance(nodes, edge[line[1]][0], edge[line[1]][1], added_centers) == (nodes-added_centers)/2);
    bool diameter_flag = ((flag0 || flag1) && groups%2 == 0);

    if(diameter_flag){
      if(edge_1g_opt(edge, nodes, lines, degree, based_nodes, based_lines, groups, line[0], added_centers))
	return;
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
	for(int i=0;i<groups;i++){
	  swap(&tmp_edge[i*2][1], &tmp_edge[i*2+1][1]);
	}
      }
      else{
	for(int i=0;i<groups;i++)
	  swap(&tmp_edge[i*2][1], &tmp_edge[i*2+1][0]);
      }
      
      assert(check_loop(groups*2, tmp_edge));
      if(!check_duplicate_edge(groups*2, tmp_edge)) continue;
      if(!check_duplicate_current_edge(lines, groups*2, line, edge, tmp_edge, groups, nodes, added_centers))
	continue;
      for(int i=0;i<groups*2;i++)
      	if(order(nodes, tmp_edge[i][0], tmp_edge[i][1], added_centers) == RIGHT)
      	  swap(&tmp_edge[i][0], &tmp_edge[i][1]); // RIGHT -> LEFT

#pragma omp parallel for
      for(int i=0;i<groups*2;i++){
	edge[line[i]][0] = tmp_edge[i][0];
        edge[line[i]][1] = tmp_edge[i][1];
      }

      break;
    }
  }
}

static bool accept(const double ASPL, const double current_ASPL, const double temp, const int nodes, const int groups,
		   const bool hill_climbing_flag, const bool detect_temp_flag, const long long i,
		   double *max_diff_energy, long long *total_accepts, long long *accepts, long long *rejects)
{
  if(ASPL <= current_ASPL){
    *accepts += 1;
    if(i > SKIP_ACCEPTS) *total_accepts +=1;
    return true;
  }
  if(hill_climbing_flag){ // Only accept when ASPL <= current_ASPL.
    *rejects += 1;
    return false;
  }

  double diff = ((current_ASPL-ASPL)*nodes*(nodes-1))/groups;

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
	     int edge[lines][2], int *diam, double *ASPL, const int cooling_cycle,
	     const int added_centers, const int based_nodes, long long *total_accepts)
{
  long long i, accepts = 0, rejects = 0;
  int best_edge[lines][2], tmp_edge[lines][2];
  edge_copy((int *)best_edge, (int *)edge, lines*2);

  // Create adjacency matrix
  int (*adjacency)[degree] = malloc(sizeof(int)*nodes*degree); // int adjacency[nodes][degree];
  create_adjacency(nodes, lines, degree, edge, adjacency);
  evaluation(nodes, based_nodes, groups, lines, degree, adjacency, diam, ASPL, added_centers);
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
      edge_copy((int *)tmp_edge, (int *)edge, lines*2);
      edge_exchange(nodes, lines, groups, degree, based_nodes, tmp_edge, added_centers, (int)i);
      create_adjacency(nodes, lines, degree, tmp_edge, adjacency);
      assert(check(nodes, based_nodes, lines, degree, groups, tmp_edge, added_centers, adjacency, (int)i));
      if(evaluation(nodes, based_nodes, groups, lines, degree, adjacency, diam, ASPL, added_centers)) break;
    }

    if(accept(*ASPL, current_ASPL, temp, nodes, groups, hill_climbing_flag,
	      detect_temp_flag, i, max_diff_energy, total_accepts, &accepts, &rejects)){
      current_ASPL  = *ASPL;
      current_diam  = *diam;
      edge_copy((int *)edge, (int *)tmp_edge, lines*2);
      if((best_diam > current_diam) || (best_diam == current_diam && best_ASPL > current_ASPL)){
	edge_copy((int *)best_edge, (int *)edge, lines*2);
	best_ASPL = current_ASPL;
	best_diam = current_diam;
      }

      if(best_diam == current_diam && best_ASPL == low_ASPL){
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
double estimated_elapse_time(const int nodes, const int based_nodes, const int lines, const int degree,
			     const int groups, int edge[lines][2], const int added_centers)
{
  int diam;    // Not use
  double ASPL; // Not use
  int (*adjacency)[degree] = malloc(sizeof(int)*nodes*degree); // int adjacency[nodes][degree];
  int (*tmp_edge)[2]       = malloc(sizeof(int)*lines*2);      // int tmp_edge[lines][2];

  timer_start(TIMER_ESTIMATED);
  for(int i=0;i<ESTIMATED_TIMES;i++){
    edge_copy((int *)tmp_edge, (int *)edge, lines*2);
    edge_exchange(nodes, lines, groups, degree, based_nodes, tmp_edge, added_centers, (int)i);
    create_adjacency(nodes, lines, degree, tmp_edge, adjacency);
    assert(check(nodes, based_nodes, lines, degree, groups, tmp_edge, added_centers, adjacency, (int)i));
    evaluation(nodes, based_nodes, groups, lines, degree, adjacency, &diam, &ASPL, added_centers);
  }  
  timer_stop(TIMER_ESTIMATED);
  
  free(tmp_edge);
  free(adjacency);

  return timer_read(TIMER_ESTIMATED)/ESTIMATED_TIMES;
}

// This function is mainly useful when groupe is 1.
void check_current_edge(const int nodes, const int degree, const int lines, const int groups,
			const int based_nodes, int edge[lines][2], const double low_ASPL, const int added_centers)
{
  int diam;    // Not use
  double ASPL;
  int (*adjacency)[degree] = malloc(sizeof(int)*nodes*degree); // int adjacency[nodes][degree];

  create_adjacency(nodes, lines, degree, edge, adjacency);
  if(! evaluation(nodes, based_nodes, groups, lines, degree, adjacency, &diam, &ASPL, added_centers))
    ERROR("The input file has a node which is never reached by another node.\n");

  if(ASPL == low_ASPL)
    END("The input file has already optimum solution.\n");

  free(adjacency);
}
