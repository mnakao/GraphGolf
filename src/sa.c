#include "common.h"

#ifndef NDEBUG
static void check_edge_restore(const int lines, const int edge[lines][2], const int prev_edge[lines][2])
{
  for(int i=0;i<lines;i++)
    if(edge[i][0] != prev_edge[i][0] || edge[i][1] != prev_edge[i][1])
      ERROR("Error in check_edge_restore");
}
#endif

static void edge_restore(const int lines, const int groups,
			 int edge[lines][2], const int restore_edge[groups*2][2],
			 const int restore_line[groups*2], const int restores)
{
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
  for(int i=0;i<restores;i++) // Value of "restores" is groups or groups * 2.
    for(int j=0;j<2;j++)
      edge[restore_line[i]][j] = restore_edge[i][j];
}

static void adjacency_restore(const int nodes, const int degree, const int groups, int adjacency[nodes][degree],
			      const int restore_adjacency[groups*2][2][3], const int restores)
{
  if(restores == groups){
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
    for(int i=0;i<groups;i++){
      for(int j=0;j<2;j++){
	int t0 = restore_adjacency[i][j][0];
	int t1 = restore_adjacency[i][j][1];
	adjacency[t0][t1] = restore_adjacency[i][j][2];
      }
    }
  }
  else
    for(int i=restores-2;i>=0;i-=2){
      int t0 = restore_adjacency[i  ][1][0];
      int t1 = restore_adjacency[i  ][1][1];
      int t2 = restore_adjacency[i+1][1][0];
      int t3 = restore_adjacency[i+1][1][1];
      swap(&adjacency[t0][t1], &adjacency[t2][t3]);
      
      t0 = restore_adjacency[i  ][0][0];
      t1 = restore_adjacency[i  ][0][1];
      t2 = restore_adjacency[i+1][0][0];
      t3 = restore_adjacency[i+1][0][1];
      swap(&adjacency[t0][t1], &adjacency[t2][t3]);
  }
}

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

static void change_adjacency_2g_opt(const int nodes, const int degree, int adjacency[nodes][degree], const int groups,
				    int tmp_edge[groups*2][2], int restore_adjacency[groups*2][2][3], const int r)
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int i=0;i<groups*2;i+=2){
    int t0, t1, t2, t3;
    if(r==0){
      t0 = tmp_edge[i][0];
      for(t1=0;t1<degree;t1++)
        if(adjacency[t0][t1] == tmp_edge[i][1])
          break;
      t2 = tmp_edge[i+1][0];
      for(t3=0;t3<degree;t3++)
        if(adjacency[t2][t3] == tmp_edge[i+1][1])
          break;
      restore_adjacency[i  ][0][0] = t0;
      restore_adjacency[i  ][0][1] = t1;
      restore_adjacency[i  ][0][2] = adjacency[t0][t1];
      restore_adjacency[i+1][0][0] = t2;
      restore_adjacency[i+1][0][1] = t3;
      restore_adjacency[i+1][0][2] = adjacency[t2][t3];
      
      t0 = tmp_edge[i+1][1];
      for(t1=0;t1<degree;t1++)
        if(adjacency[t0][t1] == tmp_edge[i+1][0])
          break;
      t2 = tmp_edge[i][1];
      for(t3=0;t3<degree;t3++)
        if(adjacency[t2][t3] == tmp_edge[i][0])
          break;
      restore_adjacency[i  ][1][0] = t0;
      restore_adjacency[i  ][1][1] = t1;
      restore_adjacency[i  ][1][2] = adjacency[t0][t1];
      restore_adjacency[i+1][1][0] = t2;
      restore_adjacency[i+1][1][1] = t3;
      restore_adjacency[i+1][1][2] = adjacency[t2][t3];
    }
    else{
      t0 = tmp_edge[i][0];
      for(t1=0;t1<degree;t1++)
        if(adjacency[t0][t1] == tmp_edge[i][1])
          break;
      t2 = tmp_edge[i+1][1];
      for(t3=0;t3<degree;t3++)
        if(adjacency[t2][t3] == tmp_edge[i+1][0])
          break;
      restore_adjacency[i  ][0][0] = t0;
      restore_adjacency[i  ][0][1] = t1;
      restore_adjacency[i  ][0][2] = adjacency[t0][t1];
      restore_adjacency[i+1][0][0] = t2;
      restore_adjacency[i+1][0][1] = t3;
      restore_adjacency[i+1][0][2] = adjacency[t2][t3];
      
      t0 = tmp_edge[i][1];
      for(t1=0;t1<degree;t1++)
        if(adjacency[t0][t1] == tmp_edge[i][0])
          break;
      t2 = tmp_edge[i+1][0];
      for(t3=0;t3<degree;t3++)
        if(adjacency[t2][t3] == tmp_edge[i+1][1])
          break;
      restore_adjacency[i  ][1][0] = t0;
      restore_adjacency[i  ][1][1] = t1;
      restore_adjacency[i  ][1][2] = adjacency[t0][t1];
      restore_adjacency[i+1][1][0] = t2;
      restore_adjacency[i+1][1][1] = t3;
      restore_adjacency[i+1][1][2] = adjacency[t2][t3];
    }
  }

  for(int i=0;i<groups*2;i+=2){
    for(int j=0;j<2;j++){
      int t0 = restore_adjacency[i  ][j][0];
      int t1 = restore_adjacency[i  ][j][1];
      int t2 = restore_adjacency[i+1][j][0];
      int t3 = restore_adjacency[i+1][j][1];
      swap(&adjacency[t0][t1], &adjacency[t2][t3]);
    }
  }
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

#ifdef _OPENMP
#pragma omp parallel for
#endif
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
  
#ifdef _OPENMP
#pragma omp parallel for
#endif
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

#ifdef _OPENMP
#pragma omp parallel for
#endif
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
			  int edge[lines][2], const int added_centers, int adjacency[nodes][degree],
			  int restore_edge[groups*2][2], int restore_adjacency[groups*2][2][3], int restore_line[groups*2],
			  int *restores, const double *node_value, const int ii)
{
  int line[groups*2], tmp_edge[groups*2][2], edge_for_adj[groups*2][2], r;
  int based_lines = lines / groups;
  int fst_edges[based_lines], snd_edges[based_lines];
  double edge_value[based_lines];
  int num_of_fsts = 0, num_of_snds = 0;

  if(node_value != NULL){
    double max = 0;
    for(int i=0;i<based_lines;i++){
      edge_value[i] = node_value[edge[i][0]%based_nodes] + node_value[edge[i][1]%based_nodes];
      //      printf("%d: %d : %d %d\n", i, (int)edge_value[i], edge[i][0], edge[i][1]);
      if(edge_value[i] > max){
	max = edge_value[i];
	fst_edges[0] = i;
	num_of_fsts  = 1;
      }
      else if(edge_value[i] == max)
	fst_edges[num_of_fsts++] = i;
    }
    //    printf("1 max = %d, num = %d\n", (int)max, num_of_fsts);

    if(num_of_fsts == 1){
      max = 0;
      for(int i=0;i<based_lines;i++){
	if(i == fst_edges[0]) continue;
	edge_value[i] = node_value[edge[i][0]%based_nodes] + node_value[edge[i][1]%based_nodes];
	if(edge_value[i] > max){
	  max = edge_value[i];
	  snd_edges[0] = i;
	  num_of_snds  = 1;
	}
	else if(edge_value[i] == max)
	  snd_edges[num_of_snds++] = i;
      }
      //      printf("2 max = %d, num = %d\n", (int)max, num_of_snds);
    }
    //    for(int i=0;i<num_of_fsts;i++)
    //      printf("%3d", fst_edges[i]);
    //    printf("\n");
    //    if(num_of_fsts == 1){
    //      for(int i=0;i<num_of_snds;i++)
    //	printf("%3d", snd_edges[i]);
    //      printf("\n--\n");
    //    }
  }

  while(1){
    while(1){
      while(1){
	if(node_value != NULL && getRandom(2) == 0){
	  if(num_of_fsts >= 2){
	    line[0] = fst_edges[getRandom(num_of_fsts)] + based_lines * getRandom(groups);
	    line[1] = fst_edges[getRandom(num_of_fsts)] + based_lines * getRandom(groups);
	  }
	  else{ // num_of_bests must be 1
	    line[0] = fst_edges[0] + based_lines * getRandom(groups);
	    line[1] = snd_edges[getRandom(num_of_snds)] + based_lines * getRandom(groups);
	  }
	}
	else{
	  line[0] = getRandom(lines);
	  line[1] = getRandom(lines);
	}

	if(line[0] != line[1]) break;
      }
      if(has_duplicated_vertex(edge[line[0]][0], edge[line[0]][1], edge[line[1]][0], edge[line[1]][1])){
        continue;
      }
      else if((line[0] - line[1]) % based_lines == 0){
	if(edge_1g_opt(edge, nodes, lines, degree, based_nodes, based_lines, groups, line[0], added_centers, adjacency,
		       restore_edge, restore_adjacency, restore_line, restores)){
	  return;
	}
	else continue;
      }
      else break;
    }

    bool flag0 = (distance(nodes, edge[line[0]][0], edge[line[0]][1], added_centers) == (nodes-added_centers)/2);
    bool flag1 = (distance(nodes, edge[line[1]][0], edge[line[1]][1], added_centers) == (nodes-added_centers)/2);
    bool diameter_flag = ((flag0 || flag1) && groups%2 == 0);

    if(diameter_flag){
      if(edge_1g_opt(edge, nodes, lines, degree, based_nodes, based_lines, groups, line[0], added_centers, adjacency,
		     restore_edge, restore_adjacency, restore_line, restores))
	
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
	  edge_for_adj[i][j] = tmp_edge[i][j] = edge[line[i]][j];

      r = getRandom(2);
      if(r == 0){
	for(int i=0;i<groups;i++)
	  swap(&tmp_edge[i*2][1], &tmp_edge[i*2+1][1]);
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

#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(int i=0;i<groups*2;i++){
	restore_line[i]    = line[i];
	restore_edge[i][0] = edge[line[i]][0];
	restore_edge[i][1] = edge[line[i]][1];
	edge[line[i]][0]   = tmp_edge[i][0];
        edge[line[i]][1]   = tmp_edge[i][1];
      }

      break;
    }
  }

  *restores = groups * 2;
  change_adjacency_2g_opt(nodes, degree, adjacency, groups, edge_for_adj, restore_adjacency, r);
}

static bool accept(const double ASPL, const double current_ASPL, const double temp, const int nodes, const int groups,
		   const bool hill_climbing_flag, const bool detect_temp_flag, const int i, double *max_diff_energy,
		   long long *total_accepts, int *accepts, int *rejects)
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
	     int edge[lines][2], int *diam, double *ASPL, const int cooling_cycle, const int added_centers,
	     const int based_nodes, long long *total_accepts, double *node_value)
{
  int restore_edge[groups*2][2], restore_adjacency[groups*2][2][3], restore_line[groups*2], restores = 0;
  int best_edge[lines][2], accepts = 0, rejects = 0;
  double tmp_value[based_nodes];
#ifndef NDEBUG
  int prev_edge[lines][2];
#endif
  long long i;
  edge_copy((int *)best_edge, (int *)edge, lines*2);

  // Create adjacency matrix
  int (*adjacency)[degree] = malloc(sizeof(int)*nodes*degree); // int adjacency[nodes][degree];
  create_adjacency(nodes, lines, degree, edge, adjacency);
  evaluation(nodes, based_nodes, groups, lines, degree, adjacency, diam, ASPL, added_centers, node_value);
  memcpy(tmp_value, node_value, sizeof(double)*based_nodes);
  
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
      edge_restore(lines, groups, edge, restore_edge, restore_line, restores);
      adjacency_restore(nodes, degree, groups, adjacency, restore_adjacency, restores);
#ifndef NDEBUG
      if(restores != 0)
	check_edge_restore(lines, edge, prev_edge);
      memcpy(prev_edge, edge, sizeof(int)*lines*2);
#endif
      edge_exchange(nodes, lines, groups, degree, based_nodes, edge, added_centers, adjacency,
		    restore_edge, restore_adjacency, restore_line, &restores, tmp_value, (int)i);
      assert(check(nodes, based_nodes, lines, degree, groups, edge, added_centers, adjacency, (int)i));
      if(evaluation(nodes, based_nodes, groups, lines, degree, adjacency, diam, ASPL, added_centers, node_value)) break;
    }

    if(accept(*ASPL, current_ASPL, temp, nodes, groups, hill_climbing_flag,
	      detect_temp_flag, i, max_diff_energy, total_accepts, &accepts, &rejects)){
      memcpy(tmp_value, node_value, sizeof(double)*based_nodes);
      current_ASPL  = *ASPL;
      current_diam  = *diam;
      restores = 0;
      if(best_ASPL > current_ASPL && best_diam >= current_diam){
	edge_copy((int *)best_edge, (int *)edge, lines*2);
	best_ASPL = current_ASPL;
	best_diam = current_diam;
      }

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
double estimated_elapse_time(const long long ncals, const int nodes, const int based_nodes, const int lines,
			     const int degree, const int groups, int edge[lines][2], const int added_centers)
{
  int diam;    // Not use
  double ASPL; // Not use
  int restore_edge[groups*2][2], restore_adjacency[groups*2][2][3], restore_line[groups*2], restores = 0;
  int (*adjacency)[degree] = malloc(sizeof(int)*nodes*degree); // int adjacency[nodes][degree];
  int (*tmp_edge)[2]       = malloc(sizeof(int)*lines*2);      // int tmp_edge[lines][2];
  
  create_adjacency(nodes, lines, degree, edge, adjacency);
  edge_copy((int *)tmp_edge, (int *)edge, lines*2);

  timer_start(TIMER_ESTIMATED);
  for(int i=0;i<ESTIMATED_TIMES;i++){
    edge_restore(lines, groups, tmp_edge, restore_edge, restore_line, restores);
    adjacency_restore(nodes, degree, groups, adjacency, restore_adjacency, restores);
    edge_exchange(nodes, lines, groups, degree, based_nodes, tmp_edge, added_centers,
		  adjacency, restore_edge, restore_adjacency, restore_line, &restores, NULL, (int)i);
    assert(check(nodes, based_nodes, lines, degree, groups, tmp_edge, added_centers, adjacency, (int)i));
    evaluation(nodes, based_nodes, groups, lines, degree, adjacency, &diam, &ASPL, added_centers, NULL);
  }
  timer_stop(TIMER_ESTIMATED);
  
  free(tmp_edge);
  free(adjacency);

  return timer_read(TIMER_ESTIMATED)/ESTIMATED_TIMES;
}

// This function is mainly useful when groupe is 1.
void check_current_edge(const int nodes, const int degree, const int lines, const int groups, const int based_nodes,
			int edge[lines][2], const double low_ASPL, const int added_centers)
{
  int diam;    // Not use
  double ASPL;
  int (*adjacency)[degree] = malloc(sizeof(int)*nodes*degree); // int adjacency[nodes][degree];

  create_adjacency(nodes, lines, degree, edge, adjacency);
  if(! evaluation(nodes, based_nodes, groups, lines, degree, adjacency, &diam, &ASPL, added_centers, NULL))
    ERROR("The input file has a node which is never reached by another node.\n");

  if(ASPL == low_ASPL)
    END("The input file has already optimum solution.\n");

  free(adjacency);
}
