#include "common.h"

static double uniform_rand()
{
  return ((double)random()+1.0)/((double)RAND_MAX+2.0);
}

static void print_result_header()
{
  PRINT_R0("   Times\tTemp\t\tCur. ASPL GAP\t\tBest ASPL GAP              ");
  PRINT_R0("Cur. Dia. GAP    Best Dia. GAP    Cur. Len. GAP    Best Len. GAP    Accept Rate\n");
}

static void print_results(const long long num, const double temp, 
			  const double current_ASPL, const double best_ASPL, const double low_ASPL,
			  const int current_diam,    const int best_diam,    const int low_diam,
			  const int current_length,  const int best_length,  const int low_length,
			  const long long accepts, const long long rejects)
{
  PRINT_R0("%8lld\t%f\t", num, temp);
  PRINT_R0("%f ( %f )   %f ( %f )     %2d ( %2d )        %2d ( %2d )        %2d ( %2d )        %2d ( %2d )         ",
	   current_ASPL,   current_ASPL-low_ASPL,     best_ASPL,   best_ASPL-low_ASPL,
	   current_diam,   current_diam-low_diam,     best_diam,   best_diam-low_diam,
	   current_length, current_length-low_length, best_length, best_length-low_length);
  if(num != 0)
    PRINT_R0("%.4f ( %lld / %lld )\n", (double)accepts/(accepts+rejects), accepts, (accepts+rejects));
  else
    PRINT_R0("-\n");
}  

static void exchange_edge(const int nodes, const int lines, const int degree, int edge[lines][2],
			  const int height, const int width, const int groups, const int low_length, 
			  const bool enable_restriction, const double max_temp, const double min_temp,
			  const double temp, const long long ii)
{
  int tmp_line[groups*2], based_lines = lines/groups, based_nodes = nodes/groups;
  assert(lines%groups == 0);
  assert(nodes%groups == 0);

  while(1){
    while(1){
      tmp_line[0] = getRandom(lines);
      tmp_line[1] = getRandom(lines);
      if(tmp_line[0] == tmp_line[1]) continue;
      else if(has_duplicated_vertex(edge[tmp_line[0]][0], edge[tmp_line[0]][1], edge[tmp_line[1]][0], edge[tmp_line[1]][1])){
        continue;
      }
      else if((tmp_line[0] - tmp_line[1]) % based_lines == 0){
        if(edge_1g_opt(edge, nodes, lines, degree, based_nodes, based_lines, height, width, groups,
		       tmp_line[0], low_length, enable_restriction, max_temp, min_temp,
		       temp, ii))
	  return;
	else
	  continue;
      }
      else
	break;
    }

    bool flag0 = ( WIDTH (edge[tmp_line[0]][0], height) + WIDTH (edge[tmp_line[0]][1], height) == (width-1) &&
		   HEIGHT(edge[tmp_line[0]][0], height) + HEIGHT(edge[tmp_line[0]][1], height) == (height-1));
    bool flag1 = ( WIDTH (edge[tmp_line[1]][0], height)	+ WIDTH (edge[tmp_line[1]][1], height) == (width-1) &&
		   HEIGHT(edge[tmp_line[1]][0], height) + HEIGHT(edge[tmp_line[1]][1], height) == (height-1));
    bool diameter_flag = ((flag0 || flag1) && groups%2==0);
    if(diameter_flag){
      if(edge_1g_opt(edge, nodes, lines, degree, based_nodes, based_lines, height, width, groups,
		     tmp_line[0], low_length, enable_restriction, max_temp, min_temp,
		     temp, ii))
        return;
      else
        continue;
    }
    
    for(int i=1;i<groups;i++){
      tmp_line[i*2  ] = tmp_line[0] + based_lines * i;
      tmp_line[i*2+1] = tmp_line[1] + based_lines * i;
      if(tmp_line[i*2  ] >= lines) tmp_line[i*2  ] -= lines;
      if(tmp_line[i*2+1] >= lines) tmp_line[i*2+1] -= lines;
    }

#ifdef _DEBUG_MSG
    printf("edge_2g_opt: ii = %lld\n", ii);
    for(int i=0;i<groups*2;i++)
      printf("line[%d] = %d\n", i, tmp_line[i]);

    for(int i=0;i<groups;i++)
      printf("Before: %d,%d-%d,%d %d,%d-%d,%d\n",
    	     WIDTH(edge[tmp_line[i*2  ]][0], height), HEIGHT(edge[tmp_line[i*2  ]][0], height),
    	     WIDTH(edge[tmp_line[i*2  ]][1], height), HEIGHT(edge[tmp_line[i*2  ]][1], height),
    	     WIDTH(edge[tmp_line[i*2+1]][0], height), HEIGHT(edge[tmp_line[i*2+1]][0], height),
	     WIDTH(edge[tmp_line[i*2+1]][1], height), HEIGHT(edge[tmp_line[i*2+1]][1], height));
#endif
    
    int r = (getRandom(2) == 0)? 1 : 0;
    int tmp_edge[groups*2][2];
    for(int i=0;i<groups;i++){
      for(int j=0;j<2;j++){
	tmp_edge[i*2  ][j] = edge[tmp_line[i*2  ]][j];
	tmp_edge[i*2+1][j] = edge[tmp_line[i*2+1]][j];
      }
      swap(&tmp_edge[i*2][1], &tmp_edge[i*2+1][r]);
    }

    bool continue_flag = false;
    if(enable_restriction){
      for(int i=0;i<2;i++){
	int w0 = WIDTH (tmp_edge[i][0], height);
	int h0 = HEIGHT(tmp_edge[i][0], height);
	int w1 = WIDTH (tmp_edge[i][1], height);
	int h1 = HEIGHT(tmp_edge[i][1], height);
	int distance = abs(w0 - w1) + abs(h0 - h1);
	double alpha = 1.0 - (max_temp-temp)/(max_temp-min_temp);
	int threshold = (int)((height+width-low_length)*alpha) + low_length;
	if(distance > threshold){
	  continue_flag = true;
	  break;
	}
      }
    }
    if(continue_flag) continue;
    
    if(!check_duplicate_tmp_edge(D_2G_OPT, groups, tmp_edge))
      continue;
    else if(!check_duplicate_current_edge(lines, (const int (*)[2])edge, groups*2, 
					  (const int (*)[2])tmp_edge, tmp_line, groups, D_2G_OPT, false))
      continue;

    for(int i=0;i<groups;i++){
      for(int j=0;j<2;j++){
	edge[tmp_line[i*2  ]][j] = tmp_edge[i*2  ][j];
	edge[tmp_line[i*2+1]][j] = tmp_edge[i*2+1][j];
      }

#ifdef _DEBUG_MSG
       printf("After : %d,%d-%d,%d %d,%d-%d,%d\n",
	      WIDTH(edge[tmp_line[i*2  ]][0], height), HEIGHT(edge[tmp_line[i*2  ]][0], height),
	      WIDTH(edge[tmp_line[i*2  ]][1], height), HEIGHT(edge[tmp_line[i*2  ]][1], height),
	      WIDTH(edge[tmp_line[i*2+1]][0], height), HEIGHT(edge[tmp_line[i*2+1]][0], height),
	      WIDTH(edge[tmp_line[i*2+1]][1], height), HEIGHT(edge[tmp_line[i*2+1]][1], height));
#endif
     }
    break;
  }
}

// When the diameter is small, the length becomes long, so the diameter isn't adopted for evaluation.
static bool accept(const double new_ASPL, const double current_ASPL,
		   const int new_total_over_length, const int current_total_over_length, const double temp,
		   const int nodes, const int degree, const bool enable_hill_climbing, const bool enable_detect_temp, 
		   double *max_diff_energy, long long *total_accepts, long long *accepts, long long *rejects,
		   const double max_temp, const double min_temp, const double weight, const long long ii)
{
  //  if(new_length < current_length){
  //    *accepts += 1;
  //    if(ii > SKIP_ACCEPTS) *total_accepts +=1;
  //    return true;
  //  }
  //  else if(new_length > current_length){
  //    *rejects += 1;
  //    return false;
  //  }

  double f = (current_ASPL-new_ASPL)*nodes*(nodes-1);
  double p = (double)(current_total_over_length - new_total_over_length) / degree * nodes;
  //   p *= (max_temp - temp) / (max_temp - min_temp);
  double diff = f + p * weight;
  if(diff >= 0){
    *accepts += 1;
    if(ii > SKIP_ACCEPTS) *total_accepts +=1;
    return true;
  }
  if(enable_hill_climbing){ // Only accept when new_ASPL <= current_ASPL.
    *rejects += 1;
    return false;
  }

  if(enable_detect_temp)
    *max_diff_energy = MAX(*max_diff_energy, -1.0 * diff);

  if(exp(diff/temp) > uniform_rand()){
    *accepts += 1;
    if(ii > SKIP_ACCEPTS) *total_accepts +=1;
    return true;
  }
  else{
    *rejects += 1;
    return false;
  }
}

static void calc_length(const int lines, int edge[lines][2], const int height,
			const int low_length, int *length, int *total_over_length)
{
  *length = 0;
  *total_over_length = 0;
  
  for(int i=0;i<lines;i++){
    int w0 = WIDTH (edge[i][0], height);
    int h0 = HEIGHT(edge[i][0], height);
    int w1 = WIDTH (edge[i][1], height);
    int h1 = HEIGHT(edge[i][1], height);
    int distance = abs(w0 - w1) + abs(h0 - h1);
    *length = MAX((*length), distance);
    if(distance > low_length)
      *total_over_length += (distance-low_length);
  }
}

long long sa(const int nodes, const int lines, double temp, const long long ncalcs,
	     const double cooling_rate, const int low_diam,  const double low_ASPL, const bool enable_bfs, 
	     const bool enable_hill_climbing, const bool enable_detect_temp,
	     double *max_diff_energy, const double max_temp, const double min_temp, int edge[lines][2],
	     int *diam, double *ASPL, const int cooling_cycle, long long *total_accepts, const int width,
	     const int based_width, const int height, const int based_height, int *length,
	     const int low_length, const double weight, const int groups, const bool enable_restriction)
{
  int degree = 2*lines/nodes, best_edge[lines][2], tmp_edge[lines][2], based_nodes = nodes/groups;
  long long ii, accepts = 0, rejects = 0;

  // Create adjacency matrix
  int (*adjacency)[degree] = malloc(sizeof(int)*nodes*degree); // int adjacency[nodes][degree];
  create_adjacency(nodes, lines, degree, (const int (*)[2])edge, adjacency);
  evaluation(nodes, degree, groups, (const int* restrict)adjacency,
	     based_nodes, height, based_height, diam, ASPL, enable_bfs);

  int tmp_total_over_length = 0;
  calc_length(lines, edge, height, low_length, length, &tmp_total_over_length);

  double current_ASPL = *ASPL,   best_ASPL   = *ASPL,   tmp_ASPL;
  int current_diam    = *diam,   best_diam   = *diam,   tmp_diam;
  int current_length  = *length, best_length = *length, tmp_length;
  int current_total_over_length = tmp_total_over_length;
  int print_interval  = (ncalcs/NUM_OF_PROGRESS == 0)? 1 : ncalcs/NUM_OF_PROGRESS;
  copy_edge((int *)best_edge, (int *)edge, lines*2);
    
  if(rank == 0 && !enable_detect_temp)
    print_result_header();

  for(ii=0;ii<ncalcs;ii++){
    if(ii % print_interval == 0 && !enable_detect_temp){
      print_results(ii, temp,
		    current_ASPL,   best_ASPL,   low_ASPL,
		    current_diam,   best_diam,   low_diam,
		    current_length, best_length, low_length,
		    accepts, rejects);
      accepts = 0;
      rejects = 0;
    }

    while(1){
      copy_edge((int *)tmp_edge, (int *)edge, lines*2);
      exchange_edge(nodes, lines, degree, tmp_edge, height, width, groups, low_length, enable_restriction, max_temp, min_temp, temp, ii);
      assert(check_loop(lines, tmp_edge));
      assert(check_duplicate_all_edge(lines, tmp_edge));
      assert(check_degree(nodes, lines, tmp_edge));
      assert(check_symmetric_edge(lines, tmp_edge, height, width, based_height, groups));
      create_adjacency(nodes, lines, degree, (const int (*)[2])tmp_edge, adjacency);
      if(evaluation(nodes, degree, groups, (const int* restrict)adjacency,
		    based_nodes, height, based_height, &tmp_diam, &tmp_ASPL, enable_bfs)){
	calc_length(lines, tmp_edge, height, low_length, &tmp_length, &tmp_total_over_length);
	break;
      }
    }

    if(accept(tmp_ASPL, current_ASPL, tmp_total_over_length, current_total_over_length,
	      temp, nodes, degree, enable_hill_climbing, enable_detect_temp,
	      max_diff_energy, total_accepts, &accepts, &rejects,
	      max_temp, min_temp, weight, ii)){
      current_ASPL   = tmp_ASPL;
      current_diam   = tmp_diam;
      current_length = tmp_length;
      current_total_over_length = tmp_total_over_length;
      copy_edge((int *)edge, (int *)tmp_edge, lines*2);
      if((low_diam <= tmp_diam) &&
	 ((best_length > tmp_length) ||
	  (best_length == tmp_length && best_diam > tmp_diam) ||
	  (best_length == tmp_length && best_diam == tmp_diam && best_ASPL > tmp_ASPL))){
	copy_edge((int *)best_edge, (int *)edge, lines*2);
	best_length = tmp_length;
	best_diam   = tmp_diam;
	best_ASPL   = tmp_ASPL;
      }

      if(best_ASPL == low_ASPL && best_length == low_length){
	if(!enable_detect_temp){
	  print_results(ii, temp, current_ASPL, best_ASPL, low_ASPL,
			current_diam, best_diam, low_diam,
			current_length, best_length, low_length,
			accepts, rejects);
	  PRINT_R0("---\nFound optimum solution.\n");
	}
	break;
      }
    }
    
    if((ii+1)%cooling_cycle == 0)
      temp *= cooling_rate;
  }

  *ASPL   = best_ASPL;
  *diam   = best_diam;
  *length = best_length;
  copy_edge((int *)edge, (int *)best_edge, lines*2);
  free(adjacency);

  return ii;
}

double estimated_elapse_time(const int nodes, const int lines, const int edge[lines][2],
			     const int height, const int width, const int based_height, const int groups,
			     const int low_length, const bool enable_bfs)
{
  int based_nodes = nodes / groups;
  int degree = 2 * lines / nodes;
  int diam;    // Not use
  double ASPL; // Not use
  int (*adjacency)[degree] = malloc(sizeof(int)*nodes*degree); // int adjacency[nodes][degree];
  int (*tmp_edge)[2]       = malloc(sizeof(int)*lines*2);      // int tmp_edge[lines][2];
  
  timer_start(TIMER_ESTIMATED);
  for(int i=0;i<ESTIMATED_TIMES;i++){
    copy_edge((int *)tmp_edge, (int *)edge, lines*2);
    exchange_edge(nodes, lines, degree, tmp_edge, height, width, groups, low_length, false,
		  (double)NOT_USED, (double)NOT_USED, (double)NOT_USED, (int)i);
    create_adjacency(nodes, lines, degree, (const int (*)[2])tmp_edge, adjacency);
    evaluation(nodes, degree, groups, (const int* restrict)adjacency,
	       based_nodes, height, based_height, &diam, &ASPL, enable_bfs);
  }
  timer_stop(TIMER_ESTIMATED);
  
  free(tmp_edge);
  free(adjacency);

  return timer_read(TIMER_ESTIMATED)/ESTIMATED_TIMES;
}

// This function is mainly useful when groupe is 1.
void check_current_edge(const int nodes, const int lines, int edge[lines][2], const double low_ASPL,
			const int groups, const int height, const int based_height, const bool enable_bfs)
{
  int based_nodes = nodes / groups;
  int degree = 2 * lines / nodes;
  int diam;    // Not use
  double ASPL;
  int (*adjacency)[degree] = malloc(sizeof(int)*nodes*degree); // int adjacency[nodes][degree];

  create_adjacency(nodes, lines, degree, (const int (*)[2])edge, adjacency);
  if(! evaluation(nodes, degree, groups, (const int* restrict)adjacency,
		  based_nodes, height, based_height, &diam, &ASPL, enable_bfs))
    ERROR("The input file has a node which is never reached by another node.\n");

  if(ASPL == low_ASPL)
    END("The input file has already optimum solution.\n");

  free(adjacency);
}
