#include "common.h"

static double uniform_rand()
{
  return ((double)random()+1.0)/((double)RAND_MAX+2.0);
}

static void print_result_header()
{
  PRINT_R0("   Times\tTemp\t\tCur. ASPL GAP\t\tBest ASPL GAP              ");
  PRINT_R0("Cur. Dia. GAP    Best Dia. GAP    Accept Rate\n");
}

static void print_results(const long long num, const double temp, const bool enable_fixed_temp,
			  const double current_ASPL, const double best_ASPL, const double low_ASPL,
			  const int current_diam,    const int best_diam,    const int low_diam,
			  const double fixed_temp,   const long long accepts, const long long rejects)
{
  if(enable_fixed_temp)
    PRINT_R0("%8lld\t%f\t", num, fixed_temp);
  else
    PRINT_R0("%8lld\t%f\t", num, temp);
  PRINT_R0("%f ( %f )   %f ( %f )     %2d ( %2d )        %2d ( %2d )         ",
	   current_ASPL,   current_ASPL-low_ASPL,     best_ASPL,   best_ASPL-low_ASPL,
	   current_diam,   current_diam-low_diam,     best_diam,   best_diam-low_diam);
  if(num != 0)
    PRINT_R0("%.4f ( %lld / %lld )\n", (double)accepts/(accepts+rejects), accepts, (accepts+rejects));
  else
    PRINT_R0("-\n");
}  

static bool edge_1g_opt(int (*edge)[2], const int nodes, const int lines, const int degree, const int based_nodes,
			const int based_lines, const int height, const int width, const int groups, const int start_line,
			const int low_length, const long long ii)
{
  assert(groups != 1);

  int tmp_line[groups], tmp_edge[groups][2], pattern;
  for(int i=0;i<groups;i++){
    int t = start_line + based_lines * i;
    tmp_line[i] = (t < lines)? t : t - lines;
  }

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

  for(int i=0;i<groups;i+=2)
    if(DISTANCE(tmp_edge[i][0], tmp_edge[i][1], height) > low_length)
      return false;

  if(!check_duplicate_current_edge(lines, (const int (*)[2])edge, groups,
                                   (const int (*)[2])tmp_edge, tmp_line,
                                   groups, D_1G_OPT, (pattern==groups)))
     return false;

  for(int i=0;i<groups;i++){
    edge[tmp_line[i]][0] = tmp_edge[i][0];
    edge[tmp_line[i]][1] = tmp_edge[i][1];
  }

  return true;
}

void exchange_edge(const int nodes, const int lines, const int degree, int edge[lines][2],
		   const int height, const int width, const int groups, const int low_length, const long long ii)
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
		       tmp_line[0], low_length, ii))
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
		     tmp_line[0], low_length, ii))
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

    int r = getRandom(2), tmp_edge[groups*2][2];
    for(int i=0;i<groups;i++){
      for(int j=0;j<2;j++){
	tmp_edge[i*2  ][j] = edge[tmp_line[i*2  ]][j];
	tmp_edge[i*2+1][j] = edge[tmp_line[i*2+1]][j];
      }
      swap(&tmp_edge[i*2][1], &tmp_edge[i*2+1][r]);
    }

    bool flag = false;
    for(int i=0;i<2;i++)
      if(DISTANCE(tmp_edge[i][0], tmp_edge[i][1], height) > low_length)
	flag = true;

    if(flag) continue;
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

static bool accept(const double new_ASPL, const double current_ASPL, const int new_diam, const int current_diam,
		   const double fixed_temp, const double temp, const int nodes, const int degree,
		   const bool enable_hill_climbing, const bool enable_detect_temp, double *max_diff_energy, long long *total_accepts,
		   long long *accepts, long long *rejects, const double max_temp, const double min_temp,
		   const int groups, const bool enable_fixed_temp, const long long ii)
{
  double diff = ((current_ASPL-new_ASPL)*nodes*(nodes-1)) / groups;

  if(enable_detect_temp){
    *max_diff_energy = MAX(*max_diff_energy, -1.0 * diff);
  }
  else{
    if(new_diam < current_diam){
      *accepts += 1;
      if(ii > SKIP_ACCEPTS) *total_accepts +=1;
      return true;
    }
    else if(new_diam > current_diam){
      *rejects += 1;
      return false;
    }
    else{ //  new_diam == current_diam
      if(diff >= 0){
	*accepts += 1;
	if(ii > SKIP_ACCEPTS) *total_accepts +=1;
	return true;
      }
      if(enable_hill_climbing){ // Only accept when new_ASPL <= current_ASPL.
	*rejects += 1;
	return false;
      }
    }
  }

  double v = (enable_fixed_temp)? exp(diff/fixed_temp) : exp(diff/temp);

  if(v > uniform_rand()){
    *accepts += 1;
    if(ii > SKIP_ACCEPTS) *total_accepts +=1;
    return true;
  }
  else{
    *rejects += 1;
    return false;
  }
}

long long sa(const int nodes, const int lines, const int degree, const int based_nodes, const long long ncalcs, const double cooling_rate,
	     const int low_diam,  const double low_ASPL, const bool enable_bfs, const bool enable_hill_climbing,
	     const bool enable_detect_temp, double *max_diff_energy, const double max_temp, const double min_temp,
	     const double fixed_temp, int edge[lines*2], int *diam, double *ASPL, const int cooling_cycle,
	     long long *total_accepts, const int width, const int based_width, const int height, const int based_height,
	     const int low_length, const int groups, const int *rotate_hash, const bool enable_fixed_temp)
{
  long long ii, accepts = 0, rejects = 0;
  double temp = max_temp;
  int (*best_edge)[2] = malloc(sizeof(int)*lines*2);
  int (*tmp_edge)[2]  = malloc(sizeof(int)*lines*2);
  int (*adjacency)[degree] = malloc(sizeof(int)*nodes*degree); // int adjacency[nodes][degree];
  
  create_adjacency(nodes, lines, degree, (const int (*)[2])edge, adjacency);
  evaluation(nodes, degree, groups, (const int* restrict)adjacency,
	     based_nodes, height, based_height, diam, ASPL, enable_bfs, rotate_hash);

  double current_ASPL = *ASPL, best_ASPL = *ASPL, tmp_ASPL;
  int current_diam    = *diam, best_diam = *diam, tmp_diam;
  int print_interval  = (ncalcs/NUM_OF_PROGRESS == 0)? 1 : ncalcs/NUM_OF_PROGRESS;
  copy_edge((int *)best_edge, (int *)edge, lines*2);
    
  if(rank == 0 && !enable_detect_temp)
    print_result_header();

  for(ii=0;ii<ncalcs;ii++){
    if(ii % print_interval == 0 && !enable_detect_temp){
      print_results(ii, temp, enable_fixed_temp, current_ASPL, best_ASPL,low_ASPL,
		    current_diam, best_diam, low_diam, fixed_temp, accepts, rejects);
      accepts = 0;
      rejects = 0;
    }

    while(1){
      copy_edge((int *)tmp_edge, (int *)edge, lines*2);
      exchange_edge(nodes, lines, degree, tmp_edge, height, width, groups, low_length, ii);
      create_adjacency(nodes, lines, degree, (const int (*)[2])tmp_edge, adjacency);
      if(evaluation(nodes, degree, groups, (const int* restrict)adjacency,
		    based_nodes, height, based_height, &tmp_diam, &tmp_ASPL, enable_bfs, rotate_hash))
	break;
    }

    if(accept(tmp_ASPL, current_ASPL, tmp_diam, current_diam, fixed_temp, temp, nodes, degree,
	      enable_hill_climbing, enable_detect_temp, max_diff_energy, total_accepts,
	      &accepts, &rejects, max_temp, min_temp, groups, enable_fixed_temp, ii)){
      current_ASPL = tmp_ASPL;
      current_diam = tmp_diam;
      copy_edge((int *)edge, (int *)tmp_edge, lines*2);
      if((best_diam > tmp_diam) ||
	 (best_diam == tmp_diam && best_ASPL > tmp_ASPL)){
	copy_edge((int *)best_edge, (int *)edge, lines*2);
	best_diam = tmp_diam;
	best_ASPL = tmp_ASPL;
      }

      if(best_diam == low_diam && best_ASPL == low_ASPL){
	if(!enable_detect_temp){
	  print_results(ii, temp, enable_fixed_temp, current_ASPL, best_ASPL, low_ASPL,
			current_diam, best_diam, low_diam, fixed_temp, accepts, rejects);
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
  copy_edge((int *)edge, (int *)best_edge, lines*2);
  
  free(best_edge);
  free(tmp_edge);
  free(adjacency);

  return ii;
}

double estimated_elapse_time(const int nodes, const int lines, const int edge[lines*2],
			     const int height, const int width, const int based_height, const int groups,
			     const int low_length, const bool enable_bfs, const int *rotate_hash)
{
  int based_nodes = nodes / groups;
  int degree = 2 * lines / nodes;
  int diam;    // Not use
  double ASPL; // Not use
  int (*adjacency)[degree] = malloc(sizeof(int)*nodes*degree); // int adjacency[nodes][degree];
  int (*tmp_edge)[2]       = malloc(sizeof(int)*lines*2);      // int tmp_edge[lines][2];
  
  timer_start(TIMER_ESTIMATED);
  for(int i=0;i<ESTIMATED_TIMES;i++){
    copy_edge((int *)tmp_edge, edge, lines*2);
    exchange_edge(nodes, lines, degree, tmp_edge, height, width, groups, low_length, i);
    create_adjacency(nodes, lines, degree, (const int (*)[2])tmp_edge, adjacency);
    evaluation(nodes, degree, groups, (const int* restrict)adjacency,
	       based_nodes, height, based_height, &diam, &ASPL, enable_bfs, rotate_hash);
  }
  timer_stop(TIMER_ESTIMATED);
  
  free(tmp_edge);
  free(adjacency);

  return timer_read(TIMER_ESTIMATED)/ESTIMATED_TIMES;
}

// This function is mainly useful when groupe is 1.
void check_current_edge(const int nodes, const int lines, const int edge[lines*2], const double low_ASPL,
			const int low_diam, const int groups, const int height, const int based_height,
			const bool enable_bfs, const int *rotate_hash)
{
  int based_nodes = nodes / groups;
  int degree = 2 * lines / nodes;
  int diam;
  double ASPL;
  int (*adjacency)[degree] = malloc(sizeof(int)*nodes*degree); // int adjacency[nodes][degree];

  create_adjacency(nodes, lines, degree, (const int (*)[2])edge, adjacency);
  if(! evaluation(nodes, degree, groups, (const int* restrict)adjacency,
		  based_nodes, height, based_height, &diam, &ASPL, enable_bfs, rotate_hash))
    ERROR("The input file has a node which is never reached by another node.\n");

  if(diam == low_diam && ASPL == low_ASPL)
    END("The input file has already optimum solution.\n");

  free(adjacency);
}
