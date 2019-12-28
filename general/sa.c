#include "common.h"

static void restore_edge(const int groups, const int kind_opt, int* restrict edge, int* restrict restored_line, const int* restrict restored_edge)
{
  if(kind_opt != D_1G_OPT && kind_opt != D_2G_OPT)
    ERROR("Wrong kind_opt: %d\n", kind_opt);

#pragma omp parallel for
  for(int i=0;i<groups*kind_opt;i++){
    edge[restored_line[i]*2  ] = restored_edge[i*2  ];
    edge[restored_line[i]*2+1] = restored_edge[i*2+1];
  }
}

static void restore_adj(const int degree, const int groups, int* restrict adj, const int kind_opt, int* restrict restored_adj_value,
			int* restrict restored_adj_idx_y, int* restrict restored_adj_idx_x)
{
  if(kind_opt != D_1G_OPT && kind_opt != D_2G_OPT)
    ERROR("Wrong kind_opt: %d\n", kind_opt);

#pragma omp parallel for
  for(int i=0;i<kind_opt*groups*2;i++){
    int y = restored_adj_idx_y[i];
    int x = restored_adj_idx_x[i];
    adj[y*degree+x] = restored_adj_value[i];
  }
}

static void copy_edge(int *restrict dst, const int *restrict src, const int n)
{
#pragma omp parallel for
  for(int i=0;i<n;i++)
    dst[i] = src[i];
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

void create_adj(const int nodes, const int lines, const int degree, 
		const int edge[lines][2], int adj[nodes][degree])
{
  int count[nodes];
  for(int i=0;i<nodes;i++)
    count[i] = 0;

  for(int i=0;i<lines;i++){
    int n1 = edge[i][0];
    int n2 = edge[i][1];
    adj[n1][count[n1]++] = n2;
    adj[n2][count[n2]++] = n1;
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
	   int edge[lines][2], const int added_centers, int* adj, const int ii)
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

  if(adj != NULL){
    int *tmp_adj = malloc(sizeof(int)*nodes*degree);
    create_adj(nodes, lines, degree, (const int (*)[2])edge, (int (*)[degree])tmp_adj);

    for(int i=0;i<nodes;i++){
      int sum[2] = {0,0};
      for(int j=0;j<degree;j++){
	sum[0] += *(adj     + i * degree + j);
	sum[1] += *(tmp_adj + i * degree + j);
      }
      if(sum[0] != sum[1]){
	PRINT_R0("[ii=%d] Error 5 %d %d\n", ii, sum[0], sum[1]);
	for(int j=0;j<degree;j++)
	  PRINT_R0("%d ", *(adj     + i * degree + j));
	PRINT_R0("\n");
	for(int j=0;j<degree;j++)
	  PRINT_R0("%d ", *(tmp_adj     + i * degree + j));
	PRINT_R0("\n");
	flag = false;
	break;
      }
    }

    for(int i=0;i<nodes;i++){
      for(int j=0;j<degree;j++){
	int tmp = *(adj + i * degree + j);
	int k;
	for(k=0;k<degree;k++)
	  if(tmp == *(tmp_adj + i * degree + k))
	    break;
	if(k == degree){
	  PRINT_R0("[ii=%d] Error 6\n", ii);
	  flag = false;
	  break;
	}
      }
    }

    for(int i=0;i<nodes;i++){
      for(int j=0;j<degree;j++){
        int tmp = *(tmp_adj + i * degree + j);
        int k;
        for(k=0;k<degree;k++)
          if(tmp == *(adj + i * degree + k))
            break;
        if(k == degree){
          PRINT_R0("[ii=%d] Error 7\n", ii);
          flag = false;
          break;
        }
      }
    }

    for(int i=0;i<nodes;i++){
      for(int j=0;j<degree;j++){
	int tmp = *(adj + i * degree + j);
	for(int k=j+1;k<degree;k++)
	  if(tmp == *(adj + i * degree + k)){
	    flag = false;
	    break;
	  }
      }
    }
    free(tmp_adj);
  }
  
  return flag;
}

bool has_duplicated_vertex(const int e00, const int e01, const int e10, const int e11)
{
  return (e00 == e10 || e01 == e11 || e00 == e11 || e01 == e10);
}

void exchange_edge_2opt(const int nodes, const int lines, const int groups, const int degree,
			const int based_nodes, int edge[lines][2], const int added_centers,
			int* restrict adj, int *kind_opt, int* restrict restored_edge, int* restrict restored_line,
			int* restrict restored_adj_value, int* restrict restored_adj_idx_y,
			int* restrict restored_adj_idx_x, const bool is_simple_graph, const int ii)
{
  int tmp_line[groups*2], tmp_edge[groups*2][2], r;
  int based_lines = lines / groups;

  while(1){
    while(1){
      while(1){
	tmp_line[0] = getRandom(lines);
	tmp_line[1] = getRandom(lines);
	if(tmp_line[0] != tmp_line[1]) break;
      }
      if(has_duplicated_vertex(edge[tmp_line[0]][0], edge[tmp_line[0]][1], edge[tmp_line[1]][0], edge[tmp_line[1]][1])){
	continue;
      }
      else if((tmp_line[0] - tmp_line[1]) % based_lines == 0){
	if(edge_1g_opt(edge, nodes, lines, degree, based_nodes, based_lines, groups, tmp_line[0], added_centers,
		       adj, kind_opt, restored_edge, restored_line, restored_adj_value, restored_adj_idx_y,
		       restored_adj_idx_x, is_simple_graph, ii))
	  return;
	else
	  continue;
      }
      else break;
    }

    bool flag0 = (distance(nodes, edge[tmp_line[0]][0], edge[tmp_line[0]][1], added_centers) == (nodes-added_centers)/2);
    bool flag1 = (distance(nodes, edge[tmp_line[1]][0], edge[tmp_line[1]][1], added_centers) == (nodes-added_centers)/2);
    bool diameter_flag = ((flag0 || flag1) && groups%2 == 0);

    if(diameter_flag){
      if(edge_1g_opt(edge, nodes, lines, degree, based_nodes, based_lines, groups, tmp_line[0], added_centers,
		     adj, kind_opt, restored_edge, restored_line, restored_adj_value, restored_adj_idx_y,
		     restored_adj_idx_x, is_simple_graph, ii))
	return;
      else
	continue;
    }

    // 2g-opt
    for(int i=1;i<groups;i++){
      int tmp0 = tmp_line[0] + based_lines * i;
      int tmp1 = tmp_line[1] + based_lines * i;
      tmp_line[0+2*i] = (tmp0 >= lines)? tmp0 - lines : tmp0;
      tmp_line[1+2*i] = (tmp1 >= lines)? tmp1 - lines : tmp1;
    }
    
    for(int i=0;i<groups*2;i++)
      for(int j=0;j<2;j++)
	tmp_edge[i][j] = edge[tmp_line[i]][j];
    
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
    if(!check_duplicate_tmp_edge(2, groups, tmp_edge))
      continue;
    else if(!check_duplicate_current_edge(lines, groups*2, tmp_line, edge, tmp_edge, groups, 2, false))
      continue;
    else
      break;
  } // end while
  
  for(int i=0;i<groups*2;i++)
    if(order(nodes, tmp_edge[i][0], tmp_edge[i][1], added_centers) == RIGHT)
      swap(&tmp_edge[i][0], &tmp_edge[i][1]); // RIGHT -> LEFT
  
  if(is_simple_graph){
    // Change a part of adj.
    int y0[groups], y1[groups], y2[groups], y3[groups];
    int x0[groups], x1[groups], x2[groups], x3[groups];
    
#pragma omp parallel for
    for(int i=0;i<groups;i++){
      y0[i] = edge[tmp_line[i*2  ]][0];
      y1[i] = edge[tmp_line[i*2  ]][1];
      y2[i] = edge[tmp_line[i*2+1]][0];
      y3[i] = edge[tmp_line[i*2+1]][1];

      for(x0[i]=0;x0[i]<degree;x0[i]++)
	if(adj[y0[i]*degree+x0[i]] == y1[i])
	  break;
      
      for(x1[i]=0;x1[i]<degree;x1[i]++)
	if(adj[y1[i]*degree+x1[i]] == y0[i])
	  break;
      
      for(x2[i]=0;x2[i]<degree;x2[i]++)
	if(adj[y2[i]*degree+x2[i]] == y3[i])
	  break;
      
      for(x3[i]=0;x3[i]<degree;x3[i]++)
	if(adj[y3[i]*degree+x3[i]] == y2[i])
	  break;
      
      if(x0[i] == degree || x1[i] == degree || x2[i] == degree || x3[i] == degree)
	ERROR("%d : %d %d %d %d\n", ii, x0[i], x1[i], x2[i], x3[i]);
    
      restored_adj_idx_y[i*4  ] = y0[i];
      restored_adj_idx_x[i*4  ] = x0[i];
      restored_adj_idx_y[i*4+1] = y1[i];
      restored_adj_idx_x[i*4+1] = x1[i];
      restored_adj_idx_y[i*4+2] = y2[i];
      restored_adj_idx_x[i*4+2] = x2[i];
      restored_adj_idx_y[i*4+3] = y3[i];
      restored_adj_idx_x[i*4+3] = x3[i];
      restored_adj_value[i*4  ] = adj[y0[i]*degree+x0[i]];
      restored_adj_value[i*4+1] = adj[y1[i]*degree+x1[i]];
      restored_adj_value[i*4+2] = adj[y2[i]*degree+x2[i]];
      restored_adj_value[i*4+3] = adj[y3[i]*degree+x3[i]];
      //
      restored_line[i*2  ] = tmp_line[i*2  ];
      restored_line[i*2+1] = tmp_line[i*2+1];
      restored_edge[i*4  ] = edge[tmp_line[i*2  ]][0];
      restored_edge[i*4+1] = edge[tmp_line[i*2  ]][1];
      restored_edge[i*4+2] = edge[tmp_line[i*2+1]][0];
      restored_edge[i*4+3] = edge[tmp_line[i*2+1]][1];
    }
    
#pragma omp parallel for
    for(int i=0;i<groups;i++){
      if(r==0){
	adj[y0[i]*degree+x0[i]] = y3[i]; adj[y1[i]*degree+x1[i]] = y2[i];
	adj[y2[i]*degree+x2[i]] = y1[i]; adj[y3[i]*degree+x3[i]] = y0[i];
      }
      else{
	adj[y0[i]*degree+x0[i]] = y2[i]; adj[y1[i]*degree+x1[i]] = y3[i];
	adj[y2[i]*degree+x2[i]] = y0[i]; adj[y3[i]*degree+x3[i]] = y1[i];
      }
    }
  }
  
#pragma omp parallel for
  for(int i=0;i<groups;i++){
    edge[tmp_line[i*2  ]][0] = tmp_edge[i*2  ][0];
    edge[tmp_line[i*2+1]][0] = tmp_edge[i*2+1][0];
    edge[tmp_line[i*2  ]][1] = tmp_edge[i*2  ][1];
    edge[tmp_line[i*2+1]][1] = tmp_edge[i*2+1][1];
  }
  
  *kind_opt = D_2G_OPT;
}

static bool accept(const int new_diam, const int current_diam, const double new_ASPL, const double current_ASPL,
		   const double temp, const int nodes, const int groups,
		   const bool hill_climbing_flag, const bool detect_temp_flag, const long long i,
		   double *max_diff_energy, long long *total_accepts, long long *accepts, long long *rejects)
{
  if(new_diam < current_diam){
    *accepts += 1;
    if(i > SKIP_ACCEPTS) *total_accepts +=1;
    return true;
  }
  else if(new_diam > current_diam){
    *rejects += 1;
    return false;
  }
  else{ //  new_diam == current_diam
    if(new_ASPL <= current_ASPL){
      *accepts += 1;
      if(i > SKIP_ACCEPTS) *total_accepts +=1;
      return true;
    }
    else if(hill_climbing_flag){ // Only accept when ASPL <= current_ASPL.
      *rejects += 1;
      return false;
    }
    
    double diff = ((current_ASPL-new_ASPL)*nodes*(nodes-1))/groups;
    
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
}

long long sa(const int nodes, const int lines, const int degree, const int groups, double temp, 
	     const long long ncalcs, const double cooling_rate,  const int low_diam,  const double low_ASPL, 
	     const bool hill_climbing_flag, const bool detect_temp_flag, double *max_diff_energy,
	     int edge[lines][2], int *diam, double *ASPL, const int cooling_cycle, const int added_centers,
	     const int based_nodes, long long *total_accepts, const bool is_simple_graph, const int algo)
{
  long long ii, accepts = 0, rejects = 0;
  int (*best_edge)[2]    = malloc(sizeof(int)*lines*2); // best_edge[lines][2]
  int (*tmp_edge)[2]     = malloc(sizeof(int)*lines*2); // tmp_edge[lines][2]
  int (*tmp_edge_nsg)[2] = malloc(sizeof(int)*lines*2); // tmp_edge_nsg[lines][2] /* nsg = not simple graph */
  int restored_adj_value[groups*4], restored_adj_idx_y[groups*4], restored_adj_idx_x[groups*4], kind_opt;
  int restored_edge[groups*4], restored_line[groups*2];
  bool restore_flag = false;
  copy_edge((int *)best_edge, (int *)edge, lines*2);
  copy_edge((int *)tmp_edge,  (int *)edge, lines*2);

  // Create adj matrix
  int *adj = malloc(sizeof(int)*nodes*degree); // int adj[nodes][degree];
  create_adj(nodes, lines, degree, (const int (*)[2])tmp_edge, (int (*)[degree])adj);
  evaluation(nodes, based_nodes, groups, lines, degree, adj, diam, ASPL, added_centers, algo);

  double current_ASPL = *ASPL;
  double best_ASPL    = *ASPL;
  int current_diam    = *diam;
  int best_diam       = *diam;
  int print_interval  = (ncalcs/NUM_OF_PROGRESS == 0)? 1 : ncalcs/NUM_OF_PROGRESS;
  if(rank == 0 && !detect_temp_flag)
    print_result_header();

  for(ii=0;ii<ncalcs;ii++){
    double tmp_ASPL;
    int tmp_diam;
    if(ii % print_interval == 0 && !detect_temp_flag){
      print_results(ii, temp, current_ASPL, best_ASPL, low_ASPL, 
		    current_diam, best_diam, low_diam, accepts, rejects);
      accepts = 0;
      rejects = 0;
    }

    while(1){
      if(is_simple_graph){
	if(restore_flag){
	  restore_adj(degree, groups, adj, kind_opt, restored_adj_value, restored_adj_idx_y, restored_adj_idx_x);
	  restore_edge(groups, kind_opt, (int *)tmp_edge, restored_line, restored_edge);
	}
      }
      else{
	copy_edge((int *)tmp_edge_nsg, (int *)tmp_edge, lines*2);
      }

      exchange_edge_2opt(nodes, lines, groups, degree, based_nodes, tmp_edge, added_centers,
			 adj, &kind_opt, restored_edge, restored_line, restored_adj_value,
      			 restored_adj_idx_y, restored_adj_idx_x, is_simple_graph, (int)ii);

      if(!is_simple_graph)
      	create_adj(nodes, lines, degree, (const int (*)[2])tmp_edge, (int (*)[degree])adj);

      assert(check(nodes, based_nodes, lines, degree, groups, tmp_edge, added_centers, adj, (int)ii));
      if(evaluation(nodes, based_nodes, groups, lines, degree, adj, &tmp_diam, &tmp_ASPL, added_centers, algo))
	break;
      else{
	if(is_simple_graph)
	  restore_flag = true;
	else
	  copy_edge((int *)tmp_edge, (int *)tmp_edge_nsg, lines*2);
      }
    }

    if(!accept(tmp_diam, current_diam, tmp_ASPL, current_ASPL, temp, nodes, groups, hill_climbing_flag,
	       detect_temp_flag, ii, max_diff_energy, total_accepts, &accepts, &rejects)){
      if(is_simple_graph)
	restore_flag = true;
      else
	copy_edge((int *)tmp_edge, (int *)tmp_edge_nsg, lines*2);
    }
    else{
      if(is_simple_graph) restore_flag = false;
      current_ASPL = tmp_ASPL;
      current_diam = tmp_diam;
      if((best_diam > current_diam) || (best_diam == current_diam && best_ASPL > current_ASPL)){
	copy_edge((int *)best_edge, (int *)tmp_edge, lines*2);
	best_ASPL = current_ASPL;
	best_diam = current_diam;
      }

      if(best_diam == current_diam && best_ASPL == low_ASPL){
	if(!detect_temp_flag){
	  print_results(ii, temp, current_ASPL, best_ASPL, low_ASPL, 
			current_diam, best_diam, low_diam, accepts, rejects);
	  PRINT_R0("---\nFound optimum solution.\n");
	}
	break;
      }
    }

    if((ii+1)%cooling_cycle == 0)
      temp *= cooling_rate;
  }

  *ASPL = best_ASPL;
  *diam = best_diam;
  copy_edge((int *)edge, (int *)best_edge, lines*2);
  free(adj);
  free(best_edge);
  free(tmp_edge);
  free(tmp_edge_nsg);
  return ii;
}

#define ESTIMATED_TIMES 5
double estimate_elapse_time(const int nodes, const int based_nodes, const int lines, const int degree,
			    const int groups, int edge[lines][2], const int added_centers,
			    const bool is_simple_graph, const int algo)
{
  int diam;    // Not use
  double ASPL; // Not use
  int *adj = malloc(sizeof(int)*nodes*degree); // int adj[nodes][degree];
  int (*tmp_edge)[2] = malloc(sizeof(int)*lines*2);      // int tmp_edge[lines][2];
  int kind_opt;
  int restored_adj_value[groups*4], restored_adj_idx_y[groups*4], restored_adj_idx_x[groups*4];
  int restored_edge[groups*4], restored_line[groups*2];

  copy_edge((int *)tmp_edge, (int *)edge, lines*2);
  create_adj(nodes, lines, degree, (const int (*)[2])tmp_edge, (int (*)[degree])adj);
  
  timer_start(TIMER_ESTIMATED);
  for(int i=0;i<ESTIMATED_TIMES;i++){
    exchange_edge_2opt(nodes, lines, groups, degree, based_nodes, tmp_edge, added_centers, adj,
		       &kind_opt, restored_edge, restored_line, restored_adj_value,
		       restored_adj_idx_y, restored_adj_idx_x, is_simple_graph, (int)i);
    if(!is_simple_graph)
      create_adj(nodes, lines, degree, (const int (*)[2])tmp_edge, (int (*)[degree])adj);
    assert(check(nodes, based_nodes, lines, degree, groups, tmp_edge, added_centers, adj, (int)i));
    evaluation(nodes, based_nodes, groups, lines, degree, adj, &diam, &ASPL, added_centers, algo);
  }  
  timer_stop(TIMER_ESTIMATED);
  
  free(tmp_edge);
  free(adj);

  return timer_read(TIMER_ESTIMATED)/ESTIMATED_TIMES;
}

// This function is mainly useful when groupe is 1.
void check_current_edge(const int nodes, const int degree, const int lines, const int groups, const int based_nodes,
			int edge[lines][2], const double low_ASPL, const int added_centers, const int algo)
{
  int diam;    // Not use
  double ASPL;
  int (*adj)[degree] = malloc(sizeof(int)*nodes*degree); // int adj[nodes][degree];

  create_adj(nodes, lines, degree, (const int (*)[2])edge, adj);
  if(! evaluation(nodes, based_nodes, groups, lines, degree, (int *)adj, &diam, &ASPL, added_centers, algo))
    ERROR("The input file has a node which is never reached by another node.\n");

  if(ASPL == low_ASPL)
    END("The input file has already optimum solution.\n");

  free(adj);
}
