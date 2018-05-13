#include "common.h"

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
  fflush(stdout);
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

static int distance(const int nodes, const int a, const int b)
{
  int v = MAX(a, b) - MIN(a, b);
  return (v < nodes/2)? v : nodes - v;
}

static bool check(const int nodes, const int lines, const int degree, const int groups, int edge[lines][2], int ii)
{
  //  verfy_regular_graph(nodes, degree, lines, edge);

  int based_lines = lines/groups;
  for(int i=0;i<based_lines;i++){
    for(int j=1;j<groups;j++){
      int k = j * based_lines + i;
      if(distance(nodes, edge[i][0], edge[i][1]) != distance(nodes, edge[k][0], edge[k][1])){
	printf("check 1: %d\n", ii);
	printf("edge[%d][0] = %d : edge[%d][1] = %d d=%d\n", i, edge[i][0], i, edge[i][1], distance(nodes, edge[i][0], edge[i][1]));
	printf("edge[%d][0] = %d : edge[%d][1] = %d d=%d\n", k, edge[k][0], k, edge[k][1], distance(nodes, edge[k][0], edge[k][1]));
	return false;
      }
    }
  }

  for(int i=0;i<based_lines;i++){
    for(int j=1;j<groups;j++){
      int k = j * based_lines + i;
      if(order(nodes, edge[i][0], edge[i][1]) != order(nodes, edge[k][0], edge[k][1])){
	printf("check 2 : %d\n", ii);
	printf("edge[%d][0] = %d : edge[%d][1] = %d %d\n", i, edge[i][0], i, edge[i][1], order(nodes, edge[i][0], edge[i][1]));
	printf("edge[%d][0] = %d : edge[%d][1] = %d %d\n", k, edge[k][0], k, edge[k][1], order(nodes, edge[k][0], edge[k][1]));
	return false;
      }
    }
  }

  for(int i=0;i<lines;i++){
    for(int j=i+1;j<lines;j++){
      if((edge[i][0] == edge[j][0] && edge[i][1] == edge[j][1]) ||
         (edge[i][0] == edge[j][1] && edge[i][1] == edge[j][0])){
        printf("check 3: %d\n", ii);
        printf("The same node conbination in the edge. %d %d\n", i, j);
	return false;
      }
    }
  }

  int based_nodes = nodes / groups;
  for(int i=0;i<based_lines;i++)
    if(order(nodes, edge[i][0], edge[i][1]) != MIDDLE)
      for(int j=0;j<groups-1;j++){
	int n = based_lines * j + i;
	int tmp0 = edge[n+based_lines][0] - edge[n][0];
	int tmp1 = edge[n+based_lines][1] - edge[n][1];
	tmp0 = (tmp0 < 0)? tmp0 + nodes : tmp0;
	tmp1 = (tmp1 < 0)? tmp1 + nodes : tmp1;
	if(tmp0 != based_nodes || tmp1 != based_nodes){
	  printf("check 4: %d\n", ii);
	  printf("The different group relationship\n");
	  for(int k=0;k<groups;k++){
	    int m = based_lines * k + i;
	    printf("edge[%2d][0], edge[%2d][1] = %2d, %2d\n", m, m, edge[m][0], edge[m][1]);
	  }
	  return false;
	}
      }

  return true;
}

static void edge_exchange(const int nodes, const int lines, const int groups, 
			  int edge[lines][2], const int ii)
{
  int line[groups*2], tmp_edge[groups*2][2];
  int based_nodes = nodes / groups;
  int based_lines = lines / groups;
  
  while(1){
    while(1){
      line[0] = getRandom(lines);
      line[1] = getRandom(lines);
      if(line[1] == line[0])	                      continue;
      else if(edge[line[0]][0] == edge[line[1]][0])   continue;
      else if(edge[line[0]][0] == edge[line[1]][1])   continue;
      else if(edge[line[0]][1] == edge[line[1]][0])   continue;
      else if(edge[line[0]][1] == edge[line[1]][1])   continue;
      else if(abs(line[0] - line[1]) % based_lines == 0){
	int start_line = line[0] % based_lines;
	if(edge_exchange_among_groups(based_nodes, based_lines, edge,
				      groups, start_line)){
	  return;
	}
	else
	  continue;
      }
      else break;
    }

    for(int i=1;i<groups;i++){
      int tmp0 = line[0] + based_lines * i;
      int tmp1 = line[1] + based_lines * i;
      line[0+2*i] = (tmp0 >= lines)? tmp0 - lines : tmp0;
      line[1+2*i] = (tmp1 >= lines)? tmp1 - lines : tmp1;
    }

    for(int i=0;i<groups*2;i++)
      for(int j=0;j<2;j++)
	tmp_edge[i][j] = edge[line[i]][j];

    bool flag0 = (distance(nodes, edge[line[0]][0], edge[line[0]][1]) == nodes/2);
    bool flag1 = (distance(nodes, edge[line[1]][0], edge[line[1]][1]) == nodes/2);
    bool double_diameter_flag = (flag0 && flag1);
    bool single_diameter_flag = (flag0 && !flag1) || (!flag0 && flag1);

    if(double_diameter_flag){
      int start_line = line[getRandom(2)] % based_lines;
      if(edge_exchange_among_groups(based_nodes, based_lines, edge,
				    groups, start_line))
	return;
      else
	continue;
    }
    else if(single_diameter_flag){
      if(flag0)
	for(int i=0;i<groups;i++)
	  for(int j=0;j<2;j++)
	    swap(&tmp_edge[i*2][j], &tmp_edge[i*2+1][j]);
      
      int pivot_lineno = getRandom(groups) * 2 + 1;
      int pivot[2] = {tmp_edge[pivot_lineno][0], tmp_edge[pivot_lineno][1]};
      int not_pivot[2] = {};

      int e0 = pivot[0] % based_nodes;
      int e1 = pivot[1] % based_nodes;
      for(int i=0;i<groups;i++){
	int n  = i * 2 + 1;
	if(pivot_lineno == n) continue;
	int e2 = tmp_edge[n][0] % based_nodes;
	int e3 = tmp_edge[n][1] % based_nodes;
	if(!has_duplicated_edge(e0, e1, e2, e3)){
	  not_pivot[0] = tmp_edge[n][0];
	  not_pivot[1] = tmp_edge[n][1];
	  break;
	}
      }

      int rand_offset = getRandom(groups/2);
      int tmp0 = pivot[0] + based_nodes * rand_offset;
      int tmp1 = pivot[1] + based_nodes * rand_offset;
      pivot[0] = (tmp0 >= nodes)? tmp0 - nodes : tmp0;
      pivot[1] = (tmp1 >= nodes)? tmp1 - nodes : tmp1;
      int pattern = getRandom(4);

      if(pattern == 0){
	swap(&tmp_edge[0][0],      &pivot[0]);
	swap(&tmp_edge[groups][0], &pivot[1]);
      }
      else if(pattern == 1){
	swap(&tmp_edge[0][0],      &pivot[1]);
	swap(&tmp_edge[groups][0], &pivot[0]);
      }
      else if(pattern == 2){
	swap(&tmp_edge[0][1],      &pivot[0]);
	swap(&tmp_edge[groups][1], &pivot[1]);
      }
      else if(pattern == 3){
	swap(&tmp_edge[0][1],      &pivot[1]);
	swap(&tmp_edge[groups][1], &pivot[0]);
      }

      for(int i=1;i<groups;i++){
	int n = i * 2;
	int tmp0 = tmp_edge[0][0] + based_nodes * i;
	int tmp1 = tmp_edge[0][1] + based_nodes * i;
	tmp_edge[n][0] = (tmp0 >= nodes)? tmp0 - nodes : tmp0;
	tmp_edge[n][1] = (tmp1 >= nodes)? tmp1 - nodes : tmp1;
      }

      for(int i=0;i<groups/2;i++){
	int n = i * 2 + 1;
	int tmp0 = pivot[0] + based_nodes * i;
	int tmp1 = pivot[1] + based_nodes * i;
	tmp_edge[n][0] = (tmp0 >= nodes)? tmp0 - nodes : tmp0;
        tmp_edge[n][1] = (tmp1 >= nodes)? tmp1 - nodes : tmp1;
      }

      for(int i=0;i<groups/2;i++){
	int n = i * 2 + 1 + groups;
	int tmp0 = not_pivot[0] + based_nodes * i;
	int tmp1 = not_pivot[1] + based_nodes * i;
	tmp_edge[n][0] = (tmp0 >= nodes)? tmp0 - nodes : tmp0;
	tmp_edge[n][1] = (tmp1 >= nodes)? tmp1 - nodes : tmp1;
      }

      // Exchange pivots
      if(pattern >= 4)
	for(int i=0;i<groups/2;i++)
	  for(int j=0;j<2;j++)
	    swap(&tmp_edge[i*2+1][j], &tmp_edge[groups+i*2+1][j]);

      if(pattern%2 == 0)
	for(int i=0;i<groups/2;i++)
	  swap(&tmp_edge[i*2+1][0], &tmp_edge[i*2+1][1]);

      if(pattern%2 == 0){
	swap(&tmp_edge[0][0], &tmp_edge[groups-1][0]);

	for(int i=0;i<groups/2-1;i++)
	  swap(&tmp_edge[i*2+1][0], &tmp_edge[groups+(i+1)*2][0]);
	
	for(int i=0;i<groups/2;i++)
	  swap(&tmp_edge[i*2+1][1], &tmp_edge[i*2+2][0]);
      } 
      else{
	swap(&tmp_edge[0][1], &tmp_edge[groups-1][0]);

	for(int i=0;i<groups/2-1;i++)
	  swap(&tmp_edge[i*2+1][0], &tmp_edge[groups+(i+1)*2][1]);
	
	for(int i=0;i<groups/2;i++)
	  swap(&tmp_edge[i*2+1][1], &tmp_edge[i*2+2][1]);
      }
    }
    else{
      if(getRandom(2) == 0){
	for(int i=0;i<groups;i++)
	  swap(&tmp_edge[i*2][1], &tmp_edge[i*2+1][1]);
      }
      else{
	for(int i=0;i<groups;i++)
	  swap(&tmp_edge[i*2][1], &tmp_edge[i*2+1][0]);
      }
    }

    if(!check_loop(groups*2, tmp_edge))            continue;
    if(!check_duplicate_edge(groups*2, tmp_edge))  continue;
    if(!check_duplicate_current_edge(lines, groups*2, line, edge, tmp_edge))
      continue;

    for(int i=0;i<groups*2;i++)
      if(order(nodes, tmp_edge[i][0], tmp_edge[i][1]) == RIGHT)
	swap(&tmp_edge[i][0], &tmp_edge[i][1]);  // RIGHT -> LEFT

    //    sort_line(groups*2,     line);  // fix me
    //    sort_tmp_edge(groups*2, tmp_edge);

    for(int i=0;i<groups*2;i++){
      edge[line[i]][0] = tmp_edge[i][0];
      edge[line[i]][1] = tmp_edge[i][1];
    }

    break;
  }
}

static bool accept(const double ASPL, const double current_ASPL, const double temp, const int nodes,
		   const bool hill_climbing_flag, const bool detect_temp_flag, double *max_diff_energy)
{
  if(ASPL <= current_ASPL) return true;
  if(hill_climbing_flag)   return false;

  double probability = (double)rand()/RAND_MAX;
  double diff = (current_ASPL - ASPL) * nodes * (nodes-1);

  if(detect_temp_flag){
    if(abs(diff) > *max_diff_energy){
      *max_diff_energy = abs(diff);
    }
  }
  if(exp(diff / temp) > probability)
    return true;
  else
    return false;
}

long long sa(const int nodes, const int lines, const int degree, const int groups, double temp, 
	     const long long ncalcs, const double cooling_rate, 
	     const int low_diam, const double low_ASPL, const bool hill_climbing_flag,
	     const bool detect_temp_flag, double *max_diff_energy,
	     int edge[lines][2], int *diam, double *ASPL, const int rank, const int size)
{
  bool complete;
  int current_edge[lines][2], best_edge[lines][2];
  size_t size_edge = sizeof(int) * lines * 2;
  long long i;
  memcpy(best_edge, edge, size_edge);

  // Create adjacency matrix
  int (*adjacency)[degree] = malloc(sizeof(int)*nodes*degree); // int adjacency[nodes][degree];
  if(rank == 0)
    create_adjacency(nodes, lines, degree, edge, adjacency);

  MPI_Bcast(adjacency, nodes*degree, MPI_INT, 0, MPI_COMM_WORLD);

  evaluation(nodes, lines, degree, adjacency, diam, ASPL, rank, size);
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
    do{
      memcpy(current_edge, edge, size_edge);

      if(rank == 0)
	edge_exchange(nodes, lines, groups, current_edge, (int)i);

      assert(check(nodes, lines, degree, groups, current_edge, (int)i));
      create_adjacency(nodes, lines, degree, current_edge, adjacency);

      MPI_Bcast(adjacency, nodes*degree, MPI_INT, 0, MPI_COMM_WORLD);
      complete = evaluation(nodes, lines, degree, adjacency, diam, ASPL, rank, size);
    } while(!complete);

    if(accept(*ASPL, current_ASPL, temp, nodes, hill_climbing_flag, detect_temp_flag, max_diff_energy)){
      current_ASPL = *ASPL;
      current_diam = *diam;

      memcpy(edge, current_edge, size_edge);
      if(best_ASPL > current_ASPL && best_diam >= current_diam)
	memcpy(best_edge, current_edge, size_edge);

      best_ASPL = MIN(best_ASPL, current_ASPL);
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
    else{
    }

    temp *= cooling_rate;
  }

  *ASPL = best_ASPL;
  *diam = best_diam;
  memcpy(edge, best_edge, size_edge);

  free(adjacency);

  return i;
}

#define ESTIMATED_TIMES 5
double estimated_elapse_time(const long long ncals, const int nodes, const int lines, const int degree,
			     int edge[lines][2], const int rank, const int size)
{
  int diam;    // Not use
  double ASPL; // Not use
  int (*adjacency)[degree] = malloc(sizeof(int)*nodes*degree); // int adjacency[nodes][degree];

  timer_start(TIMER_ESTIMATED);
  for(int i=0;i<ESTIMATED_TIMES;i++){
    if(rank == 0)
      create_adjacency(nodes, lines, degree, edge, adjacency);

    MPI_Bcast(adjacency, nodes*degree, MPI_INT, 0, MPI_COMM_WORLD);
    evaluation(nodes, lines, degree, adjacency, &diam, &ASPL, rank, size);
  }
  timer_stop(TIMER_ESTIMATED);

  free(adjacency);

  return timer_read(TIMER_ESTIMATED)/ESTIMATED_TIMES;
}

void check_current_edge(const int nodes, const int degree, const int lines, 
			int edge[lines][2], const double low_ASPL, const int rank, const int size)
{
  int diam;    // Not use
  double ASPL;
  int (*adjacency)[degree] = malloc(sizeof(int)*nodes*degree); // int adjacency[nodes][degree];
  if(rank == 0)
    create_adjacency(nodes, lines, degree, edge, adjacency);

  MPI_Bcast(adjacency, nodes*degree, MPI_INT, 0, MPI_COMM_WORLD);

  if(evaluation(nodes, lines, degree, adjacency, &diam, &ASPL, rank, size) == false){
    printf("The input file has a node which is never reached by another node.\n");
    ABORT;
  }
  if(ASPL == low_ASPL){
    printf("The input file has already optimum solution.\n");
    MPI_Finalize();
    exit(0);
  }

  free(adjacency);
}
