#include "common.h"

#define RIGHT  0
#define LEFT   1
#define MIDDLE 2

static bool has_duplicated_edge(const int e00, const int e01, const int e10, const int e11)
{
  return ((e00 == e10 && e01 == e11) || (e00 == e11 && e01 == e10));
}

static bool not_has_duplicated_edge(const int e00, const int e01, const int e10, const int e11)
{
  return !(has_duplicated_edge(e00, e01, e10, e11));
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

static void create_adjacency(const int nodes, const int lines, const int degree, 
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

static void check(const int nodes, const int lines, int edge[lines][2], int ii)
{
  int based_lines = lines/2;
  for(int i=0;i<based_lines;i++){
    int j = i + based_lines;
    if(distance(nodes, edge[i][0], edge[i][1]) != distance(nodes, edge[j][0], edge[j][1])){
      printf("check 1: %d\n", ii);
      printf("edge[%d][0] = %d : edge[%d][1] = %d d=%d\n", i, edge[i][0], i, edge[i][1], distance(nodes, edge[i][0], edge[i][1]));
      printf("edge[%d][0] = %d : edge[%d][1] = %d d=%d\n", j, edge[j][0], j, edge[j][1], distance(nodes, edge[j][0], edge[j][1]));
      exit(0);
    }
  }

  for(int i=0;i<based_lines;i++){
    int j = i + based_lines;
    if(order(nodes, edge[i][0], edge[i][1]) != order(nodes, edge[j][0], edge[j][1])){
      printf("check 2 : %d\n", ii);
      printf("edge[%d][0] = %d : edge[%d][1] = %d %d\n", i, edge[i][0], i, edge[i][1], order(nodes, edge[i][0], edge[i][1]));
      printf("edge[%d][0] = %d : edge[%d][1] = %d %d\n", j, edge[j][0], j, edge[j][1], order(nodes, edge[j][0], edge[j][1]));
      exit(0);
    }
  }

  for(int i=0;i<lines;i++){
    for(int j=i+1;j<lines;j++){
      if((edge[i][0] == edge[j][0] && edge[i][1] == edge[j][1]) ||
         (edge[i][0] == edge[j][1] && edge[i][1] == edge[j][0])){
        printf("check 3: %d\n", ii);
        printf("The same node conbination in the edge. %d %d\n", i, j);
        exit(0);
      }
    }
  }
}

static bool target_line(const int i, const int groups, const int line[2], const int opp_line[groups*2-2])
{
  if(i == line[0] || i == line[1]) 
    return false;

  for(int j=0;j<groups*2-2;j++)
    if(i == opp_line[j] || i == opp_line[j])
      return false;

  return true;
}

static void edge_exchange(const int nodes, const int lines, const int groups, 
			  int edge[lines][2], const int ii)
{
  int line[2], opp_line[groups*2-2], tmp_edge[groups*2][2];

  while(1){
    while(1){
      line[0] = rand() % lines;
      line[1] = rand() % lines;
      if(line[1] == line[0])	                         continue;
      else if(edge[line[0]][0] == edge[line[1]][0])      continue;
      else if(edge[line[0]][0] == edge[line[1]][1])      continue;
      else if(edge[line[0]][1] == edge[line[1]][0])      continue;
      else if(edge[line[0]][1] == edge[line[1]][1])      continue;
      else if(abs(line[0] - line[1]) == lines/2 && groups%2 == 0){
	if(rand()%2 == 0){
	  tmp_edge[0][0] = edge[line[0]][0]; tmp_edge[0][1] = edge[line[1]][1];
	  tmp_edge[1][0] = edge[line[1]][0]; tmp_edge[1][1] = edge[line[0]][1];
	}
	else{
	  tmp_edge[0][0] = edge[line[0]][0]; tmp_edge[0][1] = edge[line[1]][0];
          tmp_edge[1][0] = edge[line[0]][1]; tmp_edge[1][1] = edge[line[1]][1];
	}

	bool flag = true;
	for(int i=0;i<lines;i++){
	  if(i != line[0] && i != line[1]){ // Not needed
	    for(int j=0;j<2;j++){
	      if(has_duplicated_edge(edge[i][0], edge[i][1], tmp_edge[j][0], tmp_edge[j][1])){
		flag = false;
		break;
	      }
	    }
	  }
	  if(!flag) break;
	}
	if(!flag) continue;
	
	edge[line[0]][0] = tmp_edge[0][0]; edge[line[0]][1] = tmp_edge[0][1];
	edge[line[1]][0] = tmp_edge[1][0]; edge[line[1]][1] = tmp_edge[1][1];
	
	if(order(nodes, edge[line[0]][0], edge[line[0]][1]) != order(nodes, edge[line[1]][0], edge[line[1]][1]))
	  swap(&edge[line[1]][0], &edge[line[1]][1]);

	return;
      }
      else break;
    }

    for(int i=0;i<2;i++)
      opp_line[i] = (line[i] >= lines/2)? line[i]-lines/2 : line[i]+lines/2;

    bool flag0 = (distance(nodes, edge[line[0]][0], edge[line[0]][1]) == nodes/2 && groups%2 == 0);
    bool flag1 = (distance(nodes, edge[line[1]][0], edge[line[1]][1]) == nodes/2 && groups%2 == 0);
    bool double_diameter_flag = (flag0 && flag1);
    bool single_diameter_flag = (flag0 && !flag1) || (!flag0 && flag1);

    if(double_diameter_flag){
      int pattern = rand() % 4;
      if(pattern == 0){
	tmp_edge[0][0] = edge[line[0]][0];       tmp_edge[0][1] = edge[line[1]][1];
	tmp_edge[1][0] = edge[line[1]][0];       tmp_edge[1][1] = edge[line[0]][1];
	tmp_edge[2][0] = edge[opp_line[0]][0];   tmp_edge[2][1] = edge[opp_line[0]][1];
	tmp_edge[3][0] = edge[opp_line[1]][0];   tmp_edge[3][1] = edge[opp_line[1]][1];
      }
      else if(pattern == 1){
	tmp_edge[0][0] = edge[line[0]][0];       tmp_edge[0][1] = edge[line[1]][0];
	tmp_edge[1][0] = edge[line[0]][1];       tmp_edge[1][1] = edge[line[1]][1];
	tmp_edge[2][0] = edge[opp_line[0]][0];   tmp_edge[2][1] = edge[opp_line[0]][1];
	tmp_edge[3][0] = edge[opp_line[1]][0];   tmp_edge[3][1] = edge[opp_line[1]][1];
      }
      else if(pattern == 2){
	tmp_edge[0][0] = edge[line[0]][0];       tmp_edge[0][1] = edge[line[0]][1];
	tmp_edge[1][0] = edge[line[1]][0];       tmp_edge[1][1] = edge[line[1]][1];
	tmp_edge[2][0] = edge[opp_line[0]][0];   tmp_edge[2][1] = edge[opp_line[1]][1];
	tmp_edge[3][0] = edge[opp_line[1]][0];   tmp_edge[3][1] = edge[opp_line[0]][1];
      }
      else{
	tmp_edge[0][0] = edge[line[0]][0];       tmp_edge[0][1] = edge[line[0]][1];
	tmp_edge[1][0] = edge[line[1]][0];       tmp_edge[1][1] = edge[line[1]][1];
	tmp_edge[2][0] = edge[opp_line[0]][0];   tmp_edge[2][1] = edge[opp_line[1]][0];
	tmp_edge[3][0] = edge[opp_line[0]][1];   tmp_edge[3][1] = edge[opp_line[1]][1];
      }
      swap(&tmp_edge[1][0], &tmp_edge[2][0]);
      swap(&tmp_edge[1][1], &tmp_edge[2][1]);
    }
    else if(single_diameter_flag){
      if(flag0){
	swap(&line[0], &line[1]);
	swap(&opp_line[0], &opp_line[1]);
      }

      int pattern = rand() % 8;
      if(pattern == 0){
	tmp_edge[0][0] = edge[line[1]][0];     tmp_edge[0][1] = edge[line[0]][1];
	tmp_edge[1][0] = edge[line[0]][0];     tmp_edge[1][1] = edge[opp_line[0]][0];
	tmp_edge[2][0] = edge[line[1]][1];     tmp_edge[2][1] = edge[opp_line[0]][1];
	tmp_edge[3][0] = edge[opp_line[1]][0]; tmp_edge[3][1] = edge[opp_line[1]][1];
      }
      else if(pattern == 1){
	tmp_edge[0][0] = edge[line[1]][1];     tmp_edge[0][1] = edge[line[0]][1];
	tmp_edge[1][0] = edge[line[0]][0];     tmp_edge[1][1] = edge[opp_line[0]][0];
	tmp_edge[2][0] = edge[line[1]][0];     tmp_edge[2][1] = edge[opp_line[0]][1];
	tmp_edge[3][0] = edge[opp_line[1]][0]; tmp_edge[3][1] = edge[opp_line[1]][1];
      }
      else if(pattern == 2){
	tmp_edge[0][0] = edge[line[0]][0];     tmp_edge[0][1] = edge[line[1]][0];
	tmp_edge[1][0] = edge[line[0]][1];     tmp_edge[1][1] = edge[opp_line[0]][1];
	tmp_edge[2][0] = edge[opp_line[0]][0]; tmp_edge[2][1] = edge[line[1]][1];
	tmp_edge[3][0] = edge[opp_line[1]][0]; tmp_edge[3][1] = edge[opp_line[1]][1];
      }
      else if(pattern == 3){
	tmp_edge[0][0] = edge[line[0]][0];     tmp_edge[0][1] = edge[line[1]][1];
	tmp_edge[1][0] = edge[line[0]][1];     tmp_edge[1][1] = edge[opp_line[0]][1];
	tmp_edge[2][0] = edge[opp_line[0]][0]; tmp_edge[2][1] = edge[line[1]][0];
	tmp_edge[3][0] = edge[opp_line[1]][0]; tmp_edge[3][1] = edge[opp_line[1]][1];
      }
      else if(pattern == 4){
	tmp_edge[0][0] = edge[opp_line[1]][0]; tmp_edge[0][1] = edge[line[0]][1];
	tmp_edge[1][0] = edge[line[1]][0];     tmp_edge[1][1] = edge[line[1]][1];
	tmp_edge[2][0] = edge[opp_line[1]][1]; tmp_edge[2][1] = edge[opp_line[0]][1];
	tmp_edge[3][0] = edge[line[0]][0];     tmp_edge[3][1] = edge[opp_line[0]][0];
      }
      else if(pattern == 5){
	tmp_edge[0][0] = edge[opp_line[1]][1]; tmp_edge[0][1] = edge[line[0]][1];
	tmp_edge[1][0] = edge[line[1]][0];     tmp_edge[1][1] = edge[line[1]][1];
	tmp_edge[2][0] = edge[opp_line[1]][0]; tmp_edge[2][1] = edge[opp_line[0]][1];
	tmp_edge[3][0] = edge[line[0]][0];     tmp_edge[3][1] = edge[opp_line[0]][0];
      }
      else if(pattern == 6){
	tmp_edge[0][0] = edge[line[0]][0];     tmp_edge[0][1] = edge[opp_line[1]][0];
	tmp_edge[1][0] = edge[line[1]][0];     tmp_edge[1][1] = edge[line[1]][1];
	tmp_edge[2][0] = edge[opp_line[0]][0]; tmp_edge[2][1] = edge[opp_line[1]][1];
	tmp_edge[3][0] = edge[line[0]][1];     tmp_edge[3][1] = edge[opp_line[0]][1];
      }
      else if(pattern == 7){
	tmp_edge[0][0] = edge[line[0]][0];     tmp_edge[0][1] = edge[opp_line[1]][1];
	tmp_edge[1][0] = edge[line[1]][0];     tmp_edge[1][1] = edge[line[1]][1];
	tmp_edge[2][0] = edge[opp_line[0]][0]; tmp_edge[2][1] = edge[opp_line[1]][0];
	tmp_edge[3][0] = edge[line[0]][1];     tmp_edge[3][1] = edge[opp_line[0]][1];
      }
    }
    else{
      if(rand()%2 == 0){
	tmp_edge[0][0] = edge[line[0]][0];       tmp_edge[0][1] = edge[line[1]][1];
	tmp_edge[1][0] = edge[line[1]][0];       tmp_edge[1][1] = edge[line[0]][1];
	tmp_edge[2][0] = edge[opp_line[0]][0];   tmp_edge[2][1] = edge[opp_line[1]][1];
	tmp_edge[3][0] = edge[opp_line[1]][0];   tmp_edge[3][1] = edge[opp_line[0]][1];
      }
      else{
	tmp_edge[0][0] = edge[line[0]][0];       tmp_edge[0][1] = edge[line[1]][0];
	tmp_edge[1][0] = edge[line[0]][1];       tmp_edge[1][1] = edge[line[1]][1];
	tmp_edge[2][0] = edge[opp_line[0]][0];   tmp_edge[2][1] = edge[opp_line[1]][0];
	tmp_edge[3][0] = edge[opp_line[0]][1];   tmp_edge[3][1] = edge[opp_line[1]][1];
      }
    }

    // Remove loop
    bool flag = true;
    for(int i=0;i<groups*2;i++)
      flag &= (tmp_edge[i][0] != tmp_edge[i][1]);

    if(!flag) continue;

    // Remove duplicate edge in tmp_edges
    for(int i=0;i<groups*2;i++)
      for(int j=i+1;j<groups*2;j++)
	flag &= not_has_duplicated_edge(tmp_edge[i][0], tmp_edge[i][1], tmp_edge[j][0], tmp_edge[j][1]);

    if(!flag) continue;

    // Remove duplicate edge in current edges
    for(int i=0;i<lines;i++)
      if(target_line(i, groups, line, opp_line))
        for(int j=0;j<groups*2;j++)
          flag &= not_has_duplicated_edge(edge[i][0], edge[i][1], tmp_edge[j][0], tmp_edge[j][1]);

    if(!flag) continue;

    edge[line[0]][0] = tmp_edge[0][0];   edge[line[0]][1] = tmp_edge[0][1];
    edge[line[1]][0] = tmp_edge[1][0];   edge[line[1]][1] = tmp_edge[1][1];

    for(int i=0;i<groups*2-2;i++){
      edge[opp_line[i]][0] = tmp_edge[i+2][0];
      edge[opp_line[i]][1] = tmp_edge[i+2][1];
    }
    
    break;
  }

  for(int i=0;i<2;i++)
    if(order(nodes, edge[line[i]][0], edge[line[i]][1]) != order(nodes, edge[opp_line[i]][0], edge[opp_line[i]][1]))
      swap(&edge[line[i]][0], &edge[line[i]][1]);
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

      check(nodes, lines, current_edge, (int)i);
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
