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

bool has_duplicated_vertex(const int e00, const int e01, const int e10, const int e11)
{
  return (e00 == e10 || e01 == e11 || e00 == e11 || e01 == e10);
}

static void edge_exchange(const int nodes, const int lines, const int degree,
			  int edge[lines][2], int adjacency[nodes][degree], const int ii)
{
  int line[2];
  while(1){
    while(1){
      line[0] = getRandom(lines);
      line[1] = getRandom(lines);
      if(line[0] == line[1]) continue;
      if(has_duplicated_vertex(edge[line[0]][0], edge[line[0]][1], edge[line[1]][0], edge[line[1]][1]))
        continue;
      else
	break;
    }

    int tmp_edge[2][2];
    tmp_edge[0][0] = edge[line[0]][0]; tmp_edge[0][1] = edge[line[0]][1];
    tmp_edge[1][0] = edge[line[1]][0]; tmp_edge[1][1] = edge[line[1]][1];
    int r = (getRandom(2) == 0)? 1 : 0;;
    swap(&tmp_edge[0][1], &tmp_edge[1][r]);

    if(!check_duplicate_current_edge(lines, edge, line, tmp_edge))
      continue;
  
    edge[line[0]][0] = tmp_edge[0][0]; edge[line[0]][1] = tmp_edge[0][1];
    edge[line[1]][0] = tmp_edge[1][0]; edge[line[1]][1] = tmp_edge[1][1];
    break;
  }
}

static bool accept(const double ASPL, const double current_ASPL, const double temp, const int nodes, 
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
  
  double diff = (current_ASPL-ASPL)*nodes*(nodes-1);

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

long long sa(const int nodes, const int lines, double temp, const long long ncalcs,
	     const double cooling_rate, const int low_diam,  const double low_ASPL, 
	     const bool hill_climbing_flag, const bool detect_temp_flag,
	     double *max_diff_energy, int edge[lines][2], int *diam, double *ASPL,
	     const int cooling_cycle, long long *total_accepts)
{
  int degree = 2 * lines / nodes;
  int best_edge[lines][2], tmp_edge[lines][2];
  long long i, accepts = 0, rejects = 0;
  edge_copy((int *)best_edge, (int *)edge, lines*2);

  // Create adjacency matrix
  int (*adjacency)[degree] = malloc(sizeof(int)*nodes*degree); // int adjacency[nodes][degree];
  create_adjacency(nodes, lines, degree, edge, adjacency);
  evaluation(nodes, lines, degree, adjacency, diam, ASPL);
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
      edge_exchange(nodes, lines, degree, tmp_edge, adjacency, (int)i);
      create_adjacency(nodes, lines, degree, tmp_edge, adjacency);
      if(evaluation(nodes, lines, degree, adjacency, diam, ASPL)) break;
    }

    if(accept(*ASPL, current_ASPL, temp, nodes, hill_climbing_flag, detect_temp_flag,
	      i, max_diff_energy, total_accepts, &accepts, &rejects)){
      current_ASPL  = *ASPL;
      current_diam  = *diam;
      edge_copy((int *)edge, (int *)tmp_edge, lines*2);
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
double estimated_elapse_time(const int nodes, const int lines, int edge[lines][2])
{
  int degree = 2 * lines / nodes;
  int diam;    // Not use
  double ASPL; // Not use
  int (*adjacency)[degree] = malloc(sizeof(int)*nodes*degree); // int adjacency[nodes][degree];
  int (*tmp_edge)[2]       = malloc(sizeof(int)*lines*2);      // int tmp_edge[lines][2];
  
  edge_copy((int *)tmp_edge, (int *)edge, lines*2);
  create_adjacency(nodes, lines, degree, tmp_edge, adjacency);

  timer_start(TIMER_ESTIMATED);
  for(int i=0;i<ESTIMATED_TIMES;i++){
    edge_copy((int *)tmp_edge, (int *)edge, lines*2);
    edge_exchange(nodes, lines, degree, tmp_edge, adjacency, (int)i);
    evaluation(nodes, lines, degree, adjacency, &diam, &ASPL);
  }
  timer_stop(TIMER_ESTIMATED);
  
  free(tmp_edge);
  free(adjacency);

  return timer_read(TIMER_ESTIMATED)/ESTIMATED_TIMES;
}

// This function is mainly useful when groupe is 1.
void check_current_edge(const int nodes, const int lines, int edge[lines][2], const double low_ASPL)
{
  int degree = 2 * lines / nodes;
  int diam;    // Not use
  double ASPL;
  int (*adjacency)[degree] = malloc(sizeof(int)*nodes*degree); // int adjacency[nodes][degree];

  create_adjacency(nodes, lines, degree, edge, adjacency);
  if(! evaluation(nodes, lines, degree, adjacency, &diam, &ASPL))
    ERROR("The input file has a node which is never reached by another node.\n");

  if(ASPL == low_ASPL)
    END("The input file has already optimum solution.\n");

  free(adjacency);
}
