#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <stdint.h>
#include <omp.h>
#include <stdbool.h>
#include <assert.h>
#include <time.h>

#define MAX_FILENAME_LENGTH 256
#define NUM_TIMERS   2
#define TIMER_BFS    0
#define TIMER_ADJ    1
#define NOT_VISITED -1
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
static double elapsed[NUM_TIMERS], start[NUM_TIMERS];

static void verify(const int lines, const int degree, const int nodes, int edge[lines][2])
{
  if((2*lines)%degree != 0){
    fprintf(stderr,"NG. lines or n or d is invalid. lines = %d nodes = %d degree = %d\n", lines, nodes, degree);
    exit(1);
  }
  
  int n[nodes];
  for(int i=0;i<nodes;i++)
    n[i] = 0;

  for(int i=0;i<lines;i++){
    n[edge[i][0]]++;
    n[edge[i][1]]++;
  }

  for(int i=0;i<nodes;i++)
    if(degree != n[i]){
      fprintf(stderr, "NG\nNot regular graph. degree = %d n[%d] = %d\n", degree, i, n[i]);
      exit(1);
    }

  printf("Verification : OK\n");
}

// This function is inherited from "http://research.nii.ac.jp/graphgolf/py/create-random.py".
static void lower_bound_of_diam_aspl(int *low_diam, double *low_ASPL, const int nodes, const int degree)
{
  int diam = -1, n = 1, r = 1;
  double aspl = 0.0;

  while(1){
    int tmp = n + degree * pow(degree-1, r-1);
    if(tmp >= nodes)
      break;

    n = tmp;
    aspl += r * degree * pow(degree-1, r-1);
    diam = r++;
  }

  diam++;
  aspl += diam * (nodes - n);
  aspl /= (nodes - 1);

  *low_diam = diam;
  *low_ASPL = aspl;
}

static double elapsed_time()
{
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec + 1.0e-6 * t.tv_usec;
}

void timer_clear_all()
{
  for(int i=0;i<NUM_TIMERS;i++)
    elapsed[i] = 0.0;
}

void timer_clear(const int n)
{
  elapsed[n] = 0.0;
}

void timer_start(const int n)
{
  start[n] = elapsed_time();
}

void timer_stop(const int n)
{
  double now = elapsed_time();
  double t = now - start[n];
  elapsed[n] += t;
}

double timer_read(const int n)
{
  return( elapsed[n] );
}

static void clear_buffer(const int n, int buffer[n], int value)
{
#pragma omp parallel for
  for(int i=0;i<n;i++)
    buffer[i] = value;
}

static int top_down_step(const int nodes, const int num_frontier, const int snode, const int degree,
			 const int* restrict adjacency, int* restrict frontier, int* restrict next,
			 int* restrict parents, char* restrict bitmap)
{
  int count = 0;
#pragma omp parallel
  {
    int local_count = 0;
    int *local_frontier = malloc(nodes * sizeof(int));  // (num_frontier*degree)/threads * sizeof(int)
#pragma omp for nowait
     for(int i=0;i<num_frontier;i++){
       int v = frontier[i];
       for(int j=0;j<degree;j++){
         int n = *(adjacency + v * degree + j);  // adjacency[v][j];
         if(bitmap[n] == 1 || snode == n) continue;
         if(__sync_bool_compare_and_swap(&parents[n], NOT_VISITED, parents[v]+1)){
           bitmap[n] = 1;
           local_frontier[local_count++] = n;
         }
       }
     }  // end for i
#pragma omp critical
     {
       memcpy(&next[count], local_frontier, local_count*sizeof(int));
       count += local_count;
     }
     free(local_frontier);
  }
  return count;
}


static bool evaluation(const int nodes, int based_nodes, const int groups, const int lines, const int degree,
		       int adjacency[nodes][degree], int *diameter, double *ASPL, const int added_centers, double *sum,
		       bool verbose)
{
  bool reeached = true;
  int max_diam = 0;
  *sum = 0.0;

  int first_task  = based_nodes;
  int second_task = added_centers - 1;
  int whole_task  = (added_centers)? first_task+second_task : first_task;
  int per_task    = whole_task;
  int start_task  = 0;
  int *frontier = malloc(sizeof(int) * nodes);
  int *next     = malloc(sizeof(int) * nodes);
  int *parents  = malloc(sizeof(int) * nodes);
  char *bitmap  = malloc(sizeof(char) * nodes);

    //// First Search ////
  int tmp_per_task = per_task;
  tmp_per_task = MIN(tmp_per_task, first_task-start_task);
  if(tmp_per_task < 0)
    tmp_per_task = 0;

  for(int snode=0;snode<tmp_per_task;snode++){
    if(verbose) printf("%d/%d\n", snode, start_task+tmp_per_task);
    
    // Initialize
    frontier[0] = snode;
    int num_frontier = 1;
    int tmp_diam = 0;
    clear_buffer(nodes, parents, NOT_VISITED);
    memset(bitmap, 0, nodes);

    do{
      num_frontier = top_down_step(nodes, num_frontier, snode, degree, (int *)adjacency,
                                   frontier, next, parents, bitmap);

      // Swap frontier <-> next
      int *tmp = frontier;
      frontier = next;
      free(tmp);
      next = malloc(sizeof(int) * nodes);
      if(num_frontier != 0) tmp_diam++;
    } while(num_frontier);

    if(max_diam < tmp_diam)
      max_diam = tmp_diam;

    for(int i=snode+1;i<groups*based_nodes;i++){
      if(parents[i] == NOT_VISITED)  // Never visit a node
        reeached = false;

      *sum += (parents[i] + 1) * (groups - i/based_nodes);
    }

    for(int i=groups*based_nodes;i<nodes;i++){
      if(parents[i] == NOT_VISITED)  // Never visit a node
        reeached = false;

      *sum += (parents[i] + 1) * groups;
    }
  }

    //// Second Search (only if add_centers exists) ////
  if(start_task > first_task)
    start_task = based_nodes * groups + (start_task-first_task);
  else
    start_task = based_nodes * groups;

  tmp_per_task = per_task - tmp_per_task;
  for(int snode=start_task;snode<start_task+tmp_per_task;snode++){
    // Initialize
    frontier[0] = snode;
    int num_frontier = 1;
    int tmp_diam = 0;
    clear_buffer(nodes, parents, NOT_VISITED);
    memset(bitmap, 0, nodes);

    do{
      num_frontier = top_down_step(nodes, num_frontier, snode, degree, (int *)adjacency,
                                   frontier, next, parents, bitmap);

      // Swap frontier <-> next
      int *tmp = frontier;
      frontier = next;
      free(tmp);
      next = malloc(sizeof(int) * nodes);
      if(num_frontier != 0) tmp_diam++;
    } while(num_frontier);

    if(max_diam < tmp_diam)
      max_diam = tmp_diam;

    for(int i=snode+1;i<nodes;i++){
      if(parents[i] == NOT_VISITED)  // Never visit a node
        reeached = false;

      *sum += parents[i] + 1;
    }
  }

  free(frontier);
  free(next);
  free(parents);
  free(bitmap);

  if(! reeached){
    return false;
  }

  *diameter = max_diam;
  *ASPL     = *sum / ((((double)nodes-1)*nodes)/2);

  return true;
}

static void print_help(char *argv)
{
  printf("%s -f <edge_file>\n", argv);
}

static void set_args(const int argc, char **argv, char *infname, int *groups, bool *verbose)
{
  if(argc < 3)
    print_help(argv[0]);
  
  int result;
  while((result = getopt(argc,argv,"vf:g:"))!=-1){
    switch(result){
    case 'v':
      *verbose = true;
      break;
    case 'g':
      *groups = atoi(optarg);
      if(*groups < 1){
        fprintf(stderr, "-g value >= 1\n");
	exit(0);
      }
      break;
    case 'f':
      if(strlen(optarg) > MAX_FILENAME_LENGTH){
        fprintf(stderr, "Input filename is long (%s). Please change MAX_FILENAME_LENGTH.\n", optarg);
	exit(1);
      }
      strcpy(infname, optarg);
      break;
    default:
      print_help(argv[0]);
    }
  }
}

static int count_lines(const char *fname)
{
  FILE *fp = NULL;
  if((fp = fopen(fname, "r")) == NULL){
    fprintf(stderr, "File not found\n");
    exit(1);
  }
  
  int lines = 0, c;
  while((c = fgetc(fp)) != EOF)
    if(c == '\n')
      lines++;

  fclose(fp);
  return lines;
}

static void read_file(int (*edge)[2], const char *fname)
{
  FILE *fp;
  if((fp = fopen(fname, "r")) == NULL){
    fprintf(stderr, "File not found\n");
    exit(1);
  }

  int n1, n2, i = 0;
  while(fscanf(fp, "%d %d", &n1, &n2) != EOF){
    edge[i][0] = n1;
    edge[i][1] = n2;
    i++;
  }

  fclose(fp);
}

static int max_node_num(const int lines, const int edge[lines*2])
{
  int max = edge[0];
  for(int i=1;i<lines*2;i++)
    max = MAX(max, edge[i]);

  return max;
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

int main(int argc, char *argv[])
{
  char infname[MAX_FILENAME_LENGTH];
  int groups = 1, diameter = -1, added_centers = 0, low_diam;
  double ASPL = -1, sum = 0.0, low_ASPL;
  bool verbose = false;
  
  set_args(argc, argv, infname, &groups, &verbose);
  int lines      = count_lines(infname);
  int (*edge)[2] = malloc(sizeof(int)*lines*2); // int edge[lines][2];
  read_file(edge, infname);
  
  int nodes       = max_node_num(lines, (int *)edge) + 1;
  int degree      = (lines*2) / nodes;
  int based_nodes = nodes/groups;
  int (*adjacency)[degree] = malloc(sizeof(int)*nodes*degree); // int adjacency[nodes][degree];

  lower_bound_of_diam_aspl(&low_diam, &low_ASPL, nodes, degree);
  printf("Nodes Degrees Groups: %d %d %d\n", nodes, degree, groups);
  verify(lines, degree, nodes, edge);
  
  timer_clear_all();
  timer_start(TIMER_ADJ);
  create_adjacency(nodes, lines, degree, edge, adjacency);
  timer_stop(TIMER_ADJ);
  
  timer_start(TIMER_BFS);
  evaluation(nodes, based_nodes, groups, lines, degree, adjacency, &diameter, &ASPL, added_centers, &sum, verbose);
  timer_stop(TIMER_BFS);

  printf("Diameter = %d (Lowest Diameter = %d)\n", diameter, low_diam);
  printf("ASPL     = %.10f (%.0f/%.0f)\n", ASPL, sum, (double)nodes*(nodes-1)/2);
  printf("ASPL Gap = %f (%.10f - %.10f)\n", ASPL - low_ASPL, ASPL, low_ASPL);
  printf("TIME:\n");
  printf("  ADJ: %f sec.\n", timer_read(TIMER_ADJ));
  printf("  BFS: %f sec.\n", timer_read(TIMER_BFS));
  return 0;
}
