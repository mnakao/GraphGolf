#include "common.h"

static void print_help(char *argv, const int rank)
{
  END("%s -f <edge_file> [-o <output_file>] [-s <random_seed>] [-t <num_threads>] [-g <gruops>] \
       [-n <num_calculations>] [-w <max_temperature>] [-c <min_temperature>] [-d] [-a <accept_rate>] \
       [-O <optimization>] [-y] [-h]\n", argv);
}

static void set_args(const int argc, char **argv, const int rank, char *infname, char *outfname,
                     bool *outfnameflag, int *random_seed, int *thread_num, long long *ncalcs,
                     double *max_temp, bool *max_temp_flag, double *min_temp, bool *min_temp_flag,
                     bool *auto_temp_flag, double *accept_rate, bool *hill_climbing_flag,
                     bool *detect_temp_flag, int *groups, int *opt)
{
  if(argc < 3)
    print_help(argv[0], rank);

  int result;
  while((result = getopt(argc,argv,"f:o:s:t:n:w:c:g:a:O:dyph"))!=-1){
    switch(result){
    case 'f':
      if(strlen(optarg) > MAX_FILENAME_LENGTH)
        ERROR("Input filename is long (%s). Please change MAX_FILENAME_LENGTH.\n", optarg);
      strcpy(infname, optarg);
      break;
    case 'o':
      if(strlen(optarg) > MAX_FILENAME_LENGTH)
        ERROR("Output filename is long (%s). Please change MAX_FILENAME_LENGTH.\n", optarg);
      strcpy(outfname, optarg);
      *outfnameflag = true;
      break;
    case 's':
      *random_seed = atoi(optarg);
      if(*random_seed < 0)
        ERROR("-s value >= 0\n");
      break;
    case 't':
      *thread_num = atoi(optarg);
      if(*thread_num < 1)
        ERROR("-s value >= 1\n");
      break;
    case 'n':
      *ncalcs = atoll(optarg);
      if(*ncalcs <= 0)
        ERROR("-n value > 0\n");
      break;
    case 'w':
      *max_temp = atof(optarg);
      if(*max_temp <= 0)
        ERROR("-w value > 0\n");
      *max_temp_flag = true;
      break;
    case 'c':
      *min_temp = atof(optarg);
      if(*min_temp <= 0)
        ERROR("MIN value > 0\n");
      *min_temp_flag = true;
      break;
    case 'g':
      *groups = atoi(optarg);
      if(*groups < 1)
        ERROR("-g value >= 1\n");
      break;
    case 'a':
      *accept_rate = atof(optarg);
      if(*accept_rate <= 0 || *accept_rate >= 1.0)
        ERROR("0 < -a value < 1.0\n");
      *auto_temp_flag = true;
      break;
    case 'O':
      *opt = atoi(optarg);
      if(*opt != 0 && *opt != 1)
	ERROR("-O=0 or -O=1\n");
      break;
    case 'd':
      *detect_temp_flag = true;
      break;
    case 'y':
      *hill_climbing_flag = true;
      break;
    case 'h':
    default:
      print_help(argv[0], rank);
    }
  }
}

static int count_lines(const int rank, const char *fname)
{
  FILE *fp = NULL;
  if((fp = fopen(fname, "r")) == NULL)
    ERROR("File not found\n");
  
  int lines = 0, c;
  while((c = fgetc(fp)) != EOF)
    if(c == '\n')
      lines++;

  fclose(fp);
  return lines;
}

static void read_file(int (*edge)[2], const int rank, const char *fname)
{
  FILE *fp;
  if((fp = fopen(fname, "r")) == NULL)
    ERROR("File not found\n");

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

static void create_symmetric_edge(int (*edge)[2], const int based_nodes, const int based_lines,
                                  const int groups, const int degree, const int rank, const int size)
{
  for(int j=1;j<groups;j++)
    for(int i=0;i<based_lines;i++){
      edge[based_lines*j+i][0] = edge[i][0] + based_nodes * j;
      edge[based_lines*j+i][1] = edge[i][1] + based_nodes * j;
    }

  int diam;    // Not use
  double ASPL; // Not use
  int nodes = based_nodes * groups;
  int lines = based_lines * groups;
  int (*adjacency)[degree] = malloc(sizeof(int)*nodes*degree); // int adjacency[nodes][degree];
  int total_distance[based_nodes];

  while(1){
    int start_line = getRandom(based_lines);
    edge_1g_opt(edge, based_nodes, based_lines, groups, start_line);
    create_adjacency(nodes, lines, degree, edge, adjacency);
    if(evaluation(nodes, groups, lines, degree, adjacency, &diam, &ASPL, total_distance, rank, size)) break;
  }
  free(adjacency);
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

static void verfy_regular_graph(const int rank, const int n, const int d, const int lines, int edge[lines][2])
{
  PRINT_R0("Verifing a regular graph... ");

  if(n < d)
    ERROR("NG. n is too small. n = %d d = %d\n", n, d);

  if((2*lines)%d != 0)
    ERROR("NG. lines or n or d is invalid. lines = %d n = %d d = %d\n", lines, n, d);

  int nodes[n];
  for(int i=0;i<n;i++)
    nodes[i] = 0;

  for(int i=0;i<lines;i++){
    nodes[edge[i][0]]++;
    nodes[edge[i][1]]++;
  }

  for(int i=0;i<n;i++)
    if(d != nodes[i])
      ERROR("NG\nNot regular graph. d = %d nodes[%d] = %d\n", d, i, nodes[i]);

  if(!check_loop(lines, edge))
    ERROR("NG\nThe same node in the edge.\n");

  if(!check_duplicate_edge(lines, edge))
    ERROR("NG\nThe same node conbination in the edge.\n");
  
  PRINT_R0("OK\n");
}

static void output_params(const int size, const int nodes, const int degree, const int groups, const int opt,
                          const int random_seed, const int thread_num, const double max_temp, const double min_temp,
                          const double accept_rate, const long long ncalcs, const double cooling_rate,
                          const char *infname, const char *outfname, const bool outfnameflag,
                          const double average_time, const bool hill_climbing_flag, const bool auto_temp_flag)
{
  printf("---\n");
  printf("Seed: %d\n", random_seed);
  printf("Optimization: %d\n", opt);
  printf("Num. of processes: %d\n", size);
  printf("Num. of threads  : %d\n", thread_num);
  if(hill_climbing_flag == false){
    printf("Algorithm: Simulated Annealing\n");
    if(!auto_temp_flag){
      printf("   MAX Temperature: %f\n", max_temp);
      printf("   MIN Temperature: %f\n", min_temp);
      printf("   Cooling Rate: %f\n", cooling_rate);
    }
    else{
      printf("   Accept Rate: %f\n", accept_rate);
    }
  }
  else{
    printf("Algorithm: Hill climbing Method\n");
  }

  printf("Num. of Calulations: %lld\n", ncalcs);
  printf("   Average BFS time: %f sec.\n", average_time);
  printf("   Estimated elapse time: %f sec.\n", average_time * ncalcs);
  printf("Input filename: %s\n", infname);
  printf("   Num. of nodes: %d\n", nodes);
  printf("   Degree: %d\n", degree);
  printf("   Groups: %d\n", groups);
  if(outfnameflag)
    printf("Output filename: %s\n", outfname);
  printf("---\n");
}

static void output_file(FILE *fp, const int lines, int edge[lines][2])
{
  for(int i=0;i<lines;i++)
    fprintf(fp, "%d %d\n", edge[i][0], edge[i][1]);
}

int main(int argc, char *argv[])
{
  int rank, size, namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  time_t t = time(NULL);

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Get_processor_name(processor_name, &namelen);
  PRINT_R0("Run on %s\n", processor_name);
  PRINT_R0("%s---\n", ctime(&t));
  
  int diam = 0, low_diam = 0;
  double ASPL = 0, low_ASPL = 0, cooling_rate = 0;
  FILE *fp = NULL;

  // Initial parameters
  long long ncalcs = 10000;
  int random_seed = 0, thread_num = 1, groups = 1, opt = 0;
  double max_temp = 80.0, min_temp = 0.2, accept_rate = 1.0, max_diff_energy = 0;
  bool max_temp_flag = false, min_temp_flag = false, outfnameflag = false;
  bool hill_climbing_flag = false, auto_temp_flag = false, detect_temp_flag = false;
  char *infname  = malloc(MAX_FILENAME_LENGTH);
  char *outfname = malloc(MAX_FILENAME_LENGTH);

  // Set arguments
  set_args(argc, argv, rank, infname, outfname, &outfnameflag, 
	   &random_seed, &thread_num, &ncalcs, &max_temp, &max_temp_flag, 
	   &min_temp, &min_temp_flag, &auto_temp_flag, &accept_rate, 
	   &hill_climbing_flag, &detect_temp_flag, &groups, &opt);

  if((max_temp_flag && auto_temp_flag && hill_climbing_flag) || 
     (min_temp_flag && auto_temp_flag && hill_climbing_flag))
    ERROR("Two of (-w or -c), -a, and -y cannot be used.\n");
  
  if(hill_climbing_flag && detect_temp_flag)
    ERROR("Both -h and -d cannot be used.\n");

  srandom(random_seed);
  omp_set_num_threads(thread_num);

  int based_lines = count_lines(rank, infname);
  int lines       = based_lines * groups;
  int (*edge)[2]  = malloc(sizeof(int)*lines*2); // int edge[lines][2];

  read_file(edge, rank, infname);
  int based_nodes = max_node_num(based_lines, (int *)edge) + 1;
  if(based_nodes < size)
    ERROR("Number of processes is too big. (Vertexs (%d) < Processes (%d))\n", based_nodes, size);
  int nodes       = based_nodes * groups;
  int degree      = 2 * lines / nodes;

  create_symmetric_edge(edge, based_nodes, based_lines, groups, degree, rank, size);
  verfy_regular_graph(rank, nodes, degree, lines, edge);
  lower_bound_of_diam_aspl(&low_diam, &low_ASPL, nodes, degree);
  check_current_edge(nodes, degree, lines, groups, edge, low_ASPL, rank, size);
  double average_time = estimated_elapse_time(ncalcs, nodes, lines, degree, groups, edge, rank, size);

  if(hill_climbing_flag){
    max_temp = min_temp = 0.0;
    cooling_rate = 1.0;
  }
  else if(auto_temp_flag){
    double neighborhood = 2.0 / (nodes * (nodes-1));
    max_temp = min_temp = -1.0 * neighborhood / log(accept_rate) * nodes * (nodes-1);
    cooling_rate = 1.0;
  }
  else{
    cooling_rate = (max_temp != min_temp)? pow(min_temp/max_temp, 1.0/(double)ncalcs) : 1.0;
  }

  if(outfnameflag){
    struct stat stat_buf;
    if(stat(outfname, &stat_buf) == 0)
      ERROR("Output file %s exsits. \n", outfname);
    
    if((fp = fopen(outfname, "w")) == NULL)
      ERROR("Cannot open %s\n", outfname);
  }

  if(rank == 0)
    output_params(size, nodes, degree, groups, opt, random_seed, thread_num, max_temp, 
		  min_temp, accept_rate, ncalcs, cooling_rate, infname, outfname, 
		  outfnameflag, average_time, hill_climbing_flag, auto_temp_flag);
  
  // Optimization
  timer_clear_all();
  timer_start(TIMER_SA);
  long long step = sa(nodes, lines, degree, groups, max_temp, ncalcs, cooling_rate, low_diam, low_ASPL,
		      hill_climbing_flag, detect_temp_flag, &max_diff_energy, edge, &diam, &ASPL, rank, size, opt);
  timer_stop(TIMER_SA);

  if(detect_temp_flag){
    // Set max temperature to accept it 50% in maximum diff energy.
    PRINT_R0("Proposed max temperature is %f\n", -1.0 * max_diff_energy / log(0.5));
    // Set min temperature to accept it 1% in minimum diff energy.
    END("Proposed min temperature is %f\n", -1.0 * 2.0 * groups / log(0.1));
  }

  // Output results
  if(rank == 0){
    printf("---\n");
    printf("Diam. k = %d  ASPL l = %f  Diam. gap = %d  ASPL gap = %f\n",
	   diam, ASPL, diam-low_diam, ASPL-low_ASPL);
    printf("Steps: %lld  Elapse time: %f sec.\n", step, timer_read(TIMER_SA));
  }
  verfy_regular_graph(rank, nodes, degree, lines, edge);

  if(rank == 0 && outfnameflag){
    output_file(fp, lines, edge);
    fclose(fp);
  }

  MPI_Finalize();
  return 0;
}
