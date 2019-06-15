#include "common.h"

static void print_help(char *argv)
{
  END("%s -f <edge_file> [-o <output_file>] [-s <random_seed>] [-n <num_calculations>] \
[-w <max_temperature>] [-c <min_temperature>] [-C <cooling_cycle>] [-N] [-d] [-y] [-h]\n", argv);
}

static void set_args(const int argc, char **argv, char *infname, char *outfname, bool *outfnameflag,
		     int *random_seed, long long *ncalcs, double *max_temp, bool *max_temp_flag,
		     double *min_temp, bool *min_temp_flag, int *cooling_cycle,
		     bool *hill_climbing_flag, bool *detect_temp_flag, bool *verify_flag)
{
  if(argc < 3)
    print_help(argv[0]);

  int result;
  while((result = getopt(argc,argv,"f:o:s:n:w:c:C:Ndyh"))!=-1){
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
    case 'C':
      *cooling_cycle = atoi(optarg);
      if(*cooling_cycle <= 0)
	ERROR("Cooling Cycle > 0\n");
      break;
    case 'N':
      *verify_flag = false;
      break;
    case 'd':
      *detect_temp_flag = true;
      break;
    case 'y':
      *hill_climbing_flag = true;
      break;
    case 'h':
    default:
      print_help(argv[0]);
    }
  }
}

static int count_lines(const char *fname)
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

static void read_file(int (*edge)[2], const char *fname)
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

static void verfy_graph(const int nodes, const int lines, int edge[lines][2])
{
  int degree = 2 * lines / nodes;
  
  PRINT_R0("Verifing a regular graph... ");

  if((2*lines)%degree != 0)
    ERROR("NG. lines or n or d is invalid. nodes = %d degree = %d\n", nodes, degree);

  int n[nodes];
  for(int i=0;i<nodes;i++)
    n[i] = 0;

  for(int i=0;i<lines;i++){
    n[edge[i][0]]++;
    n[edge[i][1]]++;
  }

  for(int i=0;i<nodes;i++)
    if(degree != n[i])
      ERROR("NG\nNot regular graph. degree = %d n[%d] = %d\n", degree, i, n[i]);

  if(!check_loop(lines, edge))
    ERROR("NG\nThe same node in the edge.\n");

  if(!check_duplicate_all_edge(lines, edge))
    ERROR("NG\nThe same node conbination in the edge.\n");

  PRINT_R0("OK\n");
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

static void output_params(const int nodes, const int degree, const int random_seed, const double max_temp,
			  const double min_temp, const long long ncalcs, const int cooling_cycle,
			  const double cooling_rate, const char *infname, const char *outfname,
			  const bool outfnameflag, const double average_time, const bool hill_climbing_flag)
			  
{
#ifdef NDEBUG
  PRINT_R0("NO DEBUG MODE\n");
#else
  PRINT_R0("DEBUG MODE\n");
#endif
  PRINT_R0("Seed: %d\n", random_seed);
  PRINT_R0("Num. of processes: %d\n", size);
#ifdef _OPENMP
  PRINT_R0("Num. of threads  : %d\n", omp_get_max_threads());
#endif
  if(hill_climbing_flag == false){
    PRINT_R0("Algorithm: Simulated Annealing\n");
    PRINT_R0("   MAX Temperature: %f\n", max_temp);
    PRINT_R0("   MIN Temperature: %f\n", min_temp);
    PRINT_R0("   Cooling Cycle: %d\n", cooling_cycle);
    PRINT_R0("   Cooling Rate : %f\n", cooling_rate);
  }
  else{
    PRINT_R0("Algorithm: Hill climbing Method\n");
  }

  PRINT_R0("Num. of Calulations: %lld\n", ncalcs);
  PRINT_R0("   Average BFS time: %f sec.\n", average_time);
  PRINT_R0("   Estimated elapse time: %f sec.\n", average_time * ncalcs);
  PRINT_R0("Input filename: %s\n", infname);
  PRINT_R0("   Vertexes: %d\n", nodes);
  PRINT_R0("   Degree:   %d\n", degree);
  PRINT_R0("   Degree:   %d\n", degree);
  if(outfnameflag)
    PRINT_R0("Output filename: %s\n", outfname);
  PRINT_R0("---\n");
}

static void output_file(FILE *fp, const int lines, int edge[lines][2])
{
  for(int i=0;i<lines;i++)
    fprintf(fp, "%d %d\n", edge[i][0], edge[i][1]);
}

int main(int argc, char *argv[])
{
  bool max_temp_flag = false, min_temp_flag = false, outfnameflag = false, verify_flag = true;
  bool hill_climbing_flag = false, detect_temp_flag = false;
  char hostname[MPI_MAX_PROCESSOR_NAME], infname[MAX_FILENAME_LENGTH], outfname[MAX_FILENAME_LENGTH];
  int namelen, diam = 0, low_diam = 0, random_seed = 0, cooling_cycle = 1;
  long long ncalcs = 10000, num_accepts = 0;
  double ASPL = 0, low_ASPL = 0, cooling_rate = 0;
  double max_temp = 100.0, min_temp = 0.217147, max_diff_energy = 0;
  FILE *fp = NULL;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Get_processor_name(hostname, &namelen);
  PRINT_R0("Run on %s\n", hostname);
  time_t t = time(NULL);
  PRINT_R0("%s---\n", ctime(&t));

  // Set arguments
  set_args(argc, argv, infname, outfname, &outfnameflag,
	   &random_seed, &ncalcs, &max_temp, &max_temp_flag,
	   &min_temp, &min_temp_flag, &cooling_cycle,
	   &hill_climbing_flag, &detect_temp_flag, &verify_flag);

  if(hill_climbing_flag){
    if(max_temp_flag)
      ERROR("Both -y and -w cannot be used.\n");
    else if(min_temp_flag)
      ERROR("Both -y and -c cannot be used.\n");
    else if(detect_temp_flag)
      ERROR("Both -y and -d cannot be used.\n");
  }
  
  srandom(random_seed);
  int lines      = count_lines(infname);
  int (*edge)[2] = malloc(sizeof(int)*lines*2); // int edge[lines][2];
  read_file(edge, infname);
  int nodes      = max_node_num(lines, (int *)edge) + 1;
  int degree     = 2 * lines / nodes;
  if(nodes <= degree)
    ERROR("n is too small. nodes = %d degree = %d\n", nodes, degree);

  if(verify_flag)
    verfy_graph(nodes, lines, edge);

  lower_bound_of_diam_aspl(&low_diam, &low_ASPL, nodes, degree);
  check_current_edge(nodes, lines, edge, low_ASPL);
  double average_time = estimated_elapse_time(nodes, lines, edge);
  if(hill_climbing_flag){
    max_temp = min_temp = 0.0;
    cooling_rate = 1.0;
  }
  else{
    cooling_rate = (max_temp != min_temp)? pow(min_temp/max_temp, (double)cooling_cycle/ncalcs) : 1.0;
  }

  if(outfnameflag && rank == 0){
    struct stat stat_buf;
    if(stat(outfname, &stat_buf) == 0)
      ERROR("Output file %s exsits. \n", outfname);
    
    if((fp = fopen(outfname, "w")) == NULL)
      ERROR("Cannot open %s\n", outfname);
  }

  output_params(nodes, degree, random_seed, max_temp, min_temp, ncalcs, cooling_cycle,
		cooling_rate, infname, outfname, outfnameflag, average_time, hill_climbing_flag);
  // Optimization
  timer_clear_all();
  timer_start(TIMER_SA);
  long long step = sa(nodes, lines, max_temp, ncalcs, cooling_rate, low_diam, low_ASPL,
		      hill_climbing_flag, detect_temp_flag, &max_diff_energy, edge, &diam,
		      &ASPL, cooling_cycle, &num_accepts);
  timer_stop(TIMER_SA);
  
  if(detect_temp_flag){
    // Set max temperature to accept it   50% in maximum diff energy.
    PRINT_R0("Proposed max temperature is %f\n", (-1.0 * max_diff_energy) / log(0.5));
    // Set min temperature to accept it 0.01% in minimum diff energy.
    END("Proposed min temperature is %f\n", (-2.0) / log(0.0001));
  }

  // Output results
  PRINT_R0("---\n");
  PRINT_R0("Diam. k = %d  ASPL l = %f  Diam. gap = %d  ASPL gap = %f\n",
	   diam, ASPL, diam-low_diam, ASPL-low_ASPL);

  double time_sa    = timer_read(TIMER_SA);
  double time_bfs   = timer_read(TIMER_BFS);
  double time_check = timer_read(TIMER_CHECK);
  PRINT_R0("Steps: %lld  Elapse time: %f sec. (BFS: %f sec. Check: %f sec. Other: %f sec.)\n",
	   step, time_sa, time_bfs, time_check, time_sa-(time_bfs+time_check));
  if(ncalcs > SKIP_ACCEPTS)
    PRINT_R0("Accept rate: %f (= %lld/%lld)\n",
	     (double)num_accepts/(ncalcs-SKIP_ACCEPTS), num_accepts, ncalcs-SKIP_ACCEPTS);
  if(rank == 0 && outfnameflag){
    output_file(fp, lines, edge);
    fclose(fp);
  }

  if(verify_flag)
    verfy_graph(nodes, lines, edge);

  MPI_Finalize();
  return 0;
}
