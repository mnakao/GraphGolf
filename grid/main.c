#include "common.h"

static void print_help(char *argv)
{
  END("%s -f <edge_file> [-r length] [-o <output_file>] [-s <random_seed>] [-n <calculations>] [-w <max_temperature>]\
 [-c <min_temperature>] [-g <gruops>] [-C <cooling_cycle>] [-W <weight>] [-B] [-D] [-F] [-H] [-M] [-N] [-R] [-h]\n", argv);
}

static void set_args(const int argc, char **argv, char *infname, int *low_length, char *outfname, bool *enable_outfname,
		     int *random_seed, long long *ncalcs, double *max_temp, bool *enable_max_temp, double *min_temp,
		     bool *enable_min_temp, int *groups, int *cooling_cycle, double *weight, bool *enable_hill_climbing,
		     bool *enable_detect_temp, bool *enable_verify, bool *enable_bfs, bool *enable_halfway,
		     double *fixed_temp, bool *enable_fixed_temp, bool *enable_restriction)
{
  if(argc < 3)
    print_help(argv[0]);

  int result;
  while((result = getopt(argc,argv,"f:o:r:s:n:w:c:g:C:W:BDF:HMNRh"))!=-1){
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
      *enable_outfname = true;
      break;
     case 'r':
      *low_length = atoi(optarg);
      if(*low_length <= 0)
        ERROR("-r value > 0\n");
      break;
    case 's':
      *random_seed = atoi(optarg);
      if(*random_seed < 0)
        ERROR("-s value >= 0\n");
      break;
    case 'n':
      *ncalcs = atoll(optarg);
      if(*ncalcs < 0)
        ERROR("-n value >= 0\n");
      break;
    case 'w':
      *max_temp = atof(optarg);
      if(*max_temp <= 0)
        ERROR("-w value > 0\n");
      *enable_max_temp = true;
      break;
    case 'c':
      *min_temp = atof(optarg);
      if(*min_temp <= 0)
        ERROR("MIN value > 0\n");
      *enable_min_temp = true;
      break;
    case 'g':
      *groups = atoi(optarg);
      if(*groups != 1 && *groups != 2 && *groups != 4)
        ERROR("-g value == 1 or 2 or 4\n");
      break;
    case 'C':
      *cooling_cycle = atoi(optarg);
      if(*cooling_cycle <= 0)
	ERROR("Cooling Cycle > 0\n");
      break;
    case 'W':
      *weight = atof(optarg);
      break;
    case 'B':
      *enable_bfs = true;
      break;
    case 'D':
      *enable_detect_temp = true;
      break;
    case 'F':
      *fixed_temp = atof(optarg);
      if(*fixed_temp <= 0)
	ERROR("-F value > 0\n");
      *enable_fixed_temp = true;
      break;
    case 'H':
      *enable_hill_climbing = true;
      break;
    case 'M':
      *enable_halfway = true;
      break;
    case 'N':
      *enable_verify = false;
      break;
    case 'R':
      *enable_restriction = true;
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

static void read_file_lattice(int (*edge)[2], int *w, int *h, const char *fname)
{
  FILE *fp;
  if((fp = fopen(fname, "r")) == NULL){
    PRINT_R0("File not found\n");
    EXIT(1);
  }

  int n[4];
  *w = 0;
  *h = 0;
  while(fscanf(fp, "%d,%d %d,%d", &n[0], &n[1], &n[2], &n[3]) != EOF){
    *w = MAX(*w, n[0]);
    *h = MAX(*h, n[1]);
    *w = MAX(*w, n[2]);
    *h = MAX(*h, n[3]);
  }
  *w += 1;
  *h += 1;
  rewind(fp);

  int i = 0;
  while(fscanf(fp, "%d,%d %d,%d", &n[0], &n[1], &n[2], &n[3]) != EOF){
    edge[i][0] = n[0] * (*h) + n[1];
    edge[i][1] = n[2] * (*h) + n[3];
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
				  const int groups, const int degree, const int nodes, const int lines,
				  const int height, const int width, const int based_height, const bool enable_bfs)
{
  for(int i=0;i<based_lines;i++)
    for(int j=0;j<2;j++)
      edge[i][j] = WIDTH(edge[i][j], based_height) * height + HEIGHT(edge[i][j], based_height);

  if(groups == 2){
    for(int i=0;i<based_lines;i++)
      for(int j=0;j<2;j++)
        edge[based_lines+i][j] = ROTATE(edge[i][j], height, width, groups, 180);
  }
  else if(groups == 4){
    for(int i=0;i<based_lines;i++){
      for(int j=0;j<2;j++){
	edge[based_lines  +i][j] = ROTATE(edge[i][j], height, width, groups, 90);
	edge[based_lines*2+i][j] = ROTATE(edge[i][j], height, width, groups, 180);
	edge[based_lines*3+i][j] = ROTATE(edge[i][j], height, width, groups, 270);
      }
    }
  }

  // Create adjacency matrix
  int (*adjacency)[degree] = malloc(sizeof(int)*nodes*degree); // int adjacency[nodes][degree];
  int diam;    // NOT_USED 
  double ASPL; // NOT_USED
  int i = 0;
  while(1){
    int start_line = getRandom(lines);
    if(! edge_1g_opt(edge, nodes, lines, degree, based_nodes, based_lines, height, width,
		     groups, start_line, NOT_USED, false, (double)NOT_USED, (double)NOT_USED,
		     (double)NOT_USED, NOT_USED))
      continue;
    if(INITIAL_TIMES == i++)
      break;
  }

  create_adjacency(nodes, lines, degree, (const int (*)[2])edge, adjacency);
  evaluation(nodes, degree, groups, (const int* restrict)adjacency,
	     based_nodes, height, based_height, &diam, &ASPL, enable_bfs);

  assert(check_loop(lines, edge));
  assert(check_duplicate_all_edge(lines, edge));
  assert(check_degree(nodes, lines, edge));
  assert(check_symmetric_edge(lines, edge, height, width, based_height, groups));
  
  free(adjacency);
}

static void verfy_graph(const int nodes, const int lines, int edge[lines][2])
{
  PRINT_R0("Verifing a regular graph... ");
  
  int n[nodes];
  for(int i=0;i<nodes;i++)
    n[i] = 0;

  for(int i=0;i<lines;i++){
    n[edge[i][0]]++;
    n[edge[i][1]]++;
  }

  int degree = 2 * lines / nodes;
  for(int i=0;i<nodes;i++)
    if(degree != n[i])
      ERROR("NG\nNot regular graph. degree = %d n[%d] = %d.\n", degree, i, n[i]);

  if(!check_loop(lines, (const int (*)[2])edge))
    ERROR("NG\nThe same node in the edge.\n");

  if(!check_duplicate_all_edge(lines, (const int (*)[2])edge))
    ERROR("NG\nThe same node conbination in the edge.\n");

  PRINT_R0("OK\n");
}

static int dist(const int x1, const int y1, const int x2, const int y2)
{
  return(abs(x1 - x2) + abs(y1 - y2));
}

// This function is inherited from "http://research.nii.ac.jp/graphgolf/pl/lower-lattice.pl".
static void lower_bound_of_diam_aspl(int *low_diam, double *low_ASPL, const int m, const int n,
				     const int degree, const int length)
{
  int mn = m * n;
  int maxhop = MAX((m+n-2)/length,log(mn/degree)/log(degree-1)-1)+2;
  double sum = 0, current = degree;
  double moore[maxhop+1], hist[maxhop+1], mh[maxhop+1];

  for(int i=0;i<=maxhop;i++)
    moore[i] = hist[i] = 0;

  moore[0] = 1;
  moore[1] = degree + 1;
  for(int i=2;i<=maxhop;i++){
    current = current * (degree - 1);
    moore[i] = moore[i-1] + current;
    if(moore[i] > mn)
      moore[i] = mn;
  }

  for(int i=0;i<m;i++){
    for(int j=0;j<n;j++){
      for(int k=0;k<=maxhop;k++)
        hist[k]=0;

      for (int i2=0;i2<m;i2++)
        for(int j2=0;j2<n;j2++)
          hist[(dist(i,j,i2,j2)+length-1)/length]++;

      for(int k=1;k<=maxhop;k++)
        hist[k] += hist[k-1];

      for(int k=0;k<=maxhop;k++)
        mh[k] = MIN(hist[k], moore[k]);

      for(int k=1;k<=maxhop;k++){
        sum += (double)(mh[k] - mh[k-1]) * k;
      }
    }
  }

  int dboth = 0;
  for(dboth=0;;dboth++)
    if(mh[dboth] == mn)
      break;

  *low_diam = dboth;
  *low_ASPL = sum/((double)mn*(mn-1));
}

static void output_params(const int degree, const int groups, const int low_length, const int random_seed,
			  const double max_temp, const double min_temp, const long long ncalcs,
			  const int cooling_cycle, const double weight, const double cooling_rate, const char *infname,
			  const char *outfname, const bool enable_outfname, const double average_time,
			  const bool enable_hill_climbing, const int width, const int height, const bool enable_bfs,
			  const bool enable_restriction)
			  
{
#ifdef NDEBUG
  PRINT_R0("NO DEBUG MODE\n");
#else
  PRINT_R0("DEBUG MODE\n");
#endif
  PRINT_R0("Seed     : %d\n", random_seed);
  PRINT_R0("Processes: %d\n", procs);
#ifdef _OPENMP
  PRINT_R0("Threads  : %d\n", omp_get_max_threads());
#endif
  if(enable_bfs) PRINT_R0("APSP     : BFS\n");
  else           PRINT_R0("APSP     : MATRIX Opetation\n");

  if(enable_hill_climbing == false){
    if(enable_restriction)
      PRINT_R0("Algorithm: Simulated Annealing (Restricted-2opt)\n");
    else
      PRINT_R0("Algorithm: Simulated Annealing (2-opt)\n");
    PRINT_R0("   MAX Temperature: %f\n", max_temp);
    PRINT_R0("   MIN Temperature: %f\n", min_temp);
    PRINT_R0("   Cooling Cycle: %d\n", cooling_cycle);
    PRINT_R0("   Cooling Rate : %f\n", cooling_rate);
    PRINT_R0("   Weight       : %f\n", weight);
    if(groups != 1)
      PRINT_R0("   Groups       : %d\n", groups);
  }
  else{
    if(enable_restriction)
      PRINT_R0("Algorithm: Hill climbing Method (Restricted-2opt)\n");
    else
      PRINT_R0("Algorithm: Hill climbing Method (2-opt)\n");
  }

  PRINT_R0("Num. of Calulations: %lld\n", ncalcs);
  PRINT_R0("   Average APSP time    : %f sec.\n", average_time);
  PRINT_R0("   Estimated elapse time: %f sec.\n", average_time * ncalcs);
  PRINT_R0("Input filename: %s\n", infname);
  PRINT_R0("   (w x h, d, r) = (%d x %d, %d, %d)\n", width, height, degree, low_length);
  if(enable_outfname)
    PRINT_R0("Output filename: %s\n", outfname);
  PRINT_R0("---\n");
}

static void output_file(FILE *fp, const int lines, const int height, int edge[lines][2])
{
  for(int i=0;i<lines;i++)
    fprintf(fp, "%d,%d %d,%d\n", WIDTH(edge[i][0], height), HEIGHT(edge[i][0], height),
	    WIDTH(edge[i][1], height), HEIGHT(edge[i][1], height));
}

static void check_length(const int lines, const int height, const int length, int edge[lines][2])
{
  for(int i=0;i<lines;i++){
    int w0 = WIDTH (edge[i][0], height);
    int h0 = HEIGHT(edge[i][0], height);
    int w1 = WIDTH (edge[i][1], height);
    int h1 = HEIGHT(edge[i][1], height);
    int distance = abs(w0 - w1) + abs(h0 - h1);
    if(distance > length)
      printf("Over length in line %d: %d,%d %d,%d : length = %d, distance = %d\n",
	     i+1, w0, h0, w1, h1, length, distance);
  }
}

int main(int argc, char *argv[])
{
  bool enable_max_temp = false, enable_min_temp = false, enable_outfname = false, enable_verify = true;
  bool enable_hill_climbing = false, enable_detect_temp = false, enable_bfs = false, enable_halfway = false;
  bool enable_fixed_temp = false, enable_restriction = false;
  char hostname[MPI_MAX_PROCESSOR_NAME], infname[MAX_FILENAME_LENGTH], outfname[MAX_FILENAME_LENGTH];
  int namelen, diam = 0, low_diam = 0, random_seed = 0, cooling_cycle = 1;
  int based_width = 0, based_height = 0, width = 0, height = 0;
  int length = -1, low_length = NOT_DEFINED, groups = 1;
  long long ncalcs = 10000, num_accepts = 0;
  double ASPL = 0, low_ASPL = 0, cooling_rate = 0, weight = 1.0;
  double max_temp = 100.0, min_temp = 0.217147, fixed_temp = 0, max_diff_energy = 0;
  FILE *fp = NULL;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);
  MPI_Get_processor_name(hostname, &namelen);
  PRINT_R0("Run on %s\n", hostname);
  time_t t = time(NULL);
  PRINT_R0("%s---\n", ctime(&t));

  // Set arguments
  set_args(argc, argv, infname, &low_length, outfname, &enable_outfname, &random_seed,
	   &ncalcs, &max_temp, &enable_max_temp, &min_temp, &enable_min_temp, &groups, &cooling_cycle,
	   &weight, &enable_hill_climbing, &enable_detect_temp, &enable_verify, &enable_bfs, &enable_halfway,
	   &fixed_temp, &enable_fixed_temp, &enable_restriction);

  if(low_length == NOT_DEFINED)
    ERROR("Must need -r\n");
  else if(enable_hill_climbing && enable_max_temp)
    ERROR("Both -H and -w cannot be used.\n");
  else if(enable_hill_climbing && enable_min_temp)
    ERROR("Both -H and -c cannot be used.\n");
  else if(enable_hill_climbing && enable_detect_temp)
    ERROR("Both -H and -d cannot be used.\n");
  else if(max_temp == min_temp)
    ERROR("The same values in -w and -c.\n");
  else if(enable_fixed_temp && min_temp > fixed_temp)
    ERROR("The value in -F (%f) must be less than min_temp in -c (%f)\n",
	  fixed_temp, min_temp);
  
  srandom(random_seed);
  int based_lines = count_lines(infname);
  int lines       = (enable_halfway)? based_lines : based_lines * groups;
  int (*edge)[2]  = malloc(sizeof(int)*lines*2); // int edge[lines][2];
  read_file_lattice(edge, &based_width, &based_height, infname);
  if(groups == 4 && (based_width != based_height))
    ERROR("When g = 4, width(%d) must be equal to height(%d).\n", based_width, based_height);
  
  int based_nodes = max_node_num(based_lines, (int *)edge) + 1;
  if(enable_halfway){
    if(based_nodes%groups != 0)
      ERROR("based_nodes must be divisible by groups\n");

    based_nodes /= groups;
    based_lines /= groups;
    
    if(groups == 2){
      based_height /= 2;
    }
    else if(groups == 4){
      based_width  /= 2;
      based_height /= 2;
    }
  }

  if(groups == 1){
    height = based_height;
    width  = based_width;
  }
  else if(groups == 2){
    height = based_height * 2;
    width  = based_width;
  }
  else{ // groups == 4
    height = based_height * 2;
    width  = based_width  * 2;
  }
  
  int nodes  = based_nodes * groups;
  int degree = 2 * lines / nodes;
  if(nodes <= degree)
    ERROR("n is too small. nodes = %d degree = %d\n", nodes, degree);
  else if(based_width*based_height != based_nodes)
    ERROR("Not grid graph (width %d x height %d != nodes %d).\n", based_width, based_height, based_nodes);

  if(!enable_halfway && groups != 1)
    create_symmetric_edge(edge, based_nodes, based_lines, groups, degree, nodes, lines,
			  height, width, based_height, enable_bfs);

  if(enable_verify)
    verfy_graph(nodes, lines, edge);

  lower_bound_of_diam_aspl(&low_diam, &low_ASPL, width, height, degree, low_length);
  check_current_edge(nodes, lines, edge, low_ASPL, groups, height, based_height, enable_bfs);
  double average_time = estimated_elapse_time(nodes, lines, (const int (*)[2])edge,
					      height, width, based_height, groups,
					      low_length, enable_bfs);
  if(enable_hill_climbing){
    fixed_temp = max_temp = min_temp = 0.0;
    cooling_rate = 1.0;
  }
  else if(enable_fixed_temp){
    cooling_rate = 1.0;
  }
  else{
    fixed_temp = max_temp;
    cooling_rate = (max_temp != min_temp)? pow(min_temp/max_temp, (double)cooling_cycle/ncalcs) : 1.0;
  }

  if(enable_outfname && rank == 0){
    struct stat stat_buf;
    if(stat(outfname, &stat_buf) == 0)
      ERROR("Output file %s exsits. \n", outfname);
    
    if((fp = fopen(outfname, "w")) == NULL)
      ERROR("Cannot open %s\n", outfname);
  }

  output_params(degree, groups, low_length, random_seed,
		max_temp, min_temp, ncalcs,
		cooling_cycle, weight, cooling_rate, infname,
		outfname, enable_outfname, average_time,
		enable_hill_climbing, width, height, enable_bfs,
		enable_restriction);

  // Optimization
  timer_clear_all();
  timer_start(TIMER_SA);
  long long step = sa(nodes, lines, fixed_temp, ncalcs,
		      cooling_rate, low_diam, low_ASPL, enable_bfs, 
		      enable_hill_climbing, enable_detect_temp,
		      &max_diff_energy, max_temp, min_temp, edge,
		      &diam, &ASPL, cooling_cycle, &num_accepts, width,
		      based_width, height, based_height, &length,
		      low_length, weight, groups, enable_restriction);
  timer_stop(TIMER_SA);
  
  if(enable_detect_temp){
    // Set max temperature to accept it   50% in maximum diff energy.
    PRINT_R0("Proposed max temperature is %f\n", (-1.0 * max_diff_energy) / log(0.5));
    // Set min temperature to accept it 0.01% in minimum diff energy.
    END("Proposed min temperature is %f\n", (-2.0) / log(0.0001));
  }

  // Output results
  PRINT_R0("---\n");
  PRINT_R0("Diam. k = %d  ASPL l = %f  Diam. gap = %d  ASPL gap = %f Length Gap = %d\n",
	   diam, ASPL, diam-low_diam, ASPL-low_ASPL, length-low_length);

  double time_sa    = timer_read(TIMER_SA);
  double time_apsp  = timer_read(TIMER_APSP);
  double time_check = timer_read(TIMER_CHECK);
  PRINT_R0("Steps: %lld  Elapse time: %f sec. (APSP: %f sec. Check: %f sec. Other: %f sec.)\n",
	   step, time_sa, time_apsp, time_check, time_sa-(time_apsp+time_check));
  if(ncalcs > SKIP_ACCEPTS)
    PRINT_R0("Accept rate: %f (= %lld/%lld)\n",
	     (double)num_accepts/(ncalcs-SKIP_ACCEPTS), num_accepts, ncalcs-SKIP_ACCEPTS);
  if(rank == 0 && enable_outfname){
    output_file(fp, lines, height, edge);
    fclose(fp);
  }

  check_length(lines, height, low_length, edge);
  if(enable_verify)
    verfy_graph(nodes, lines, edge);

  MPI_Finalize();
  return 0;
}
