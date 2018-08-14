#include "common.h"

static void print_help(char *argv, const int rank)
{
  END("%s -f <edge_file> [-o <output_file>] [-s <random_seed>] [-t <num_threads>] [-g <gruops>] \
[-n <num_calculations>] [-w <max_temperature>] [-c <min_temperature>] [-C <cooling_cycle>] [-d] [-a <accept_rate>] \
[-O <optimization>] [-p <add lines for center>] [-H] [-y] [-h]\n", argv);
}

static void set_args(const int argc, char **argv, const int rank, char *infname, char *outfname,
                     bool *outfnameflag, int *random_seed, int *thread_num, long long *ncalcs,
                     double *max_temp, bool *max_temp_flag, double *min_temp, bool *min_temp_flag,
                     bool *auto_temp_flag, double *accept_rate, int *cooling_cycle, 
		     bool *hill_climbing_flag, bool *detect_temp_flag, int *groups, int *opt,
		     bool *center_flag, int *add_degree_to_center, bool *halfway_flag)
{
  if(argc < 3)
    print_help(argv[0], rank);

  int result;
  while((result = getopt(argc,argv,"f:o:s:t:n:w:c:C:g:a:O:p:dHyh"))!=-1){
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
    case 'C':
      *cooling_cycle = atoi(optarg);
      if(*cooling_cycle <= 0)
	ERROR("Cooling Cycle > 0\n");
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
    case 'p':
      *center_flag = true;
      *add_degree_to_center = atoi(optarg);
      if(*add_degree_to_center <= 0)
	ERROR("-p > 0\n");
      break;
    case 'H':
      *halfway_flag = true;
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

static void create_symmetric_edge_with_center(int (*edge)[2], const int based_nodes, const int based_lines,
					      const int groups, const int degree, const int rank, const int size,
					      const int original_based_lines, const int add_degree_to_center, int const nodes)
{
  if(add_degree_to_center > based_nodes)
    ERROR("Number of degree or nodes is invalid\n");
  
  if(((based_nodes-add_degree_to_center)*groups)%2 == 1)
    ERROR("Number of degree or nodes is invalid under handshaking lemma\n");
  
  int center_vertex  = based_nodes * groups;
  int connect_vertex = 0; // getRandom(based_nodes);
  if(groups%2 == 1){
    for(int i=0;i<based_lines-original_based_lines;i++){
      int line = original_based_lines + i;
      edge[line][0] = i;
      edge[line][1] = (connect_vertex != i)? (nodes-i-1) : center_vertex;
    }
    //      printf("degree         = %d\n", degree);
    //      printf("nodes          = %d\n", nodes);
    //      printf("lines          = %d\n", lines);
    //      printf("based_lines    = %d\n", based_lines);
    //      printf("center_vertex  = %d\n", center_vertex);
    //      printf("connect_vertex = %d\n", connect_vertex);
    
    for(int j=1;j<groups;j++){
      for(int i=0;i<based_lines;i++){
	edge[based_lines*j+i][0] = edge[i][0] + based_nodes * j;
	if(edge[i][1] == center_vertex)
	  edge[based_lines*j+i][1] = center_vertex;
	else{
	  int next = edge[i][1] + based_nodes * j;
	  edge[based_lines*j+i][1] = (next < nodes)? next : next-nodes+1;
	}
      }
    }
  }
}

static void create_symmetric_edge(int (*edge)[2], const int based_nodes, const int based_lines,
                                  const int groups, const int degree, const int rank, const int size,
				  const int opt, const int center_flag)
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
    int start_line = getRandom(based_lines*groups);
    edge_1g_opt(edge, nodes, based_nodes, based_lines, groups, start_line, center_flag);
    create_adjacency(nodes, lines, degree, edge, adjacency);
    if(evaluation(nodes, based_nodes, groups, lines, degree, adjacency, &diam, &ASPL, total_distance, rank, size, opt, center_flag)) break;
  }
  free(adjacency);
}

static void verfy_graph(const int rank, const int nodes, const int based_nodes, const int degree,
			const int groups, const int lines, int edge[lines][2], bool center_flag,
			const int add_degree_to_center)

{
  PRINT_R0("Verifing a regular graph... ");

  if(nodes < degree)
    ERROR("NG. n is too small. nodes = %d degree = %d\n", nodes, degree);

  if((2*lines)%degree != 0)
    ERROR("NG. lines or n or d is invalid. lines = %d nodes = %d degree = %d\n", lines, nodes, degree);

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

  if(!check_duplicate_edge(lines, edge))
    ERROR("NG\nThe same node conbination in the edge.\n");

  if(!check(rank, nodes, based_nodes, lines, degree, groups, edge, center_flag, add_degree_to_center, 0))
    ERROR("NG\nNot symmetric graph.\n");
  
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

static void output_params(const int size, const int nodes, const int degree, const int groups, const int opt,
                          const int random_seed, const int thread_num, const double max_temp, const double min_temp,
                          const double accept_rate, const long long ncalcs, const int cooling_cycle,
			  const double cooling_rate, const char *infname, const char *outfname, const bool outfnameflag,
                          const double average_time, const bool hill_climbing_flag, const bool auto_temp_flag, const bool center_flag)
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
      printf("   Cooling Cycle: %d\n", cooling_cycle);
      printf("   Cooling Rate : %f\n", cooling_rate);
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
  if(center_flag) printf("   Center Point: Exist\n");
  else            printf("   Center Point: None\n");
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

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Get_processor_name(processor_name, &namelen);
  PRINT_R0("Run on %s\n", processor_name);
  time_t t = time(NULL);
  PRINT_R0("%s---\n", ctime(&t));
  
  int diam = 0, low_diam = 0;
  double ASPL = 0, low_ASPL = 0, cooling_rate = 0;
  FILE *fp = NULL;

  // Initial parameters
  long long ncalcs = 10000, num_accepts = 0;
  int random_seed = 0, thread_num = 1, groups = 1, opt = 0, cooling_cycle = 1;
  int add_degree_to_center = -1;
  double max_temp = 80.0, min_temp = 0.2, accept_rate = 1.0, max_diff_energy = 0;
  bool max_temp_flag = false, min_temp_flag = false, outfnameflag = false, center_flag = false;
  bool hill_climbing_flag = false, auto_temp_flag = false, detect_temp_flag = false;
  bool halfway_flag = false;
  char *infname  = malloc(MAX_FILENAME_LENGTH);
  char *outfname = malloc(MAX_FILENAME_LENGTH);

  // Set arguments
  set_args(argc, argv, rank, infname, outfname, &outfnameflag, 
	   &random_seed, &thread_num, &ncalcs, &max_temp, &max_temp_flag, 
	   &min_temp, &min_temp_flag, &auto_temp_flag, &accept_rate, &cooling_cycle, 
	   &hill_climbing_flag, &detect_temp_flag, &groups, &opt,
	   &center_flag, &add_degree_to_center, &halfway_flag);

  if((max_temp_flag && auto_temp_flag && hill_climbing_flag) || 
     (min_temp_flag && auto_temp_flag && hill_climbing_flag))
    ERROR("Two of (-w or -c), -a, and -y cannot be used.\n");
  
  if(hill_climbing_flag && detect_temp_flag)
    ERROR("Both -h and -d cannot be used.\n");

  srandom(random_seed);
  omp_set_num_threads(thread_num);

  int based_lines = count_lines(rank, infname);
  int lines       = (halfway_flag)? based_lines : based_lines * groups;
  int (*edge)[2]  = malloc(sizeof(int)*lines*2); // int edge[lines][2];

  read_file(edge, rank, infname);
  int based_nodes = max_node_num(based_lines, (int *)edge) + 1;
  if(halfway_flag){
    if(based_nodes%groups != 0) ERROR("based_nodes must be divisible by groups\n");
    based_nodes /= groups;
  }

  if(based_nodes < size)
    ERROR("Number of processes is too big. (Vertexs (%d) < Processes (%d))\n", based_nodes, size);
  
  int nodes  = based_nodes * groups;
  int degree = 2 * lines / nodes;

  if(center_flag){
    if(add_degree_to_center != 1 || groups%2 == 0 || halfway_flag)
      ERROR("Not implemented yet1\n");

    degree += 1;
    nodes  += 1;
    lines  += ((based_nodes-add_degree_to_center)/2+add_degree_to_center)*groups;
    int original_based_lines = based_lines;
    based_lines = lines / groups;
    int (*tmp_edge)[2] = malloc(sizeof(int)*lines*2);
    memcpy(tmp_edge, edge, sizeof(int)*original_based_lines*2);
    free(edge);
    edge = tmp_edge;
    create_symmetric_edge_with_center(edge, based_nodes, based_lines, groups, degree, rank, size,
				      original_based_lines,add_degree_to_center, nodes);
    //    for(int i=0;i<lines;i++)
    //      printf("%d %d\n", edge[i][0], edge[i][1]);
  }
  else{
    if(!halfway_flag)
      create_symmetric_edge(edge, based_nodes, based_lines, groups, degree, rank, size, opt, center_flag);
  }
  
  verfy_graph(rank, nodes, based_nodes, degree, groups, lines, edge, center_flag, add_degree_to_center);
  lower_bound_of_diam_aspl(&low_diam, &low_ASPL, nodes, degree);
  check_current_edge(nodes, degree, lines, groups, based_nodes, edge, low_ASPL, rank, size, center_flag);
  double average_time = estimated_elapse_time(ncalcs, nodes, based_nodes, lines, degree, groups,
					      edge, rank, size, opt, center_flag, add_degree_to_center);

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
    cooling_rate = (max_temp != min_temp)? pow(min_temp/max_temp, 1.0/((double)ncalcs/cooling_cycle)) : 1.0;
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
		  min_temp, accept_rate, ncalcs, cooling_cycle, cooling_rate, infname, outfname, 
		  outfnameflag, average_time, hill_climbing_flag, auto_temp_flag, center_flag);
  
  // Optimization
  timer_clear_all();
  timer_start(TIMER_SA);
  long long step = sa(nodes, lines, degree, groups, max_temp, ncalcs, cooling_rate, low_diam, low_ASPL,
		      hill_climbing_flag, detect_temp_flag, &max_diff_energy, edge, &diam, &ASPL, rank, 
		      size, opt, cooling_cycle, center_flag, add_degree_to_center, based_nodes, &num_accepts);
  timer_stop(TIMER_SA);
  
  if(detect_temp_flag){
    // Set max temperature to accept it 50% in maximum diff energy.
    PRINT_R0("Proposed max temperature is %f\n", (-1.0 * max_diff_energy) / log(0.5));
    // Set min temperature to accept it  5% in minimum diff energy.
    END("Proposed min temperature is %f\n", (-2.0) / log(0.05));
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

  verfy_graph(rank, nodes, based_nodes, degree, groups, lines, edge, center_flag, add_degree_to_center);

  MPI_Finalize();
  return 0;
}
