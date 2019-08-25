#include "common.h"

static void print_help(char *argv)
{
  END("%s -f <edge_file> [-o <output_file>] [-s <random_seed>] \
[-n <num_calculations>] [-w <max_temperature>] [-c <min_temperature>] [-C <cooling_cycle>] [-N] [-d] \
[-g <gruops>] [-v <add vertexes>] [-e <add edges>] [-H] [-y] [-B] [-L] [-h]\n", argv);
}

static void set_args(const int argc, char **argv, char *infname, char *outfname, bool *outfnameflag,
		     int *random_seed, long long *ncalcs, double *max_temp,
		     bool *max_temp_flag, double *min_temp, bool *min_temp_flag, 
		     int *cooling_cycle, bool *hill_climbing_flag, bool *detect_temp_flag, int *groups,
		     int *added_centers, int *added_edges_to_center, bool *halfway_flag, bool *verify_flag,
		     bool *enable_bfs, bool *enable_low_mem)
{
  if(argc < 3)
    print_help(argv[0]);

  int result;
  while((result = getopt(argc,argv,"f:o:s:n:w:c:C:g:v:e:dHNyBLh"))!=-1){
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
        ERROR("-n value >= 1\n");
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
	ERROR("Cooling Cycle >= 1\n");
      break;
    case 'g':
      *groups = atoi(optarg);
      if(*groups < 1)
        ERROR("-g value >= 1\n");
      break;
    case 'd':
      *detect_temp_flag = true;
      break;
    case 'y':
      *hill_climbing_flag = true;
      break;
    case 'v':
      *added_centers = atoi(optarg);
      if(*added_centers < 0)
	ERROR("-v value >= 0\n");
      break;
    case 'e':
      *added_edges_to_center = atoi(optarg);
      if(*added_edges_to_center <= 0)
	ERROR("-e value >= 1\n");
      break;
    case 'H':
      *halfway_flag = true;
      break;
    case 'N':
      *verify_flag = false;
      break;
    case 'B':
      *enable_bfs = true;
      break;
    case 'L':
      *enable_low_mem = true;
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

static void create_symmetric_edge_with_vertexes(int (*edge)[2], const int based_nodes, const int based_lines,
						const int org_based_lines, const int groups, const int degree,
						const int added_centers, const int added_edges_to_center,
						int const nodes)
{
  if((based_nodes-added_centers*added_edges_to_center)%2 != 0)
    ERROR("Sorry. Not implemented. based_nodes (%d) - -v value (%d) * -e value (%d) must be even number\n",
	  based_nodes, added_centers, added_edges_to_center);

  int org_nodes = nodes - added_centers;
  for(int i=0;i<added_centers;i++){
    for(int j=0;j<added_edges_to_center;j++){
      int line = org_based_lines + i * added_edges_to_center + j;
      edge[line][0] = i * added_edges_to_center + j;
      edge[line][1] = org_nodes + i;
    }
  }

  int j = 1;
  for(int i=added_centers*added_edges_to_center;i<based_lines-org_based_lines;i++){
    int line = org_based_lines + i;
    edge[line][0] = i;
    edge[line][1] = 2 * based_nodes - j++;
  }

  for(int j=1;j<groups;j++){
    for(int i=0;i<org_based_lines;i++){
      edge[based_lines*j+i][0] = edge[i][0] + based_nodes * j;
      edge[based_lines*j+i][1] = edge[i][1] + based_nodes * j;
    }
    for(int i=org_based_lines;i<org_based_lines+added_centers*added_edges_to_center;i++){
       edge[based_lines*j+i][0] = edge[i][0] + based_nodes * j;
       edge[based_lines*j+i][1] = edge[i][1];
    }
     for(int i=org_based_lines+added_centers*added_edges_to_center;i<based_lines;i++){
       edge[based_lines*j+i][0] = edge[i][0] + based_nodes * j;
       int tmp = edge[i][1] + based_nodes * j;
       edge[based_lines*j+i][1] = (tmp < org_nodes)? tmp : tmp - org_nodes;
     }
  }
}

static void create_symmetric_edge(int (*edge)[2], const int based_nodes, const int based_lines, const int groups,
				  const int degree, const int nodes, const int lines, const int algo)
{
  for(int j=1;j<groups;j++)
    for(int i=0;i<based_lines;i++){
      edge[based_lines*j+i][0] = edge[i][0] + based_nodes * j;
      edge[based_lines*j+i][1] = edge[i][1] + based_nodes * j;
    }

  int diam;    // Not use
  double ASPL; // Not use
  int (*adj)[degree] = malloc(sizeof(int)*nodes*degree); // int adj[nodes][degree];
  int kind_opt;
  int restored_adj_value[groups*4], restored_adj_idx_y[groups*4], restored_adj_idx_x[groups*4];
  int restored_edge[groups*4], restored_line[groups*2];

  create_adj(nodes, lines, degree, (const int (*)[2])edge, adj);
  while(1){
    int start_line = getRandom(lines);
    edge_1g_opt(edge, nodes, lines, degree, based_nodes, based_lines, groups, start_line, 0, (int *)adj,
		&kind_opt, restored_edge, restored_line, restored_adj_value, restored_adj_idx_y, restored_adj_idx_x, NOT_USED);
    if(evaluation(nodes, based_nodes, groups, lines, degree, adj, &diam, &ASPL, 0, algo))
      break;
    // NOT_RESTORE
  }

  free(adj);
}

static void verfy_graph(const int nodes, const int based_nodes, const int degree, const int groups,
			const int lines, int edge[lines][2], const int added_centers)

{
  PRINT_R0("Verifing a regular graph... ");

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

  if(!check(nodes, based_nodes, lines, degree, groups, edge, added_centers, NULL, 0))
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

static void output_params(const int nodes, const int degree, const int groups, const int random_seed,
			  const double max_temp, const double min_temp, const long long ncalcs, const int cooling_cycle,
			  const double cooling_rate, const char *infname, const char *outfname,
			  const bool outfnameflag, const double average_time, const bool hill_climbing_flag,
			  const int added_centers, const int added_edges_to_center, const int algo)
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
  if(algo == BFS)            PRINT_R0("APSP: BFS\n");
  else if(algo == MATRIX_OP) PRINT_R0("APSP: MATRIX Opetation\n");
  else                       PRINT_R0("APSP: MATRIX Operation (LOW MEMORY)\n");

  PRINT_R0("Num. of Calulations: %lld\n", ncalcs);
  PRINT_R0("   Average APSP time: %f sec.\n", average_time);
  PRINT_R0("   Estimated elapse time: %f sec.\n", average_time * ncalcs);
  PRINT_R0("Input filename: %s\n", infname);
  PRINT_R0("   Vertexes: %d\n", nodes);
  PRINT_R0("   Degree:   %d\n", degree);
  PRINT_R0("   Groups:   %d\n", groups);
  if(added_centers){
    PRINT_R0("   Center points:    %d\n", added_centers);
    PRINT_R0("   Edges to centers: %d\n", added_edges_to_center);
  }
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
  bool hill_climbing_flag = false, detect_temp_flag = false, halfway_flag = false;
  bool enable_bfs = false, enable_low_mem = false;
  char hostname[MPI_MAX_PROCESSOR_NAME], infname[MAX_FILENAME_LENGTH], outfname[MAX_FILENAME_LENGTH];
  int namelen, diam = 0, low_diam = 0, random_seed = 0, algo;
  int groups = 1, cooling_cycle = 1, added_centers = 0, added_edges_to_center = 1;
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
  set_args(argc, argv, infname, outfname, &outfnameflag, &random_seed, &ncalcs,
	   &max_temp, &max_temp_flag, &min_temp, &min_temp_flag, &cooling_cycle,
	   &hill_climbing_flag, &detect_temp_flag, &groups, &added_centers,
	   &added_edges_to_center, &halfway_flag, &verify_flag, &enable_bfs, &enable_low_mem);

  if(hill_climbing_flag && max_temp_flag)
    ERROR("Both -y and -w cannot be used.\n");
  else if(hill_climbing_flag && min_temp_flag)
    ERROR("Both -y and -c cannot be used.\n");
  else if(hill_climbing_flag && detect_temp_flag)
    ERROR("Both -y and -d cannot be used.\n");
  else if(detect_temp_flag && groups != 1)
    ERROR("When using -d option, -g must be 1.\n");
  else if(!added_centers && added_edges_to_center != 1)
    ERROR("-e option must be with -v\n");
  
  srandom(random_seed);
  int based_lines = count_lines(infname);
  int lines       = (halfway_flag)? based_lines : based_lines * groups;
  int (*edge)[2]  = malloc(sizeof(int)*lines*2); // int edge[lines][2];
  read_file(edge, infname);
  int based_nodes = max_node_num(based_lines, (int *)edge) + 1;

  if(halfway_flag){
    if(added_centers) based_nodes -= added_centers;
    if(based_nodes%groups != 0)
      ERROR("based_nodes must be divisible by groups\n");
    else
      based_nodes /= groups;
  }

  int nodes  = based_nodes * groups;
  int degree = 2 * lines / nodes;
  if(nodes <= degree)
    ERROR("n is too small. nodes = %d degree = %d\n", nodes, degree);
    
  if(!halfway_flag && added_centers){
    degree += 1;
    nodes  += added_centers;
    
    if(added_centers > based_nodes)
      ERROR("-v value (%d) is invalid\n", added_centers);
    else if((degree*nodes)%2 != 0)
      ERROR("Number of degrees or number of nodes is invalid %d %d\n", degree, nodes);
    else if(groups*added_edges_to_center != degree)
      ERROR("-g value (%d) * -e value (%d) == degree + 1 (%d)\n", groups, added_edges_to_center, degree);
    else if((nodes-added_centers)/groups%2 != 0)
      ERROR("Number of nodes, groups, or added_centers is invalid 1\n");
    else if(((nodes-added_centers)/groups - degree/groups*added_centers)%2 != 0)
      ERROR("Number of nodes, groups, or added_centers is invalid 2\n");
    else if(((nodes-added_centers)/groups - degree/groups*added_centers)%2 < 0)
      ERROR("Number of nodes, groups, or added_centers is invalid 3\n");
    
    lines = (nodes*degree)/2;
    int org_based_lines = based_lines;
    based_lines = lines / groups;
    void *tmp = realloc(edge, lines*2*sizeof(int));
    if(tmp == NULL)
      ERROR("realloc is failed\n");
    else
      edge = tmp;
    
    create_symmetric_edge_with_vertexes(edge, based_nodes, based_lines, org_based_lines, groups,
					degree, added_centers, added_edges_to_center, nodes);
  }
  else if(halfway_flag && added_centers){
    nodes += added_centers;
  }

  if(enable_bfs) algo = BFS;
  else if(enable_low_mem) algo = MATRIX_OP_LOW_MEM;
  else{
    unsigned int elements = (nodes+(UINT64_BITS-1))/UINT64_BITS;
    double s = (double)nodes/groups * elements * sizeof(uint64_t);
    algo = (s <= (double)MATRIX_OP_THRESHOLD)? MATRIX_OP : MATRIX_OP_LOW_MEM;
  }
  if(added_centers != 0) algo = BFS; // fix me

  if(!halfway_flag && !added_centers)
     create_symmetric_edge(edge, based_nodes, based_lines, groups, degree, nodes, lines, algo);
    
  if(verify_flag)
    verfy_graph(nodes, based_nodes, degree, groups, lines, edge, added_centers);

  lower_bound_of_diam_aspl(&low_diam, &low_ASPL, nodes, degree);
  if(groups == 1)
    check_current_edge(nodes, degree, lines, groups, based_nodes, edge, low_ASPL, added_centers, algo);
  double average_time = estimate_elapse_time(nodes, based_nodes, lines, degree,
					     groups, edge, added_centers, algo);
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

  output_params(nodes, degree, groups, random_seed, max_temp, min_temp, 
		ncalcs, cooling_cycle, cooling_rate, infname, outfname, outfnameflag,
		average_time, hill_climbing_flag, added_centers, added_edges_to_center, algo);
  // Optimization
  timer_clear_all();
  timer_start(TIMER_SA);
  long long step = sa(nodes, lines, degree, groups, max_temp, ncalcs, cooling_rate, low_diam, low_ASPL,
		      hill_climbing_flag, detect_temp_flag, &max_diff_energy, edge, &diam, &ASPL, 
		      cooling_cycle, added_centers, based_nodes, &num_accepts, algo);
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
  double time_apsp  = timer_read(TIMER_APSP);
  double time_check = timer_read(TIMER_CHECK);
  PRINT_R0("Steps: %lld  Elapse time: %f sec. (APSP: %f sec. Check: %f sec. Other: %f sec.)\n",
	   step, time_sa, time_apsp, time_check, time_sa-(time_apsp+time_check));
  if(ncalcs > SKIP_ACCEPTS)
    PRINT_R0("Accept rate: %f (= %lld/%lld)\n",
	     (double)num_accepts/(ncalcs-SKIP_ACCEPTS), num_accepts, ncalcs-SKIP_ACCEPTS);
  if(rank == 0 && outfnameflag){
    output_file(fp, lines, edge);
    fclose(fp);
  }

  if(verify_flag)
    verfy_graph(nodes, based_nodes, degree, groups, lines, edge, added_centers);

  MPI_Finalize();
  return 0;
}
