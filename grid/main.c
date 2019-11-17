#include "common.h"

static void print_help(char *argv)
{
  END("%s [-f edge_file] [-W width] [-H height] [-D degree] [-R length] [-o output_file] [-s random_seed]\
 [-n calculations] [-w max_temperature] [-c min_temperature] [-g groups] [-C cooling_cycle] [-B] [-d]\
 [-F fixed_temperature] [-Y] [-M] [-h]\n", argv);
}

static void set_args(const int argc, char **argv, char *infname, int *low_length, char *outfname, 
		     int *random_seed, long long *ncalcs, double *max_temp, double *min_temp, int *groups,
		     int *cooling_cycle, bool *enable_hill_climbing, bool *enable_detect_temp, bool *enable_bfs,
		     bool *enable_halfway, double *fixed_temp, int *width, int *height, int *degree)
{
  if(argc < 3)
    print_help(argv[0]);

  int result;
  while((result = getopt(argc,argv,"f:W:H:D:R:o:s:n:w:c:g:C:BdF:YMh"))!=-1){
    switch(result){
    case 'f':
      if(strlen(optarg) > MAX_FILENAME_LENGTH)
        ERROR("Input filename is long (%s). Please change MAX_FILENAME_LENGTH.\n", optarg);
      strcpy(infname, optarg);
      break;
    case 'W':
      *width = atoi(optarg);
      if(*width <= 0)
        ERROR("-W value > 0\n");
      break;
    case 'H':
      *height = atoi(optarg);
      if(*height <= 0)
        ERROR("-H value > 0\n");
      break;
    case 'D':
      *degree = atoi(optarg);
      if(*degree <= 0)
        ERROR("-D value > 0\n");
      break;
    case 'R':
      *low_length = atoi(optarg);
      if(*low_length <= 0)
        ERROR("-R value > 0\n");
      break;
    case 'o':
      if(strlen(optarg) > MAX_FILENAME_LENGTH)
        ERROR("Output filename is long (%s). Please change MAX_FILENAME_LENGTH.\n", optarg);
      strcpy(outfname, optarg);
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
      break;
    case 'c':
      *min_temp = atof(optarg);
      if(*min_temp <= 0)
        ERROR("-c value > 0\n");
      break;
    case 'g':
      *groups = atoi(optarg);
      if(*groups != 1 && *groups != 2 && *groups != 4)
        ERROR("-g value == 1 or 2 or 4\n");
      break;
    case 'C':
      *cooling_cycle = atoi(optarg);
      if(*cooling_cycle <= 0)
	ERROR("-C value > 0\n");
      break;
    case 'B':
      *enable_bfs = true;
      break;
    case 'd':
      *enable_detect_temp = true;
      break;
    case 'F':
      *fixed_temp = atof(optarg);
      if(*fixed_temp <= 0)
	ERROR("-F value > 0\n");
      break;
    case 'Y':
      *enable_hill_climbing = true;
      break;
    case 'M':
      *enable_halfway = true;
      break;
    case 'h':
    default:
      print_help(argv[0]);
    }
  }
}

static bool confirm_dist(const int v, const int w, const int height, const int low_length)
{
  int w0 = WIDTH (v, height);
  int h0 = HEIGHT(v, height);
  int w1 = WIDTH (w, height);
  int h1 = HEIGHT(w, height);
  int distance = abs(w0 - w1) + abs(h0 - h1);
  return (distance <= low_length);
}

static void simple_exchange_edge(const int height, const int low_length, const int lines, int* edge)
{
  while(1){
    int e1, e2, new_e1_v, new_e1_w, new_e2_v, new_e2_w;
    do{
      e1 = random() % lines;
      e2 = random() % lines;
    } while( e1 == e2 );
    int e1_v = edge[e1*2]; int e1_w = edge[e1*2+1];
    int e2_v = edge[e2*2]; int e2_w = edge[e2*2+1];
    if(confirm_dist(e1_v, e2_v, height, low_length) && confirm_dist(e1_w, e2_w, height, low_length)){
      new_e1_v = e1_v;  new_e1_w = e2_v;
      new_e2_v = e1_w;  new_e2_w = e2_w;
    }
    else if(confirm_dist(e1_v, e2_w, height, low_length) && confirm_dist(e1_w, e2_v, height, low_length)){
      new_e1_v = e1_v;  new_e1_w = e2_w;
      new_e2_v = e1_w;  new_e2_w = e2_v;
    }
    else{
      continue;
    }
    edge[2*e1] = new_e1_v;  edge[2*e1+1] = new_e1_w;
    edge[2*e2] = new_e2_v;  edge[2*e2+1] = new_e2_w;
    break;
  }
}

#ifdef _OPENMP
static int top_down_step(const int nodes, const int num_frontier, const int degree,
			 const int* restrict adjacency, int* restrict frontier,
			 int* restrict next, char* restrict bitmap)
{
  int count = 0;
  int local_frontier[nodes];
#pragma omp parallel private(local_frontier)
  {
    int local_count = 0;
#pragma omp for nowait
     for(int i=0;i<num_frontier;i++){
       int v = frontier[i];
       for(int j=0;j<degree;j++){
         int n = *(adjacency + v * degree + j);  // adjacency[v][j];
	 if(bitmap[n] == NOT_VISITED){
           bitmap[n] = VISITED;
           local_frontier[local_count++] = n;
         }
       }
     }  // end for i
#pragma omp critical
     {
       memcpy(&next[count], local_frontier, local_count*sizeof(int));
       count += local_count;
     }
  }
  return count;
}
#else
static int top_down_step(const int nodes, const int num_frontier, const int degree,
			 const int* restrict adjacency, int* restrict frontier,
                         int* restrict next, char* restrict bitmap)
{
  int count = 0;
  for(int i=0;i<num_frontier;i++){
    int v = frontier[i];
    for(int j=0;j<degree;j++){
      int n = *(adjacency + v * degree + j);  // int n = adjacency[v][j];
      if(bitmap[n] == NOT_VISITED){
        bitmap[n] = VISITED;
        next[count++] = n;
      }
    }
  }

  return count;
}
#endif

static int simple_bfs(const int nodes, const int degree, int *adjacency)
{
  char *bitmap  = malloc(sizeof(char) * nodes);
  int *frontier = malloc(sizeof(int)  * nodes);
  int *next     = malloc(sizeof(int)  * nodes);
  int num_frontier = 1, root = 0, num = 0;
  
  for(int i=0;i<nodes;i++)
    bitmap[i] = NOT_VISITED;
  
  frontier[0]  = root;
  bitmap[root] = VISITED;

  while(1){
    num_frontier = top_down_step(nodes, num_frontier, degree,
				 adjacency, frontier, next, bitmap);
    if(num_frontier == 0) break;
    
    int *tmp = frontier;
    frontier = next;
    next     = tmp;
  }

  for(int i=0;i<nodes;i++)
    if(bitmap[i] == NOT_VISITED)
      num++;

  free(bitmap);
  free(frontier);
  free(next);

  return num;
}

// Inherited from http://research.nii.ac.jp/graphgolf/c/create-lattice.c
static void create_lattice(const int nodes, const int lines, const int width, const int height,
			   const int degree, const int low_length, int edge[lines*2])
{
  int i = 0;
  for(int x=0;x<width/2;x++){
    for(int y=0;y<height;y++){
      for(int k=0;k<degree;k++){
        edge[i*2]   = y + 2 * x * height;
        edge[i*2+1] = edge[2*i] + height;
        i++;
      }
    }
  }
  
  if(width%2 == 1){
    for(int y=0;y<height/2;y++){
      for(int k=0;k<degree;k++){
        edge[i*2]   = (width - 1) * height + 2 * y;
        edge[i*2+1] = edge[i*2] + 1;
        i++;
      }
    }

    /* add self-loop */
    if(height%2 == 1){
      for(int k=0;k<degree/2;k++){
        edge[i*2] = edge[i*2+1] = nodes - 1;
        i++;
      }
    }
  }
  
  int *tmp_edge = malloc(lines*2*sizeof(int));
  int (*adjacency)[degree] = malloc(sizeof(int)*nodes*degree); // int adjacency[nodes][degree];
  create_adjacency(nodes, lines, degree, (const int (*)[2])edge, adjacency);
  int min_num = simple_bfs(nodes, degree, (int *)adjacency);

  while(1){
    memcpy(tmp_edge, edge, sizeof(int)*lines*2);
    simple_exchange_edge(height, low_length, lines, tmp_edge);
    create_adjacency(nodes, lines, degree, (const int (*)[2])tmp_edge, adjacency);
    int tmp_num = simple_bfs(nodes, degree, (int *)adjacency);
    if(tmp_num == 0){
      memcpy(edge, tmp_edge, sizeof(int)*lines*2);
      break;
    }
    else{
      if(tmp_num <= min_num){
	min_num = tmp_num;
	memcpy(edge, tmp_edge, sizeof(int)*lines*2);
      }
    }
  }
  free(tmp_edge);
  free(adjacency);
  //  for(int i=0;i<lines;i++)
  //    printf("%d,%d %d,%d\n", WIDTH(edge[i*2], height), HEIGHT(edge[i*2], height),
  //	   WIDTH(edge[i*2+1], height), HEIGHT(edge[i*2+1], height));
  //EXIT(0);
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

static void read_file_lattice(int *edge, int *w, int *h, const char *fname)
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
    edge[i*2  ] = n[0] * (*h) + n[1];
    edge[i+2+1] = n[2] * (*h) + n[3];
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

static void create_symmetric_edge(int *edge, const int based_nodes, const int based_lines,
				  const int groups, const int degree, const int nodes, const int lines,
				  const int height, const int width, const int based_height,
				  const int low_length)
{
  for(int i=0;i<based_lines;i++)
    for(int j=0;j<2;j++)
      edge[i*2+j] = WIDTH(edge[i*2+j], based_height) * height + HEIGHT(edge[i*2+j], based_height);

  if(groups == 2){
    for(int i=0;i<based_lines;i++)
      for(int j=0;j<2;j++)
        edge[(based_lines+i)*2+j] = ROTATE(edge[i*2+j], height, width, groups, 180);
  }
  else if(groups == 4){
    for(int i=0;i<based_lines;i++){
      for(int j=0;j<2;j++){
	edge[(based_lines  +i)*2+j] = ROTATE(edge[i*2+j], height, width, groups, 90);
	edge[(based_lines*2+i)*2+j] = ROTATE(edge[i*2+j], height, width, groups, 180);
	edge[(based_lines*3+i)*2+j] = ROTATE(edge[i*2+j], height, width, groups, 270);
      }
    }
  }

  int *tmp_edge = malloc(lines*2*sizeof(int));
  int (*adjacency)[degree] = malloc(sizeof(int)*nodes*degree); // int adjacency[nodes][degree];
  create_adjacency(nodes, lines, degree, (const int (*)[2])edge, adjacency);
  int min_num = simple_bfs(nodes, degree, (int *)adjacency);

  while(1){
    memcpy(tmp_edge, edge, sizeof(int)*lines*2);
    exchange_edge(nodes, lines, degree, (int (*)[2])tmp_edge, height, width, groups, low_length, 0);
    create_adjacency(nodes, lines, degree, (const int (*)[2])tmp_edge, adjacency);
    int tmp_num = simple_bfs(nodes, degree, (int *)adjacency);
    if(tmp_num == 0){
      memcpy(edge, tmp_edge, sizeof(int)*lines*2);
      break;
    }
    else{
      if(tmp_num <= min_num){
	min_num = tmp_num;
	memcpy(edge, tmp_edge, sizeof(int)*lines*2);
      }
    }
  }
  free(tmp_edge);
  free(adjacency);
}

static void verfy_graph(const int nodes, const int lines, const int edge[lines*2],
			const int height, const int low_length)
{
  PRINT_R0("Verifing a regular graph... ");
  
  int n[nodes];
  for(int i=0;i<nodes;i++)
    n[i] = 0;

  for(int i=0;i<lines;i++){
    n[edge[i*2  ]]++;
    n[edge[i*2+1]]++;
  }

  int degree = 2 * lines / nodes;
  for(int i=0;i<nodes;i++)
    if(degree != n[i])
      ERROR("NG\nNot regular graph. degree = %d n[%d] = %d.\n", degree, i, n[i]);

  for(int i=0;i<lines;i++){
    int w0 = WIDTH (edge[i*2  ], height);
    int h0 = HEIGHT(edge[i*2  ], height);
    int w1 = WIDTH (edge[i*2+1], height);
    int h1 = HEIGHT(edge[i*2+1], height);
    int distance = abs(w0 - w1) + abs(h0 - h1);
    if(distance > low_length)
      ERROR("Over length in line %d: %d,%d %d,%d : length = %d, distance = %d\n",
	    i+1, w0, h0, w1, h1, low_length, distance);
  }
  
  PRINT_R0("OK\n");
}

static int dist(const int x1, const int y1, const int x2, const int y2)
{
  return(abs(x1 - x2) + abs(y1 - y2));
}

static void lower_bound_of_diam_aspl(int *low_diam, double *low_ASPL, const int m, const int n,
                                     const int degree, const int length)
{
  int moore[m*n], hist[m*n], mh[m*n];
  int mn = m * n, current = degree, ii;
  double sum = 0;

  moore[0] = 1;
  moore[1] = degree + 1;
  for(ii=2;;ii++){
    current = current * (degree - 1);
    moore[ii] = moore[ii-1] + current;
    if(moore[ii] >= mn){
      moore[ii] = mn;
      break;
    }
  }

  int maxhop = MAX((m+n-2+(length-1))/length, ii);
  for(int i=ii+1;i<=maxhop;i++)
    moore[i] = mn;

  for(int i=0;i<m;i++){
    for(int j=0;j<n;j++){
      for(int k=0;k<=maxhop;k++)
        hist[k] = 0;

      for (int i2=0;i2<m;i2++)
        for(int j2=0;j2<n;j2++)
          hist[(dist(i,j,i2,j2)+length-1)/length]++;

      for(int k=1;k<=maxhop;k++)
        hist[k] += hist[k-1];

      for(int k=0;k<=maxhop;k++)
        mh[k] = MIN(hist[k], moore[k]);

      for(int k=1;k<=maxhop;k++)
        sum += (double)(mh[k] - mh[k-1]) * k;
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
			  const int cooling_cycle, const double cooling_rate, const char *infname,
			  const char *outfname, const double average_time, const bool enable_hill_climbing,
			  const int width, const int height, const bool enable_bfs, const bool enable_fixed_temp,
			  const double fixed_temp)
			  
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

  if(enable_hill_climbing)
    PRINT_R0("Algorithm: Hill climbing Method\n");
  else{
    if(enable_fixed_temp)
      PRINT_R0("Algorithm: Fixed Temperature Simulated Annealing : %f\n", fixed_temp);
    else
      PRINT_R0("Algorithm: Simulated Annealing\n");
    
    PRINT_R0("   MAX Temperature: %f\n", max_temp);
    PRINT_R0("   MIN Temperature: %f\n", min_temp);
    PRINT_R0("   Cooling Cycle: %d\n", cooling_cycle);
    PRINT_R0("   Cooling Rate : %f\n", cooling_rate);
  }

  if(groups != 1)
    PRINT_R0("   Groups       : %d\n", groups);
  PRINT_R0("Num. of Calulations: %lld\n", ncalcs);
  PRINT_R0("   Average APSP time    : %f sec.\n", average_time);
  PRINT_R0("   Estimated elapse time: %f sec.\n", average_time * ncalcs);
  if(infname[0] != NOT_C_DEFINED)
    PRINT_R0("Input filename: %s\n", infname);
  PRINT_R0("   (w x h, d, r) = (%d x %d, %d, %d)\n", width, height, degree, low_length);
  if(outfname[0] != NOT_C_DEFINED)
    PRINT_R0("Output filename: %s\n", outfname);
  PRINT_R0("---\n");
}

static void output_file(FILE *fp, const int lines, const int height, const int edge[lines*2])
{
  for(int i=0;i<lines;i++)
    fprintf(fp, "%d,%d %d,%d\n", WIDTH(edge[i*2], height), HEIGHT(edge[i*2], height),
	    WIDTH(edge[i*2+1], height), HEIGHT(edge[i*2+1], height));
}

int main(int argc, char *argv[])
{
  bool enable_hill_climbing = false, enable_detect_temp = false, enable_bfs = false, enable_halfway = false;
  char hostname[MPI_MAX_PROCESSOR_NAME];
  char infname[MAX_FILENAME_LENGTH] = {NOT_C_DEFINED}, outfname[MAX_FILENAME_LENGTH] = {NOT_C_DEFINED};
  int random_seed = 0, cooling_cycle = 1, groups = 1;
  int namelen, based_lines, lines, based_width, based_height, based_nodes, nodes, *edge;
  int diam  = NOT_N_DEFINED, degree = NOT_N_DEFINED, low_diam   = NOT_N_DEFINED;
  int width = NOT_N_DEFINED, height = NOT_N_DEFINED, low_length = NOT_N_DEFINED;
  long long ncalcs = DEFAULT_NCALCS, num_accepts = 0;
  double ASPL     = NOT_N_DEFINED, low_ASPL = NOT_N_DEFINED, cooling_rate = NOT_N_DEFINED, max_diff_energy = NOT_N_DEFINED;
  double max_temp = NOT_N_DEFINED, min_temp = NOT_N_DEFINED, fixed_temp   = NOT_N_DEFINED;
  FILE *fp = NULL;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);
  MPI_Get_processor_name(hostname, &namelen);
  PRINT_R0("Run on %s\n", hostname);
  time_t t = time(NULL);
  PRINT_R0("%s---\n", ctime(&t));

  // Set arguments
  set_args(argc, argv, infname, &low_length, outfname, &random_seed, &ncalcs, &max_temp,
	   &min_temp, &groups, &cooling_cycle, &enable_hill_climbing, &enable_detect_temp,
	   &enable_bfs, &enable_halfway, &fixed_temp, &width, &height, &degree);

  // Set other arguments
  bool enable_max_temp   = (max_temp    != NOT_N_DEFINED);
  bool enable_min_temp   = (min_temp    != NOT_N_DEFINED);
  bool enable_fixed_temp = (fixed_temp  != NOT_N_DEFINED);
  bool enable_infname    = (infname[0]  != NOT_C_DEFINED);
  bool enable_outfname   = (outfname[0] != NOT_C_DEFINED);
  bool enable_whd        = (width != NOT_N_DEFINED && height != NOT_N_DEFINED && degree != NOT_N_DEFINED);
  
  // Check arguments
  if(low_length == NOT_N_DEFINED)                         ERROR("Must need -R\n");
  else if(enable_hill_climbing && enable_max_temp)        ERROR("Both -Y and -w cannot be used.\n");
  else if(enable_hill_climbing && enable_min_temp)        ERROR("Both -Y and -c cannot be used.\n");
  else if(enable_hill_climbing && enable_detect_temp)     ERROR("Both -Y and -d cannot be used.\n");
  else if(!enable_infname && !enable_whd)                 ERROR("Must set -f or \"-W and -H and -D\"\n");
  else if(enable_halfway && !enable_infname)              ERROR("Must set both -M and -f\n");
  if(!enable_max_temp) max_temp = 100.0;
  if(!enable_min_temp) min_temp = 0.217147;
  if(max_temp == min_temp)                                ERROR("The same values in -w and -c.\n");
  if(enable_detect_temp) ncalcs = DEFAULT_DETECT_NCALS;
  
  srandom(random_seed);
  if(enable_infname){
    based_lines = count_lines(infname);
    lines       = (enable_halfway)? based_lines : based_lines * groups;
    edge        = malloc(sizeof(int)*lines*2); // int edge[lines][2];
    read_file_lattice(edge, &based_width, &based_height, infname);
    based_nodes = max_node_num(based_lines, (int *)edge) + 1;
    
    if(enable_halfway){
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

    nodes  = based_nodes * groups;
    degree = 2 * lines / nodes;
  }
  else{
    nodes       = width * height;
    based_nodes = nodes / groups;
    lines       = nodes * degree / 2;
    based_lines = lines / groups;
    edge        = malloc(sizeof(int)*lines*2); // int edge[lines][2];
    
    if(groups == 1){
      based_width  = width;
      based_height = height;
    }
    else if(groups == 2){
      based_width  = width;
      based_height = height/2;
    }
    else{ // groups == 4
      based_width  = width/2;
      based_height = height/2;
    }
  }
  
  if(groups == 4 && (based_width != based_height))
    ERROR("When g = 4, width(%d) must be equal to height(%d).\n", based_width, based_height);
  else if(groups == 4 && width%2 != 0 && height%2 != 0)
    ERROR("When g = 4, width(%d) and height(%d) are divisible by 2.\n", width, height);
  else if(groups == 2 && height%2 != 0)
    ERROR("When g = 2, height(%d) os divisible by 2.\n", height);
  else if(nodes%groups != 0)
    ERROR("nodes(%d) must be divisible by groups(%d)\n", nodes, groups);
  else if(lines%groups != 0)
    ERROR("(nodes*degree/2) must be divisible by groups(%d)\n", groups);
  else if(based_width*based_height != based_nodes)
    ERROR("Not grid graph (width %d x height %d != nodes %d).\n", based_width, based_height, based_nodes);

  if(!enable_infname)
    create_lattice(based_nodes, based_lines, based_width, based_height, degree, low_length, edge);
  
  int *rotate_hash = malloc(nodes * sizeof(int));
  create_rotate_hash(nodes, height, width, groups, rotate_hash);
  
  if(!enable_halfway && groups != 1)
    create_symmetric_edge(edge, based_nodes, based_lines, groups, degree, nodes,
			  lines, height, width, based_height, low_length);

  verfy_graph(nodes, lines, edge, height, low_length);
  lower_bound_of_diam_aspl(&low_diam, &low_ASPL, width, height, degree, low_length);
  check_current_edge(nodes, lines, edge, low_ASPL, low_diam, groups, height, based_height, enable_bfs, rotate_hash);
  double average_time = estimated_elapse_time(nodes, lines, edge, height, width, based_height, groups,
					      low_length, enable_bfs, rotate_hash);
  if(enable_hill_climbing){
    fixed_temp = max_temp = min_temp = 0.0;
    cooling_rate = 1.0;
  }
  else{
    cooling_rate = pow(min_temp/max_temp, (double)cooling_cycle/ncalcs);
  }

  if(enable_outfname && rank == 0){
    struct stat stat_buf;
    if(stat(outfname, &stat_buf) == 0)
      ERROR("Output file %s exsits. \n", outfname);
    
    if((fp = fopen(outfname, "w")) == NULL)
      ERROR("Cannot open %s\n", outfname);
  }

  output_params(degree, groups, low_length, random_seed, max_temp, min_temp, ncalcs,
		cooling_cycle, cooling_rate, infname, outfname, average_time,
		enable_hill_climbing, width, height, enable_bfs, enable_fixed_temp, fixed_temp);

  // Optimization
  timer_clear_all();
  timer_start(TIMER_SA);
  long long step = sa(nodes, lines, degree, based_nodes, ncalcs, cooling_rate, low_diam, low_ASPL, enable_bfs,
		      enable_hill_climbing, enable_detect_temp, &max_diff_energy, max_temp,
		      min_temp, fixed_temp, edge, &diam, &ASPL, cooling_cycle, &num_accepts, width,
		      based_width, height, based_height, low_length, groups, rotate_hash, enable_fixed_temp);
  timer_stop(TIMER_SA);
  
  if(enable_detect_temp){
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
  if(rank == 0 && enable_outfname){
    output_file(fp, lines, height, edge);
    fclose(fp);
  }

  verfy_graph(nodes, lines, edge, height, low_length);

  MPI_Finalize();
  return 0;
}
