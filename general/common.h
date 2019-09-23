#ifndef _COMMON_H
#define _COMMON_H

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#include <unistd.h>
#include <sys/stat.h>
#include <stdint.h>
#include <stdbool.h>
#include <assert.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifndef _KCOMPUTER
#include <nmmintrin.h>
#endif

int rank, procs, threads;
#define VISITED     1
#define NOT_VISITED 0

#define RIGHT  0
#define LEFT   1
#define MIDDLE 2
#define D_1G_OPT 1
#define D_2G_OPT 2
#define NOT_USED -1
#define R_DYNAMIC -1
#define ENABLE_CHECK  true
#define DISABLE_CHECK false
#define NO_CHANGE     -1

#define NUM_TIMERS          4
#define TIMER_SA            0
#define TIMER_ESTIMATED     1
#define TIMER_APSP          2
#define TIMER_CHECK         3

#define MAX_FILENAME_LENGTH 255
#define NUM_OF_PROGRESS     100
#define SKIP_ACCEPTS        10000

#define PRINT_ADJ(a)  do{ print_adj(nodes,  degree, a);}while(0);
#define PRINT_EDGE(e) do{ print_edge(nodes, degree, e);}while(0)

#ifdef _KCOMPUTER
#define POPCNT(a) __builtin_popcountll(a)
#else
#define POPCNT(a) _mm_popcnt_u64(a)
#endif

#define BFS                  0
#define MATRIX_OP            1
#define MATRIX_OP_MEM_SAVING 2
#define MATRIX_OP_THRESHOLD  2147483648
#define UINT64_BITS          64
#define CHUNK                64 /* (multiple of sizeof(uint64_t)*8 for AVX-512) */

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define ABORT()         do{MPI_Abort(MPI_COMM_WORLD, 1); exit(1);}while(0)
#define PRINT_R0(...)   do{if(rank==0) printf(__VA_ARGS__);}while(0)
#define END(...)        do{if(rank==0) printf(__VA_ARGS__); MPI_Finalize(); exit(0);}while(0)
#define ERROR(...)      do{if(rank==0) printf(__VA_ARGS__); ABORT();}while(0)
#define EXIT(r)         do{MPI_Finalize(); exit(r);}while(0)

extern void printb(uint64_t v);
extern void print_adj(const int nodes, const int degree, const int adj[nodes][degree]);
extern void print_edge(const int nodes, const int degree, const int edge[nodes*degree/2][2]);
extern void swap(int *a, int *b);
extern int order(int nodes, const int a, const int b, const int added_centers);
extern long long sa(const int nodes, const int lines, const int degree, const int groups,
		    double temp, const long long ncalcs, const double cooling_rate, const int low_diam, const double low_ASPL,
		    const bool hill_climbing_flag, const bool detect_temp_flag, double *max_diff_energy, int edge[lines][2], int *diameter, double *ASPL,
		    const int cooling_cyclie, const int added_centers, const int k_opt, const int based_nodes, long long *num_accepts, const int algo);
extern void check_current_edge(const int nodes, const int degree, const int lines, const int groups,
			       const int based_nodes, int edge[lines][2], const double low_ASPL, const int added_centers, const int algo);
extern double estimate_elapse_time(const int nodes, const int based_nodes, const int lines, const int degree,
				   const int groups, int edge[lines][2], const int add_degree_to_center, const int k_opt, const int algo);
extern bool edge_1g_opt(int (*edge)[2], const int nodes, const int lines, const int degree, const int based_nodes, const int based_lines, const int groups,
			const int start_line, const int add_centers, int* restrict adj,  int *kind_opt, int* restrict restored_adj_edge, int* restrict restored_adj_line,
			int* restrict restored_adj_val, int* restrict restored_adj_idx_y, int* restrict restored_adj_idx_x, const bool enable_check, const int ii);
extern bool has_duplicated_edge(const int e00, const int e01, const int e10, const int e11);
extern bool check_loop(const int lines, int (*edge)[2]);
extern bool check_duplicate_tmp_edge(const int k_opt, const int lines, int (*edge)[2]);
extern bool check_duplicate_current_edge(const int lines, const int tmp_lines, const int new_line[tmp_lines], int (*edge)[2],
					 int tmp_edge[tmp_lines][2], const int original_groups, const int g_opt, const bool is_center);
extern void create_adj(const int nodes, const int lines, const int degree,
                             const int edge[lines][2], int adj[nodes][degree]);
extern int getRandom(const int max);
extern void timer_clear_all();
extern void timer_clear(const int n);
extern void timer_start(const int n);
extern void timer_stop(const int n);
extern double timer_read(const int n);
extern bool evaluation(const int nodes, const int based_nodes, const int groups, const int lines, const int degree, 
		       int* restrict adj, int* restrict diameter, double* restrict ASPL, const int added_centers, const int algo);
extern int distance(int nodes, const int a, const int b, const int added_centers);
extern bool check(const int nodes, const int based_nodes, const int lines, const int degree, const int groups,
		  int edge[lines][2], const int add_degree_to_center, int* adj, const int ii);
extern void clear_buffer(int *buffer, const int n);
extern void clear_buffers(uint64_t* restrict A, uint64_t* restrict B, const int s);
#endif
