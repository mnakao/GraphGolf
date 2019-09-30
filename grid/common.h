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
#include <nmmintrin.h>
#ifdef _OPENMP
#include <omp.h>
#endif

int rank, procs, threads;
#define VISITED         1
#define NOT_VISITED     0
#define NOT_DEFINED    -1
#define NUM_TIMERS      4
#define TIMER_SA        0
#define TIMER_ESTIMATED 1
#define TIMER_APSP      2
#define TIMER_CHECK     3
#define NOT_USED       -1
#define D_2G_OPT        2
#define D_1G_OPT        1

#ifdef _KCOMPUTER
#define POPCNT(a) __builtin_popcountll(a)
#else
#define POPCNT(a) _mm_popcnt_u64(a)
#endif
#define UINT64_BITS         64
#define CHUNK               64 /* (multiple of sizeof(uint64_t)*8 for AVX-512) */

#define MAX_FILENAME_LENGTH 255
#define NUM_OF_PROGRESS     100
#define SKIP_ACCEPTS        10000
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define ABORT()        do{MPI_Abort(MPI_COMM_WORLD, 1); exit(1);}while(0)
#define PRINT_R0(...)  do{if(rank==0) printf(__VA_ARGS__);}while(0)
#define END(...)       do{if(rank==0) printf(__VA_ARGS__); MPI_Finalize(); exit(0);}while(0)
#define ERROR(...)     do{if(rank==0) printf(__VA_ARGS__); ABORT();}while(0)
#define EXIT(r)        do{MPI_Finalize(); exit(r);}while(0)
#define OUTPUT_EDGE()  do{output_edge(lines, edge, height);}while(0)

extern void swap(int *a, int *b);
extern long long sa(const int nodes, const int lines, double temp, const long long ncalcs,
		    const double cooling_rate, const int low_diam, const double low_ASPL, const bool enable_bfs, 
		    const bool hill_climbing_flag, const bool detect_temp_flag, double *max_diff_energy,
		    const double max_temp, const double min_temp, int edge[lines][2], int *diam, double *ASPL,
		    const int cooling_cyclie, long long *num_accepts, const int width, const int based_width,
		    const int height, const int based_height, int *length, const int low_length, const double weight, const int groups);
extern void check_current_edge(const int nodes, const int lines, int edge[lines][2], const double low_ASPL, const bool enable_bfs);
extern double estimated_elapse_time(const int nodes, const int lines, const int edge[lines][2], const int height, const int width,
				    const int groups, const bool enable_bfs);
extern bool has_duplicated_edge(const int e00, const int e01, const int e10, const int e11);
extern bool has_duplicated_vertex(const int e00, const int e01, const int e10, const int e11);
extern bool check_loop(const int lines, const int edge[lines][2]);
extern bool check_duplicate_all_edge(const int lines, const int edge[lines][2]);
extern bool check_duplicate_current_edge(const int lines, const int edge[lines][2], const int tmp_lines,
					 const int tmp_edge[tmp_lines][2], const int tmp_line[2],
					 const int groups, const int g_opt, const bool is_center);
extern void create_adjacency(const int nodes, const int lines, const int degree,
                             const int edge[lines][2], int adjacency[nodes][degree]);
extern int getRandom(const int max);
extern void timer_clear_all();
extern void timer_clear(const int n);
extern void timer_start(const int n);
extern void timer_stop(const int n);
extern double timer_read(const int n);
extern bool evaluation(const int nodes, const int lines, const int degree, 
		       const int* restrict adjacency, int *diameter, double *ASPL, const bool enable_bfs);
extern void copy_edge(int *restrict buf1, const int *restrict buf2, const int n);
extern bool edge_1g_opt(int (*edge)[2], const int nodes, const int lines, const int degree, const int based_nodes,
			const int based_lines, const int height, const int groups, const int start_line, const long long ii);
extern void output_edge(const int lines, const int edge[lines][2], const int height);
extern bool check_symmetric_edge(const int lines, const int edge[lines][2], const int height, const int width,
				 const int based_height, const int groups);
extern bool check_vector(const int groups, const int lines, const int height, const int edge[lines][2]);
extern int WIDTH (const int v, const int height);
extern int HEIGHT(const int v, const int height);
extern int ROTATE(const int v, const int height, const int degree);
#endif
