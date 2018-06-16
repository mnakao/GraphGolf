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
#include <omp.h>
#include <stdbool.h>
#include <assert.h>

#define NOT_VISITED -1

#define RIGHT  0
#define LEFT   1
#define MIDDLE 2

#define NUM_TIMERS      2
#define TIMER_SA        0
#define TIMER_ESTIMATED 1
#define MAX_FILENAME_LENGTH 255
#define NUM_OF_PROGRESS 100
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define ABORT() {MPI_Abort(MPI_COMM_WORLD, 1); exit(1);}
#define PRINT_R0(...)  {if(rank==0) printf(__VA_ARGS__);}
#define END(...) {if(rank==0) printf(__VA_ARGS__); MPI_Finalize(); exit(0);}
#define ERROR(...) {if(rank==0) printf(__VA_ARGS__); ABORT();}

extern void swap(int *a, int *b);
extern int order(const int nodes, const int a, const int b);
extern long long sa(const int nodes, const int lines, const int degree, const int groups,  double temp, 
		    const long long ncalcs, const double cooling_rate, const int low_diam, 
		    const double low_ASPL, const bool hill_climbing_flag, 
		    const bool detect_temp_flag, double *max_diff_energy, int edge[lines][2],
		    int *diameter, double *ASPL, const int rank, const int size);
extern void check_current_edge(const int nodes, const int degree, const int lines, const int groups,
			       int edge[lines][2], const double low_ASPL, const int rank, const int size);
extern double estimated_elapse_time(const long long ncals, const int nodes, const int lines, const int degree,
				    const int groups, int edge[lines][2], const int rank, const int size);
extern bool edge_exchange_among_groups(int (*edge)[2], const int based_nodes, const int based_lines,
				       const int groups, const int start_line);
extern bool has_duplicated_edge(const int e00, const int e01, const int e10, const int e11);
extern bool check_loop(const int lines, int (*edge)[2]);
extern bool check_duplicate_edge(const int lines, int (*edge)[2]);
extern bool check_duplicate_current_edge(const int lines, const int groups, const int line[groups],
					 int (*edge)[2], int tmp_edge[groups][2]);
extern void create_adjacency(const int nodes, const int lines, const int degree,
                             int edge[lines][2], int adjacency[nodes][degree]);
extern int getRandom(const int max);
extern void timer_clear_all();
extern void timer_clear(const int n);
extern void timer_start(const int n);
extern void timer_stop(const int n);
extern double timer_read(const int n);
extern bool evaluation(const int nodes, const int groups, const int lines, const int degree, 
		       int adjacency[nodes][degree], int *diameter, double *ASPL, const int rank, const int size);
#endif
