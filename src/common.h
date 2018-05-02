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

#define NOT_VISITED -1

#define NUM_TIMERS      2
#define TIMER_SA        0
#define TIMER_ESTIMATED 1
#define MAX_FILENAME_LENGTH 255
#define NUM_OF_PROGRESS 100
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define ABORT {MPI_Abort(MPI_COMM_WORLD, 1); exit(1);}

extern void swap(int *a, int *b);
extern int order(const int nodes, const int a, const int b);
extern long long sa(const int nodes, const int lines, const int degree, double temp, 
		    const long long ncalcs, const double cooling_rate, const int low_diam, 
		    const double low_ASPL, const bool hill_climbing_flag, 
		    const bool detect_temp_flag, double *max_diff_energy, int edge[lines][2],
		    int *diameter, double *ASPL, const int rank, const int size);
extern void check_current_edge(const int nodes, const int degree, const int lines, int edge[lines][2], 
			       const double low_ASPL, const int rank, const int size);
extern double estimated_elapse_time(const long long ncals, const int nodes, const int lines, const int degree,
				    int edge[lines][2], const int rank, const int size);
extern void timer_clear_all();
extern void timer_clear(const int n);
extern void timer_start(const int n);
extern void timer_stop(const int n);
extern double timer_read(const int n);

extern bool evaluation(const int nodes, const int lines, const int degree, int adjacency[nodes][degree], 
		       int *diameter, double *ASPL, const int rank, const int size);
#endif
