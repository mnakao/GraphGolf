#include <sys/time.h>
#include "common.h"
static double elapsed[NUM_TIMERS], start[NUM_TIMERS];

static double elapsed_time()
{
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec + 1.0e-6 * t.tv_usec;
}

void timer_clear_all()
{
  for(int i=0;i<NUM_TIMERS;i++)
    elapsed[i] = 0.0;
}

void timer_clear(const int n)
{
  elapsed[n] = 0.0;
}

void timer_start(const int n)
{
  start[n] = elapsed_time();
}

void timer_stop(const int n)
{
  double now = elapsed_time();
  double t = now - start[n];
  elapsed[n] += t;
}

double timer_read(const int n)
{
  return( elapsed[n] );
}
