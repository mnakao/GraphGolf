#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define DIMS 2
#define N    16

int main()
{
  int lines = 0, d[DIMS];

#if DIMS == 2
  fprintf(stderr, "%d x %d, degree = %d\n", N, N, 4);
  
  int n[N][N], edge[(N*N*4)/2][2];
  for(int j=0;j<N;j++)
    for(int i=0;i<N;i++)
      n[j][i] = j * N + i;

  for(int j=0;j<N;j++){
    for(int i=0;i<N;i++){
      int c = n[j][i];
      d[1] = c / N;
      d[0] = c % N;
      if(d[0] != N - 1){ // Add right line
      	edge[lines][0] = c;
      	edge[lines][1] = c + 1;
      	lines++;
      }
      if(d[1] != N - 1){ // Add down line
	edge[lines][0] = c;
	edge[lines][1] = c + N;
	lines++;
      }
    }
  }

#elif DIMS == 3
  fprintf(stderr, "%d x %d x %d, degree = %d\n", N, N, N, 6);

  int n[N][N][N], edge[(N*N*N*6)/2][2];
  for(int k=0;k<N;k++)
    for(int j=0;j<N;j++)
      for(int i=0;i<N;i++)
	n[k][j][i] = k * (N * N) + j * N + i;

  for(int k=0;k<N;k++){
    for(int j=0;j<N;j++){
      for(int i=0;i<N;i++){
	int c = n[k][j][i];
	d[2] = c / (N*N);
	d[1] = (c - d[2]*(N*N)) / N;
	d[0] = c % N;
	if(d[0] != N - 1){
	  edge[lines][0] = c;
	  edge[lines][1] = c + 1;
	  lines++;
	}
	if(d[1] != N - 1){
	  edge[lines][0] = c;
	  edge[lines][1] = c + N;
	  lines++;
	}
	if(d[2] != N - 1){
          edge[lines][0] = c;
          edge[lines][1] = c + N * N;
          lines++;
        }
      }
    }
  }

#elif DIMS == 4
  fprintf(stderr, "%d x %d x %d x %d, degree = %d\n", N, N, N, N, 8);
  
  int n[N][N][N][N], edge[(N*N*N*N*8)/2][2];
  for(int m=0;m<N;m++)
    for(int k=0;k<N;k++)
      for(int j=0;j<N;j++)
	for(int i=0;i<N;i++)
	  n[m][k][j][i] = m * (N * N * N) + k * (N * N) + j * N + i;
  
  for(int m=0;m<N;m++){
    for(int k=0;k<N;k++){
      for(int j=0;j<N;j++){
	for(int i=0;i<N;i++){
	  int c = n[m][k][j][i];
	  d[3] = c / (N*N*N);
	  d[2] = (c - d[3]*(N*N*N)) / (N*N);
	  d[1] = (c - d[3]*(N*N*N) - d[2]*(N*N)) / N;
	  d[0] = c % N;
	  if(d[0] != N - 1){
	    edge[lines][0] = c;
	    edge[lines][1] = c + 1;
	    lines++;
	  }
	  if(d[1] != N - 1){
	    edge[lines][0] = c;
	    edge[lines][1] = c + N;
	    lines++;
	  }
	  if(d[2] != N - 1){
	    edge[lines][0] = c;
	    edge[lines][1] = c + N * N;
	    lines++;
	  }
	  if(d[3] != N - 1){
            edge[lines][0] = c;
            edge[lines][1] = c + N * N * N;
            lines++;
          }
        }
      }
    }
  }

#else
  printf("DIMS must be 2 or 3 or 4\n");
  exit(0);
#endif

  for(int i=0;i<lines;i++)
    printf("%d %d\n", edge[i][0], edge[i][1]);
  return 0;
}
