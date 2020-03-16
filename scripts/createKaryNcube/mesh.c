#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define DIMS 4
#define N1   10
#define N2   10
#define N3   10
#define N4   10

int main()
{
  int lines = 0, d[DIMS];

#if DIMS == 2
  fprintf(stderr, "%d x %d, degree = %d\n", N2, N1, 4);
  
  int n[N2][N1], edge[(N2*N1*4)/2][2];
  for(int j=0;j<N2;j++)
    for(int i=0;i<N1;i++)
      n[j][i] = j * N1 + i;

  for(int j=0;j<N2;j++){
    for(int i=0;i<N1;i++){
      int c = n[j][i];
      d[1] = c / N1;
      d[0] = c % N1;
      if(d[0] != N1 - 1){ // Add right line
      	edge[lines][0] = c;
      	edge[lines][1] = c + 1;
      	lines++;
      }
      if(d[1] != N2 - 1){ // Add down line
	edge[lines][0] = c;
	edge[lines][1] = c + N1;
	lines++;
      }
    }
  }

#elif DIMS == 3
  fprintf(stderr, "%d x %d x %d, degree = %d\n", N3, N2, N1, 6);

  int n[N3][N2][N1], edge[(N3*N2*N1*6)/2][2];
  for(int k=0;k<N3;k++)
    for(int j=0;j<N2;j++)
      for(int i=0;i<N1;i++)
	n[k][j][i] = k * (N2 * N1) + j * N1 + i;

  for(int k=0;k<N3;k++){
    for(int j=0;j<N2;j++){
      for(int i=0;i<N1;i++){
	int c = n[k][j][i];
	d[2] = c / (N2*N1);
	d[1] = (c - d[2]*(N2*N1)) / N1;
	d[0] = c % N1;
	if(d[0] != N1 - 1){
	  edge[lines][0] = c;
	  edge[lines][1] = c + 1;
	  lines++;
	}
	if(d[1] != N2 - 1){
	  edge[lines][0] = c;
	  edge[lines][1] = c + N1;
	  lines++;
	}
	if(d[2] != N3 - 1){
          edge[lines][0] = c;
          edge[lines][1] = c + N2 * N1;
          lines++;
        }
      }
    }
  }

#elif DIMS == 4
  fprintf(stderr, "%d x %d x %d x %d, degree = %d\n", N4, N3, N2, N1, 8);
  
  int n[N4][N3][N2][N1], edge[(N4*N3*N2*N1*8)/2][2];
  for(int m=0;m<N4;m++)
    for(int k=0;k<N3;k++)
      for(int j=0;j<N2;j++)
	for(int i=0;i<N1;i++)
	  n[m][k][j][i] = m * (N3 * N2 * N1) + k * (N2 * N1) + j * N1 + i;
  
  for(int m=0;m<N4;m++){
    for(int k=0;k<N3;k++){
      for(int j=0;j<N2;j++){
	for(int i=0;i<N1;i++){
	  int c = n[m][k][j][i];
	  d[3] = c / (N3*N2*N1);
	  d[2] = (c - d[3]*(N3*N2*N1)) / (N2*N1);
	  d[1] = (c - d[3]*(N3*N2*N1) - d[2]*(N2*N1)) / N1;
	  d[0] = c % N1;
	  if(d[0] != N1 - 1){
	    edge[lines][0] = c;
	    edge[lines][1] = c + 1;
	    lines++;
	  }
	  if(d[1] != N2 - 1){
	    edge[lines][0] = c;
	    edge[lines][1] = c + N1;
	    lines++;
	  }
	  if(d[2] != N3 - 1){
	    edge[lines][0] = c;
	    edge[lines][1] = c + N2 * N1;
	    lines++;
	  }
	  if(d[3] != N4 - 1){
            edge[lines][0] = c;
            edge[lines][1] = c + N3 * N2 * N1;
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
