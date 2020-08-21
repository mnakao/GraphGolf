#include <stdio.h>
#include <stdlib.h>
#define MAX_G 20

int main(int argc, char *argv[])
{
  if(argc == 1 || argc != 4){
    printf("./%s nodes degree centers\n", argv[0]);
    return 1;
  }

  int nodes   = atoi(argv[1]);
  int degree  = atoi(argv[2]);
  int centers = atoi(argv[3]);
  int g[MAX_G], num = 0;

  printf("nodes, degree, centers = %d %d %d\n", nodes, degree, centers);
  for(int i=2;i<=degree;i++){
    if(degree%i == 0){
      if(MAX_G == num + 1){
	printf("MAX_G is too small\n");
	return 1;
      }
      g[num++] = i;
    }
  }

  for(int i=0;i<num;i++)
    if( (nodes - centers)%g[i] == 0 )
      if( ((nodes - centers)/g[i])%2 == 0 || (degree-1)%2 == 0 )
	if( (nodes-centers-centers*degree) % (2 * g[i]) == 0 )
	  printf("g = %d\n", g[i]);
  
  return 0;
}
