#include <stdio.h>
#define WIDTH  100
#define HEIGHT 100
#define DEGREE 4
#define LENGTH 2

int main()
{
  int lines = (WIDTH*HEIGHT*DEGREE)/2, eid = 0;
  int edge[lines][2];
  
  for(int j=0;j<HEIGHT;j++){
    for(int i=0;i<WIDTH;i++){
      if(j != HEIGHT-1){
	int n = i * HEIGHT + j;
	edge[eid][0] = n;
	edge[eid][1] = n + 1;
	eid++;
      }
    }
  }

  for(int j=0;j<HEIGHT;j++){
    for(int i=0;i<WIDTH;i++){
      if(i != WIDTH-1){
        int n = i * HEIGHT + j;
        edge[eid][0] = n;
      	edge[eid][1] = n + HEIGHT;
      	eid++;
      }
    }
  }

  for(int j=0;j<HEIGHT/2;j++){
    int n = j;
    edge[eid][0] = n * 2;
    edge[eid][1] = edge[eid][0] + 1;
    eid++;
  }

  for(int i=0;i<WIDTH/2;i++){
    int n = i * HEIGHT;
    edge[eid][0] = n * 2;
    edge[eid][1] = edge[eid][0] + HEIGHT;
    eid++;
  }

  for(int j=0;j<HEIGHT/2;j++){
    int n = (WIDTH-1) * HEIGHT;
    edge[eid][0] = n + j * 2;
    edge[eid][1] = edge[eid][0] + 1;
    eid++;
  }

  for(int i=0;i<WIDTH/2;i++){
    int n = (HEIGHT-1) + HEIGHT * i * 2;
    edge[eid][0] = n;
    edge[eid][1] = edge[eid][0] + HEIGHT;
    eid++;
  }
  
  // output
  for(int i=0;i<eid;i++)
    printf("%d,%d %d,%d\n",
	   edge[i][0]%HEIGHT, edge[i][0]/HEIGHT, edge[i][1]%HEIGHT, edge[i][1]/HEIGHT);
  
  return 0;
}
