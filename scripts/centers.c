#include <stdio.h>

int main()
{
  //  int n = 3019, d = 30;
  //  int n = 4855, d = 30;
  //  int n=1726, d=30;
  int n=4856, d=15;
#define N 3
  int b_g[N] = {3, 5, 15};
  //  #define N 7
  //    int b_g[N] = {2, 3, 5, 6, 10, 15, 30};
  
  for(int i=0;i<N;i++){
    int g = b_g[i];
    for(int c=1;c<n;c++){
      if((n-c)%g != 0) continue;
      if(((n-c)/g) % 2 != 0) continue;
      if((((n-c)/g) - d/g * c)% 2 != 0) continue;
      if((((n-c)/g) - d/g * c) < 0) continue;
      int based_nodes = (n-c)/g;
      //      printf("g = %d c = %d : make K DATA=\"..\/data\/n%dd%d.random.edges\" G=%d V=%d E=%d S=$i\n", g, c, based_nodes, d-1, g, c, d/g);
      //printf("make K DATA=\"..\\/data\\/n%dd%d.random.edges\" G=%d V=%d E=%d S=$i\n", based_nodes, d-1, g, c, d/g);
      //      printf("make cygnus DATA=\"..\\/data\\/n%dd%d.random.edges\" G=%d V=%d E=%d N=10000000 T=2034.200008\n", based_nodes, d-1, g, c, d/g);
      printf("python ./scripts/create_random.py %d %d\n", based_nodes, d-1);
    }
  }
  return 0;
}
