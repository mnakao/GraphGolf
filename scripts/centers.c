#include <stdio.h>

int main()
{
    int n = 3019, d = 30;
    //     int n = 4855, d = 30;
  #define N 7
    int b_g[N] = {2, 3, 5, 6, 10, 15, 30};
  
  //   int n = 50, d = 4;
  //#define N 2
  //   int b_g[N] = {2, 4};

  for(int i=0;i<N;i++){
    int g = b_g[i];
    for(int c=1;c<n;c++){
      if((n-c)%g != 0) continue;
      if(((n-c)/g) % 2 != 0) continue;
      if((((n-c)/g) - d/g * c)% 2 != 0) continue;
      if((((n-c)/g) - d/g * c) < 0) continue;
      printf("g = %d c = %d\n", g, c);
      
    }
  }
  return 0;
}
