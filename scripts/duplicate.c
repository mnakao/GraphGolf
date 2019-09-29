#include <stdio.h>
#include <stdlib.h>

int main()
{
  int nodes  = 1728;
  int degree = 30;  // degree must be even number
  
  for(int n=0;n<nodes;n++){
    for(int d=0;d<degree/2;d++){
      int n_plus_1 = n + 1;
      if(n_plus_1 == nodes) n_plus_1 = 0;
      printf("%d %d\n", n, n_plus_1);
    }
  }
  return 0;
}
