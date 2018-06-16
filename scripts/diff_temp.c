#include <stdio.h>
#include <math.h>

int main()
{
  int n = 50;
  double max_temp = 7.0, min_temp = 2.5;
  double diff = pow(min_temp/max_temp, 1.0/(double)(n-1));
  for(int i=0;i<n;i++)
    printf("%f ", max_temp * pow(diff, i));
  return 0;
}
