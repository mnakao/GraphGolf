#include <stdio.h>
#include <math.h>

int main()
{
  int n = 50;
  //  double max_temp = 5425.11, min_temp = 0.2171475;
  double max_temp = 6316.118889, min_temp = 0.2171475;
  double diff = pow(min_temp/max_temp, 1.0/(double)(n-1));
  for(int i=0;i<n;i++)
    printf("%f ", max_temp * pow(diff, i));
  return 0;
}
