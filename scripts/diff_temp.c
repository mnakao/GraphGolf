#include <stdio.h>
#include <math.h>

int main()
{
  int n = 100;
  //  double max_temp = 5425.11, min_temp = 0.2171475;
    double max_temp = 523.70, min_temp = 0.2171475;
    //  double max_temp = 452.43 , min_temp = 0.2171475;
  double diff = pow(min_temp/max_temp, 1.0/(double)(n-1));
  for(int i=0;i<n;i++)
    printf("%f ", max_temp * pow(diff, i));
  return 0;
}
