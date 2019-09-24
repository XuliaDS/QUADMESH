#include<math.h>

#include <stdio.h>
#define PI 3.1415926535897931159979635

#define ERFC(a,b,x) (0.5 * ( 1 - erf(3.554147 * (x -(0.5 * (a + b))) / (0.5 * (a - b)))))
#define ANGCUT 2.9

int main (void) 
{
  int i, N = 10000;
  double theta, f;
  for ( i = 0 ; i < N; i++) {
   theta = (double)i/(double)N * 3.0 * PI/ 2.0;
   f = ERFC(ANGCUT, PI, theta);
   printf(" THETA %lf  ERFC %lf \n", theta, f);
  }
  return 0;
}

  

