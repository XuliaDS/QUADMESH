#include<math.h>

#include <stdio.h>
int main (int argc, char *argv[]) 
{
  int i, c8, l;
  double c1, c2, c3, c4, c5, c6, c7, ma, Ma;
  FILE *fil;
  
  if (argc != 2) { 
      printf(" NEED INPUT FILE! \n");
      return 0;
  }
  fil = fopen(argv[1], "r");
  if (fil == NULL ) {
      printf(" NULLL FILE! \n");
      return 0;
  }    
  l = 0;
  ma = 180.0;
  Ma = 0.0;
  while(getc(fil)!=EOF)
  {
    i = fscanf(fil, "%lf %lf %lf %lf %lf %lf %lf %d", &c1, &c2, &c3, &c4, &c5, &c6, &c7, &c8);
    if ( i !=  8) continue;
    printf(" L %d --> %lf %lf %lf %lf %lf %lf %lf %d\n",
     l,c1, c2, c3, c4, c5, c6, c7, c8);
    
    if ( c4 < ma ) {
        ma = c4;
        printf(" UPDATE   %lf \n", ma);
    }
    if ( c4 > Ma ) Ma = c5;
    l++;
  }
  printf("LINES READ %d MIN MAX ANGLES %lf  %lf  \n", l, ma, Ma);
  return 0;
}

  

