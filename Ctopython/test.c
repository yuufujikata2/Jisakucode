#include<stdio.h>

int add(int* x, int* y)
{
  int result;

  result = *x + *y;
  *x = 10;
  return result ;
}
  
