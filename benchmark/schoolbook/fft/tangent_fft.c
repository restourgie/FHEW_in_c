#include <complex.h>
#include <stdio.h>
#include <math.h>
#include "../mul.h"

void base_change(double complex *x,int lo,int m, int n,int n_new)
{ int j = 0;
  double temp;
  for (int i = lo; i <= lo+m; ++i)
  { 
    temp = s(n_new,j)/s(n,j);
    printf("temp = %f, x[%d] = %f + I %f\n",temp,i,creal(x[i]),cimag(x[i]));
    x[i] = x[i] * temp;
    ++j;
  }
}

void cost_4_twist(double complex *cplx_x,int n,int m,int lo)
{
  int j = 0;
  double temp;
  for (int i = lo; i <= lo+m; ++i)
  {
    temp = s(n/4,j)/s(n,j);
    printf("temp = %f, x[%d] = %f + I %f\n",temp,i,creal(cplx_x[i]),cimag(cplx_x[i]));
    cplx_x[i] = cplx_x[i] * W(n,j)*temp;
    ++j;
  }
}

void cost_4_untwist(double complex *cplx_x,int n,int m,int lo)
{
  int j = 0;
  double temp;
  for (int i = lo; i <= lo+m; ++i)
  {
    temp = s(n/4,j)/s(n,j);
    cplx_x[i] = cplx_x[i] * conj(W(n,j))*temp;
    ++j;
  }
}

void tangent_8(double complex *x,int n,int lo)
{ 
  double complex temp;
  if(n == 2){
    temp = x[lo];
    x[lo] = temp + x[lo+1];
    x[lo+1] = temp - x[lo+1];
    print_complex(x,8);
  }
  else if(n>=8)
  {
    int m = n/2;
    //Go from (x^8n +1 to x^4n -1 and x^4n +1)
    for(int i=lo; i < lo+m;++i){
      temp = x[i];
      x[i] = temp + x[i+m];
      x[i+m] = temp - x[i+m];
    }
    //NOW take the left side
    //Go from (x^4n-1) to (x^2n-1) and (x^2n+1)
    int right_lo = lo + m;
    m = m/2;
    for(int i=lo; i < lo+m;++i){
      temp = x[i];
      x[i] = temp + x[i+m];
      x[i+m] = temp - x[i+m];
    }
    //NOW take the right side
    //Go from (x^4n+1) to (x^2n-i) and (x^2n+i)
    for (int i = right_lo; i < right_lo+m; ++i)
    {
      // printf("i = %d, m = %d\n",i,m );
      temp = x[i];
      x[i] = temp + I * x[i+m];
      x[i+m] = temp - I * x[i+m];
    }
    //BEFORE CHANGING M SET OTHER LO VALUES FOR LEFT_RIGHT AND RIGHT_RIGHT BRANCH
    int left_right_lo = lo + m;
    int right_right_lo = right_lo + m;
    m = m/2;
    //DO COST-4 TWIST
    cost_4_twist(x,n,m,right_lo);
    tangent_8(x,n/4,right_lo);
    //DO COST-4 UNTWIST
    cost_4_untwist(x,n,m,right_right_lo);
    tangent_8(x,n/4,right_right_lo);
    //BACK TO THE LEFT
    //NOW TO A BASE TWIST GO FROM X^2n-1 BASE S(8n,k) to X^2n-1 BASE S(2n,k)
    base_change(x,lo,m,n,n/4);
    tangent_8(x,n/4,lo);

    //FINISH THE LEFT_RIGHT BRANCH
    //ALSO BASE CHANGE FOR MID GO FROM X^2n+1 BASE S(8n,k) to X^2n+1 BASE S(4n,k)
    base_change(x,left_right_lo,m,n,n/2);
    //GO FROM (X^2n+1) TO (X^n-i) AND (x^n+i)
    for (int i = left_right_lo; i < left_right_lo+m; ++i)
    {
      // printf("i = %d, m = %d\n",i,m );
      temp = x[i];
      x[i] = temp + I * x[i+m];
      x[i+m] = temp - I * x[i+m];
    }
    m=m/2;
    cost_4_twist(x,n/2,m,left_right_lo);
    tangent_8(x,m,left_right_lo);
    cost_4_untwist(x,n/2,m,left_right_lo+m);
    tangent_8(x,m,left_right_lo+m);
  }
}
