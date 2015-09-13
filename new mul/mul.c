#include <complex.h>
#include "mul.h"
#include <stdio.h>
#include <math.h>

#define W(N,k) (cexp(2.0 * M_PI * I * (double)k / (double) N))


/******************************************************************
*
*	SUPPORT CODE
*
******************************************************************/
void print_complex(const double complex *a, int N){
    for(int i=0;i<N;i++)
      printf("cplxpoly[%d] = %f + i * %f\n",i,creal(a[i]),cimag(a[i]));
    printf("\n");
}

/******************************************************************
*
*	CONVERSION
*
******************************************************************/
void to_complex(const ring_t *x, double complex *cplx_x)
{
  for(int i=0;i<CPLXDIM;++i)
    cplx_x[i] = x->v[i] + I*x->v[i+CPLXDIM];
}

void to_real(const double complex *cplx_x, ring_t *x)
{
  for(int i=0;i<CPLXDIM;++i){
    x->v[i] = round(creal(cplx_x[i]));
    x->v[i+CPLXDIM] = round(cimag(cplx_x[i]));
  }
}

void twist(double complex *cplx_x,int n,int m,int lo)
{
  // printf("n = %d, m = %d, lo = %d\n",n,m,lo );
  int j = 1;
  for (int i = lo+1; i < lo+m; ++i)
  {
    // printf("i = %d, j = %d\n",i,j);
    cplx_x[i] = cplx_x[i] * W(n,j);
    ++j;
  }
}

void untwist(double complex *cplx_x,int n,int m,int lo)
{
  // printf("n = %d, m = %d, lo = %d\n",n,m,lo );
  int j = 1;
  for (int i = lo+1; i < lo+m; ++i)
  {
    // printf("i = %d, j = %d\n",i,j);
    cplx_x[i] = cplx_x[i] * conj(W(n,j));
    ++j;
  }
}

/******************************************************************
*
* SMART COMPLEX MULTIPLICATION
*
******************************************************************/

void inverse_phi(double complex *x,int n,int lo,double complex root)
{	
	if(n > 1){
		int m = n/2;
		inverse_phi(x,m,lo,csqrt(root));
		inverse_phi(x,m,lo+m,csqrt(-root));
		double complex temp;
		for(int i=lo;i<m+lo;++i){
			// printf("\nN is at the moment: %d\n",n);
			// printf("lo = %d, i = %d, m = %d\n",lo,i,m);
			// printf("x[i]= ( %f + I * %f) x[i+m] = ( %f + I * %f)\n",creal(x[i]),cimag(x[i]),creal(x[i+m]),cimag(x[i+m]) );
			temp = x[i] + x[i+m];
			x[i+m] =  (x[i] - x[i+m]) * conj(root);

			x[i] = temp;
			// printf("x[i]= ( %f + I * %f) x[i+m] = ( %f + I * %f)\n",creal(x[i]),cimag(x[i]),creal(x[i+m]),cimag(x[i+m]) );
		}
	}
}

void recursive_phi(double complex *x,int n,int lo,double complex root)
{	
  	if(n > 1){
    double complex temp;
    int m = n/2;
    // printf("\nN is at the moment: %d\n",n);
    // printf("Root: ( %f + I * %f)\n",creal(root),cimag(root));
    for(int i=lo; i < lo+m;++i){
      temp = root * x[i+m];
      // printf("lo = %d, i = %d, m = %d\n",lo,i,m);
      // printf("temp: ( %f + I * %f)\n",creal(temp),cimag(temp) );
	      //phiprime
      x[i+m] = x[i] - temp;
	      //phi
      x[i] = x[i] + temp;
    }
    // print_complex(x);
    recursive_phi(x,m,lo,csqrt(root));
    recursive_phi(x,m,lo + m,csqrt(-root));
  }
}


void smart_complex_mul(ring_t *r, const ring_t *x, const ring_t *y)
{
	double complex cplx_x[CPLXDIM];
	double complex cplx_y[CPLXDIM];
	double complex cplx_res[CPLXDIM];

	to_complex(x,cplx_x);
	to_complex(y,cplx_y);



	double complex root = I;
	root = csqrt(root);
	
	recursive_phi(cplx_x,CPLXDIM,0,root);
	recursive_phi(cplx_y,CPLXDIM,0,root);


	for (int i = 0; i < CPLXDIM; ++i)
	{
		cplx_res[i] = (cplx_x[i] * cplx_y[i])/CPLXDIM;
	}

	inverse_phi(cplx_res,CPLXDIM,0,root);
  printf("\n\n**************SMART COMPLEX MUL RESULT**************\n");
  print_complex(cplx_res,CPLXDIM);
	to_real(cplx_res,r);
}

/******************************************************************
*
* Normal FFT MULTIPLICATION
*
******************************************************************/
void recursive_FFT(double complex *x,int n,int lo,double complex root)
{ 

  if(n > 1){
    double complex temp;
    int m = n/2;
    // printf("\nN is at the moment: %d\n",n);
    // printf("Root: ( %f + I * %f)\n",creal(root),cimag(root));
    for(int i=lo; i < lo+m;++i){
      temp = root * x[i+m];
      // printf("lo = %d, i = %d, m = %d\n",lo,i,m);
      // printf("temp: ( %f + I * %f)\n",creal(temp),cimag(temp) );
        //phiprime
      x[i+m] = x[i] - temp;
        //phi
      x[i] = x[i] + temp;
    }
    // print_complex(x,REALDIM);
    recursive_FFT(x,m,lo,csqrt(root));
    if(root == 1){
      root = I;
      recursive_FFT(x,m,lo+m,root);
    }
    else  
      recursive_FFT(x,m,lo + m,csqrt(-root));
  }
}

void inverse_FFT(double complex *x,int n,int lo,double complex root)
{ 
  if(n > 1){
    int m = n/2;
    inverse_FFT(x,m,lo,csqrt(root));
    if(root == 1){
      inverse_FFT(x,m,lo+m,I);
    }
    else
      inverse_FFT(x,m,lo+m,csqrt(-root));

    // printf("\nN is at the moment: %d\n",n);
    // printf("lo = %d, m = %d\n",lo,m);
    // printf("Root: ( %f + I * %f)\n",creal(root),cimag(root));
    double complex temp;
    for(int i=lo;i<m+lo;++i){
      // printf("x[i]= ( %f + I * %f) x[i+m] = ( %f + I * %f)\n",creal(x[i]),cimag(x[i]),creal(x[i+m]),cimag(x[i+m]) );
      temp = x[i] + x[i+m];
      x[i+m] =  ((x[i] - x[i+m]) * conj(root));

      x[i] = temp;
      // printf("x[i]= ( %f + I * %f) x[i+m] = ( %f + I * %f)\n",creal(x[i]),cimag(x[i]),creal(x[i+m]),cimag(x[i+m]) );
    }
    //print_complex(x,REALDIM);
  }
}

//DOES CALCULATES X * Y = R mod (x^n - 1)
void normal_FFT_mul(ring_t *r, const ring_t *x, const ring_t *y)
{
  double complex cplx_x[REALDIM];
  double complex cplx_y[REALDIM];
  double complex cplx_res[REALDIM];

  for (int i = 0; i < REALDIM; ++i)
  {
    cplx_x[i] = x->v[i];
    cplx_y[i] = y->v[i];
  }

  double complex root = 1;
  
  // printf("\n\n**************APPLY FFT**************\n");
  recursive_FFT(cplx_x,REALDIM,0,root);
  // print_complex(cplx_x,REALDIM);



  //printf("\n\n**************NORMAL FFT Y**************\n");
  recursive_FFT(cplx_y,REALDIM,0,root);
  //print_complex(cplx_y,CPLXDIM);


  for (int i = 0; i < REALDIM; ++i)
  {
    cplx_res[i] = (cplx_x[i] * cplx_y[i])/REALDIM;
  }

  // printf("\n\n**************NORMAL FFT MUL RESULT**************\n");
  inverse_FFT(cplx_res,REALDIM,0,root);
  
  // print_complex(cplx_res,REALDIM);
  for (int i = 0; i < REALDIM; ++i)
  {
    r->v[i] = round(creal(cplx_res[i]));
  }
  
}

/******************************************************************
*
* TWISTED FFT MULTIPLICATION
*
******************************************************************/
void twisted_recursive_FFT(double complex *x,int n,int lo)
{
  if(n > 1)
  {
    int m = n/2;
    double complex temp;

    for(int i=lo; i < lo+m;++i){
      temp = x[i];
      x[i] = temp + x[i+m];
      x[i+m] = temp - x[i+m];
    }
    // printf("n = %d\n",n );
    // printf("\n\n**************RECURSION**************\n");
    twisted_recursive_FFT(x,m,lo);
    // printf("\n\n**************TWISTING**************\n");
    // print_complex(x,REALDIM);
    twist(x,n,m,lo+m);
    // print_complex(x,REALDIM);
    twisted_recursive_FFT(x,m,lo+m);
  }
}

void twisted_inverse_FFT(double complex *x,int n,int lo)
{
  if(n > 1)
  {
    int m = n/2;
    double complex temp;

    twisted_inverse_FFT(x,m,lo);
    twisted_inverse_FFT(x,m,lo+m);
    untwist(x,n,m,lo+m);

    for (int i = lo; i < lo+m; ++i)
    {
      temp = x[i];
      x[i] = temp + x[i+m];
      x[i+m] = temp - x[i+m];
    }

  }
}

void twisted_FFT_mul(ring_t *r, const ring_t *x, const ring_t *y)
{
  double complex cplx_x[REALDIM];
  double complex cplx_y[REALDIM];
  double complex cplx_res[REALDIM];

  for (int i = 0; i < REALDIM; ++i)
  {
    cplx_x[i] = x->v[i];
    cplx_y[i] = y->v[i];
  }
  // printf("\n\n**************APPLY TWISTED FFT**************\n");
  // print_complex(cplx_x,REALDIM);
  twisted_recursive_FFT(cplx_x,REALDIM,0);
  // print_complex(cplx_x,REALDIM);
  twisted_recursive_FFT(cplx_y,REALDIM,0);

  for (int i = 0; i < REALDIM; ++i)
  {
    cplx_res[i] = (cplx_x[i] * cplx_y[i])/REALDIM;
  }
  // printf("\n\n**************TWISTED FFT MUL RESULT Before inverse**************\n");
  // print_complex(cplx_res,REALDIM);
  // printf("\n\n**************TWISTED FFT MUL RESULT**************\n");
  twisted_inverse_FFT(cplx_res,REALDIM,0);
  // print_complex(cplx_res,REALDIM);

  for (int i = 0; i < REALDIM; ++i)
  {
    r->v[i] = round(creal(cplx_res[i]));
  }

}

/******************************************************************
*
* SPLIT RADIX FFT MULTIPLICATION
*
******************************************************************/
void split_radix_recursive(double complex *x,int n,int lo)
{
  double complex temp;
  if(n == 2){
    temp = x[lo];
    x[lo] = temp + x[lo+1];
    x[lo+1] = temp - x[lo+1];
  }
  else if(n > 2){
    int m = n/2;
    //Go from (x^4n +1 to x^2n -1 and x^2n +1)
    for(int i=lo; i < lo+m;++i){
      temp = x[i];
      x[i] = temp + x[i+m];
      x[i+m] = temp - x[i+m];
    }
    //Do recursive step for (x^2n -1)
    split_radix_recursive(x,m,lo);

    lo = lo+m;
    //temp_n = m;
    m = m/2;

    //Go from (x^2n +1 to x^n -i and x^n +i)
    for (int i = lo; i < lo+m; ++i)
    {
      temp = x[i];
      x[i] = temp + I * x[i+m];
      x[i+m] = temp - I * x[i+m];
    }
    twist(x,n,m,lo);
    split_radix_recursive(x,m,lo);
    untwist(x,n,m,lo+m);
    split_radix_recursive(x,m,lo+m);
  }

}

void split_radix_inverse(double complex *x,int n,int lo)
{
  // printf("N = %d\n",n);
  // print_complex(x,REALDIM);
  double complex temp;
  if(n == 2){
    temp = x[lo];
    x[lo] = temp + x[lo+1];
    x[lo+1] = temp - x[lo+1];
  }
  else if(n > 2){
    // printf("n = %d lo = %d\n",n,lo );
    int m = n/4;
    lo = lo+n/2;
    // printf("m = %d lo = %d\n",m,lo );
    split_radix_inverse(x,m,lo+m);
    twist(x,n,m,lo+m);
    split_radix_inverse(x,m,lo);
    untwist(x,n,m,lo);
    for (int i = lo; i < lo+m; ++i)
    {
      temp = x[i];
      x[i] = temp + x[i+m];
      x[i+m] = (temp - x[i+m])*-I;
    }
    m = m*2;
    lo = lo -m;
    // printf("m = %d lo = %d\n",m,lo );
    split_radix_inverse(x,m,lo);
    for(int i=lo; i < lo+m;++i){
      temp = x[i];
      x[i] = temp + x[i+m];
      x[i+m] = temp - x[i+m];
    }

  }

}

void split_radix_FFT_mul(ring_t *r, const ring_t *x, const ring_t *y)
{
  double complex cplx_x[REALDIM];
  double complex cplx_y[REALDIM];
  double complex cplx_res[REALDIM];

  for (int i = 0; i < REALDIM; ++i)
  {
    cplx_x[i] = x->v[i];
    cplx_y[i] = y->v[i];
  }
  // printf("\n\n**************APPLY split-radix FFT**************\n");
  // print_complex(cplx_x,REALDIM);
  split_radix_recursive(cplx_x,REALDIM,0);

  // print_complex(cplx_x,REALDIM);
  split_radix_recursive(cplx_y,REALDIM,0);

  for (int i = 0; i < REALDIM; ++i)
  {
    cplx_res[i] = (cplx_x[i] * cplx_y[i])/REALDIM;
  }
  // printf("\n\n**************TWISTED FFT MUL RESULT Before inverse**************\n");
  // print_complex(cplx_res,REALDIM);
  // printf("\n\n**************split-radix FFT MUL RESULT**************\n");
  split_radix_inverse(cplx_res,REALDIM,0);
  // print_complex(cplx_res,REALDIM);

  for (int i = 0; i < REALDIM; ++i)
  {
    r->v[i] = round(creal(cplx_res[i]));
  }

}

/******************************************************************
*
*	NAIVE SCHOOLBOOK MULTIPLICATION
*
******************************************************************/
/* Very simple schoolbook multiplication. Works. */
void naive_real_mul(ring_t *r, const ring_t *x, const ring_t *y)
{
  int i,j;
  for(i=0;i<REALDIM;i++)
    r->v[i] = 0;

  for(i=0;i<REALDIM;i++)
  {
    for(j=0;j<REALDIM;j++)
    {
      if(i+j < REALDIM)
        r->v[i+j] += x->v[i] * y->v[j];
      else
        r->v[i+j-REALDIM] -= x->v[i] * y->v[j];
    }
  }
}

/******************************************************************
*
* NAIVE CYCLIC SCHOOLBOOK MULTIPLICATION
*
******************************************************************/
/* Very simple schoolbook multiplication. Works. */
void naive_cyclic_real_mul(ring_t *r, const ring_t *x, const ring_t *y)
{
  int i,j;
  for(i=0;i<REALDIM;i++)
    r->v[i] = 0;

  for(i=0;i<REALDIM;i++)
  {
    for(j=0;j<REALDIM;j++)
    {
      if(i+j < REALDIM)
        r->v[i+j] += x->v[i] * y->v[j];
      else
        r->v[i+j-REALDIM] += x->v[i] * y->v[j];
    }
  }
}

/******************************************************************
*
* COMPLEX MULTIPLICATION
*
******************************************************************/
void naive_complex_mul(ring_t *r, const ring_t *x, const ring_t *y)
{ 
  double complex cplx_x[CPLXDIM];
  double complex cplx_y[CPLXDIM];
  double complex cplx_res[CPLXDIM];

  to_complex(x,cplx_x);
  to_complex(y,cplx_y);

  double complex t;
  double complex big[REALDIM];

  for(int i=0;i<REALDIM;++i)
    big[i] = 0;

  for(int i =0;i<CPLXDIM;++i)
    for(int j=0;j<CPLXDIM;++j)
      big[i+j] += cplx_x[i] * cplx_y[j];    

  for(int i=CPLXDIM;i<REALDIM;++i){
    // printf("%f + i%f\n",creal(big[i+512]),cimag(big[i+512]));
    t = big[i] * (0+I*1);
    // printf("%f + i%f\n",creal(t),cimag(t));
    cplx_res[i-CPLXDIM] = big[i-CPLXDIM] + t;
  }
  printf("\n\n**************NAIVE COMPLEX CALC**************\n");
  print_complex(cplx_res,CPLXDIM);
  to_real(cplx_res,r);   
}
