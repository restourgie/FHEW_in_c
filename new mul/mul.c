#include <complex.h>
#include "mul.h"
#include <stdio.h>
#include <math.h>

#define M_PI 3.14159265358979323846

/******************************************************************
*
*	SUPPORT CODE
*
******************************************************************/
void print_complex(const double complex *a){
    for(int i=0;i<CPLXDIM;i++)
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

/******************************************************************
*
*	COMPLEX MULTIPLICATION
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
    // printf("\nNAIVE COMPLEX CALC\n");
    // print_complex(cplx_res);
    to_real(cplx_res,r);   
}

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

  	double complex root = 0+I*1;
  	root = csqrt(root);
  	
  	// printf("\n*********************** POLY X start ***********************\n");	
  	// print_complex(cplx_x);
  	recursive_phi(cplx_x,CPLXDIM,0,root);
  	// print_complex(cplx_x);
  	// printf("\n*********************** POLY X end ***********************\n");
  	// printf("\n***********************INVERSE POLY X start***********************\n");
  	// inverse_phi(cplx_x,CPLXDIM,0,root,1);
  	// printf("\n***********************INVERSE POLY X end***********************\n");
  	// print_complex(cplx_x);

  	// printf("\nSMART COMPLEX CALC POLY Y\n");
  	//print_complex(cplx_y);
  	recursive_phi(cplx_y,CPLXDIM,0,root);
  	// print_complex(cplx_y);


  	
 // 	double complex test[CPLXDIM];
 // 	test[0] = -266339.000000 + I * 110623.000000;
	// test[1] = -225798.000000 + I * 50944.000000;
	// test[2] = -195563.000000 + I * 71298.000000;
	// test[3] = -181704.000000 + I * 136320.000000;
	// test[4] = -146400.000000 + I * 116352.000000;
	// test[5] = -128248.000000 + I * 160754.000000;
	// test[6] = -109350.000000 + I * 278158.000000;
	// test[7] = -24202.000000 + I * 301782.000000;

	// printf("\nRESULT TEST\n");
	// print_complex(test);
	// printf("\nFFTING THE RESULT \n");
	// recursive_phi(test,CPLXDIM,0,root);
	// print_complex(test);
	// printf("\nBACKWARDS\n");
	// inverse_phi(test,CPLXDIM,0,root);
	// print_complex(test);


  	for (int i = 0; i < CPLXDIM; ++i)
  	{
  		cplx_res[i] = (cplx_x[i] * cplx_y[i])/CPLXDIM;
  	}

  	// printf("\nSMART COMPLEX CALC RESULT\n");
  	// print_complex(cplx_res);
  	inverse_phi(cplx_res,CPLXDIM,0,root);
  	// print_complex(cplx_res);

  	to_real(cplx_res,r);

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