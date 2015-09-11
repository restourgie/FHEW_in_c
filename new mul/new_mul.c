#include <stdio.h>
#include <math.h>


/******************************************************************
*
*	FFT PHI
*
******************************************************************/
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

void iterative_phi(double complex *x)
{
  //variables for the fft
  unsigned long mmax,m,j,istep,i;


}



/******************************************************************
*
* SMART NEW MULTIPLICATION
*
******************************************************************/
void smart_mul(ring_t *r, const ring_t *x, const ring_t *y)
{
	double cplx_x[REALDIM];
	double cplx_y[REALDIM];
	double cplx_r[REALDIM];

	for (int i = 0; i < REALDIM; ++i)
	{
		cplx_x = x.v[i];
		cplx_y = y.v[i];
	}


}