#include <stdio.h> 
#include <math.h> 
#include <stdlib.h> 
#include <complex.h> 

#define CPLXDIM 8
#define REALDIM (2*CPLXDIM)

typedef double complex data_t;

#define W(N,k) (cexp(-2.0f * M_PI * I * (double)k / (double) N))

void print_complex(const data_t *a){
    for(int i=0;i<CPLXDIM;i++)
      printf("cplxpoly[%d] = %f + i * %f\n",i,creal(a[i]),cimag(a[i]));
    printf("\n");
}

void ditfft2(data_t *in,data_t *out,int stride,int N){
	// print_complex(out);
	if(N == 2){
		out[0] = in[0] + in[stride];
		out[N/2] = in[0] - in[stride];
	}
	else{
		ditfft2(in,out,stride << 1,N >> 1);
		ditfft2(in+stride,out+N/2,stride <<1, N>>1);

		{ /* k=0 -> no mult */
			data_t Ek = out[0];
			data_t Ok = out[N/2];
			out[0] = Ek + Ok;
			out[N/2] = Ek - Ok;
		}

		int k;
		for(k=1;k<N/2;k++){
			data_t Ek = out[k];
			data_t Ok = out[(k+N/2)];
			out[k] = 		Ek + W(N,k) * Ok;
			out[(k+N/2)] = 	Ek - W(N,k) *Ok;
		}
	}
}

int main()
{	data_t in[CPLXDIM];
	data_t out[CPLXDIM];

	for (int i = 0; i < CPLXDIM; ++i)
		{	
			out[i] =0;
			if(i < CPLXDIM/2)
				in[i] = 1;
			else
				in[i] = 0;
		}

	print_complex(in);	
	ditfft2(in,out,1,CPLXDIM);
	print_complex(out);

	return 0;
}