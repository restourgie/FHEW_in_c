#include <stdio.h> 
#include <math.h> 
#include <stdlib.h> 
#include <complex.h> 

typedef complex float data_t; 
#define W(N,k) (cexp(-2.0f * M_PI * I * (float)(k) / (float)(N)))
 

float s(int n, int k) 
{
	if (n <= 4) 
		return 1.0f;

	int k4 = k % (n/4);

	if (k4 <= n/8) 
		return (s(n/4,k4) * cosf(2.0f * M_PI * (float)k4 / (float)n));
	
	return (s(n/4,k4) * sinf(2.0f * M_PI * (float)k4 / (float)n));
}

void tangentfft8(data_t *base, int TN, data_t *in, data_t *out, int stride, int N) 
{
	if(N == 1) 
	{
		if(in < base) 
			in += TN;
		out[0] = in[0];
	}
	else if(N == 2) 
	{
		data_t *i0 = in, *i1 = in + stride;
		if(i0 < base) 
			i0 += TN;
		if(i1 < base) 
			i1 += TN;
		out[0] = *i0 + *i1;
		out[N/2] = *i0 - *i1;
	}
	else if(N == 4) 
	{
		tangentfft8(base, TN, in, out, stride << 1, N >> 1);
		tangentfft8(base, TN, in+stride, out+2, stride << 1, N >> 1);

		data_t temp1 = out[0] + out[2];
		data_t temp2 = out[0] - out[2];
		out[0] = temp1;
		out[2] = temp2;
		temp1 = out[1] - I*out[3];
		temp2 = out[1] + I*out[3];
		out[1] = temp1;
		out[3] = temp2;
	}
	else
	{
		tangentfft8(base, TN, in, out, stride << 2, N >> 2);
		tangentfft8(base, TN, in+(stride*2), out+2*N/8, stride << 3, N >> 3);
		tangentfft8(base, TN, in-(stride*2), out+3*N/8, stride << 3, N >> 3);
		tangentfft8(base, TN, in+(stride), out+4*N/8, stride << 2, N >> 2);
		tangentfft8(base, TN, in-(stride), out+6*N/8, stride << 2, N >> 2);
		int k;
		for(k=0;k<N/8;k++) 
		{
			float s4 = s(N/4,k)/s(N,k);
			float s4_n8 = s(N/4,k+N/8)/s(N,k+N/8);

			float s2 = s(N/2,k)/s(N,k);
			float s2_n8 = s(N/2,k+N/8)/s(N,k+N/8);
			data_t w0 = W(N,k)*s4;
			data_t w1 = W(N,k+N/8)*s4_n8;
			data_t w2 = W(N,2*k)*s(N/8,k)/s(N/2,k);

			data_t zk_p = w0 * out[k+4*N/8];
			data_t zk_n = conj(w0) * out[k+6*N/8];
			data_t zk2_p = w1 * out[k+5*N/8];
			data_t zk2_n = conj(w1) * out[k+7*N/8];
			data_t uk = out[k] * s4;
			data_t uk2 = out[k+N/8] * s4_n8;
			data_t yk_p = w2 * out[k+2*N/8];
			data_t yk_n = conj(w2) * out[k+3*N/8];

			data_t y0 = (yk_p + yk_n)*s2;
			data_t y1 = (yk_p - yk_n)*I*s2_n2;

			out[k] = uk + y0 + (zk_p + zk_n);
			out[k+4*N/8] = uk + y0 - (zk_p + zk_n);
			out[k+2*N/8] = uk - y0 - I*(zk_p - zk_n);
			out[k+6*N/8] = uk - y0 + I*(zk_p - zk_n);
			out[k+1*N/8] = uk2 - y1 + (zk2_p + zk2_n);
			out[k+3*N/8] = uk2 + y1 - I*(zk2_p - zk2_n);
			out[k+5*N/8] = uk2 - y1 - (zk2_p + zk2_n);
			out[k+7*N/8] = uk2 + y1 + I*(zk2_p - zk2_n);
		}
	}

}

void tangentfft4(data_t *base, int TN, data_t *in, data_t *out, int stride, int N)
{
	if(N == 1) 
	{
		if(in < base) 
			in += TN;
		out[0] = in[0];
	}
	else if(N == 2)
	{
		data_t *i0 = in, *i1 = in + stride;
		if(i0 < base) 
			i0 += TN;
		if(i1 < base) 
			i1 += TN;
		out[0] = *i0 + *i1;
		out[N/2] = *i0 - *i1;
	}
	else
	{
		tangentfft4(base, TN, in, out, stride << 1, N >> 1);
		tangentfft8(base, TN, in+stride, out+N/2, stride << 2, N >> 2);
		tangentfft8(base, TN, in-stride, out+3*N/4, stride << 2, N >> 2);

		{
			data_t Uk = out[0];
			data_t Zk = out[0+N/2];
			data_t Uk2 = out[0+N/4];
			data_t Zdk = out[0+3*N/4];
			out[0] = Uk + (Zk + Zdk);
			out[0+N/2] = Uk - (Zk + Zdk);
			out[0+N/4] = Uk2 - I*(Zk - Zdk);
			out[0+3*N/4] = Uk2 + I*(Zk - Zdk);
		}
		int k;
		for(k=1;k<N/4;k++) 
		{
			data_t Uk = out[k];
			data_t Zk = out[k+N/2];
			data_t Uk2 = out[k+N/4];
			data_t Zdk = out[k+3*N/4];
			data_t w = W(N,k)*s(N/4,k);
			out[k] = Uk + (w*Zk + conj(w)*Zdk);
			out[k+N/2] = Uk - (w*Zk + conj(w)*Zdk);
			out[k+N/4] = Uk2 - I*(w*Zk - conj(w)*Zdk);
			out[k+3*N/4] = Uk2 + I*(w*Zk - conj(w)*Zdk);
		}
	}
}

int main(){

	data_t test;
	int N = 8;
	test[0] = 1 + I*0;
	test[1] = 1 + I*0;
	test[2] = 1 + I*0;
	test[3] = 1 + I*0;
	test[4] = 0 + I*0;


	return 0;
}