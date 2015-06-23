#include <stdio.h>
#include <math.h>
#include <complex.h>
 
double PI;
typedef double complex cplx;
 
void _fft(cplx buf[], cplx out[], int n, int step)
{
	if (step < n) {
		_fft(out, buf, n, step * 2);
		_fft(out + step, buf + step, n, step * 2);
 
		for (int i = 0; i < n; i += 2 * step) {
			cplx t = cexp(-I * PI * i / n) * out[i + step];
			buf[i / 2]     = out[i] + t;
			buf[(i + n)/2] = out[i] - t;
		}
	}
}
 
void fft(cplx buf[], int n)
{
	cplx out[n];
	for (int i = 0; i < n; i++) out[i] = buf[i];
 
	_fft(buf, out, n, 1);
}
 
 
void show(const char * s, cplx buf[]) {
	printf("%s\n", s);
	for (int i = 0; i < 8; i++)
		if (!cimag(buf[i]))
			printf("index i: %d (%g)\n",i, creal(buf[i]));
		else
			printf("index i: %d (%g, %g)\n",i, creal(buf[i]), cimag(buf[i]));
}
 
int main()
{
	PI = atan2(1, 1) * 4;
	cplx buf[]  = {2,4,6,2,0,0,0,0};
	cplx buf2[] = {1,2,8,5,0,0,0,0};

 
	show("Data: ", buf);
	fft(buf, 8);
	show("\nFFT : ", buf);
 
	return 0;
}