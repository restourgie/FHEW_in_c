#include <stdio.h>
#include <math.h>
#include <complex.h>
 
#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#define N 8

double complex we[N/2];

void wait_on_enter()
{
    printf("Press ENTER to continue\n");
    getchar();
}

void setup()
{
	for (int i = 0; i < (N / 2); i++){
	    we[i] = cos(2.0 * M_PI * i / N) + I*sin(2.0 * M_PI * i / N); ;
	}
}

// void BitInvert(double complex a[], int n)
// {
//     int i, j, mv = n/2;
//     int k, rev = 0;
//     double complex b;
//     for (i = 1; i < n; i++) // run tru all the indexes from 1 to n
//     {
//         k = i;
//         mv = n / 2;
//         rev = 0;
//         while (k > 0) // invert the actual index
//         {
//             if ((k % 2) > 0)
//             	rev = rev + mv;
//             k = k / 2;
//             mv = mv / 2; 
//         }
//         {  // switch the actual sample and the bitinverted one
//             if (i < rev)
//             {
//                 b = a[rev];
//                 a[rev] = a[i];
//                 a[i] = b;
//             }
//         }
//     }
// }

void BitInvert(double complex data[]){
  int i,mv,k,rev;
  double complex temp;

  for(i = 1; i<(N);i++){//run through all index 1 to N
    k = i;
    mv = (N)/2;
    rev = 0;
    while(k > 0){//Invert the index
      if ((k % 2) > 0)
              rev = rev + mv;
            k = k / 2;
            mv = mv / 2;
    }
    if(i < rev){
      temp = data[rev];
      data[rev] = data[i];
      data[i] = temp;
    }
  }
}

void CalcFFT(double complex data[], int n){
  BitInvert(data);
  
  int i,k,m;
  double complex h;
  k = 1;
  while(k <= n/2){
  	m =0;

  	while(m <= (n-2*k)){
  		for(i = m; i<m+k;i++){
  			// printf("\ni:%d m:%d k:%d\n",i,m,k);
  			h = data[i+k] * we[((i-m)*N / k/ 2)];
  			data[i+k] = data[i] - h;
  			data[i] += h;
  		}
  		m += 2 * k;
  	}
  	k = k*2;
  }
  
}

//Ring_FFT => complex_double[513] => double[513][2]
//Ring_ModQ => ZmodQ[1024] => int32_t[1024] 
void FFTforward(int val[],double complex data[]) {
    // double complex data[N];
    for(int k=0;k<N;++k){
      data[k] = val[k] + 0.0*I;
    }
    CalcFFT(data,N);
}

void FFTbackward(int val[],double complex data[]){
	CalcFFT(data,-1);
	for(int k=0;k<N;++k){
		val[k] = creal(data[k])/N;
	}
}
 
int main()
{	
	setup();
	int buf[] = {2, 4, 6, 2, 0, 0, 0, 0};
	double complex data[N];
	int buf2[] = {1,2,8,5,0,0,0,0};
	double complex data2[N];
 	FFTforward(buf,data);

 	for(int i=0;i<N;++i){
 		printf("index i: %d (%f,%f)\n",i,creal(data[i]),cimag(data[i]));
 	}

 	// FFTforward(buf2,data2);
 	// for(int i=0;i<N;++i){
 	// 	data[i] = data[i] * data2[i];
 	// }
 // 	for(int i=0;i<N/2;++i){
 // 		data[2*i] = 0.0 + I*0.0;
 // 		data[i*2+1] = data[i*2+1] * data2[i*2+1];	
	// }
	// data[N-1] = 0.0 + I*0.0;
	// FFTbackward(buf,data);
	// for(int i=0;i<N;++i){
 // 		printf("index i:%d value: %d\n",i,buf[i]);
	// }

	
 
	return 0;
}