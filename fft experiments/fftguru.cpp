#include <cstdio> 
#include <cstdlib> // required by DevC++
#include <iostream> // for I/O and calculation
#include <cmath>

using namespace std; 

//****************************** fft_recur ****************************
void fft_recur(long N,double complex data[]) {
	long h, i1, j = 0, k, i2 = N/2, l1, l2 = 1;
	double c = -1.0, s = 0.0, t1, t2, u1, u2;
	double complex temp;
	for (h = 0; h < N-2; h++) // ***** bit reverse starts here ***
	{ 
		if (h < j) 
		{
			temp = data[h];
			data[h] = data[j];
			data[j] = temp;
		}
		k = i2;
		while (k <= j) 
		{
	 		j = j - k; 
	 		k = k/2;
	 	}
		j = j + k;
	} //****** bit reverse done ******

	for (k = 1; k < N; k = k*2) // ****** butterfly calculations start here
	{ 
 		l1 = l2; 
 		l2 = l2*2;
		u1 = 1.0; 
		u2 = 0.0;
		for (j = 0; j < l1; j++) 
		{
			for (h = j; h < N; h = h + l2) 
			{
				i1 = h + l1;
				t2 = (r[i1] - i[i1])*u2 ;
				t1 = t2 + r[i1]*(u1 - u2) ;
				t2 = t2 + i[i1]*(u1 + u2) ;
				r[i1] = r[h] - t1;
				i[i1] = i[h] - t2;
				r[h] = r[h] + t1;
				i[h] = i[h] + t2;
			} // end for over h
			t1 = u1 * c - u2 * s;
			u2 = u1 * s + u2 * c;
			u1 = t1;
 		} // end for over j
	 s = - sqrt((1.0 - c) / 2.0);
	 c = sqrt((1.0 + c) / 2.0);
	} // end for over k
} // end fft

int main (int nNumberofArgs, char* pszArgs[ ] ){ // arguments needed for Dev C++ I/O
	const long N = 32; double r[N] = {0.}, i[N] = {0.} ; // REPLACE WITH YOUR NUMBERS
	long n; 
	double twopi = 6.2831853071795865, t;
	for (n = 0; n < N; n++) { 	// generate saw tooth test data - DELETE IF USING OTHER DATA
		t = twopi*n/N; 			//unconventional way of doing things - see Chapter 6
		r[n]=1.1+sin(t)+.5*sin(2.*t)+(1./3.)*sin(3.*t)+.25*sin(4.*t)+.2*sin(5.*t)+(1./6.)*sin(6.*t)+(1./7.)*sin(7.*t);
	}
	fft_recur(N, r, i);
	cout<<"\nfft_recur results for N = " << N;
	cout<<"\n n r[n] i[n] amplitude phase(radians)\n";
	double amp, phase;
	for (n = 0; n < N; n++) {
		 r[n] = r[n]/N ; i[n] = i[n]/N ; // scale outputs
		 amp = sqrt( r[n]*r[n] + i[n]*i[n] ) ;
		 if (r[n]==0 && i[n]==0) 
		 	phase = 0;
		 else 
		 	phase = atan2(i[n],r[n]);
		 printf("%2d\t%9.4f\t%9.4f\t%9.4f\t%9.4f\n", n, r[n], i[n], amp, phase);
	} //end for loop over n
	system ("PAUSE");
	return 0;
} // end main