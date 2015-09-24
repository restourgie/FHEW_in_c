#include <stdio.h>
#include <math.h>
#include "../support.h"
#include "split_radix_non_rec.h"
#include "split_radix_fast.h"

#define calc_cos(N,k) (cos(2.0 * M_PI * (double)k / (double) N))
#define calc_sin(N,k) (sin(2.0 * M_PI * (double)k / (double) N))

void new_twist(cplx *cplx_x,int n,int m,int lo)
{
  double a,b,c,d;
  int j = 1;
  for (int i = lo+1; i < lo+m; ++i)
  {	
  	a = cplx_x->real[i];
  	b = cplx_x->imag[i];	
  	c = calc_cos(n,j);
  	d = calc_sin(n,j);
  	 
    //(a + ib) * (c + id) = (ac - bd) + i(ad+bc)
  	cplx_x->real[i] = (a*c) - (b*d);
  	cplx_x->imag[i] = (a*d) + (b*c);
    ++j;
  }
}

void new_untwist(cplx *cplx_x,int n,int m,int lo)
{
  double a,b,c,d;
  int j = 1;
  for (int i = lo+1; i < lo+m; ++i)
  {
    a = cplx_x->real[i];
  	b = cplx_x->imag[i];	
  	c = calc_cos(n,j);
  	d = -calc_sin(n,j);
  	 
    //(a + ib) * (c + id) = (ac - bd) + i(ad+bc)
  	cplx_x->real[i] = (a*c) - (b*d);
  	cplx_x->imag[i] = (a*d) + (b*c);
    ++j;
  }
}

/******************************************************************
*
* SPLIT RADIX SPECIFIC FOR 4 SIZE
*
******************************************************************/
void split_radix_4(cplx *x,int lo)
{
	double temp_real,temp_im;
	int m = 2;
	//Go from (x^4 -1) to (x^2 -1) and (x^2 +1)
	for(int i=lo; i < lo+m;++i)
	{
	  temp_real = x->real[i];
	  temp_im = x->imag[i];

	  x->real[i] = temp_real + x->real[i+m];
	  x->imag[i] = temp_im + x->imag[i+m];

	  x->real[i+m] = temp_real - x->real[i+m];
	  x->imag[i+m] = temp_im - x->imag[i+m];
	}
	//BRANCH TAKE THE LEFT 2 ARRAYPARTS
	temp_real = x->real[lo];
	temp_im = x->imag[lo];

	x->real[lo] = temp_real + x->real[lo+1];
	x->imag[lo] = temp_im + x->imag[lo+1];

	x->real[lo+1] = temp_real - x->real[lo+1];
	x->imag[lo+1] = temp_im - x->imag[lo+1];
	
	lo += m;
	//Go from (x^2 +1) to (x -i) and (x +i)
	temp_real = x->real[lo];
	temp_im = x->imag[lo];

	//(a + ib) + i(c + id) = (a + ib) + (-d + ic) = (a-d + i(b+c))
	x->real[lo] = temp_real - x->imag[lo+1];
	x->imag[lo] = temp_im + x->real[lo+1];

	//(a + ib) - i(c + id) = (a + ib) - (-d + ic) = (a+d + i(b-c))
	temp_real = temp_real + x->imag[lo+1];
	x->imag[lo+1] = temp_im - x->real[lo+1];
	x->real[lo+1] = temp_real;
}

/******************************************************************
*
* SPLIT RADIX SPECIFIC FOR 8 SIZE
*
******************************************************************/
void split_radix_8(cplx *x,int lo)
{
	double temp_real,temp_im;
	int m = 4;
	//Go from (x^8 -1) to (x^4 -1) and (x^4 +1)
	for(int i=lo; i < lo+m;++i)
	{
	  temp_real = x->real[i];
	  temp_im = x->imag[i];

	  x->real[i] = temp_real + x->real[i+m];
	  x->imag[i] = temp_im + x->imag[i+m];

	  x->real[i+m] = temp_real - x->real[i+m];
	  x->imag[i+m] = temp_im - x->imag[i+m];
	}
	//BRANCH TAKE THE LEFT 4 ARRAYPARTS
	split_radix_4(x,lo);
	lo += m;
	m = 2;
	//Go from (x^4 +1) to (x^2 -i) and (x^2 +i)
	for (int i = lo; i < lo+m; ++i)
	{
		temp_real = x->real[i];
		temp_im = x->imag[i];

		//(a + ib) + i(c + id) = (a + ib) + (-d + ic) = (a-d + i(b+c))
		x->real[i] = temp_real - x->imag[i+m];
		x->imag[i] = temp_im + x->real[i+m];

		//(a + ib) - i(c + id) = (a + ib) - (-d + ic) = (a+d + i(b-c))
		temp_real = temp_real + x->imag[i+m];
		x->imag[i+m] = temp_im - x->real[i+m];
		x->real[i+m] = temp_real;
	}
	//left branch of right branch
	new_twist(x,8,m,lo);
	temp_real = x->real[lo];
	temp_im = x->imag[lo];

	x->real[lo] = temp_real + x->real[lo+1];
	x->imag[lo] = temp_im + x->imag[lo+1];

	x->real[lo+1] = temp_real - x->real[lo+1];
	x->imag[lo+1] = temp_im - x->imag[lo+1];
	//Right branch of right branch
	new_untwist(x,8,m,lo+m);
	lo +=m;
	temp_real = x->real[lo];
	temp_im = x->imag[lo];

	x->real[lo] = temp_real + x->real[lo+1];
	x->imag[lo] = temp_im + x->imag[lo+1];

	x->real[lo+1] = temp_real - x->real[lo+1];
	x->imag[lo+1] = temp_im - x->imag[lo+1];
}

/******************************************************************
*
* SPLIT RADIX SPECIFIC FOR 16 SIZE
*
******************************************************************/
void split_radix_16(cplx *x,int lo)
{
	double temp_real,temp_im;
	int m = 8;
	//Go from (x^16 -1) to (x^8 -1) and (x^8 +1)
	for(int i=lo; i < lo+m;++i)
	{
	  temp_real = x->real[i];
	  temp_im = x->imag[i];

	  x->real[i] = temp_real + x->real[i+m];
	  x->imag[i] = temp_im + x->imag[i+m];

	  x->real[i+m] = temp_real - x->real[i+m];
	  x->imag[i+m] = temp_im - x->imag[i+m];
	}
	//BRANCH TAKE THE LEFT 8 ARRAYPARTS
	split_radix_8(x,lo);
	lo += m;
	m = 4;
	//Go from (x^8 +1) to (x^4 -i) and (x^4 +i)
	for (int i = lo; i < lo+m; ++i)
	{
		temp_real = x->real[i];
		temp_im = x->imag[i];

		//(a + ib) + i(c + id) = (a + ib) + (-d + ic) = (a-d + i(b+c))
		x->real[i] = temp_real - x->imag[i+m];
		x->imag[i] = temp_im + x->real[i+m];

		//(a + ib) - i(c + id) = (a + ib) - (-d + ic) = (a+d + i(b-c))
		temp_real = temp_real + x->imag[i+m];
		x->imag[i+m] = temp_im - x->real[i+m];
		x->real[i+m] = temp_real;
	}
	//left branch of right branch
	new_twist(x,16,m,lo);
	split_radix_4(x,lo);
	//Right branch of right branch
	new_untwist(x,16,m,lo+m);
	split_radix_4(x,lo+m);
}

/******************************************************************
*
* SPLIT RADIX SPECIFIC FOR 32 SIZE
*
******************************************************************/
void split_radix_32(cplx *x,int lo)
{
	double temp_real,temp_im;
	int m = 16;
	//Go from (x^32 -1) to (x^16 -1) and (x^16 +1)
	for(int i=lo; i < lo+m;++i)
	{
	  temp_real = x->real[i];
	  temp_im = x->imag[i];

	  x->real[i] = temp_real + x->real[i+m];
	  x->imag[i] = temp_im + x->imag[i+m];

	  x->real[i+m] = temp_real - x->real[i+m];
	  x->imag[i+m] = temp_im - x->imag[i+m];
	}
	//BRANCH TAKE THE LEFT 16 ARRAYPARTS
	split_radix_16(x,lo);
	lo += m;
	m = 8;
	//Go from (x^16 +1) to (x^8 -i) and (x^8 +i)
	for (int i = lo; i < lo+m; ++i)
	{
		temp_real = x->real[i];
		temp_im = x->imag[i];

		//(a + ib) + i(c + id) = (a + ib) + (-d + ic) = (a-d + i(b+c))
		x->real[i] = temp_real - x->imag[i+m];
		x->imag[i] = temp_im + x->real[i+m];

		//(a + ib) - i(c + id) = (a + ib) - (-d + ic) = (a+d + i(b-c))
		temp_real = temp_real + x->imag[i+m];
		x->imag[i+m] = temp_im - x->real[i+m];
		x->real[i+m] = temp_real;
	}
	//left branch of right branch
	new_twist(x,32,m,lo);
	split_radix_8(x,lo);
	//Right branch of right branch
	new_untwist(x,32,m,lo+m);
	split_radix_8(x,lo+m);
}

/******************************************************************
*
* SPLIT RADIX SPECIFIC FOR 64 SIZE
*
******************************************************************/
void split_radix_64(cplx *x,int lo)
{
	double temp_real,temp_im;
	int m = 32;
	//Go from (x^64 -1) to (x^32 -1) and (x^32 +1)
	for(int i=lo; i < lo+m;++i)
	{
	  temp_real = x->real[i];
	  temp_im = x->imag[i];

	  x->real[i] = temp_real + x->real[i+m];
	  x->imag[i] = temp_im + x->imag[i+m];

	  x->real[i+m] = temp_real - x->real[i+m];
	  x->imag[i+m] = temp_im - x->imag[i+m];
	}
	//BRANCH TAKE THE LEFT 32 ARRAYPARTS
	split_radix_32(x,lo);
	lo += m;
	m = 16;
	//Go from (x^32 +1) to (x^16 -i) and (x^16 +i)
	for (int i = lo; i < lo+m; ++i)
	{
		temp_real = x->real[i];
		temp_im = x->imag[i];

		//(a + ib) + i(c + id) = (a + ib) + (-d + ic) = (a-d + i(b+c))
		x->real[i] = temp_real - x->imag[i+m];
		x->imag[i] = temp_im + x->real[i+m];

		//(a + ib) - i(c + id) = (a + ib) - (-d + ic) = (a+d + i(b-c))
		temp_real = temp_real + x->imag[i+m];
		x->imag[i+m] = temp_im - x->real[i+m];
		x->real[i+m] = temp_real;
	}
	//left branch of right branch
	new_twist(x,64,m,lo);
	split_radix_16(x,lo);
	//Right branch of right branch
	new_untwist(x,64,m,lo+m);
	split_radix_16(x,lo+m);
}

/******************************************************************
*
* SPLIT RADIX SPECIFIC FOR 128 SIZE
*
******************************************************************/
void split_radix_128(cplx *x,int lo)
{
	double temp_real,temp_im;
	int m = 64;
	//Go from (x^128 -1) to (x^64 -1) and (x^64 +1)
	for(int i=lo; i < lo+m;++i)
	{
	  temp_real = x->real[i];
	  temp_im = x->imag[i];

	  x->real[i] = temp_real + x->real[i+m];
	  x->imag[i] = temp_im + x->imag[i+m];

	  x->real[i+m] = temp_real - x->real[i+m];
	  x->imag[i+m] = temp_im - x->imag[i+m];
	}
	//BRANCH TAKE THE LEFT 64 ARRAYPARTS
	split_radix_64(x,lo);
	lo += m;
	m = 32;
	//Go from (x^64 +1) to (x^32 -i) and (x^32 +i)
	for (int i = lo; i < lo+m; ++i)
	{
		temp_real = x->real[i];
		temp_im = x->imag[i];

		//(a + ib) + i(c + id) = (a + ib) + (-d + ic) = (a-d + i(b+c))
		x->real[i] = temp_real - x->imag[i+m];
		x->imag[i] = temp_im + x->real[i+m];

		//(a + ib) - i(c + id) = (a + ib) - (-d + ic) = (a+d + i(b-c))
		temp_real = temp_real + x->imag[i+m];
		x->imag[i+m] = temp_im - x->real[i+m];
		x->real[i+m] = temp_real;
	}
	//left branch of right branch
	new_twist(x,128,m,lo);
	split_radix_32(x,lo);
	//Right branch of right branch
	new_untwist(x,128,m,lo+m);
	split_radix_32(x,lo+m);
}

/******************************************************************
*
* SPLIT RADIX SPECIFIC FOR 256 SIZE
*
******************************************************************/
void split_radix_256(cplx *x)
{
	double temp_real,temp_im;
	int m = 128;
	//Go from (x^256 -1) to (x^128 -1) and (x^128 +1)
	for(int i=0; i < 256;++i)
	{
	  temp_real = x->real[i];
	  temp_im = x->imag[i];

	  x->real[i] = temp_real + x->real[i+m];
	  x->imag[i] = temp_im + x->imag[i+m];

	  x->real[i+m] = temp_real - x->real[i+m];
	  x->imag[i+m] = temp_im - x->imag[i+m];
	}
	//BRANCH TAKE THE LEFT 128 ARRAYPARTS
	split_radix_128(x,0);

	m = 64;
	//Go from (x^128 +1) to (x^64 -i) and (x^64 +i)
	for (int i = 128; i < 128+m; ++i)
	{
		temp_real = x->real[i];
		temp_im = x->imag[i];

		//(a + ib) + i(c + id) = (a + ib) + (-d + ic) = (a-d + i(b+c))
		x->real[i] = temp_real - x->imag[i+m];
		x->imag[i] = temp_im + x->real[i+m];

		//(a + ib) - i(c + id) = (a + ib) - (-d + ic) = (a+d + i(b-c))
		temp_real = temp_real + x->imag[i+m];
		x->imag[i+m] = temp_im - x->real[i+m];
		x->real[i+m] = temp_real;
	}
	//left branch of right branch
	new_twist(x,256,m,128);
	split_radix_64(x,128);
	//Right branch of right branch
	new_untwist(x,256,m,128+m);
	split_radix_64(x,128+m);
}

/******************************************************************
*
* SPLIT RADIX SPECIFIC FOR 512 SIZE
*
******************************************************************/
void split_radix_512(cplx *x)
{
	double temp_real,temp_im;
	int m = 256;
	//Go from (x^512 -1) to (x^256 -1) and (x^256 +1)
	for(int i=0; i < 512;++i)
	{
	  temp_real = x->real[i];
	  temp_im = x->imag[i];

	  x->real[i] = temp_real + x->real[i+m];
	  x->imag[i] = temp_im + x->imag[i+m];

	  x->real[i+m] = temp_real - x->real[i+m];
	  x->imag[i+m] = temp_im - x->imag[i+m];
	}
	//BRANCH TAKE THE LEFT 256 ARRAYPARTS
	split_radix_256(x);

	m = 128;
	//Go from (x^256 +1) to (x^128 -i) and (x^128 +i)
	for (int i = 256; i < 256+m; ++i)
	{
		temp_real = x->real[i];
		temp_im = x->imag[i];

		//(a + ib) + i(c + id) = (a + ib) + (-d + ic) = (a-d + i(b+c))
		x->real[i] = temp_real - x->imag[i+m];
		x->imag[i] = temp_im + x->real[i+m];

		//(a + ib) - i(c + id) = (a + ib) - (-d + ic) = (a+d + i(b-c))
		temp_real = temp_real + x->imag[i+m];
		x->imag[i+m] = temp_im - x->real[i+m];
		x->real[i+m] = temp_real;
	}
	//left branch of right branch
	new_twist(x,512,m,256);
	split_radix_128(x,256);
	//Right branch of right branch
	new_untwist(x,512,m,256+m);
	split_radix_128(x,385);
}
