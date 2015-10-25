#include <complex.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <immintrin.h>
#include "../mul.h"
#include "negacyclic.h"


double ***wortel;

void print_complx(const cplx_ptr *x,int N){
  for (int i = 0; i < N; ++i)
    printf("cplxpoly[%d] = %f + i * %f\n",i,x->real[i],x->imag[i]);
  printf("\n");
}

void print256_num(__m256d var) 
{
    double *v64val = (double*) &var;
    printf("%f %f %f %f\n", v64val[0], v64val[1],v64val[2],v64val[3]);
}

void init_wortel(int n,int lo,int level,double complex root)
{	
	if(n > 1){
		int m = n/2;
		wortel[0][level][lo/n] = creal(root);
		wortel[1][level][lo/n] = cimag(root);
		wortel[2][level][lo/n] = -cimag(root);
		--level;
		init_wortel(m,lo,level,csqrt(root));
		init_wortel(m,lo+m,level,csqrt(-root));
	}
}


void init_negacyc()
{	
	wortel = malloc(3*sizeof(**wortel));
	int loga = log2(CPLXDIM);
	wortel[0] = malloc(loga*sizeof(*wortel));
	wortel[1] = malloc(loga*sizeof(*wortel));
	wortel[2] = malloc(loga*sizeof(*wortel));
	int j = CPLXDIM/2;
	for (int i = 0; i < loga; ++i)
	{
		wortel[0][i] = malloc(j*sizeof(double));
		wortel[1][i] = malloc(j*sizeof(double));
		wortel[2][i] = malloc(j*sizeof(double));
		// printf("i = %d j = %d\n",i,j);
		j = j>>1;
	}
	init_wortel(CPLXDIM,0,(loga-1),csqrt(I));
}

/******************************************************************
*
* SMART COMPLEX MULTIPLICATION
*
******************************************************************/
void inverse_phi(cplx_ptr *x)
{	
	__m256d real_x,imag_x,real_y,imag_y,imag_twid,real_twid,temp_real,temp_imag,sub_real,sub_imag;
  	int lo;
	// printf("hi\n");
	for(lo = 0;lo < CPLXDIM;lo +=8)
  	{
		real_x = _mm256_load_pd(x->real+lo);
		imag_x = _mm256_load_pd(x->imag+lo);
		real_y = _mm256_load_pd(x->real+lo+4);
		imag_y = _mm256_load_pd(x->imag+lo+4);
		//LAYER 2
		//we need to collect all the left and the right parts
		//real part
		sub_real = _mm256_unpacklo_pd(real_x,real_y);
		sub_imag = _mm256_unpackhi_pd(real_x,real_y);

		real_y = _mm256_add_pd(sub_real,sub_imag);
		real_x = _mm256_sub_pd(sub_real,sub_imag);

		//imag part
		sub_real = _mm256_unpacklo_pd(imag_x,imag_y);
		sub_imag = _mm256_unpackhi_pd(imag_x,imag_y);

		imag_y = _mm256_add_pd(sub_real,sub_imag);
		imag_x = _mm256_sub_pd(sub_real,sub_imag);

		//Now do complex inverse computation with roots of unity
		real_twid = _mm256_setr_pd(wortel[0][0][lo/2],wortel[0][0][(lo+4)/2],wortel[0][0][(lo+2)/2],wortel[0][0][(lo+6)/2]);
    	imag_twid = _mm256_setr_pd(wortel[2][0][lo/2],wortel[2][0][(lo+4)/2],wortel[2][0][(lo+2)/2],wortel[2][0][(lo+6)/2]);
	    temp_real = _mm256_mul_pd(imag_x,imag_twid);
    	temp_imag = _mm256_mul_pd(imag_x,real_twid);

		//TEMP_imag = ad + bc
		imag_x = _mm256_fmadd_pd(real_x,imag_twid,temp_imag);
		//TEMP_real = ac - bd
		real_x = _mm256_fmsub_pd(real_x,real_twid,temp_real);

		//LAYER 4
		//real part
		sub_real = _mm256_permute2f128_pd(real_y,real_x,0x20);
		sub_imag = _mm256_permute2f128_pd(real_y,real_x,0x31);

		real_y = _mm256_add_pd(sub_real,sub_imag);
		real_x = _mm256_sub_pd(sub_real,sub_imag);
		//imag part
		sub_real = _mm256_permute2f128_pd(imag_y,imag_x,0x20);
		sub_imag = _mm256_permute2f128_pd(imag_y,imag_x,0x31);

		imag_y = _mm256_add_pd(sub_real,sub_imag);
		imag_x = _mm256_sub_pd(sub_real,sub_imag);
		//mult
		real_twid = _mm256_setr_pd(wortel[0][1][lo/4],wortel[0][1][(lo+4)/4],wortel[0][1][lo/4],wortel[0][1][(lo+4)/4]);
	    imag_twid = _mm256_setr_pd(wortel[2][1][lo/4],wortel[2][1][(lo+4)/4],wortel[2][1][lo/4],wortel[2][1][(lo+4)/4]);
	    temp_real = _mm256_mul_pd(imag_x,imag_twid);
    	temp_imag = _mm256_mul_pd(imag_x,real_twid);
		//TEMP_imag = ad + bc
		imag_x = _mm256_fmadd_pd(real_x,imag_twid,temp_imag);
		//TEMP_real = ac - bd
		real_x = _mm256_fmsub_pd(real_x,real_twid,temp_real);
		
		//LAYER 8
		//real part
		temp_real = _mm256_unpacklo_pd(real_y,real_x);
		temp_imag = _mm256_unpackhi_pd(real_y,real_x);

		real_x = _mm256_add_pd(temp_real,temp_imag);
		real_y = _mm256_sub_pd(temp_real,temp_imag);
		//imag part
		sub_real  = _mm256_unpacklo_pd(imag_y,imag_x);
		sub_imag  = _mm256_unpackhi_pd(imag_y,imag_x);

		imag_x = _mm256_add_pd(sub_real,sub_imag);
		imag_y = _mm256_sub_pd(sub_real,sub_imag);
		//mult
	    real_twid = _mm256_set1_pd(wortel[0][2][lo/8]);
    	imag_twid = _mm256_set1_pd(wortel[2][2][lo/8]);
	    temp_real = _mm256_mul_pd(imag_y,imag_twid);
    	temp_imag = _mm256_mul_pd(imag_y,real_twid);

		//TEMP_imag = ad + bc
		imag_y = _mm256_fmadd_pd(real_y,imag_twid,temp_imag);
		//TEMP_real = ac - bd
		real_y = _mm256_fmsub_pd(real_y,real_twid,temp_real);

		real_x =_mm256_permute4x64_pd(real_x,0xd8);
		real_y =_mm256_permute4x64_pd(real_y,0xd8);
		imag_x =_mm256_permute4x64_pd(imag_x,0xd8);
		imag_y =_mm256_permute4x64_pd(imag_y,0xd8);

		_mm256_store_pd(x->real+lo,real_x);
		_mm256_store_pd(x->imag+lo,imag_x);
		_mm256_store_pd(x->real+lo+4,real_y);
		_mm256_store_pd(x->imag+lo+4,imag_y);
  	}
  
	__m256d v0_r,v8_r,v16_r,v24_r,v32_r,v40_r,v48_r,v56_r,v0_i,v8_i,v16_i,v24_i,v32_i,v40_i,v48_i,v56_i;
	__m256d real_twiddle,imag_twiddle;
	//Make sure n = 64
	for(lo = 0;lo < CPLXDIM;lo +=64)
	{	

		for(int offset = 0; offset < 8; offset+=4)
		{	
			//LOAD A LOT of STUF
	  		//FIRST ROUND WE NEED 0..3,8..11,16..19,24..27 LEFT SIDE (Real and Imag)
	  		//32..35,40..43,48..51,56..59 (Real and Imag)
		  	//REAL PART
		  	v0_r  = _mm256_load_pd(x->real+lo+offset);
		  	v8_r  = _mm256_load_pd(x->real+lo+offset+8);
		  	v16_r = _mm256_load_pd(x->real+lo+offset+16);
		  	v24_r = _mm256_load_pd(x->real+lo+offset+24);
		  	v32_r = _mm256_load_pd(x->real+lo+offset+32);
		  	v40_r = _mm256_load_pd(x->real+lo+offset+40);
		  	v48_r = _mm256_load_pd(x->real+lo+offset+48);
		  	v56_r = _mm256_load_pd(x->real+lo+offset+56);
		  	//IMAG PART
			v0_i  = _mm256_load_pd(x->imag+lo+offset);
		  	v8_i  = _mm256_load_pd(x->imag+lo+offset+8);
		  	v16_i = _mm256_load_pd(x->imag+lo+offset+16);
		  	v24_i = _mm256_load_pd(x->imag+lo+offset+24);
		  	v32_i = _mm256_load_pd(x->imag+lo+offset+32);
		  	v40_i = _mm256_load_pd(x->imag+lo+offset+40);
		  	v48_i = _mm256_load_pd(x->imag+lo+offset+48);
		  	v56_i = _mm256_load_pd(x->imag+lo+offset+56);

		  	//WE START WITH LAYER 2
  		  	real_twiddle = _mm256_set1_pd(wortel[0][3][(lo)/16]);
			imag_twiddle = _mm256_set1_pd(wortel[2][3][(lo)/16]);
		  	//TWIDDLE 0 - 8
		    temp_real = _mm256_sub_pd(v0_r,v8_r);
	  		temp_imag = _mm256_sub_pd(v0_i,v8_i);

		    v0_r = _mm256_add_pd(v0_r,v8_r);
	  		v0_i = _mm256_add_pd(v0_i,v8_i);

		    sub_real = _mm256_mul_pd(temp_imag,imag_twiddle);
		    temp_imag = _mm256_mul_pd(temp_imag,real_twiddle);

			v8_i = _mm256_fmadd_pd(temp_real,imag_twiddle,temp_imag);
		    v8_r = _mm256_fmsub_pd(temp_real,real_twiddle,sub_real);

  		  	real_twiddle = _mm256_set1_pd(wortel[0][3][(lo+16)/16]);
			imag_twiddle = _mm256_set1_pd(wortel[2][3][(lo+16)/16]);
		  	//TWIDDLE 16 - 24
		    temp_real = _mm256_sub_pd(v16_r,v24_r);
	  		temp_imag = _mm256_sub_pd(v16_i,v24_i);

		    v16_r = _mm256_add_pd(v16_r,v24_r);
	  		v16_i = _mm256_add_pd(v16_i,v24_i);

		    sub_real = _mm256_mul_pd(temp_imag,imag_twiddle);
		    temp_imag = _mm256_mul_pd(temp_imag,real_twiddle);

			v24_i = _mm256_fmadd_pd(temp_real,imag_twiddle,temp_imag);
		    v24_r = _mm256_fmsub_pd(temp_real,real_twiddle,sub_real);

  		  	real_twiddle = _mm256_set1_pd(wortel[0][3][(lo+32)/16]);
			imag_twiddle = _mm256_set1_pd(wortel[2][3][(lo+32)/16]);
		  	//TWIDDLE 32 - 40
		    temp_real = _mm256_sub_pd(v32_r,v40_r);
	  		temp_imag = _mm256_sub_pd(v32_i,v40_i);

		    v32_r = _mm256_add_pd(v32_r,v40_r);
	  		v32_i = _mm256_add_pd(v32_i,v40_i);

		    sub_real = _mm256_mul_pd(temp_imag,imag_twiddle);
		    temp_imag = _mm256_mul_pd(temp_imag,real_twiddle);

			v40_i = _mm256_fmadd_pd(temp_real,imag_twiddle,temp_imag);
		    v40_r = _mm256_fmsub_pd(temp_real,real_twiddle,sub_real);

  		  	
  		  	real_twiddle = _mm256_set1_pd(wortel[0][3][(lo+48)/16]);
			imag_twiddle = _mm256_set1_pd(wortel[2][3][(lo+48)/16]);
		  	//TWIDDLE 48 - 56
		    temp_real = _mm256_sub_pd(v48_r,v56_r);
	  		temp_imag = _mm256_sub_pd(v48_i,v56_i);

		    v48_r = _mm256_add_pd(v48_r,v56_r);
	  		v48_i = _mm256_add_pd(v48_i,v56_i);

		    sub_real = _mm256_mul_pd(temp_imag,imag_twiddle);
		    temp_imag = _mm256_mul_pd(temp_imag,real_twiddle);

			v56_i = _mm256_fmadd_pd(temp_real,imag_twiddle,temp_imag);
		    v56_r = _mm256_fmsub_pd(temp_real,real_twiddle,sub_real);

		    //START WITH LAYER 32
		    //LEFT SIDE
  		  	real_twiddle = _mm256_set1_pd(wortel[0][4][(lo)/32]);
			imag_twiddle = _mm256_set1_pd(wortel[2][4][(lo)/32]);
		  	//TWIDDLE 0 - 16
		    temp_real = _mm256_sub_pd(v0_r,v16_r);
	  		temp_imag = _mm256_sub_pd(v0_i,v16_i);

		    v0_r = _mm256_add_pd(v0_r,v16_r);
	  		v0_i = _mm256_add_pd(v0_i,v16_i);

		    sub_real = _mm256_mul_pd(temp_imag,imag_twiddle);
		    temp_imag = _mm256_mul_pd(temp_imag,real_twiddle);

			v16_i = _mm256_fmadd_pd(temp_real,imag_twiddle,temp_imag);
		    v16_r = _mm256_fmsub_pd(temp_real,real_twiddle,sub_real);
		  	//TWIDDLE 8 - 24
		    temp_real = _mm256_sub_pd(v8_r,v24_r);
	  		temp_imag = _mm256_sub_pd(v8_i,v24_i);

		    v8_r = _mm256_add_pd(v8_r,v24_r);
	  		v8_i = _mm256_add_pd(v8_i,v24_i);

		    sub_real = _mm256_mul_pd(temp_imag,imag_twiddle);
		    temp_imag = _mm256_mul_pd(temp_imag,real_twiddle);

			v24_i = _mm256_fmadd_pd(temp_real,imag_twiddle,temp_imag);
		    v24_r = _mm256_fmsub_pd(temp_real,real_twiddle,sub_real);

		    //RIGHT SIDE
  		  	real_twiddle = _mm256_set1_pd(wortel[0][4][(lo+32)/32]);
			imag_twiddle = _mm256_set1_pd(wortel[2][4][(lo+32)/32]);
		  	//TWIDDLE 32 - 48
		    temp_real = _mm256_sub_pd(v32_r,v48_r);
	  		temp_imag = _mm256_sub_pd(v32_i,v48_i);

		    v32_r = _mm256_add_pd(v32_r,v48_r);
	  		v32_i = _mm256_add_pd(v32_i,v48_i);

		    sub_real = _mm256_mul_pd(temp_imag,imag_twiddle);
		    temp_imag = _mm256_mul_pd(temp_imag,real_twiddle);

			v48_i = _mm256_fmadd_pd(temp_real,imag_twiddle,temp_imag);
		    v48_r = _mm256_fmsub_pd(temp_real,real_twiddle,sub_real);
		  	//TWIDDLE 40 - 56
		    temp_real = _mm256_sub_pd(v40_r,v56_r);
	  		temp_imag = _mm256_sub_pd(v40_i,v56_i);

		    v40_r = _mm256_add_pd(v40_r,v56_r);
	  		v40_i = _mm256_add_pd(v40_i,v56_i);

		    sub_real = _mm256_mul_pd(temp_imag,imag_twiddle);
		    temp_imag = _mm256_mul_pd(temp_imag,real_twiddle);

			v56_i = _mm256_fmadd_pd(temp_real,imag_twiddle,temp_imag);
		    v56_r = _mm256_fmsub_pd(temp_real,real_twiddle,sub_real);

		    //START WITH LAYER 64
  		  	real_twiddle = _mm256_set1_pd(wortel[0][5][(lo)/64]);
			imag_twiddle = _mm256_set1_pd(wortel[2][5][(lo)/64]);
		  	//TWIDDLE 0 - 32
		    temp_real = _mm256_sub_pd(v0_r,v32_r);
	  		temp_imag = _mm256_sub_pd(v0_i,v32_i);

		    v0_r = _mm256_add_pd(v0_r,v32_r);
	  		v0_i = _mm256_add_pd(v0_i,v32_i);

		    sub_real = _mm256_mul_pd(temp_imag,imag_twiddle);
		    temp_imag = _mm256_mul_pd(temp_imag,real_twiddle);

			v32_i = _mm256_fmadd_pd(temp_real,imag_twiddle,temp_imag);
		    v32_r = _mm256_fmsub_pd(temp_real,real_twiddle,sub_real);
		  	//TWIDDLE 8 - 40
		    temp_real = _mm256_sub_pd(v8_r,v40_r);
	  		temp_imag = _mm256_sub_pd(v8_i,v40_i);

		    v8_r = _mm256_add_pd(v8_r,v40_r);
	  		v8_i = _mm256_add_pd(v8_i,v40_i);

		    sub_real = _mm256_mul_pd(temp_imag,imag_twiddle);
		    temp_imag = _mm256_mul_pd(temp_imag,real_twiddle);

			v40_i = _mm256_fmadd_pd(temp_real,imag_twiddle,temp_imag);
		    v40_r = _mm256_fmsub_pd(temp_real,real_twiddle,sub_real);
		  	//TWIDDLE 16 - 48
		    temp_real = _mm256_sub_pd(v16_r,v48_r);
	  		temp_imag = _mm256_sub_pd(v16_i,v48_i);

		    v16_r = _mm256_add_pd(v16_r,v48_r);
	  		v16_i = _mm256_add_pd(v16_i,v48_i);

		    sub_real = _mm256_mul_pd(temp_imag,imag_twiddle);
		    temp_imag = _mm256_mul_pd(temp_imag,real_twiddle);

			v48_i = _mm256_fmadd_pd(temp_real,imag_twiddle,temp_imag);
		    v48_r = _mm256_fmsub_pd(temp_real,real_twiddle,sub_real);
		  	//TWIDDLE 24 - 56
		    temp_real = _mm256_sub_pd(v24_r,v56_r);
	  		temp_imag = _mm256_sub_pd(v24_i,v56_i);

		    v24_r = _mm256_add_pd(v24_r,v56_r);
	  		v24_i = _mm256_add_pd(v24_i,v56_i);

		    sub_real = _mm256_mul_pd(temp_imag,imag_twiddle);
		    temp_imag = _mm256_mul_pd(temp_imag,real_twiddle);

			v56_i = _mm256_fmadd_pd(temp_real,imag_twiddle,temp_imag);
		    v56_r = _mm256_fmsub_pd(temp_real,real_twiddle,sub_real);

			//STORE ALL RESULTS
			_mm256_store_pd(x->real+lo+offset,v0_r);
			_mm256_store_pd(x->real+lo+offset+8,v8_r);
			_mm256_store_pd(x->real+lo+offset+16,v16_r);
			_mm256_store_pd(x->real+lo+offset+24,v24_r);
			_mm256_store_pd(x->real+lo+offset+32,v32_r);
			_mm256_store_pd(x->real+lo+offset+40,v40_r);
			_mm256_store_pd(x->real+lo+offset+48,v48_r);
			_mm256_store_pd(x->real+lo+offset+56,v56_r);

			_mm256_store_pd(x->imag+lo+offset,v0_i);
			_mm256_store_pd(x->imag+lo+offset+8,v8_i);
			_mm256_store_pd(x->imag+lo+offset+16,v16_i);
			_mm256_store_pd(x->imag+lo+offset+24,v24_i);
			_mm256_store_pd(x->imag+lo+offset+32,v32_i);
			_mm256_store_pd(x->imag+lo+offset+40,v40_i);
			_mm256_store_pd(x->imag+lo+offset+48,v48_i);
			_mm256_store_pd(x->imag+lo+offset+56,v56_i);
		}
	}

	//init twiddle for 512
	real_twid = _mm256_set1_pd(wortel[0][8][0]);
	imag_twid = _mm256_set1_pd(wortel[2][8][0]);
	lo =0;
	//FUSING 512-256-128
  	for(int offset = 0; offset < 64; offset+=4)
  	{	
		//LOAD A LOT of STUF
	  	//FIRST ROUND WE NEED 0..3,64..67,128..131,192..195 LEFT SIDE (Real and Imag)
	  	//256..259,320..323,384..387,448..451 (Real and Imag)

	  	//REAL PART
	  	v0_r  = _mm256_load_pd(x->real+offset);
	  	v8_r  = _mm256_load_pd(x->real+offset+64);
	  	v16_r = _mm256_load_pd(x->real+offset+128);
	  	v24_r = _mm256_load_pd(x->real+offset+192);
	  	v32_r = _mm256_load_pd(x->real+offset+256);
	  	v40_r = _mm256_load_pd(x->real+offset+320);
	  	v48_r = _mm256_load_pd(x->real+offset+384);
	  	v56_r = _mm256_load_pd(x->real+offset+448);
	  	//IMAG PART
		v0_i  = _mm256_load_pd(x->imag+offset);
	  	v8_i  = _mm256_load_pd(x->imag+offset+64);
	  	v16_i = _mm256_load_pd(x->imag+offset+128);
	  	v24_i = _mm256_load_pd(x->imag+offset+192);
	  	v32_i = _mm256_load_pd(x->imag+offset+256);
	  	v40_i = _mm256_load_pd(x->imag+offset+320);
	  	v48_i = _mm256_load_pd(x->imag+offset+384);
	  	v56_i = _mm256_load_pd(x->imag+offset+448);

		//WE START WITH LAYER 128
		real_twiddle = _mm256_set1_pd(wortel[0][6][0]);
		imag_twiddle = _mm256_set1_pd(wortel[2][6][0]);
		//TWIDDLE 0 - 64
		temp_real = _mm256_sub_pd(v0_r,v8_r);
		temp_imag = _mm256_sub_pd(v0_i,v8_i);

		v0_r = _mm256_add_pd(v0_r,v8_r);
		v0_i = _mm256_add_pd(v0_i,v8_i);

		sub_real = _mm256_mul_pd(temp_imag,imag_twiddle);
		temp_imag = _mm256_mul_pd(temp_imag,real_twiddle);

		v8_i = _mm256_fmadd_pd(temp_real,imag_twiddle,temp_imag);
		v8_r = _mm256_fmsub_pd(temp_real,real_twiddle,sub_real);

		real_twiddle = _mm256_set1_pd(wortel[0][6][1]);
		imag_twiddle = _mm256_set1_pd(wortel[2][6][1]);
		//TWIDDLE 128 - 192
		temp_real = _mm256_sub_pd(v16_r,v24_r);
		temp_imag = _mm256_sub_pd(v16_i,v24_i);

		v16_r = _mm256_add_pd(v16_r,v24_r);
		v16_i = _mm256_add_pd(v16_i,v24_i);

		sub_real = _mm256_mul_pd(temp_imag,imag_twiddle);
		temp_imag = _mm256_mul_pd(temp_imag,real_twiddle);

		v24_i = _mm256_fmadd_pd(temp_real,imag_twiddle,temp_imag);
		v24_r = _mm256_fmsub_pd(temp_real,real_twiddle,sub_real);

		real_twiddle = _mm256_set1_pd(wortel[0][6][2]);
		imag_twiddle = _mm256_set1_pd(wortel[2][6][2]);
		//TWIDDLE 256 - 320
		temp_real = _mm256_sub_pd(v32_r,v40_r);
		temp_imag = _mm256_sub_pd(v32_i,v40_i);

		v32_r = _mm256_add_pd(v32_r,v40_r);
		v32_i = _mm256_add_pd(v32_i,v40_i);

		sub_real = _mm256_mul_pd(temp_imag,imag_twiddle);
		temp_imag = _mm256_mul_pd(temp_imag,real_twiddle);

		v40_i = _mm256_fmadd_pd(temp_real,imag_twiddle,temp_imag);
		v40_r = _mm256_fmsub_pd(temp_real,real_twiddle,sub_real);


		real_twiddle = _mm256_set1_pd(wortel[0][6][3]);
		imag_twiddle = _mm256_set1_pd(wortel[2][6][3]);
		//TWIDDLE 384 - 448
		temp_real = _mm256_sub_pd(v48_r,v56_r);
		temp_imag = _mm256_sub_pd(v48_i,v56_i);

		v48_r = _mm256_add_pd(v48_r,v56_r);
		v48_i = _mm256_add_pd(v48_i,v56_i);

		sub_real = _mm256_mul_pd(temp_imag,imag_twiddle);
		temp_imag = _mm256_mul_pd(temp_imag,real_twiddle);

		v56_i = _mm256_fmadd_pd(temp_real,imag_twiddle,temp_imag);
		v56_r = _mm256_fmsub_pd(temp_real,real_twiddle,sub_real);

		//START WITH LAYER 256
		//LEFT SIDE
		real_twiddle = _mm256_set1_pd(wortel[0][7][0]);
		imag_twiddle = _mm256_set1_pd(wortel[2][7][0]);
		//TWIDDLE 0 - 128
		temp_real = _mm256_sub_pd(v0_r,v16_r);
		temp_imag = _mm256_sub_pd(v0_i,v16_i);

		v0_r = _mm256_add_pd(v0_r,v16_r);
		v0_i = _mm256_add_pd(v0_i,v16_i);

		sub_real = _mm256_mul_pd(temp_imag,imag_twiddle);
		temp_imag = _mm256_mul_pd(temp_imag,real_twiddle);

		v16_i = _mm256_fmadd_pd(temp_real,imag_twiddle,temp_imag);
		v16_r = _mm256_fmsub_pd(temp_real,real_twiddle,sub_real);
		//TWIDDLE 64 - 192
		temp_real = _mm256_sub_pd(v8_r,v24_r);
		temp_imag = _mm256_sub_pd(v8_i,v24_i);

		v8_r = _mm256_add_pd(v8_r,v24_r);
		v8_i = _mm256_add_pd(v8_i,v24_i);

		sub_real = _mm256_mul_pd(temp_imag,imag_twiddle);
		temp_imag = _mm256_mul_pd(temp_imag,real_twiddle);

		v24_i = _mm256_fmadd_pd(temp_real,imag_twiddle,temp_imag);
		v24_r = _mm256_fmsub_pd(temp_real,real_twiddle,sub_real);

		//RIGHT SIDE
		real_twiddle = _mm256_set1_pd(wortel[0][7][1]);
		imag_twiddle = _mm256_set1_pd(wortel[2][7][1]);
		//TWIDDLE 256 - 384
		temp_real = _mm256_sub_pd(v32_r,v48_r);
		temp_imag = _mm256_sub_pd(v32_i,v48_i);

		v32_r = _mm256_add_pd(v32_r,v48_r);
		v32_i = _mm256_add_pd(v32_i,v48_i);

		sub_real = _mm256_mul_pd(temp_imag,imag_twiddle);
		temp_imag = _mm256_mul_pd(temp_imag,real_twiddle);

		v48_i = _mm256_fmadd_pd(temp_real,imag_twiddle,temp_imag);
		v48_r = _mm256_fmsub_pd(temp_real,real_twiddle,sub_real);
		//TWIDDLE 320 - 448
		temp_real = _mm256_sub_pd(v40_r,v56_r);
		temp_imag = _mm256_sub_pd(v40_i,v56_i);

		v40_r = _mm256_add_pd(v40_r,v56_r);
		v40_i = _mm256_add_pd(v40_i,v56_i);

		sub_real = _mm256_mul_pd(temp_imag,imag_twiddle);
		temp_imag = _mm256_mul_pd(temp_imag,real_twiddle);

		v56_i = _mm256_fmadd_pd(temp_real,imag_twiddle,temp_imag);
		v56_r = _mm256_fmsub_pd(temp_real,real_twiddle,sub_real);


		//START 512
		//TWIDDLE 0 - 256
		temp_real = _mm256_sub_pd(v0_r,v32_r);
		temp_imag = _mm256_sub_pd(v0_i,v32_i);

		v0_r = _mm256_add_pd(v0_r,v32_r);
		v0_i = _mm256_add_pd(v0_i,v32_i);

		sub_real = _mm256_mul_pd(temp_imag,imag_twid);
		temp_imag = _mm256_mul_pd(temp_imag,real_twid);

		v32_i = _mm256_fmadd_pd(temp_real,imag_twid,temp_imag);
		v32_r = _mm256_fmsub_pd(temp_real,real_twid,sub_real);
		//TWIDDLE 64 - 320
		temp_real = _mm256_sub_pd(v8_r,v40_r);
		temp_imag = _mm256_sub_pd(v8_i,v40_i);

		v8_r = _mm256_add_pd(v8_r,v40_r);
		v8_i = _mm256_add_pd(v8_i,v40_i);

		sub_real = _mm256_mul_pd(temp_imag,imag_twid);
		temp_imag = _mm256_mul_pd(temp_imag,real_twid);

		v40_i = _mm256_fmadd_pd(temp_real,imag_twid,temp_imag);
		v40_r = _mm256_fmsub_pd(temp_real,real_twid,sub_real);
		//TWIDDLE 128 - 384
		temp_real = _mm256_sub_pd(v16_r,v48_r);
		temp_imag = _mm256_sub_pd(v16_i,v48_i);

		v16_r = _mm256_add_pd(v16_r,v48_r);
		v16_i = _mm256_add_pd(v16_i,v48_i);

		sub_real = _mm256_mul_pd(temp_imag,imag_twid);
		temp_imag = _mm256_mul_pd(temp_imag,real_twid);

		v48_i = _mm256_fmadd_pd(temp_real,imag_twid,temp_imag);
		v48_r = _mm256_fmsub_pd(temp_real,real_twid,sub_real);
		//TWIDDLE 192 - 448
		temp_real = _mm256_sub_pd(v24_r,v56_r);
		temp_imag = _mm256_sub_pd(v24_i,v56_i);

		v24_r = _mm256_add_pd(v24_r,v56_r);
		v24_i = _mm256_add_pd(v24_i,v56_i);

		sub_real = _mm256_mul_pd(temp_imag,imag_twid);
		temp_imag = _mm256_mul_pd(temp_imag,real_twid);

		v56_i = _mm256_fmadd_pd(temp_real,imag_twid,temp_imag);
		v56_r = _mm256_fmsub_pd(temp_real,real_twid,sub_real);

	  		//STORE ALL RESULTS
		_mm256_store_pd(x->real+offset,v0_r);
		_mm256_store_pd(x->real+offset+64,v8_r);
		_mm256_store_pd(x->real+offset+128,v16_r);
		_mm256_store_pd(x->real+offset+192,v24_r);
		_mm256_store_pd(x->real+offset+256,v32_r);
		_mm256_store_pd(x->real+offset+320,v40_r);
		_mm256_store_pd(x->real+offset+384,v48_r);
		_mm256_store_pd(x->real+offset+448,v56_r);

		_mm256_store_pd(x->imag+offset,v0_i);
		_mm256_store_pd(x->imag+offset+64,v8_i);
		_mm256_store_pd(x->imag+offset+128,v16_i);
		_mm256_store_pd(x->imag+offset+192,v24_i);
		_mm256_store_pd(x->imag+offset+256,v32_i);
		_mm256_store_pd(x->imag+offset+320,v40_i);
		_mm256_store_pd(x->imag+offset+384,v48_i);
		_mm256_store_pd(x->imag+offset+448,v56_i);
	}

}

void iterative_phi(cplx_ptr *x)
{
  __m256d imag_twid,real_twid,temp_real,temp_imag,sub_real,sub_imag;
  __m256d real_twiddle,imag_twiddle;
  int lo =0;
  __m256d v0_r,v8_r,v16_r,v24_r,v32_r,v40_r,v48_r,v56_r,v0_i,v8_i,v16_i,v24_i,v32_i,v40_i,v48_i,v56_i;
  //init twiddle for 512
  real_twid = _mm256_set1_pd(wortel[0][8][0]);
  imag_twid = _mm256_set1_pd(wortel[1][8][0]);
  //FUSING 512-256-128
  for(int offset = 0; offset < 64; offset+=4)
  {	
	//LOAD A LOT of STUF
  	//FIRST ROUND WE NEED 0..3,64..67,128..131,192..195 LEFT SIDE (Real and Imag)
  	//256..259,320..323,384..387,448..451 (Real and Imag)

  	//REAL PART
  	v0_r  = _mm256_load_pd(x->real+offset);
  	v8_r  = _mm256_load_pd(x->real+offset+64);
  	v16_r = _mm256_load_pd(x->real+offset+128);
  	v24_r = _mm256_load_pd(x->real+offset+192);
  	v32_r = _mm256_load_pd(x->real+offset+256);
  	v40_r = _mm256_load_pd(x->real+offset+320);
  	v48_r = _mm256_load_pd(x->real+offset+384);
  	v56_r = _mm256_load_pd(x->real+offset+448);
  	//IMAG PART
	v0_i  = _mm256_load_pd(x->imag+offset);
  	v8_i  = _mm256_load_pd(x->imag+offset+64);
  	v16_i = _mm256_load_pd(x->imag+offset+128);
  	v24_i = _mm256_load_pd(x->imag+offset+192);
  	v32_i = _mm256_load_pd(x->imag+offset+256);
  	v40_i = _mm256_load_pd(x->imag+offset+320);
  	v48_i = _mm256_load_pd(x->imag+offset+384);
  	v56_i = _mm256_load_pd(x->imag+offset+448);
	  	
  	//START TWIDDLE
  	//WE ARE NOW IN 512 SO EVERYTHING BETWEEN 256 and 512 needs to be multiplied by root of unity
  	//TWIDDLE 0 - 256
	temp_imag = _mm256_mul_pd(v32_i,imag_twid);	
	temp_real = _mm256_mul_pd(v32_r,imag_twid);

	sub_real = _mm256_sub_pd(temp_real,temp_imag);
	temp_imag = _mm256_add_pd(temp_real,temp_imag);

	v32_r = _mm256_sub_pd(v0_r,sub_real);
	v32_i = _mm256_sub_pd(v0_i,temp_imag);

	v0_r = _mm256_add_pd(v0_r,sub_real);
	v0_i = _mm256_add_pd(v0_i,temp_imag);
	//TWIDDLE 64-320
	temp_imag = _mm256_mul_pd(v40_i,imag_twid);	
	temp_real = _mm256_mul_pd(v40_r,imag_twid);

	sub_real = _mm256_sub_pd(temp_real,temp_imag);
	temp_imag = _mm256_add_pd(temp_real,temp_imag);

	v40_r = _mm256_sub_pd(v8_r,sub_real);
	v40_i = _mm256_sub_pd(v8_i,temp_imag);

	v8_r = _mm256_add_pd(v8_r,sub_real);
	v8_i = _mm256_add_pd(v8_i,temp_imag);
	//TWIDDLE 128-384
	temp_imag = _mm256_mul_pd(v48_i,imag_twid);	
	temp_real = _mm256_mul_pd(v48_r,imag_twid);

	sub_real = _mm256_sub_pd(temp_real,temp_imag);
	temp_imag = _mm256_add_pd(temp_real,temp_imag);

	v48_r = _mm256_sub_pd(v16_r,sub_real);
	v48_i = _mm256_sub_pd(v16_i,temp_imag);

	v16_r = _mm256_add_pd(v16_r,sub_real);
	v16_i = _mm256_add_pd(v16_i,temp_imag);
	//TWIDDLE 192-448
	temp_imag = _mm256_mul_pd(v56_i,imag_twid);	
	temp_real = _mm256_mul_pd(v56_r,imag_twid);

	sub_real = _mm256_sub_pd(temp_real,temp_imag);
	temp_imag = _mm256_add_pd(temp_real,temp_imag);

	v56_r = _mm256_sub_pd(v24_r,sub_real);
	v56_i = _mm256_sub_pd(v24_i,temp_imag);

	v24_r = _mm256_add_pd(v24_r,sub_real);
	v24_i = _mm256_add_pd(v24_i,temp_imag);

	//NOW WE START WITH 256
	//FIRST LEFT SIDE
  	real_twiddle = _mm256_set1_pd(wortel[0][7][0]);
	imag_twiddle = _mm256_set1_pd(wortel[1][7][0]);
	//TWIDDLE 0 - 128
	temp_real = _mm256_mul_pd(v16_i,imag_twiddle);
	temp_imag = _mm256_mul_pd(v16_i,real_twiddle);

	temp_real = _mm256_fmsub_pd(v16_r,real_twiddle,temp_real);
	temp_imag = _mm256_fmadd_pd(v16_r,imag_twiddle,temp_imag);

	v16_r = _mm256_sub_pd(v0_r,temp_real);
	v16_i = _mm256_sub_pd(v0_i,temp_imag);

	v0_r = _mm256_add_pd(v0_r,temp_real);
	v0_i = _mm256_add_pd(v0_i,temp_imag);
	//TWIDDLE 64 - 192
	temp_real = _mm256_mul_pd(v24_i,imag_twiddle);
	temp_imag = _mm256_mul_pd(v24_i,real_twiddle);

	temp_real = _mm256_fmsub_pd(v24_r,real_twiddle,temp_real);
	temp_imag = _mm256_fmadd_pd(v24_r,imag_twiddle,temp_imag);

	v24_r = _mm256_sub_pd(v8_r,temp_real);
	v24_i = _mm256_sub_pd(v8_i,temp_imag);

	v8_r = _mm256_add_pd(v8_r,temp_real);
	v8_i = _mm256_add_pd(v8_i,temp_imag);
	//NOW WE DO THE RIGHT SIDE
  	real_twiddle = _mm256_set1_pd(wortel[0][7][1]);
	imag_twiddle = _mm256_set1_pd(wortel[1][7][1]);
	//TWIDDLE 256 - 384
	temp_real = _mm256_mul_pd(v48_i,imag_twiddle);
	temp_imag = _mm256_mul_pd(v48_i,real_twiddle);

	temp_real = _mm256_fmsub_pd(v48_r,real_twiddle,temp_real);
	temp_imag = _mm256_fmadd_pd(v48_r,imag_twiddle,temp_imag);

	v48_r = _mm256_sub_pd(v32_r,temp_real);
	v48_i = _mm256_sub_pd(v32_i,temp_imag);

	v32_r = _mm256_add_pd(v32_r,temp_real);
	v32_i = _mm256_add_pd(v32_i,temp_imag);
	//TWIDDLE 320 - 448
	temp_real = _mm256_mul_pd(v56_i,imag_twiddle);
	temp_imag = _mm256_mul_pd(v56_i,real_twiddle);

	temp_real = _mm256_fmsub_pd(v56_r,real_twiddle,temp_real);
	temp_imag = _mm256_fmadd_pd(v56_r,imag_twiddle,temp_imag);

	v56_r = _mm256_sub_pd(v40_r,temp_real);
	v56_i = _mm256_sub_pd(v40_i,temp_imag);

	v40_r = _mm256_add_pd(v40_r,temp_real);
	v40_i = _mm256_add_pd(v40_i,temp_imag);

	//START LAYER 128
	//WE NEED TO LOAD INDIVIDUAL TWIDDLES PER TWIDDLE
  	real_twiddle = _mm256_set1_pd(wortel[0][6][0]);
	imag_twiddle = _mm256_set1_pd(wortel[1][6][0]);
	//TWIDDLE 0 - 64
	temp_real = _mm256_mul_pd(v8_i,imag_twiddle);
	temp_imag = _mm256_mul_pd(v8_i,real_twiddle);

	temp_real = _mm256_fmsub_pd(v8_r,real_twiddle,temp_real);
	temp_imag = _mm256_fmadd_pd(v8_r,imag_twiddle,temp_imag);

	v8_r = _mm256_sub_pd(v0_r,temp_real);
	v8_i = _mm256_sub_pd(v0_i,temp_imag);

	v0_r = _mm256_add_pd(v0_r,temp_real);
	v0_i = _mm256_add_pd(v0_i,temp_imag);
	//WE NEED TO LOAD INDIVIDUAL TWIDDLES PER TWIDDLE
  	real_twiddle = _mm256_set1_pd(wortel[0][6][1]);
	imag_twiddle = _mm256_set1_pd(wortel[1][6][1]);
	//TWIDDLE 128 - 192
	temp_real = _mm256_mul_pd(v24_i,imag_twiddle);
	temp_imag = _mm256_mul_pd(v24_i,real_twiddle);

	temp_real = _mm256_fmsub_pd(v24_r,real_twiddle,temp_real);
	temp_imag = _mm256_fmadd_pd(v24_r,imag_twiddle,temp_imag);

	v24_r = _mm256_sub_pd(v16_r,temp_real);
	v24_i = _mm256_sub_pd(v16_i,temp_imag);

	v16_r = _mm256_add_pd(v16_r,temp_real);
	v16_i = _mm256_add_pd(v16_i,temp_imag);
	//WE NEED TO LOAD INDIVIDUAL TWIDDLES PER TWIDDLE
  	real_twiddle = _mm256_set1_pd(wortel[0][6][2]);
	imag_twiddle = _mm256_set1_pd(wortel[1][6][2]);
	//TWIDDLE 256 - 320
	temp_real = _mm256_mul_pd(v40_i,imag_twiddle);
	temp_imag = _mm256_mul_pd(v40_i,real_twiddle);

	temp_real = _mm256_fmsub_pd(v40_r,real_twiddle,temp_real);
	temp_imag = _mm256_fmadd_pd(v40_r,imag_twiddle,temp_imag);

	v40_r = _mm256_sub_pd(v32_r,temp_real);
	v40_i = _mm256_sub_pd(v32_i,temp_imag);

	v32_r = _mm256_add_pd(v32_r,temp_real);
	v32_i = _mm256_add_pd(v32_i,temp_imag);
	//WE NEED TO LOAD INDIVIDUAL TWIDDLES PER TWIDDLE
  	real_twiddle = _mm256_set1_pd(wortel[0][6][3]);
	imag_twiddle = _mm256_set1_pd(wortel[1][6][3]);
	//TWIDDLE 384 - 448
	temp_real = _mm256_mul_pd(v56_i,imag_twiddle);
	temp_imag = _mm256_mul_pd(v56_i,real_twiddle);

	temp_real = _mm256_fmsub_pd(v56_r,real_twiddle,temp_real);
	temp_imag = _mm256_fmadd_pd(v56_r,imag_twiddle,temp_imag);

	v56_r = _mm256_sub_pd(v48_r,temp_real);
	v56_i = _mm256_sub_pd(v48_i,temp_imag);

	v48_r = _mm256_add_pd(v48_r,temp_real);
	v48_i = _mm256_add_pd(v48_i,temp_imag);
	
	//STORE ALL RESULTS
	_mm256_store_pd(x->real+offset,v0_r);
	_mm256_store_pd(x->real+offset+64,v8_r);
	_mm256_store_pd(x->real+offset+128,v16_r);
	_mm256_store_pd(x->real+offset+192,v24_r);
	_mm256_store_pd(x->real+offset+256,v32_r);
	_mm256_store_pd(x->real+offset+320,v40_r);
	_mm256_store_pd(x->real+offset+384,v48_r);
	_mm256_store_pd(x->real+offset+448,v56_r);

	_mm256_store_pd(x->imag+offset,v0_i);
	_mm256_store_pd(x->imag+offset+64,v8_i);
	_mm256_store_pd(x->imag+offset+128,v16_i);
	_mm256_store_pd(x->imag+offset+192,v24_i);
	_mm256_store_pd(x->imag+offset+256,v32_i);
	_mm256_store_pd(x->imag+offset+320,v40_i);
	_mm256_store_pd(x->imag+offset+384,v48_i);
	_mm256_store_pd(x->imag+offset+448,v56_i);
 }

  //Make sure n = 64
  for(lo = 0;lo < CPLXDIM;lo +=64)
  {	
  	real_twid = _mm256_set1_pd(wortel[0][5][lo/64]);
	imag_twid = _mm256_set1_pd(wortel[1][5][lo/64]);
  	for(int offset = 0; offset < 8; offset+=4)
  	{	
  		//LOAD A LOT of STUF
	  	//FIRST ROUND WE NEED 0..3,8..11,16..19,24..27 LEFT SIDE (Real and Imag)
	  	//32..35,40..43,48..51,56..59 (Real and Imag)
	  	//REAL PART
	  	v0_r  = _mm256_load_pd(x->real+lo+offset);
	  	v8_r  = _mm256_load_pd(x->real+lo+offset+8);
	  	v16_r = _mm256_load_pd(x->real+lo+offset+16);
	  	v24_r = _mm256_load_pd(x->real+lo+offset+24);
	  	v32_r = _mm256_load_pd(x->real+lo+offset+32);
	  	v40_r = _mm256_load_pd(x->real+lo+offset+40);
	  	v48_r = _mm256_load_pd(x->real+lo+offset+48);
	  	v56_r = _mm256_load_pd(x->real+lo+offset+56);
	  	//IMAG PART
		v0_i  = _mm256_load_pd(x->imag+lo+offset);
	  	v8_i  = _mm256_load_pd(x->imag+lo+offset+8);
	  	v16_i = _mm256_load_pd(x->imag+lo+offset+16);
	  	v24_i = _mm256_load_pd(x->imag+lo+offset+24);
	  	v32_i = _mm256_load_pd(x->imag+lo+offset+32);
	  	v40_i = _mm256_load_pd(x->imag+lo+offset+40);
	  	v48_i = _mm256_load_pd(x->imag+lo+offset+48);
	  	v56_i = _mm256_load_pd(x->imag+lo+offset+56);

	  	
	  	//START TWIDDLE
	  	//WE ARE NOW IN 64 SO EVERYTHING BETWEEN 32 and 64 needs to be multiplied by root of unity
	  	//TWIDDLE 0 - 32
		temp_real = _mm256_mul_pd(v32_i,imag_twid);
		temp_imag = _mm256_mul_pd(v32_i,real_twid);

		temp_real = _mm256_fmsub_pd(v32_r,real_twid,temp_real);
		temp_imag = _mm256_fmadd_pd(v32_r,imag_twid,temp_imag);

		v32_r = _mm256_sub_pd(v0_r,temp_real);
		v32_i = _mm256_sub_pd(v0_i,temp_imag);

		v0_r = _mm256_add_pd(v0_r,temp_real);
		v0_i = _mm256_add_pd(v0_i,temp_imag);
		//TWIDDLE 8-40
		temp_real = _mm256_mul_pd(v40_i,imag_twid);
		temp_imag = _mm256_mul_pd(v40_i,real_twid);

		temp_real = _mm256_fmsub_pd(v40_r,real_twid,temp_real);
		temp_imag = _mm256_fmadd_pd(v40_r,imag_twid,temp_imag);

		v40_r = _mm256_sub_pd(v8_r,temp_real);
		v40_i = _mm256_sub_pd(v8_i,temp_imag);

		v8_r = _mm256_add_pd(v8_r,temp_real);
		v8_i = _mm256_add_pd(v8_i,temp_imag);
		//TWIDDLE 16-48
		temp_real = _mm256_mul_pd(v48_i,imag_twid);
		temp_imag = _mm256_mul_pd(v48_i,real_twid);

		temp_real = _mm256_fmsub_pd(v48_r,real_twid,temp_real);
		temp_imag = _mm256_fmadd_pd(v48_r,imag_twid,temp_imag);

		v48_r = _mm256_sub_pd(v16_r,temp_real);
		v48_i = _mm256_sub_pd(v16_i,temp_imag);

		v16_r = _mm256_add_pd(v16_r,temp_real);
		v16_i = _mm256_add_pd(v16_i,temp_imag);
		//TWIDDLE 24-56
		temp_real = _mm256_mul_pd(v56_i,imag_twid);
		temp_imag = _mm256_mul_pd(v56_i,real_twid);

		temp_real = _mm256_fmsub_pd(v56_r,real_twid,temp_real);
		temp_imag = _mm256_fmadd_pd(v56_r,imag_twid,temp_imag);

		v56_r = _mm256_sub_pd(v24_r,temp_real);
		v56_i = _mm256_sub_pd(v24_i,temp_imag);

		v24_r = _mm256_add_pd(v24_r,temp_real);
		v24_i = _mm256_add_pd(v24_i,temp_imag);

		//NOW WE START WITH 32
		//FIRST LEFT SIDE
	  	real_twiddle = _mm256_set1_pd(wortel[0][4][lo/32]);
		imag_twiddle = _mm256_set1_pd(wortel[1][4][lo/32]);
		//TWIDDLE 0 - 16
		temp_real = _mm256_mul_pd(v16_i,imag_twiddle);
		temp_imag = _mm256_mul_pd(v16_i,real_twiddle);

		temp_real = _mm256_fmsub_pd(v16_r,real_twiddle,temp_real);
		temp_imag = _mm256_fmadd_pd(v16_r,imag_twiddle,temp_imag);

		v16_r = _mm256_sub_pd(v0_r,temp_real);
		v16_i = _mm256_sub_pd(v0_i,temp_imag);

		v0_r = _mm256_add_pd(v0_r,temp_real);
		v0_i = _mm256_add_pd(v0_i,temp_imag);
		//TWIDDLE 8 - 24
		temp_real = _mm256_mul_pd(v24_i,imag_twiddle);
		temp_imag = _mm256_mul_pd(v24_i,real_twiddle);

		temp_real = _mm256_fmsub_pd(v24_r,real_twiddle,temp_real);
		temp_imag = _mm256_fmadd_pd(v24_r,imag_twiddle,temp_imag);

		v24_r = _mm256_sub_pd(v8_r,temp_real);
		v24_i = _mm256_sub_pd(v8_i,temp_imag);

		v8_r = _mm256_add_pd(v8_r,temp_real);
		v8_i = _mm256_add_pd(v8_i,temp_imag);
		//NOW WE DO THE RIGHT SIDE
	  	real_twiddle = _mm256_set1_pd(wortel[0][4][(lo+32)/32]);
		imag_twiddle = _mm256_set1_pd(wortel[1][4][(lo+32)/32]);
		//TWIDDLE 32 - 48
		temp_real = _mm256_mul_pd(v48_i,imag_twiddle);
		temp_imag = _mm256_mul_pd(v48_i,real_twiddle);

		temp_real = _mm256_fmsub_pd(v48_r,real_twiddle,temp_real);
		temp_imag = _mm256_fmadd_pd(v48_r,imag_twiddle,temp_imag);

		v48_r = _mm256_sub_pd(v32_r,temp_real);
		v48_i = _mm256_sub_pd(v32_i,temp_imag);

		v32_r = _mm256_add_pd(v32_r,temp_real);
		v32_i = _mm256_add_pd(v32_i,temp_imag);
		//TWIDDLE 40 - 56
		temp_real = _mm256_mul_pd(v56_i,imag_twiddle);
		temp_imag = _mm256_mul_pd(v56_i,real_twiddle);

		temp_real = _mm256_fmsub_pd(v56_r,real_twiddle,temp_real);
		temp_imag = _mm256_fmadd_pd(v56_r,imag_twiddle,temp_imag);

		v56_r = _mm256_sub_pd(v40_r,temp_real);
		v56_i = _mm256_sub_pd(v40_i,temp_imag);

		v40_r = _mm256_add_pd(v40_r,temp_real);
		v40_i = _mm256_add_pd(v40_i,temp_imag);

		//START LAYER 16
		//WE NEED TO LOAD INDIVIDUAL TWIDDLES PER TWIDDLE
	  	real_twiddle = _mm256_set1_pd(wortel[0][3][(lo)/16]);
		imag_twiddle = _mm256_set1_pd(wortel[1][3][(lo)/16]);
		//TWIDDLE 0 - 8
		temp_real = _mm256_mul_pd(v8_i,imag_twiddle);
		temp_imag = _mm256_mul_pd(v8_i,real_twiddle);

		temp_real = _mm256_fmsub_pd(v8_r,real_twiddle,temp_real);
		temp_imag = _mm256_fmadd_pd(v8_r,imag_twiddle,temp_imag);

		v8_r = _mm256_sub_pd(v0_r,temp_real);
		v8_i = _mm256_sub_pd(v0_i,temp_imag);

		v0_r = _mm256_add_pd(v0_r,temp_real);
		v0_i = _mm256_add_pd(v0_i,temp_imag);
		//WE NEED TO LOAD INDIVIDUAL TWIDDLES PER TWIDDLE
	  	real_twiddle = _mm256_set1_pd(wortel[0][3][(lo+16)/16]);
		imag_twiddle = _mm256_set1_pd(wortel[1][3][(lo+16)/16]);
		//TWIDDLE 16 - 24
		temp_real = _mm256_mul_pd(v24_i,imag_twiddle);
		temp_imag = _mm256_mul_pd(v24_i,real_twiddle);

		temp_real = _mm256_fmsub_pd(v24_r,real_twiddle,temp_real);
		temp_imag = _mm256_fmadd_pd(v24_r,imag_twiddle,temp_imag);

		v24_r = _mm256_sub_pd(v16_r,temp_real);
		v24_i = _mm256_sub_pd(v16_i,temp_imag);

		v16_r = _mm256_add_pd(v16_r,temp_real);
		v16_i = _mm256_add_pd(v16_i,temp_imag);
		//WE NEED TO LOAD INDIVIDUAL TWIDDLES PER TWIDDLE
	  	real_twiddle = _mm256_set1_pd(wortel[0][3][(lo+32)/16]);
		imag_twiddle = _mm256_set1_pd(wortel[1][3][(lo+32)/16]);
		//TWIDDLE 32 - 40
		temp_real = _mm256_mul_pd(v40_i,imag_twiddle);
		temp_imag = _mm256_mul_pd(v40_i,real_twiddle);

		temp_real = _mm256_fmsub_pd(v40_r,real_twiddle,temp_real);
		temp_imag = _mm256_fmadd_pd(v40_r,imag_twiddle,temp_imag);

		v40_r = _mm256_sub_pd(v32_r,temp_real);
		v40_i = _mm256_sub_pd(v32_i,temp_imag);

		v32_r = _mm256_add_pd(v32_r,temp_real);
		v32_i = _mm256_add_pd(v32_i,temp_imag);
		//WE NEED TO LOAD INDIVIDUAL TWIDDLES PER TWIDDLE
	  	real_twiddle = _mm256_set1_pd(wortel[0][3][(lo+48)/16]);
		imag_twiddle = _mm256_set1_pd(wortel[1][3][(lo+48)/16]);
		//TWIDDLE 48 - 56
		temp_real = _mm256_mul_pd(v56_i,imag_twiddle);
		temp_imag = _mm256_mul_pd(v56_i,real_twiddle);

		temp_real = _mm256_fmsub_pd(v56_r,real_twiddle,temp_real);
		temp_imag = _mm256_fmadd_pd(v56_r,imag_twiddle,temp_imag);

		v56_r = _mm256_sub_pd(v48_r,temp_real);
		v56_i = _mm256_sub_pd(v48_i,temp_imag);

		v48_r = _mm256_add_pd(v48_r,temp_real);
		v48_i = _mm256_add_pd(v48_i,temp_imag);

		//STORE ALL RESULTS
		_mm256_store_pd(x->real+lo+offset,v0_r);
		_mm256_store_pd(x->real+lo+offset+8,v8_r);
		_mm256_store_pd(x->real+lo+offset+16,v16_r);
		_mm256_store_pd(x->real+lo+offset+24,v24_r);
		_mm256_store_pd(x->real+lo+offset+32,v32_r);
		_mm256_store_pd(x->real+lo+offset+40,v40_r);
		_mm256_store_pd(x->real+lo+offset+48,v48_r);
		_mm256_store_pd(x->real+lo+offset+56,v56_r);

		_mm256_store_pd(x->imag+lo+offset,v0_i);
		_mm256_store_pd(x->imag+lo+offset+8,v8_i);
		_mm256_store_pd(x->imag+lo+offset+16,v16_i);
		_mm256_store_pd(x->imag+lo+offset+24,v24_i);
		_mm256_store_pd(x->imag+lo+offset+32,v32_i);
		_mm256_store_pd(x->imag+lo+offset+40,v40_i);
		_mm256_store_pd(x->imag+lo+offset+48,v48_i);
		_mm256_store_pd(x->imag+lo+offset+56,v56_i);
  	}
  }
  for(lo=0;lo < CPLXDIM;lo +=8)
  {	
    __m256d real_x,imag_x,real_y,imag_y;
    real_twid = _mm256_set1_pd(wortel[0][2][lo/8]);
    imag_twid = _mm256_set1_pd(wortel[1][2][lo/8]);
	//(a + ib) * (c + id) = (ac - bd) + i(ad+bc)
	real_x = _mm256_load_pd(x->real+lo);
	imag_x = _mm256_load_pd(x->imag+lo);
	real_y = _mm256_load_pd(x->real+lo+4);
	imag_y = _mm256_load_pd(x->imag+lo+4);
	//TEMP_real = bd
	temp_real = _mm256_mul_pd(imag_y,imag_twid);
	//TEMP_imag = bc
	temp_imag = _mm256_mul_pd(imag_y,real_twid);

	//TEMP_real = ac - bd
	temp_real = _mm256_fmsub_pd(real_y,real_twid,temp_real);
	//TEMP_imag = ad + bc
	temp_imag = _mm256_fmadd_pd(real_y,imag_twid,temp_imag);

	real_y = _mm256_sub_pd(real_x,temp_real);
	imag_y = _mm256_sub_pd(imag_x,temp_imag);

	real_x = _mm256_add_pd(real_x,temp_real);
	imag_x = _mm256_add_pd(imag_x,temp_imag);

	//START LAYER 4
	//We know that real_x = |a|b|c|d| imag_x = |ai|bi|ci|di| real_y = |e|f|g|h| imag_y = |ei|fi|gi|hi| 
	//We need to twidle (c+ci),(d+di) and (g+gi),(h+hi)

	sub_real = _mm256_permute2f128_pd(real_x,real_y,0x31);
	sub_imag = _mm256_permute2f128_pd(imag_x,imag_y,0x31);
	real_twid = _mm256_setr_pd(wortel[0][1][lo/4],wortel[0][1][lo/4],wortel[0][1][(lo+4)/4],wortel[0][1][(lo+4)/4]);
    imag_twid = _mm256_setr_pd(wortel[1][1][lo/4],wortel[1][1][lo/4],wortel[1][1][(lo+4)/4],wortel[1][1][(lo+4)/4]);

    temp_real = _mm256_mul_pd(sub_imag,imag_twid);
    temp_imag = _mm256_mul_pd(sub_imag,real_twid);
	//TEMP_real = ac - bd
	temp_real = _mm256_fmsub_pd(sub_real,real_twid,temp_real);
	//TEMP_imag = ad + bc
	temp_imag = _mm256_fmadd_pd(sub_real,imag_twid,temp_imag);
	
	//REAL PART
	//get abef
	sub_real = _mm256_permute2f128_pd(real_x,real_y,0x20);
	//abef  
	//cdgh-
	sub_imag = _mm256_sub_pd(sub_real,temp_real);
	//abef
	//cdgh+
	sub_real = _mm256_add_pd(sub_real,temp_real);

	// //NEEDED TO COMPLETE LAYER 4
	// real_x = _mm256_permute2f128_pd(sub_imag,sub_real,0x02);
	// real_y = _mm256_permute2f128_pd(sub_imag,sub_real,0x13);
	//PREPARE REALS FOR LAYER 2
	//STORE ALL FACTORS THAT NEED TO BE MULT WITH ROOTS OF UNITY IN REAL_X
	real_x = _mm256_unpackhi_pd(sub_real,sub_imag);
	real_y = _mm256_unpacklo_pd(sub_real,sub_imag);


	//IMAG PART
	//get ai bi ei fi
	sub_real = _mm256_permute2f128_pd(imag_x,imag_y,0x20);
	//abef  
	//cdgh-
	sub_imag = _mm256_sub_pd(sub_real,temp_imag);
	//abef
	//cdgh+
	sub_real = _mm256_add_pd(sub_real,temp_imag);

	//NEEDED TO COMPLETE LAYER 4
	// imag_x = _mm256_permute2f128_pd(sub_imag,sub_real,0x02);
	// imag_y = _mm256_permute2f128_pd(sub_imag,sub_real,0x13);
	//PREPARE IMAGS FOR LAYER 2
	//STORE ALL FACTORS THAT NEED TO BE MULT WITH ROOTS OF UNITY IN IMAG_X
	//STORE ALL NON TWIDLE IMAGS IN IMAG_Y 
	imag_x = _mm256_unpackhi_pd(sub_real,sub_imag);
	imag_y = _mm256_unpacklo_pd(sub_real,sub_imag);

	// // //START LAYER 2!!
	real_twid = _mm256_setr_pd(wortel[0][0][lo/2],wortel[0][0][(lo+2)/2],wortel[0][0][(lo+4)/2],wortel[0][0][(lo+6)/2]);
    imag_twid = _mm256_setr_pd(wortel[1][0][lo/2],wortel[1][0][(lo+2)/2],wortel[1][0][(lo+4)/2],wortel[1][0][(lo+6)/2]);

    temp_real = _mm256_mul_pd(imag_x,imag_twid);
    temp_imag = _mm256_mul_pd(imag_x,real_twid);

	//TEMP_real = ac - bd
	temp_real = _mm256_fmsub_pd(real_x,real_twid,temp_real);
	//TEMP_imag = ad + bc
	temp_imag = _mm256_fmadd_pd(real_x,imag_twid,temp_imag);

	sub_real = _mm256_sub_pd(real_y,temp_real);
	sub_imag = _mm256_add_pd(real_y,temp_real); 

	temp_real = _mm256_unpacklo_pd(sub_imag,sub_real);
	sub_real  = _mm256_unpackhi_pd(sub_imag,sub_real);

	real_x = _mm256_permute2f128_pd(temp_real,sub_real,0x20);
	real_y = _mm256_permute2f128_pd(temp_real,sub_real,0x31);

	sub_real = _mm256_sub_pd(imag_y,temp_imag);
	sub_imag = _mm256_add_pd(imag_y,temp_imag);

	temp_real = _mm256_unpacklo_pd(sub_imag,sub_real);
	sub_real  = _mm256_unpackhi_pd(sub_imag,sub_real);

	imag_x = _mm256_permute2f128_pd(temp_real,sub_real,0x20);
	imag_y = _mm256_permute2f128_pd(temp_real,sub_real,0x31);

	_mm256_store_pd(x->real+lo,real_x);
	_mm256_store_pd(x->imag+lo,imag_x);
	_mm256_store_pd(x->real+lo+4,real_y);
	_mm256_store_pd(x->imag+lo+4,imag_y);	
  }
}

void phi_forward(cplx_ptr *x,const ring_t *ring)
{
  int j = CPLXDIM;
  for (int i = 0; i < CPLXDIM; ++i)
  {
    x->real[i] = ring->v[i];
    x->imag[i] = ring->v[j];
    ++j;
  }
  iterative_phi(x);
}

void phi_backward(cplx_ptr *x, ring_t *ring)
{
  inverse_phi(x);

  int j = CPLXDIM;
  for (int i = 0; i < CPLXDIM; ++i)
  {
    ring->v[i] = x->real[i];
    ring->v[j] = x->imag[i];
    ++j;
  }
}