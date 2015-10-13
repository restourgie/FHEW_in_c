#include <complex.h>
#include <stdio.h>
#include <math.h>
#include <immintrin.h>
#include "../mul.h"

#define NRLAYERS 8
int layer_size[NRLAYERS] = {1,1,3,5,11,21,43,85};

int locations[NRLAYERS][85] = 
{
	{0},
	{0},
	{0, 256, 384},
	{0, 128, 192, 256, 384},
	{0,64, 96, 128, 192, 256, 320, 352, 384, 448, 480},
	{0, 32, 48,64, 96, 128, 160, 176, 192, 224, 240,  256, 288, 304, 320, 352, 384, 416, 432, 448, 480},
	{0, 16, 24, 32, 48,64, 80, 88, 96, 112, 120, 128, 144, 152, 160, 176, 192, 208, 216, 224, 240, 256, 
	272, 280, 288, 304, 320, 336, 344, 352, 368, 376, 384, 400, 408, 416, 432, 448, 464, 472, 480, 496, 504},
	{0, 8, 12, 16, 24, 32, 40, 44, 48, 56, 60, 64, 72, 76, 80, 88, 96, 104, 108, 112, 120, 128, 136, 140, 
	144, 152, 160, 168, 172, 176, 184, 188, 192, 200, 204, 208, 216, 224, 232, 236, 240, 248, 252, 256, 264, 268, 272, 280, 
	288, 296, 300, 304, 312, 316, 320, 328, 332, 336, 344, 352, 360, 364, 368, 376, 384, 392, 396, 400, 408, 416, 424, 428, 
	432, 440, 444, 448, 456, 460, 464, 472, 480, 488, 492, 496, 504}
};

double **LUT1,**LUT2,**LUT3;

/*
* This function initialises all the sizes and locations of the loop
*
* This function initializes the lookup tables(LUTs) 
* We calculate all the roots of unity for 8,16,32,64,128,256 and 512
* We also calculate the first 512 roots of unity for 2048 
* This is needed for Bernsteins trick (x^512-i) equals (x^512-1) when twisted
* with the roots of 2048
*/
void init_vctr()
{
	int size = 8, j=8;
	LUT1 = malloc(sizeof *LUT1 * size);
	LUT2 = malloc(sizeof *LUT2 * size);
	LUT3 = malloc(sizeof *LUT3 * size);
	for (int i = 0; i < size-1; ++i)
	{
		posix_memalign((void**)&LUT1[i],32, j * sizeof(double));
		posix_memalign((void**)&LUT2[i],32, j * sizeof(double));
		posix_memalign((void**)&LUT3[i],32, j * sizeof(double));
		for (int root = 0; root < j; ++root)
		{
		  LUT1[i][root] = calc_cos(j,root);
		  LUT2[i][root] = calc_sin(j,root);
		  LUT3[i][root] = -LUT2[i][root];
		}
		j = j<<1;
	}
	posix_memalign((void**)&LUT1[7],32, 512 * sizeof(double));
	posix_memalign((void**)&LUT2[7],32, 512 * sizeof(double));
	posix_memalign((void**)&LUT3[7],32, 512 * sizeof(double));

	for (int i = 0; i < CPLXDIM; ++i)
	{
	  LUT1[7][i] = calc_cos(ROOTDIM,i);
	  LUT2[7][i] = calc_sin(ROOTDIM,i);
	  LUT3[7][i] = -LUT2[7][i]; 
	}
}

void destruct_vctr()
{
	for (int i = 7; i >= 0; --i)
	{
		free(LUT1[i]);
		free(LUT2[i]);
		free(LUT3[i]);
	}
	free(LUT1);
	free(LUT2);
	free(LUT3);
}

/******************************************************************
*
* LOOKUPTABLES FOR VECTOR TWIST
*
******************************************************************/
void vec_twist(cplx_ptr *cplx_x,int n,int m,int lo)
{
  __m256d real_x,imag_x,real_tbl,imag_tbl,imag_temp,real_temp;
  int j = 0, scale;
  if(n==ROOTDIM)
    scale = 7;
  else
    scale = log2(n)-3;
  for (int i = lo; i < lo+m; i+=4)
  { 
    // for (int bla = j; bla < j+4; ++bla)
    // {
    //   printf("root= %f + i * %f\n",table[0][scale*bla],table[1][scale*bla]);
    // }
    real_x = _mm256_load_pd(cplx_x->real+i);
    imag_x = _mm256_load_pd(cplx_x->imag+i);
    real_tbl = _mm256_load_pd(&LUT1[scale][j]);
    imag_tbl = _mm256_load_pd(&LUT2[scale][j]);
     
    //(a + ib) * (c + id) = (ac - bd) + i(ad+bc)
    //real_temp = bd
    real_temp = _mm256_mul_pd(imag_x,imag_tbl);
    //imag_temp = ad
    imag_temp = _mm256_mul_pd(real_x,imag_tbl);
    //REPLACED FOR COMMENTED SECTION
    //real_x = ac
    real_x = _mm256_mul_pd(real_x,real_tbl);
    //imag_x = bc
    imag_x = _mm256_mul_pd(imag_x,real_tbl);
    //real_x = ac - bd => real_x - real_temp
    real_x = _mm256_sub_pd(real_x,real_temp);
    //imag_x = ad + bc => imag_temp + imag_x
    imag_x = _mm256_add_pd(imag_x,imag_temp);
    //THESE ARE NOT WORKING 
    // real_x = _mm256_fmsub_pd(real_x,real_tbl,real_temp);
    // imag_x = _mm256_fmadd_pd(imag_x,real_tbl,imag_temp);
    _mm256_store_pd(cplx_x->real+i,real_x);
    _mm256_store_pd(cplx_x->imag+i,imag_x);
    j+=4;
  }
}

void vec_untwist(cplx_ptr *cplx_x,int n,int m,int lo)
{
  __m256d real_x,imag_x,real_tbl,imag_tbl,imag_temp,real_temp;
  int j = 0, scale;
  if(n==ROOTDIM)
    scale = 7;
  else
    scale = log2(n)-3;
  for (int i = lo; i < lo+m; i+=4)
  { 
    // for (int bla = j; bla < j+4; ++bla)
    // {
    //   printf("root= %f + i * %f\n",table[0][scale*bla],table[1][scale*bla]);
    // }
    real_x = _mm256_load_pd(cplx_x->real+i);
    imag_x = _mm256_load_pd(cplx_x->imag+i);
    real_tbl = _mm256_load_pd(&LUT1[scale][j]);
    imag_tbl = _mm256_load_pd(&LUT3[scale][j]);
     
    //(a + ib) * (c + id) = (ac - bd) + i(ad+bc)
    //real_temp = bd
    real_temp = _mm256_mul_pd(imag_x,imag_tbl);
    //imag_temp = ad
    imag_temp = _mm256_mul_pd(real_x,imag_tbl);
    //REPLACED FOR COMMENTED SECTION
    //real_x = ac
    real_x = _mm256_mul_pd(real_x,real_tbl);
    //imag_x = bc
    imag_x = _mm256_mul_pd(imag_x,real_tbl);

    //real_x = ac - bd => real_x - real_temp
    real_x = _mm256_sub_pd(real_x,real_temp);
    //imag_x = ad + bc => imag_temp + imag_x
    imag_x = _mm256_add_pd(imag_x,imag_temp);
    //THESE ARE NOT WORKING 
    // real_x = _mm256_fmsub_pd(real_x,real_tbl,real_temp);
    // imag_x = _mm256_fmadd_pd(imag_x,real_tbl,imag_temp);
    _mm256_store_pd(cplx_x->real+i,real_x);
    _mm256_store_pd(cplx_x->imag+i,imag_x);
    j+=4;
  }
}

void sr_vector_nonrec(cplx_ptr *x,int n)
{
	int new_n =n,lo,m;
	double temp_real,temp_im;
  	__m256d real_x,imag_x,real_y,imag_y,imag_temp,real_temp;
  	__m128d real_x_small,imag_x_small,imag_y_small,real_y_small,imag_temp_small,real_temp_small;
	//DO LAYERS 512,256,128,64,32 and 16
	for (int iter=0; iter<NRLAYERS-2; ++iter)
	{
		//This forloop wil loop for every array location defined in the current layer.
		for (int count = 0; count < layer_size[iter]; ++count)
		{
			lo = locations[iter][count];
			m = new_n/2;
			//Go from (x^4n +1 to x^2n -1 and x^2n +1)
		    for(int i=lo; i < lo+m;i+=4){
		      //LOAD REAL,LOAD IMAG i AND LOAD REAL,LOAD IMAG i + m
		      real_x = _mm256_load_pd(x->real+i);
		      imag_x = _mm256_load_pd(x->imag+i);
		      real_y = _mm256_load_pd(x->real+i+m);
		      imag_y = _mm256_load_pd(x->imag+i+m);

		      //TEMP IS X - Y 
		      real_temp = _mm256_sub_pd(real_x,real_y);
		      imag_temp = _mm256_sub_pd(imag_x,imag_y);

		      //X = X+Y
		      real_x = _mm256_add_pd(real_x,real_y);
		      imag_x = _mm256_add_pd(imag_x,imag_y);

		      _mm256_store_pd(x->real+i,real_x);
		      _mm256_store_pd(x->imag+i,imag_x);
		      _mm256_store_pd(x->real+i+m,real_temp);
		      _mm256_store_pd(x->imag+i+m,imag_temp);
		    }

		    lo = lo+m;
		    m = m/2;

		    //Go from (x^2n +1 to x^n -i and x^n +i)
		    for (int i = lo; i < lo+m; i+=4)
		    {
		      //LOAD REAL,LOAD IMAG i AND LOAD REAL,LOAD IMAG i + m
		      real_x = _mm256_load_pd(x->real+i);
		      imag_x = _mm256_load_pd(x->imag+i);
		      real_y = _mm256_load_pd(x->real+i+m);
		      imag_y = _mm256_load_pd(x->imag+i+m);   

		      //(a + ib) - i(c + id) = (a + ib) - (-d + ic) = (a+d + i(b-c))
		      real_temp = _mm256_add_pd(real_x,imag_y);
		      imag_temp = _mm256_sub_pd(imag_x,real_y);

		      //(a + ib) + i(c + id) = (a + ib) + (-d + ic) = (a-d + i(b+c))
		      real_x = _mm256_sub_pd(real_x,imag_y);
		      imag_x = _mm256_add_pd(imag_x,real_y);
		      
		      _mm256_store_pd(x->real+i,real_x);
		      _mm256_store_pd(x->imag+i,imag_x);
		      _mm256_store_pd(x->real+i+m,real_temp);
		      _mm256_store_pd(x->imag+i+m,imag_temp);
		    }
		    vec_twist(x,new_n,m,lo);
		    vec_untwist(x,new_n,m,lo+m);
		}
		new_n = new_n>>1;
	}

	//DO EVERYTHING FOR 8
	for (int count = 0; count < layer_size[6]; ++count)
	{
		lo = locations[6][count];
		m = 4;
	    //Go from (x^8+1) to (x^4-1) and (x^4+1)
	    //LOAD REAL,LOAD IMAG i AND LOAD REAL,LOAD IMAG i + m
	    real_x = _mm256_load_pd(x->real+lo);
	    imag_x = _mm256_load_pd(x->imag+lo);
	    real_y = _mm256_load_pd(x->real+lo+m);
	    imag_y = _mm256_load_pd(x->imag+lo+m);

	    //TEMP IS X - Y 
	    real_temp = _mm256_sub_pd(real_x,real_y);
	    imag_temp = _mm256_sub_pd(imag_x,imag_y);

	    //X = X+Y
	    real_x = _mm256_add_pd(real_x,real_y);
	    imag_x = _mm256_add_pd(imag_x,imag_y);

	    _mm256_store_pd(x->real+lo,real_x);
	    _mm256_store_pd(x->imag+lo,imag_x);
	    _mm256_store_pd(x->real+lo+m,real_temp);
	    _mm256_store_pd(x->imag+lo+m,imag_temp);

	    lo = lo+m;
	    m = m/2;
	    //Go from (x^4-1) to (x^2-i) and (x^2+i)
	    //LOAD REAL,LOAD IMAG i AND LOAD REAL,LOAD IMAG i + m
	    real_x_small = _mm_load_pd(x->real+lo);
	    imag_x_small = _mm_load_pd(x->imag+lo);
	    real_y_small = _mm_load_pd(x->real+lo+m);
	    imag_y_small = _mm_load_pd(x->imag+lo+m);   

	    //(a + ib) - i(c + id) = (a + ib) - (-d + ic) = (a+d + i(b-c))
	    real_temp_small = _mm_add_pd(real_x_small,imag_y_small);
	    imag_temp_small = _mm_sub_pd(imag_x_small,real_y_small);

	    //(a + ib) + i(c + id) = (a + ib) + (-d + ic) = (a-d + i(b+c))
	    real_x_small = _mm_sub_pd(real_x_small,imag_y_small);
	    imag_x_small = _mm_add_pd(imag_x_small,real_y_small);
	    
	    _mm_store_pd(x->real+lo,real_x_small);
	    _mm_store_pd(x->imag+lo,imag_x_small);
	    _mm_store_pd(x->real+lo+m,real_temp_small);
	    _mm_store_pd(x->imag+lo+m,imag_temp_small);
	    //WE ARE GOING TO DO A MANUAL TWIST SINCE 
	    //(x->real[lo+1],x->imag[lo+1]) AND (x->real[lo+m+1],x->imag[lo+m+1])
	    //ONLY NEED TO BE TWISTED
	    int other_lo = lo+m;
	    //(a + ib) * (c + id) = (ac - bd) + i(ad+bc)
	    temp_real = (x->real[lo+1] * LUT1[0][1]) - (x->imag[lo+1] * LUT2[0][1]);//ac -bd
	    x->imag[lo+1] = (x->real[lo+1] * LUT2[0][1]) + (x->imag[lo+1] * LUT1[0][1]);
	    x->real[lo+1] = temp_real;
	    //OTHER TWIST
	    temp_real = (x->real[other_lo+1] * LUT1[0][1]) - (x->imag[other_lo+1] * LUT3[0][1]);//ac -bd
	    x->imag[other_lo+1] = (x->real[other_lo+1] * LUT3[0][1]) + (x->imag[other_lo+1] * LUT1[0][1]);
	    x->real[other_lo+1] = temp_real;

	    //WE CAN NOW FINISH THE FFT BY ADDING AND SUBTRACTING
	    //LEFT PART
	    temp_real = x->real[lo];
	    temp_im = x->imag[lo];
	    x->real[lo] = temp_real + x->real[lo+1];
	    x->imag[lo] = temp_im + x->imag[lo+1];
	    x->real[lo+1] = temp_real - x->real[lo+1];
	    x->imag[lo+1] = temp_im - x->imag[lo+1];
	    //RIGHT PART
	    temp_real = x->real[other_lo];
	    temp_im = x->imag[other_lo];
	    x->real[other_lo] = temp_real + x->real[other_lo+1];
	    x->imag[other_lo] = temp_im + x->imag[other_lo+1];
	    x->real[other_lo+1] = temp_real - x->real[other_lo+1];
	    x->imag[other_lo+1] = temp_im - x->imag[other_lo+1];
	}

	//DO EVERYTHING FOR 4
	for (int count = 0; count < layer_size[7]; ++count)
	{
		lo = locations[7][count];
		    // printf("ENTERING\n");
	    int m = 2;
	    //Go from (x^4-1) to (x^2-1) and (x^2+1)
	    real_x_small = _mm_load_pd(x->real+lo);
	    imag_x_small = _mm_load_pd(x->imag+lo);
	    real_y_small = _mm_load_pd(x->real+lo+m);
	    imag_y_small = _mm_load_pd(x->imag+lo+m);

	    real_temp_small = _mm_sub_pd(real_x_small,real_y_small);
	    imag_temp_small = _mm_sub_pd(imag_x_small,imag_y_small);
	    real_x_small = _mm_add_pd(real_x_small,real_y_small);
	    imag_x_small = _mm_add_pd(imag_x_small,imag_y_small);

	    _mm_store_pd(x->real+lo,real_x_small);
	    _mm_store_pd(x->imag+lo,imag_x_small);
	    _mm_store_pd(x->real+lo+m,real_temp_small);
	    _mm_store_pd(x->imag+lo+m,imag_temp_small);
	    //DO LEFT BRANCH GO FROM (x^2-1) to (x-1) and (x+1)
	    temp_real = x->real[lo];
	    temp_im = x->imag[lo];
	    x->real[lo] = temp_real + x->real[lo+1];
	    x->imag[lo] = temp_im + x->imag[lo+1];
	    x->real[lo+1] = temp_real - x->real[lo+1];
	    x->imag[lo+1] = temp_im - x->imag[lo+1];
	    //DO RIGHT BRANCH GO FROM (x^2+1) to (x-i) and (x+i) 
	    lo = lo+m;
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
}

void sr_vector_nonrec_inverse(cplx_ptr *x)
{
	int new_n =16,lo,m;
	double temp_real,temp_im;
  __m256d real_x,imag_x,real_y,imag_y,imag_temp,real_temp;
  __m128d real_x_small,imag_x_small,imag_y_small,real_y_small,imag_temp_small,real_temp_small;

  // FIRST WE WILL DO ALL LAYERS OF 4
  for (int count = 0; count < layer_size[7]; ++count)
	{
		lo = locations[7][count];
		m = 2;
	    //We will first do the LEFT branch
	    temp_real = x->real[lo];
	    temp_im = x->imag[lo];
	    x->real[lo] = temp_real + x->real[lo+1];
	    x->imag[lo] = temp_im + x->imag[lo+1];
	    x->real[lo+1] = temp_real - x->real[lo+1];
	    x->imag[lo+1] = temp_im - x->imag[lo+1];
	    //Now we do the RIGHT branch
	    lo = lo+m;
	    temp_real = x->real[lo];
	    temp_im = x->imag[lo];
	    //x->real[lo] = a-d, x->imag[lo] = b+c, x->real[lo+1] = a+d, x->imag[lo+1] = b-c
	    //2a = (a-d)+(a+d), 2b = (b+c)+(b-c), 2c = (b+c)-(b-c), 2d = (a+d)-(a-d) 
	    x->real[lo]   = temp_real     + x->real[lo+1];
	    x->imag[lo]   = temp_im       + x->imag[lo+1];
	    temp_real     = x->real[lo+1] - temp_real;
	    x->real[lo+1] = temp_im       - x->imag[lo+1];
	    x->imag[lo+1] = temp_real;
	    //NOW WE WILL TIE EVERYTHING TOGETHER WE WILL GO FROM (x^2-1) and (x^2+1) TO (x^4-1)
	    lo = lo-m;
	    real_x_small = _mm_load_pd(x->real+lo);
	    imag_x_small = _mm_load_pd(x->imag+lo);
	    real_y_small = _mm_load_pd(x->real+lo+m);
	    imag_y_small = _mm_load_pd(x->imag+lo+m);

	    real_temp_small = _mm_sub_pd(real_x_small,real_y_small);
	    imag_temp_small = _mm_sub_pd(imag_x_small,imag_y_small);
	    real_x_small = _mm_add_pd(real_x_small,real_y_small);
	    imag_x_small = _mm_add_pd(imag_x_small,imag_y_small);

	    _mm_store_pd(x->real+lo,real_x_small);
	    _mm_store_pd(x->imag+lo,imag_x_small);
	    _mm_store_pd(x->real+lo+m,real_temp_small);
	    _mm_store_pd(x->imag+lo+m,imag_temp_small);
	}
	//NOW WE WILL DO ALL LAYERS OF 8
	for (int count = 0; count < layer_size[6]; ++count)
	{
		lo = locations[6][count];
		m   = 2;
	    lo  = lo+4;
	    int other_lo = lo+m; 
	    //WE CAN NOW START THE FFT BY TIEING TOGETHER (x-1) and (x+1) to (x^2-1)
	    //LEFT PART
	    temp_real = x->real[lo];
	    temp_im = x->imag[lo];
	    x->real[lo] = temp_real + x->real[lo+1];
	    x->imag[lo] = temp_im + x->imag[lo+1];
	    x->real[lo+1] = temp_real - x->real[lo+1];
	    x->imag[lo+1] = temp_im - x->imag[lo+1];
	    //RIGHT PART
	    temp_real = x->real[other_lo];
	    temp_im = x->imag[other_lo];
	    x->real[other_lo] = temp_real + x->real[other_lo+1];
	    x->imag[other_lo] = temp_im + x->imag[other_lo+1];
	    x->real[other_lo+1] = temp_real - x->real[other_lo+1];
	    x->imag[other_lo+1] = temp_im - x->imag[other_lo+1];
	    //NOW WE NEED TO DO A MANUAL TWIST TO GO FROM (x^2-1) to (x^2-i) and (x^2-1) to (x^2+i)
	    //WE DO AN INVERSE TWIST ON THE LEFT PART TO GO FROM (x^2-1) to (x^2-i)
	    temp_real = (x->real[lo+1] * LUT1[0][1]) - (x->imag[lo+1] * LUT3[0][1]);//ac -bd
	    x->imag[lo+1] = (x->real[lo+1] * LUT3[0][1]) + (x->imag[lo+1] * LUT1[0][1]);
	    x->real[lo+1] = temp_real;
	    //WE DO AN FORWARD TWIST ON THE RIGHT PART TO GO FROM (x^2-1) to (x^2+i)
	    temp_real = (x->real[other_lo+1] * LUT1[0][1]) - (x->imag[other_lo+1] * LUT2[0][1]);//ac -bd
	    x->imag[other_lo+1] = (x->real[other_lo+1] * LUT2[0][1]) + (x->imag[other_lo+1] * LUT1[0][1]);
	    x->real[other_lo+1] = temp_real;
	    //WE NOW NEED TO TIE (x^2-i) and (x^2+i) TOGETHER TO FORM (X^4+1)
	    //real_temp = a-d, imag_temp = b+c, real_y = a+d, imag_y = b-c
	    //2a = (a-d)+(a+d), 2b = (b+c)+(b-c), 2c = (b+c)-(b-c), 2d = (a+d)-(a-d) 
	    //LOAD REAL,LOAD IMAG i AND LOAD REAL,LOAD IMAG i + m
	    real_temp_small = _mm_load_pd(x->real+lo);
	    imag_temp_small = _mm_load_pd(x->imag+lo);
	    real_y_small = _mm_load_pd(x->real+lo+m);
	    imag_y_small = _mm_load_pd(x->imag+lo+m);   

	    real_x_small = _mm_add_pd(real_temp_small,real_y_small);
	    imag_x_small = _mm_add_pd(imag_temp_small,imag_y_small);

	    real_temp_small = _mm_sub_pd(real_y_small,real_temp_small);
	    real_y_small = _mm_sub_pd(imag_temp_small,imag_y_small);
	    
	    _mm_store_pd(x->real+lo,real_x_small);
	    _mm_store_pd(x->imag+lo,imag_x_small);
	    _mm_store_pd(x->real+lo+m,real_y_small);
	    _mm_store_pd(x->imag+lo+m,real_temp_small);
	    //NOW WE HAVE THE ENTIRE RIGHT BRANCH (x^4+1)
	    //RECURSIVE CALL TO GET LEFT BRANCH (x^4-1)
	    m = 4;
	    lo = lo-m;
	    //NOW WE NEED TO TIE TOGETHER (x^4-1) and (x^4+1) to form (x^8-1)
	    //LOAD REAL,LOAD IMAG i AND LOAD REAL,LOAD IMAG i + m
	    real_x = _mm256_load_pd(x->real+lo);
	    imag_x = _mm256_load_pd(x->imag+lo);
	    real_y = _mm256_load_pd(x->real+lo+m);
	    imag_y = _mm256_load_pd(x->imag+lo+m);

	    //TEMP IS X - Y 
	    real_temp = _mm256_sub_pd(real_x,real_y);
	    imag_temp = _mm256_sub_pd(imag_x,imag_y);

	    //X = X+Y
	    real_x = _mm256_add_pd(real_x,real_y);
	    imag_x = _mm256_add_pd(imag_x,imag_y);

	    _mm256_store_pd(x->real+lo,real_x);
	    _mm256_store_pd(x->imag+lo,imag_x);
	    _mm256_store_pd(x->real+lo+m,real_temp);
	    _mm256_store_pd(x->imag+lo+m,imag_temp);
	}

	//NOW DO THE OTHER LAYERS
	//DO LAYERS 512,256,128,64,32 and 16
	for (int iter=NRLAYERS-3; iter>-1; --iter)
	{	
		//This forloop wil loop for every array location defined in the current layer.
		for (int count = 0; count < layer_size[iter]; ++count)
		{
			lo = locations[iter][count];
			//WHEN n > 8 we have the general problem, in order to get (x^(4n)-1) we need to collect
		    //The left branch (x^(2n)-1) and the two right branches (x^n-1) and (x^n-1)
		    //We need to twist the right branches to (x^n-i) and (x^n+i) in order to fold them back to (x^(2n)+1)
		    //Afterwards we can fold (x^(2n)-1) and (x^(2n)+1) to get (x^4n-1)
		    m = new_n/4;
		    lo = lo+new_n/2;
		    //Twist right branch to get (x^n+i)
		    vec_twist(x,new_n,m,lo+m);
		    //unTwist left branch to get (x^n-i)
		    vec_untwist(x,new_n,m,lo);
		    // printf("BEFORE STITCHING\n");
		    // print_cplx(*x,CPLXDIM);
		    //TIE (x^n-i) and (x^n+i) back to (x^(2n)+1)
		    for (int i = lo; i < lo+m; i+=4)
		    {
		      //real_temp = a-d, imag_temp = b+c, real_y = a+d, imag_y = b-c
		      //2a = (a-d)+(a+d), 2b = (b+c)+(b-c), 2c = (b+c)-(b-c), 2d = (a+d)-(a-d) 
		      //LOAD REAL,LOAD IMAG i AND LOAD REAL,LOAD IMAG i + m
		      real_temp = _mm256_load_pd(x->real+i);
		      imag_temp = _mm256_load_pd(x->imag+i);
		      real_y = _mm256_load_pd(x->real+i+m);
		      imag_y = _mm256_load_pd(x->imag+i+m);   

		      real_x = _mm256_add_pd(real_temp,real_y);
		      imag_x = _mm256_add_pd(imag_temp,imag_y);

		      real_temp = _mm256_sub_pd(real_y,real_temp);
		      real_y = _mm256_sub_pd(imag_temp,imag_y);
		      
		      _mm256_store_pd(x->real+i,real_x);
		      _mm256_store_pd(x->imag+i,imag_x);
		      _mm256_store_pd(x->real+i+m,real_y);
		      _mm256_store_pd(x->imag+i+m,real_temp);
		    }
		    // printf("AFTER STITCHING\n");
		    // print_cplx(*x,CPLXDIM);
		    m = m*2;
		    lo = lo -m;
		    //TIE (x^(2n)-1) and (x^(2n)+1) to get (x^4n-1)
	    	for(int i=lo; i < lo+m;i+=4)
	    	{
		      //LOAD REAL,LOAD IMAG i AND LOAD REAL,LOAD IMAG i + m
		      real_x = _mm256_load_pd(x->real+i);
		      imag_x = _mm256_load_pd(x->imag+i);
		      real_y = _mm256_load_pd(x->real+i+m);
		      imag_y = _mm256_load_pd(x->imag+i+m);

		      //TEMP IS X - Y 
		      real_temp = _mm256_sub_pd(real_x,real_y);
		      imag_temp = _mm256_sub_pd(imag_x,imag_y);

		      //X = X+Y
		      real_x = _mm256_add_pd(real_x,real_y);
		      imag_x = _mm256_add_pd(imag_x,imag_y);

		      _mm256_store_pd(x->real+i,real_x);
		      _mm256_store_pd(x->imag+i,imag_x);
		      _mm256_store_pd(x->real+i+m,real_temp);
		      _mm256_store_pd(x->imag+i+m,imag_temp);
		  	}

		}
		new_n = new_n<<1;
	}

}

void fft_vector_nonrec_forward(cplx_ptr *x){
  vec_twist(x,ROOTDIM,CPLXDIM,0);
  sr_vector_nonrec(x,CPLXDIM);
}

void fft_vector_nonrec_backward(cplx_ptr *cplx_x,ring_t *res)
{
  sr_vector_nonrec_inverse(cplx_x);
  vec_untwist(cplx_x,ROOTDIM,CPLXDIM,0);
 
  int j = CPLXDIM;
  for (int i = 0; i < CPLXDIM; ++i)
  {
    res->v[i] = cplx_x->real[i]/CPLXDIM;
    res->v[j] = cplx_x->imag[i]/CPLXDIM;
    ++j; 
  } 
}