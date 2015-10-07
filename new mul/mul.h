#ifndef MUL_H
#define MUL_H

#define CPLXDIM 512
#define REALDIM (2*CPLXDIM)
#define ROOTDIM (2*REALDIM)

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#include <stdint.h>
#include <complex.h>

typedef struct {
  uint32_t v[REALDIM];
} ring_t;

typedef struct {
	double real[CPLXDIM];
	double imag[CPLXDIM];
} cplx;

typedef struct
{	
	double *real;
	double *imag;
}cplx_ptr;

double **LUT1,**LUT2,**LUT3;

void print_complex(const double complex *a, int N);
void print_double(const cplx *x,int N);
void to_complex(const ring_t *x, double complex *cplx_x);
void to_real(const double complex *cplx_x, ring_t *x);
void twist(double complex *cplx_x,int n,int m,int lo);
void untwist(double complex *cplx_x,int n,int m,int lo);
void table_twist(cplx *cplx_x,int n,int m,int lo);
void table_untwist(cplx *cplx_x,int n,int m,int lo);
void vector_twist(cplx_ptr *cplx_x,int n,int m,int lo);
void vector_untwist(cplx_ptr *cplx_x,int n,int m,int lo);

void init_table();
void init_table_vctr();
void naive_real_mul(ring_t *r, const ring_t *x, const ring_t *y);
void naive_complex_mul(ring_t *r, const ring_t *x, const ring_t *y);
void split_radix_mul(ring_t *r, const ring_t *x, const ring_t *y);
void normal_fft_mul(ring_t *r, const ring_t *x, const ring_t *y);
void sr_precomp_mul(ring_t *r, const ring_t *x, const ring_t *y);
void sr_vector_mul(ring_t *r, const ring_t *x, const ring_t *y);

#endif