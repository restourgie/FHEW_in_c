#ifndef SUPPORT_H
#define SUPPORT_H

#include <stdint.h>
#include <complex.h>

#define CPLXDIM 8
#define REALDIM (2*CPLXDIM)

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

typedef struct {
  uint32_t v[REALDIM];
} ring_t;

void print_complex(const double complex *a, int N);
void print_double(const double *a, int N);
void to_complex(const ring_t *x, double complex *cplx_x);
void to_real(const double complex *cplx_x, ring_t *x);
void twist(double complex *cplx_x,int n,int m,int lo);
void untwist(double complex *cplx_x,int n,int m,int lo);

#endif