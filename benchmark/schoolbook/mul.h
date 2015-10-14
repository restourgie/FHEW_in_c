#ifndef MUL_H
#define MUL_H

#define CPLXDIM 512
#define REALDIM (2*CPLXDIM)
#define ZEROPAD (2*REALDIM)

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#define W(N,k) (cexp(2.0 * M_PI * I * (double)k / (double) N))

#include <stdint.h>
#include <complex.h>

typedef struct {
  uint32_t v[REALDIM];
} ring_t;

void print_complex(const double complex *a, int N);

void tangent_mul(ring_t *r, const ring_t *x, const ring_t *y);
void naive_real_mul(ring_t *r, const ring_t *x, const ring_t *y);
void split_radix_mul(ring_t *r, const ring_t *x, const ring_t *y);
void normal_fft_mul(ring_t *r, const ring_t *x, const ring_t *y);
void twisted_mul(ring_t *r, const ring_t *x, const ring_t *y);
#endif