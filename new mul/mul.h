#ifndef MUL_H
#define MUL_H

#include <stdint.h>

#define CPLXDIM 512
#define REALDIM (2*CPLXDIM)

typedef struct {
  uint32_t v[REALDIM];
} ring_t;

void naive_cyclic_real_mul(ring_t *r, const ring_t *x, const ring_t *y);
void naive_real_mul(ring_t *r, const ring_t *x, const ring_t *y);
void naive_complex_mul(ring_t *r, const ring_t *x, const ring_t *y);
void smart_complex_mul(ring_t *r, const ring_t *x, const ring_t *y);
void twisted_FFT_mul(ring_t *r, const ring_t *x, const ring_t *y);
void split_radix_FFT_mul(ring_t *r, const ring_t *x, const ring_t *y);
void normal_FFT_mul(ring_t *r, const ring_t *x, const ring_t *y);

#endif