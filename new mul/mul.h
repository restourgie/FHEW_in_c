#ifndef MUL_H
#define MUL_H

#include <stdint.h>

#define CPLXDIM 8
#define REALDIM (2*CPLXDIM)

typedef struct {
  uint32_t v[REALDIM];
} ring_t;

void naive_real_mul(ring_t *r, const ring_t *x, const ring_t *y);
void naive_complex_mul(ring_t *r, const ring_t *x, const ring_t *y);
void smart_complex_mul(ring_t *r, const ring_t *x, const ring_t *y);
void fft_mul(ring_t *r, const ring_t *x, const ring_t *y);

#endif