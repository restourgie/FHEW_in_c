#ifndef RING_H
#define RING_H

#include <stdint.h>

typedef struct {
  uint32_t v[1024];
} ring_t;

#define CPLXDIM 8
#define REALDIM (2*CPLXDIM)

void ring_mul(ring_t *r, const ring_t *x, const ring_t *y);
void exact_mul(ring_t *r, const ring_t *x, const ring_t *y);
void complex_mul(ring_t *r, const ring_t *x, const ring_t *y);
void ring_print(const ring_t *x);

#endif
