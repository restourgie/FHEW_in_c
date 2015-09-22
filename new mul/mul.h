#ifndef MUL_H
#define MUL_H

#include "support.h"

void naive_real_mul(ring_t *r, const ring_t *x, const ring_t *y);
void naive_complex_mul(ring_t *r, const ring_t *x, const ring_t *y);
void split_radix_mul(ring_t *r, const ring_t *x, const ring_t *y);
void split_radix_fast_mul(ring_t *r, const ring_t *x, const ring_t *y);

#endif