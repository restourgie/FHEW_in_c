#ifndef SR_PRECOMP_H
#define SR_PRECOMP_H

void fft_precompsr_forward(cplx *x);
void fft_precompsr_backward(cplx *x,ring_t *r);
void init_table();
void destruct_table();

#endif