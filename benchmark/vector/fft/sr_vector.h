#ifndef SR_VECTOR_H
#define SR_VECTOR_H

void init_table_vctr();
void destruct_table_vctr();
void fft_vector_forward(cplx_ptr *x,const ring_t *ring);
void fft_vector_backward(cplx_ptr *cplx_x,ring_t *res);

#endif