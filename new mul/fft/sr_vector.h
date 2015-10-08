#ifndef SR_VECTOR_H
#define SR_VECTOR_H

// void sr_vector(cplx_ptr *x,int n,int lo);
// void sr_vector_inverse(cplx_ptr *x,int n,int lo);
// void vector_twist(cplx_ptr *cplx_x,int n,int m,int lo);
// void vector_untwist(cplx_ptr *cplx_x,int n,int m,int lo);
void init_table_vctr();
void fft_vector_forward(cplx_ptr *cplx_x, const ring_t *x);
void fft_vector_backward(cplx_ptr *cplx_x,ring_t *x);

#endif