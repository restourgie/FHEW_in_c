#ifndef NEGACYCLIC_H
#define NEGACYCLIC_H

void init_negacyc();
void phi_forward(cplx_ptr *x,const ring_t *ring);
void phi_backward(cplx_ptr *x, ring_t *ring);

#endif