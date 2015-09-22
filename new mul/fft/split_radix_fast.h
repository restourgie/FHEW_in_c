#ifndef SPLIT_RADIX_FAST_H
#define SPLIT_RADIX_FAST_H

void split_radix_fast(cplx *x,int n,int lo);
void split_radix_fast_inverse(cplx *x,int n,int lo);
void fast_twist(cplx *cplx_x,int n,int m,int lo);
void fast_untwist(cplx *cplx_x,int n,int m,int lo);

#endif