#ifndef SPLIT_RADIX_FAST_H
#define SPLIT_RADIX_FAST_H

void split_radix_fast(double complex *x,int n,int lo);
void split_radix_fast_inverse(double complex *x,int n,int lo);
void fast_twist(double *cplx_x,int n,int m,int lo);

#endif