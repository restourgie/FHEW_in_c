#ifndef SPLIT_RADIX_FFT_H
#define SPLIT_RADIX_FFT_H

void split_radix_recursive(double complex *x,int n,int lo);
void split_radix_recursive_inverse(double complex *x,int n,int lo);
void twist(double complex *cplx_x,int n,int m,int lo);
void untwist(double complex *cplx_x,int n,int m,int lo);

#endif