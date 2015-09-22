#ifndef TWISTED_FFT_H
#define TWISTED_FFT_H

void twisted_recursive_FFT(double complex *x,int n,int lo);
void twisted_recursive_inverse_FFT(double complex *x,int n,int lo);

#endif