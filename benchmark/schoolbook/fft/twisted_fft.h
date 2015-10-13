#ifndef TWISTED_FFT_H
#define TWISTED_FFT_H

void twisted_recursive_inverse(double complex *x,int n,int lo);
void twisted_recursive(double complex *x,int n,int lo);

#endif