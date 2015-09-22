#ifndef NORMAL_FFT_H
#define NORMAL_FFT_H

void recursive_FFT(double complex *x,int n,int lo,double complex root);
void inverse_FFT(double complex *x,int n,int lo,double complex root);

#endif