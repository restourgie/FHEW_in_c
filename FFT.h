#ifndef FFT_H
#define FFT_H

#include "params.h"

void FFTsetup();
void FFTforward(Ring_FFT res, Ring_ModQ val);
void FFTbackward(Ring_ModQ res, Ring_FFT val);

#endif