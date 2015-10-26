#ifndef MUL_H
#define MUL_H

#define CPLXDIM 512
#define REALDIM (2*CPLXDIM)
#define ROOTDIM (2*REALDIM)

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#include <stdint.h>
#include <complex.h>

typedef struct {
  uint32_t v[REALDIM];
} ring_t;

typedef struct
{	
	double *real;
	double *imag;
}cplx_ptr;

#endif