#ifndef DISTRIB_NEW_H
#define DISTRIB_NEW_H

#include "params.h"

typedef struct {
  double std_dev;
  int max; // Size of table. Set to 0 to ignore table and use rejection sampling
  int offset;
  const float* table; // CDF of Gaussian of standard deviation std_dev centered around offset
} Distrib;


int Sample(const Distrib Chi);

// Distribution of std dev 1.2
extern const float Chi1_Table[23];

extern const Distrib Chi1;

extern const float NoTable[1];

extern const Distrib Chi2;

extern const Distrib Chi3;

extern const float Binary_Table[3];

extern const Distrib Chi_Binary;

#endif