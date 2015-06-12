#ifndef DISTRIB_NEW_H
#define DISTRIB_NEW_H

#include "params.h"

typedef struct {
  double std_dev;
  int max; // Size of table. Set to 0 to ignore table and use rejection sampling
  int offset;
  const float* table; // CDF of Gaussian of standard deviation std_dev centered around offset
} Distrib;


int Sample_1(const Distrib Chi);
int Sample_3(const Distrib Chi);
int random_int();
int Sample_2(const Distrib Chi);

/*********TEMP***********/
// void file_setup();
// void close_files();
// int function1();
// int function2();
// int function3();
// int function4();
// int function5();
// int function6();
/*********TEMP***********/

// Distribution of std dev 1.2
extern const float Chi1_Table[23];

extern const Distrib Chi1;

extern const float NoTable[1];

extern const Distrib Chi2;

extern const Distrib Chi3;

extern const float Binary_Table[3];

extern const Distrib Chi_Binary;

#endif