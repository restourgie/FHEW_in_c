#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

const float Chi1_Table[23] = {
1.12011750313263e-14, 2.38717233762211e-12, 3.04966971020178e-10,
2.34394541319773e-8, 1.08538196465647e-6, 0.0000303513786306856,
0.000514575939439740, 0.00532464699317562, 0.0340111330223921,
0.136723892128727, 0.357520614142345, 0.642479385857655,
0.863276107871273, 0.965988866977608, 0.994675353006824,
0.999485424060560, 0.999969648621369, 0.999998914618035,
0.999999976560546, 0.999999999695033, 0.999999999997613,
0.999999999999989, 1.00000000000000};

const Distrib Chi1 = {
  1.4,
  23,
  11,
  Chi1_Table
};

const float NoTable[1] = {1};

const Distrib Chi2 = {
  (double) (1 << 17),
  0,
  0,
  NoTable
};

const Distrib Chi3 = {
  (double) 6,
  0,
  0,
  NoTable
};

const float Binary_Table[3] = {
 .25,
 .75,
1.0 };

const Distrib Chi_Binary = {
  0,
  3,
  1,
  Binary_Table
};

static int fd = -1;

void randombytes(unsigned char *x,unsigned long long xlen)
{
  int i;

  if (fd == -1) {
    for (;;) {
      fd = open("/dev/urandom",O_RDONLY);
      if (fd != -1) break;
      sleep(1);
    }
  }

  while (xlen > 0) {
    if (xlen < 1048576) i = xlen; else i = 1048576;

    i = read(fd,x,i);
    if (i < 1) {
      sleep(1);
      continue;
    }

    x += i;
    xlen -= i;
  }
}
//25% chance on -1 50% on 0 and 25% on 1
int some_numbers(const Distrib Chi){
 
  unsigned char x[1];
  randombytes(x,1);
  double r = x[1]/255;
  for (int i = 0; i < Chi.max; ++i) 
      if (r <= Chi.table[i]) 
        return i - 1;
}