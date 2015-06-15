#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>

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

int toInt(unsigned char* bytes) {
    return (int)(((unsigned char)bytes[3] << 24) |
                 ((unsigned char)bytes[2] << 16) |
                 ((unsigned char)bytes[1] << 8) |
                 (unsigned char)bytes[0]);
}

int random_int(){
  int length = 4;
  unsigned char x[length];
  randombytes(x,length);
  int var = toInt(x);
  return abs(var);
}

int main(){
  for(int i=0; i < 100;++i)
    for(int j=0; j<1024;++j){
       printf("%d\n",random_int());
    }

}