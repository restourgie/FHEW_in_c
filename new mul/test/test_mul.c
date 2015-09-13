#include <stdio.h>
#include "../mul.h"
#include <stdbool.h>

#define NTESTS 1000
#define THRESHOLD 2

uint32_t myabs(uint32_t a)
{
  return (a>(1UL<<31)?-a:a);
}

void rand_test(){
  FILE *urandom = fopen("/dev/urandom", "r");
  ring_t r,re,x,y;
  int n,i;
  int success = 0;
  
  for(n=0;n<NTESTS;n++)
  { 

    for(i=0;i<REALDIM;++i){
      x.v[i] = fgetc(urandom);
      y.v[i] = fgetc(urandom);
    }


    // x.v[0] = 58;
    // x.v[1] =37;
    // x.v[2] =238;
    // x.v[3] =155;
    // x.v[4] =1;
    // x.v[5] =153;
    // x.v[6] =239;
    // x.v[7] =69;
    // x.v[8] =239;
    // x.v[9] =152;
    // x.v[10] =60;
    // x.v[11] =43;
    // x.v[12] =148;
    // x.v[13] = 235;
    // x.v[14] = 149;
    // x.v[15] = 77;

    // y.v[0] = 58;
    // y.v[1] =37;
    // y.v[2] =238;
    // y.v[3] =155;
    // y.v[4] =1;
    // y.v[5] =153;
    // y.v[6] =239;
    // y.v[7] =69;
    // y.v[8] =239;
    // y.v[9] =152;
    // y.v[10] =60;
    // y.v[11] =43;
    // y.v[12] =148;
    // y.v[13] = 235;
    // y.v[14] = 149;
    // y.v[15] = 77;

    // normal_FFT_mul(&re,&x,&y);
    split_radix_FFT_mul(&re,&x,&y);
    // twisted_FFT_mul(&r,&x,&y);
    naive_cyclic_real_mul(&r,&x,&y);
    // smart_complex_mul(&r,&x,&y);
    // naive_complex_mul(&r,&x,&y);
    // naive_real_mul(&r,&x,&y);
    bool error = false;
    
    for(i=0;i<REALDIM;i++)
    {

      if(myabs(r.v[i] - re.v[i]) > THRESHOLD){
        printf("school: %u \n", re.v[i]);
        printf("mine: %u \n",r.v[i]);
        printf("difference: %u\n\n",myabs(r.v[i] - re.v[i]));
        error = true;
      }
      
    }
    if(!error)
      success++;

  }
  printf("Ammount of successful multiplications: %d\n", success);

  fclose(urandom);

}

int main()
{
  rand_test();

  return 0;
}