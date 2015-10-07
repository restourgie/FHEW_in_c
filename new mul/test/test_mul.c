#include <stdio.h>
#include "../mul.h"
// #include "../support.h"
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
  init_table();
  init_table_vctr();

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
    // x.v[8] = 66;
    // x.v[9] =336;
    // x.v[10] =188;
    // x.v[11] =277;
    // x.v[12] =427;
    // x.v[13] =485;
    // x.v[14] =456;
    // x.v[15] =74;
    // x.v[16] =239;
    // x.v[17] =152;
    // x.v[18] =60;
    // x.v[19] =43;
    // x.v[20] =148;
    // x.v[21] = 235;
    // x.v[22] = 149;
    // x.v[23] = 77;
    // x.v[24] =40;
    // x.v[25] =437;
    // x.v[26] =367;
    // x.v[27] =223;
    // x.v[28] =396;
    // x.v[29] = 105;
    // x.v[30] = 324;
    // x.v[31] = 285;

    // y.v[0] = 58;
    // y.v[1] =37;
    // y.v[2] =238;
    // y.v[3] =155;
    // y.v[4] =1;
    // y.v[5] =153;
    // y.v[6] =239;
    // y.v[7] =69;
    // y.v[8] = 66;
    // y.v[9] =336;
    // y.v[10] =188;
    // y.v[11] =277;
    // y.v[12] =427;
    // y.v[13] =485;
    // y.v[14] =456;
    // y.v[15] =74;
    // y.v[16] =239;
    // y.v[17] =152;
    // y.v[18] =60;
    // y.v[19] =43;
    // y.v[20] =148;
    // y.v[21] = 235;
    // y.v[22] = 149;
    // y.v[23] = 77;
    // y.v[24] =40;
    // y.v[25] =437;
    // y.v[26] =367;
    // y.v[27] =223;
    // y.v[28] =396;
    // y.v[29] = 105;
    // y.v[30] = 324;
    // y.v[31] = 285;

    // normal_fft_mul(&r,&x,&y);
    split_radix_mul(&re,&x,&y);
    sr_vector_mul(&r,&x,&y);
    //sr_precomp_mul(&r,&x,&y);
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