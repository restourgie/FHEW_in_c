#include <stdio.h>
#include "../ring.h"
#include <stdbool.h>
#include <assert.h>

#define NTESTS 1
#define THRESHOLD 2

void wait_on_enter(){
  printf("\nPress ENTER to continue\n");
  getchar();
}

uint32_t myabs(uint32_t a)
{
  return (a>(1UL<<31)?-a:a);
}

void ones_test(){
  ring_t r,x,y;

  for(int i=0;i<REALDIM;++i){
    x.v[i] = 0;
    y.v[i] = 0;
  }

  for(int i=0;i<REALDIM;++i){
    x.v[i] = 1;
    for(int j=0;j<REALDIM;++j){
      y.v[j] = 1;
      ring_mul(&r,&x,&y);
      if((i+j) < REALDIM){
        assert(r.v[i+j] == 1);
        // printf("%u\n",r.v[i+j]);
      }
      else
      {
        assert(r.v[i+j-REALDIM] == (unsigned int) -1);
        // printf("%u\n",r.v[i+j-1024]);
      }
      y.v[j] = 0;
    }
    x.v[i] = 0;
  }
}

void rand_test(){
  FILE *urandom = fopen("/dev/urandom", "r");
  ring_t r,re,x,y;
  int n,i;
  int success = 0;
  
  for(n=0;n<NTESTS;n++)
  { 
    // fread(x.v,sizeof(uint32_t),REALDIM,urandom);
    // fread(y.v,sizeof(uint32_t),REALDIM,urandom);
    // for(i=0;i<REALDIM;i++){
    // //   // printf("y.v[%d] = %d\n",i,y.v[i]);
    // //   //y.v[i] &= 0x7ff; gets 957 of 1000 success
    //   y.v[i] &= 0x3ff;
    // //   // printf("new value y.v[%d] = %d\n",i,y.v[i]);
    // //   //wait_on_enter();
    // }

    for(i=0;i<REALDIM;++i){
      x.v[i] = fgetc(urandom);
      y.v[i] = fgetc(urandom);
    }

    // ring_mul(&r,&x,&y);
    complex_mul(&r,&x,&y);
    exact_mul(&re,&x,&y);
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
  //ones_test();

  return 0;
}
