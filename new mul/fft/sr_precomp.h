#ifndef SR_PRECOMP_H
#define SR_PRECOMP_H

void sr_precomp(cplx *x,int n,int lo);
void sr_precomp_inverse(cplx *x,int n,int lo);
void table_twist(cplx *cplx_x,int n,int m,int lo);
void table_untwist(cplx *cplx_x,int n,int m,int lo);
void init_table();

#endif