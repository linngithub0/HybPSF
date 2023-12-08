#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"
#include <string.h>
void  gaus_estimate(float **star0,int Ng1,int Ng2,double *mean,
         double *sigma)
{
  int i,j,k,ic=0,nmax;
  double *tmp,r2,r2l,det1,det2,meant,sigmat,Nh;
  double sum;
  Nh=Ng1/2-0.5;
  nmax=Ng1*Ng2; 
  tmp=dvector(0,nmax);
  r2l=Nh*Nh;
  ic=0;
  meant=0;
  for(i=0;i<Ng1;i++){
    det1=i-Nh;
    det1=det1*det1;
    for(j=0;j<Ng2;j++){
      det2=j-Nh;
      det2=det2*det2;
      r2=det1+det2;
      if(r2>=r2l && star0[i][j]!=0){
      //if(r2>=r2l){
  tmp[ic]=star0[i][j];
  meant+=tmp[ic];
  ic++;
      }
    }
  }
  //printf("ng1,ng2=%d %d ,iccc=%d\n",Ng1,Ng2,ic);
    
  if(ic<16){
    printf("to few data to estimate the gausion noise in star\n");
    getchar();  
  }
  meant/=ic;
  sum=0;
  for(i=0;i<ic;i++){
    r2=tmp[i]-meant;
    sum+=r2*r2;
  }
  sum/=(ic-1.);
  *sigma=sqrt(sum);
  *mean=meant;
  //printf("mean=%e,sigma=%e\n",mean,sigma);
  //  return sigma;
  free_dvector(tmp,0,nmax);
    
}
