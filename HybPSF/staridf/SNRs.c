#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"
float SNRs_estimate(float **image,int nx,int ny,float gain){
  void  gaus_estimate(float **star0,int Ng1,int Ng2,double *mean,
         double *sigma);
  int i,j,k;
  double mean,sigma,snrs,det1,det2;
  double sum1,sum2;
  float **obj,**error;
  error=matrix(0,nx,0,ny);
  obj=matrix(0,nx,0,ny);
  for(i=0;i<nx;i++)for(j=0;j<ny;j++){
    obj[i][j]=image[i][j];
  }
  gaus_estimate(obj,nx,ny,&mean,&sigma);
  //printf("mean=%f,sigma=%f\t",mean,sigma );
  for(i=nx/2-1;i<nx-nx/2;i++){
    for(j=ny/2-1;j<ny-ny/2;j++){
      obj[i][j]-=mean;
      det1=sigma*sigma*gain*gain;
      det2=gain*fabs(obj[i][j]);
      if(det2<det1){error[i][j]=sigma;}
      else{error[i][j]=sqrt(det1+det2)/gain;}
    }
  }
  sum1=sum2=0.;
  for(i=nx/2-1;i<nx-nx/2;i++){
    for(j=ny/2-1;j<ny-ny/2;j++){
      sum1+=obj[i][j];
      sum2+=error[i][j]*error[i][j];
    }
  }//printf("sum=%f,sum2=%f\n",sum1,sqrt(sum2) );
  snrs=(float)sum1/sqrt(sum2);
  return(snrs);
}


void  gaus_estimate(float **star0,int Ng1,int Ng2,double *mean,
         double *sigma)
{
  int i,j,k,ic=0,nmax;
  double *tmp,r2,r2l,det1,det2,meant,Nh;
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
      if(r2>=r2l ){
  tmp[ic]=(double)star0[i][j];
  meant+=tmp[ic];
  ic++;
      }
    }
  }
  //printf("ng1,ng2=%d %d ,iccc=%d\n",Ng1,Ng2,ic);
    free_dvector(tmp,0,nmax);
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
  //ooprintf("mean=%e,sigma=%e\n",mean,sigma);
  //  return sigma;
    
}