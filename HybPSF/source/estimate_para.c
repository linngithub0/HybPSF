#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"
#include <string.h>

void estimate_para(double **star,double **center,int Nstar,int Ng, int nbound,
       double *para)
{
  void moffat_fit(double **star,int m,int n,double *xcen,
      double *para,double *gamma);
  void moffat_fit_help(double **star,int m,int n,double *xcen,
      double *para,double *gamma);
  double k_sigma_1D(double *data,int np,double *dm, int *indx);
  void sort(int n, float arr[]);
  int Nstart,ic,i,j,k,Nstartt, *indx;
  double tmp,sigmam,**moffat_para,sigmag,xcens[2],**stard,dms[2],gamma[2];
  double paratmp[5],rdm,betam,sum[5];
  float **paratmps;
  FILE *fp;
  Nstart=Nstar;
  if(Nstart<20)Nstart=Nstar;
  moffat_para=dmatrix(0,5,0,Nstar);
  paratmps=matrix(0,3,0,Nstart);
  indx=ivector(0,Nstar);
  stard=dmatrix(0,Ng,0,Ng);
  sigmam=0; 
  Nstartt=0;
  for(ic=0;ic<Nstart;ic++){
    k=0;
    for(i=0;i<Ng;i++){
      for(j=0;j<Ng;j++){
        stard[i][j]=star[ic][k];
  xcens[0]=center[ic][0]-nbound;
  xcens[1]=center[ic][1]-nbound;
        k++;
      }
    }
    moffat_fit(stard,Ng,Ng,xcens,paratmp,gamma);
    if(paratmp[0]>0 && paratmp[1]>0 &&paratmp[2]>0){
      moffat_para[0][Nstartt]=paratmp[0];
      moffat_para[1][Nstartt]=paratmp[1];
      moffat_para[2][Nstartt]=paratmp[2];
      indx[Nstartt]=Nstartt;
      Nstartt++;
    }
  }
  if(Nstartt<=3){
    printf("small number\n");
    sum[0]=sum[1]=sum[2]=sum[3]=0;
    for(ic=0;ic<Nstart;ic++){
      k=0;
      for(i=0;i<Ng;i++){
        for(j=0;j<Ng;j++){
          stard[i][j]=star[ic][k];
          xcens[0]=center[ic][0]-nbound;
          xcens[1]=center[ic][1]-nbound;
          k++;
        }
      }
      moffat_fit_help(stard,Ng,Ng,xcens,paratmp,gamma);//printf("sigma=%f\n", paratmp[0]);
      paratmps[0][ic+1]=paratmp[0];
      paratmps[1][ic+1]=paratmp[1];
      paratmps[2][ic+1]=paratmp[2];
    }
    sort(Nstart, paratmps[0]);
    sort(Nstart, paratmps[1]);
    sort(Nstart, paratmps[2]);//for(ic=1;ic<=Nstart;ic++)printf("paratmps=%f\n",paratmps[2][ic] );
    para[0]=paratmps[0][Nstart/2];
    para[1]=paratmps[1][Nstart/2];
    para[2]=paratmps[2][Nstart/2];
  }
  else{
    printf("normal number\n");
  k_sigma_1D(moffat_para[1],Nstartt,dms,indx);
  Nstartt=indx[Nstartt];
  for(i=0;i<Nstartt;i++){
    j=indx[i];
    moffat_para[0][i]=moffat_para[0][j];    
    moffat_para[1][i]=moffat_para[1][j];   
    moffat_para[2][i]=moffat_para[2][j];
    indx[i]=i;
  }
  k_sigma_1D(moffat_para[2],Nstartt,dms,indx);
  Nstartt=indx[Nstartt];
  for(i=0;i<Nstartt;i++){
    j=indx[i];
    moffat_para[0][i]=moffat_para[0][j];    
    moffat_para[1][i]=moffat_para[1][j];     
    moffat_para[2][i]=moffat_para[2][j];
    indx[i]=i;
  }
  
  k_sigma_1D(moffat_para[0],Nstartt,dms,indx);
  Nstartt=indx[Nstartt];
  sigmam=rdm=betam=0;
  for(i=0;i<Nstartt;i++){
    j=indx[i];
    sigmam+=moffat_para[0][j];    
    rdm   +=moffat_para[1][j];       
    betam +=moffat_para[2][j];
  }
  
  //printf("mean=%e sigma=%e\n",dm[0],dm[1]);
  para[0]=sigmam/Nstartt;
  para[1]=rdm/Nstartt;
  para[2]=betam/Nstartt;
  }
  if(isnormal(para[0])==0 || para[1]<=0 || para[2]<=0){
    printf("usual case\n");
    float **stack;stack=matrix(0,Ng,0,Ng);
    stack[i][j]=0.;
    for(ic=0;ic<Nstar;ic++){
      k=0;
      for(i=0;i<Ng;i++)for(j=0;j<Ng;j++){
        stack[i][j]+=star[ic][k];k++;
      }
    }
    for(i=0;i<Ng;i++)for(j=0;j<Ng;j++){
      stack[i][j]/=Nstar;
    }
    void size(float **I,int Ng,double *cen,double rg,double *se);
    double cen[2],paras[5];
    cen[0]=cen[1]=(Ng+1)/2.;
    size(stack,Ng,cen,3.5,paras);
    para[0]=paras[0];
    free_matrix(stack,0,Ng,0,Ng);
  }
  
  free_dmatrix(moffat_para,0,5,0,Nstar);
  free_dmatrix(stard,0,Ng,0,Ng);
  free_ivector(indx,0,Nstar);
  free_matrix(paratmps,0,3,0,Nstart);

}

//rg=1.25*para[0]; ierror=1; double gamma[2],gammaerr[2]; noise=matrix(0,Nstar,0,Mp); nbound=0; double para[5]; gamma[0],gamma[1]是输出的椭率

////////////////////////

//se[0] 为输出 大小



