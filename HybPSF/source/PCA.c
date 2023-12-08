#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"
#include <string.h>
void PCA(float **star,int Nstar,int Mp,double **basef,double **coff)
{
  
  void svdcmp(double **a, int m, int n, double *w, double **v);
  void indexx(int n, double arr[], int *indx);
  double **a,**v,*w,*check;
  int i,j,k,*indx;
  if(Mp>=Nstar){
    a=dmatrix(1,Mp,1,Nstar);
    v=dmatrix(1,Nstar,1,Nstar);
    w=dvector(1,Nstar);
    indx=ivector(1,Nstar);
    for(i=1;i<=Nstar;i++){for(j=1;j<=Mp;j++)a[j][i]=star[i-1][j-1];}
    svdcmp(a,Mp,Nstar,w,v);   
    indexx(Nstar,w,indx);    
    for(j=1;j<=Nstar;j++){
      k=indx[Nstar-j+1];
      for(i=1;i<=Mp;i++)basef[j-1][i-1]=a[i][k];
    }
    for(i=1;i<=Nstar;i++){
      for(j=1;j<=Nstar;j++){
  k=indx[Nstar-j+1];
  coff[i-1][j-1]= w[k]*v[i][k];
      }
    }
    //for(j=1;j<=Nstar;j++){for(i=1;i<=Mp;i++)basef[j-1][i-1]=a[i][j];}
    //for(i=1;i<=Nstar;i++){for(k=1;k<=Nstar;k++)coff[i-1][k-1]= w[k]*v[i][k];}
    
    free_dmatrix(a,1,Mp,1,Nstar);
    free_dmatrix(v,1,Nstar,1,Nstar);
    free_dvector(w,1,Nstar);
    free_ivector(indx,1,Nstar);
  }

  else{
    a=dmatrix(1,Nstar,1,Mp);
    v=dmatrix(1,Mp,1,Mp);
    w=dvector(1,Mp);
    indx=ivector(1,Mp);
    for(i=1;i<=Nstar;i++){for(j=1;j<=Mp;j++)a[i][j]=star[i-1][j-1];}
    svdcmp(a,Nstar,Mp,w,v);   
    indexx(Mp,w,indx);  
    for(j=1;j<=Mp;j++){    
      k=indx[Mp-j+1];
      for(i=1;i<=Mp;i++)basef[j-1][i-1]=v[i][k];
    }
    for(i=1;i<=Nstar;i++){
      for(j=1;j<=Mp;j++){
  k=indx[Mp-j+1];
  coff[i-1][j-1]= w[k]*a[i][k];
      }
    }  
    free_dmatrix(a,1,Nstar,1,Mp);
    free_dmatrix(v,1,Mp,1,Mp);
    free_dvector(w,1,Mp);
    free_ivector(indx,1,Mp);
  }
}
void PCA_err(float **err,int Nstar,int Mp,int nbf,double **basef,double **cofferr)
{
  float gasdev(long *iseed);
  long iseed=-1;
  int nc=100,i,j,k,ns;
  double sum,*mean,**coeffs,*start,*errt,det;
  FILE *fp;
  start=dvector(0,Mp);
  mean=dvector(0,nbf);
  coeffs=dmatrix(0,nc,0,nbf);
  for(ns=0;ns<Nstar;ns++){
    for(k=0;k<nbf;k++)mean[k]=0;
    for(k=0;k<nc;k++){
      for(j=0;j<Mp;j++)start[j]=gasdev(&iseed)*err[ns][j];
      for(i=0;i<nbf;i++){
  coeffs[k][i]=0;
  for(j=0;j<Mp;j++)coeffs[k][i]+=start[j]*basef[i][j];
  mean[i]+=coeffs[k][i];
      }
    }   
    for(k=0;k<nbf;k++)mean[k]/=nc;
    for(i=0;i<nbf;i++){
      sum=0;
      for(k=0;k<nc;k++){
  det=coeffs[k][i]-mean[i];
  sum+=det*det;
      }
      cofferr[ns][i]=sqrt(sum/(nc-1.));
    }
  }
  
  free_dvector(start,0,Mp);
  free_dvector(mean,0,nbf);
  free_dmatrix(coeffs,0,nc,0,nbf);
}
