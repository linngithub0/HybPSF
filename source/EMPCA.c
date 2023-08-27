#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"
void EMPCA(double **star,int Nstar,int Mp,double **basisf,double **coeff,double **coefferr,
     double **weight,int nbf)
{
  void EC(double *xj,double **phi,double *weight,int Mp, int nbf,double *coeff);
  void PCA(double **star,int Nstar,int Mp,double **basef,double **coff);
  float ran1(long *idum);
  double *phi,*xj,**xij,*cj,absphi,phi0,*phit,detphi,ccj,sum,chi20=0,chi2,detchi2;
  int i,j,k,niter,nmax=1000,k0;
  long iseed=-1;
  phi0=1./Mp;
  phi=dvector(0,Mp);  
  phit=dvector(0,Mp);
  xij=dmatrix(0,Nstar,0,Mp);
  cj=dvector(0,Nstar);printf("here Nstar=%d,Mp=%d\n",Nstar,Mp);
  for(i=0;i<Mp;i++){
    phi[i]=phit[i]=phi0;
    for(j=0;j<Nstar;j++)xij[j][i]=star[j][i];
  }

  /*int write_fits_3D(char *argv,float ***stars,int *dim);
  int dim[3],ng=sqrt(Mp);printf("ng=%d\n",ng);
  float ***write;write=f3tensor(0,Nstar,0,ng,0,ng);
  for(k=0;k<Nstar;k++){
    for(i=0;i<ng;i++)for(j=0;j<ng;j++){
      write[k][i][j]=(float)weight[k][i*ng+j];
    }
  }
  dim[0]=dim[1]=ng;dim[2]=Nstar;
  write_fits_3D("fits/PCA.fits",write,dim);*/

  printf("PCA\n");
  for(i=0;i<nbf;i++){
    for(j=0;j<Mp;j++){
      basisf[i][j]=1.+ran1(&iseed);
    }
  }
  niter=0;
  do{
    for(i=0;i<Mp;i++)for(j=0;j<Nstar;j++)xij[j][i]=star[j][i];
    
    //E-step
    for(i=0;i<Nstar;i++){
      EC(xij[i],basisf,weight[i],Mp,nbf,coeff[i]);
    }
    //Mstep;
    for(k=0;k<nbf;k++){
      for(j=0;j<Mp;j++){
  basisf[k][j]=0;ccj=0;
  for(i=0;i<Nstar;i++){
    basisf[k][j]+=weight[i][j]*coeff[i][k]*xij[i][j];
    ccj+=weight[i][j]*coeff[i][k]*coeff[i][k];
  }
  basisf[k][j]/=ccj;  
      }
      for(j=0;j<Mp;j++)for(i=0;i<Nstar;i++)xij[i][j]-=basisf[k][j]*coeff[i][k]; 
    }
 
    //re-orthognal
    for(k0=0;k0<nbf-1;k0++){
      absphi=0;
      for(j=0;j<Mp;j++)absphi+=basisf[k0][j]*basisf[k0][j]; //normalization
      absphi=sqrt(absphi);
      for(j=0;j<Mp;j++)basisf[k0][j]/=absphi;
      for(k=k0+1;k<nbf;k++){
  sum=0;
  for(j=0;j<Mp;j++)sum+=basisf[k0][j]*basisf[k][j];
  for(j=0;j<Mp;j++)basisf[k][j]-=basisf[k0][j]*sum;
  
      }
    }
    //renormalization the last one
    for(k=0;k<nbf;k++){
    absphi=0;
    for(j=0;j<Mp;j++)absphi+=basisf[k][j]*basisf[k][j];
    absphi=sqrt(absphi);
    for(j=0;j<Mp;j++)basisf[k][j]/=absphi;
    }
    //calculate chi^2
    chi2=0;
    for(j=0;j<Mp;j++){
      for(i=0;i<Nstar;i++){
  xij[i][j]=0;
  for(k=0;k<nbf;k++)xij[i][j]+=basisf[k][j]*coeff[i][k];
  chi2+=(xij[i][j]-star[i][j])*(xij[i][j]-star[i][j])*weight[i][j];
      }
    }
    detchi2=chi20/chi2-1;
    detchi2=fabs(detchi2);
    chi20=chi2;
    niter++;
    printf("niter=%d chi2=%e detchi2=%e\n",niter,detchi2,chi2);
    
    //}while(niter<300);
  }while (niter<nmax && detchi2>1.e-14);
  free_dvector(phi,0,Mp); 
  free_dvector(phit,0,Mp);
  free_dmatrix(xij,0,Nstar,0,Mp);
  free_dvector(cj,0,Nstar);
}


