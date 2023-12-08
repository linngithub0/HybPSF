#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"
void Rebasis(double **basisf0,int nbf0,int Mp,double **basisf,double **coeff,int *nbf)
{
  double *phi,*xj,**xij,*cj,absphi,phi0,*phit,detphi,*weight,**basisft,sum,det;
  int i,j,k,niter,nmax=5000;
  phi0=1./Mp;
  phi=dvector(0,Mp);  
  phit=dvector(0,Mp);
  xij=dmatrix(0,nbf0,0,Mp);
  basisft=dmatrix(0,nbf0,0,Mp);
  cj=dvector(0,nbf0);
  weight=dvector(0,nbf0);
  for(i=0;i<Mp;i++)for(j=0;j<nbf0;j++){xij[j][i]=basisf0[j][i];basisft[j][i]=0;}
  for(j=0;j<nbf0;j++)weight[j]=1;
  weight[0]=nbf0/2;
  k=0;
  do{
    for(i=0;i<Mp;i++)phi[i]=phit[i]=basisf0[k][i];
    niter=0;
    do{
      //E-step
      for(i=0;i<nbf0;i++){
	cj[i]=0;
	for(j=0;j<Mp;j++)cj[i]+=xij[i][j]*phi[j];
      }
      //Mstep;
      absphi=0;
      for(j=0;j<Mp;j++){
	phi[j]=0;
	for(i=0;i<nbf0;i++)phi[j]+=weight[i]*cj[i]*xij[i][j];
	absphi+=phi[j]*phi[j];
      }
      //Renormalise
      absphi=sqrt(absphi);
      detphi=0;
      for(j=0;j<Mp;j++){
	phi[j]/=absphi;
	detphi+=(phi[j]-phit[j])*(phi[j]-phit[j]);
	phit[j]=phi[j];
      }
      detphi=sqrt(detphi);
      niter++;
      
    }while(niter<nmax & detphi>1.e-14);
    //printf("nbf0=%d k=%d niter=%d\n",nbf0,k,niter);
    det=0;
    for(i=0;i<nbf0;i++){
      coeff[i][k]=cj[i];
      sum=0;
      for(j=0;j<Mp;j++){
	xij[i][j]-=cj[i]*phi[j];
	basisft[i][j]+=cj[i]*phi[j];
	sum+=basisft[i][j]*basisf0[i][j];
      }
      //printf("%e ",sum);
      sum=fabs(1-sum);
      if(det<sum)det=sum;
    }
    //printf("\n");
    for(j=0;j<Mp;j++)basisf[k][j]=phi[j];
    weight[k]=1;
    weight[k+1]=nbf0/2;
    k++;
  }while(k<nbf0 & det>0.0000001);
  *nbf=k;
  free_dmatrix(basisft,0,nbf0,0,Mp);
  free_dvector(weight,0,nbf0);
  free_dvector(phi,0,Mp); 
  free_dvector(phit,0,Mp);
  free_dmatrix(xij,0,nbf0,0,Mp);
  free_dvector(cj,0,nbf0);
}

