#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"

void krig_prepare(float **x,float **y,float **error,double ***vo,double **voy,
      int **indx,double *alpha,int npt,int nbf)
                  
{
  double Powvgram(double **x,double *y,int npt,int ndim,double beta,double nug);
  void krig(double **x,double *y,int npt,int ndim,double alph,double beta,
            double nug,double **vo,double *voy, int *indx,double *err, int fit);
  double **xtmp,**ytmp,**err,beta=1.5,nug=0;
  int i,j,ndim=2,fitsign=1;
  xtmp=dmatrix(0,npt-1,0,1);
  ytmp=dmatrix(0,nbf-1,0,npt);
  err=dmatrix(0,nbf-1,0,npt);
  for(i=0;i<npt;i++){
    for(j=0;j<2;j++)xtmp[i][j]=(double)x[i][j];
    for(j=0;j<nbf;j++){ytmp[j][i]=(double)y[i][j];err[j][i]=(double)error[i][j];}
  }
  for(i=0;i<nbf;i++){
    alpha[i]=Powvgram(xtmp,ytmp[i],npt,ndim,beta,nug);
    krig(xtmp,ytmp[i],npt,ndim,alpha[i],beta,nug,vo[i],voy[i],indx[i],
	 err[i],fitsign);
  }
  
  free_dmatrix(xtmp,0,npt-1,0,1);
  free_dmatrix(ytmp,0,nbf-1,0,npt);
  free_dmatrix(err,0,nbf-1,0,npt);

}

void creat_star(float **basisf,float **x,double ***vo,double **voy,int **indx,
		double *alpha,int npt, int nbf,int ng,int Mp,
		float *xt,float **out_star)
{
  void krig_interp(double **x,double **vo,double *voy,int *indx,float *xt,
                   int npt,int ndim,double alph,double beta,double nug,
                   double *out,int esign);
  double **xtmp,*coff,beta=1.5,out[2],nug=0,sum;
  int i,j,k=0,kt,ndim=2,esign=0;
  if(ng*ng!=Mp){
    printf("The grid in basis fucntions are different from that in star\n");
    getchar();
  }
  coff=dvector(0,nbf-1);
  xtmp=dmatrix(0,npt-1,0,1);
  for(i=0;i<npt;i++){
    for(j=0;j<2;j++)xtmp[i][j]=x[i][j];
  }
  for(i=0;i<nbf;i++){
    krig_interp(xtmp,vo[i],voy[i],indx[i],xt,npt,ndim,alpha[i],beta,nug,
		out,esign);
    coff[i]=out[0];
  }
  kt=0;
  sum=0;
  for(i=0;i<ng;i++){
    for(j=0;j<ng;j++){
      out_star[i][j]=0;
      for(k=0;k<nbf;k++)out_star[i][j] += basisf[k][kt]*coff[k];
      kt++;
      sum +=out_star[i][j];
    }
  }
  //printf("sum=%e\n",sum);
  for(i=0;i<ng;i++){
    for(j=0;j<ng;j++){
      out_star[i][j]/=sum;
    }
  }
  free_dvector(coff,0,nbf-1);
  free_dmatrix(xtmp,0,npt-1,0,1);
  

}
