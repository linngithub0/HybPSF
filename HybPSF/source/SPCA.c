#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"
#include <string.h>
void SPCA(double **images,double **weight,double **basisf,int Mp,int Nstar,int nbf,
		  double **compbf,double **coeff,int Nb){
  void EC(double *xj,double **phi,double *weight,int Mp, int nbf,double *coeff);
  void SPrepare_Dmatrix(double **weight,double **basisf,int Mp,int nimg,int nbf,
		       double ***sigmaka);
  void Screat_Dmatrix(double ***sigmaka,int Mp,int nimg,int nbf,double *coeff,
		     double **a,double *b,double **weight,double **images,double **basisf);
  void lubksb(double **a, int n, int *indx, double b[]);
  void ludcmp(double **a, int n, int *indx, double *d);
  void Rebasis(double **basisf0,int nbf0,int Mp,double **basisf,double **coeff,int *nbf);
  float ran1(long *idum);
  
  int i,j,k,ipc,nbfs;
  long iseed=-1;
  double **xij,**coff;
  double **a,*b,d,*compbf0,*D,sum,sum1,sum2=0,detbf,sumt,**basis;
  double ***sigmaka,**sigmakb,**rbasisf, **rcoeff,detchi2,absphi;
  int *indx,niter=0,nmax=1000,rnbf,k0;
  double **star;

  xij=dmatrix(0,Nstar,0,Mp);
  a=dmatrix(0,nbf,0,nbf);
  b=dvector(0,nbf);
  D=dvector(0,nbf);
  indx=ivector(0,nbf);
  sigmaka=d3tensor(0,Nstar,0,nbf,0,nbf); 
  sigmakb=dmatrix(0,Nstar,0,nbf);
  coff=dmatrix(0,Nstar,0,Nb);
  rbasisf=dmatrix(0,nbf,0,Mp);
  rcoeff=dmatrix(0,nbf,0,nbf);
  star=dmatrix(0,Nstar,0,Mp);
  basis=dmatrix(0,nbf,0,Mp);


  int write_fits_3D(char *argv,float ***stars,int *dim);
  int dim[3],ng=sqrt(Mp);printf("ng=%d\n",ng);
  float ***write;write=f3tensor(0,Nstar,0,ng,0,ng);
  for(k=0;k<Nstar;k++){
    for(i=0;i<ng;i++)for(j=0;j<ng;j++){
      write[k][i][j]=(float)weight[k][i*ng+j];
    }
  }
  dim[0]=dim[1]=ng;dim[2]=Nstar;
  //write_fits_3D("fits/PCAs.fits",write,dim);
  
  for(i=0;i<Nb;i++){
    for(j=0;j<Mp;j++){
      compbf[i][j]=1.+ran1(&iseed);
    }
  }

  for(i=0;i<nbf;i++){
    for(k=0;k<Mp;k++)basis[i][k]=basisf[i][k];
  }
  nbfs=nbf;
  Rebasis(basis,nbfs,Mp,rbasisf,rcoeff,&rnbf);
  printf("nbfs=%d rnbf=%d\n",nbfs,rnbf);
  SPrepare_Dmatrix(weight,rbasisf,Mp,Nstar,rnbf,sigmaka);
  
  
  niter=0;
  do{
    //E-step
    for(i=0;i<Nstar;i++){
      for(j=0;j<Mp;j++){
	star[i][j]=xij[i][j]=images[i][j];
      }
    }
    for(i=0;i<Nstar;i++){
      EC(xij[i],compbf,weight[i],Mp,Nb,coff[i]);
    }
    //printf("coeff\n");
    for(i=0;i<Nstar;i++){
      for(j=0;j<Nb;j++)coeff[j][i]=coff[i][j];
    }

    for(ipc=0;ipc<Nb;ipc++){

      Screat_Dmatrix(sigmaka,Mp,Nstar,rnbf,coeff[ipc],a,b,weight,star,rbasisf);
      //printf(" created matrix\n");  
      ludcmp(a, rnbf, indx, &d);
      lubksb(a, rnbf, indx, b);
      sum=0;
      for(k=0;k<Mp;k++){
	      compbf[ipc][k]=0;
	      for(i=0;i<rnbf;i++)compbf[ipc][k]+=rbasisf[i][k]*b[i];
      }

      for(i=0;i<Nstar;i++){
	      for(k=0;k<Mp;k++)star[i][k]-=compbf[ipc][k]*coeff[ipc][i];
      }

    }
    //re-orthognal
    for(k0=0;k0<Nb-1;k0++){
      absphi=0;
      for(j=0;j<Mp;j++)absphi+=compbf[k0][j]*compbf[k0][j]; //normalization
      absphi=sqrt(absphi);
      for(j=0;j<Mp;j++)compbf[k0][j]/=absphi;
      for(k=k0+1;k<Nb;k++){
	    sum=0;
	    for(j=0;j<Mp;j++)sum+=compbf[k0][j]*compbf[k][j];
	    for(j=0;j<Mp;j++)compbf[k][j]-=compbf[k0][j]*sum;
	
      }
    }
    //renormalization the last one
    for(i=0;i<Nb;i++){
    absphi=0;
    for(j=0;j<Mp;j++)absphi+=compbf[i][j]*compbf[i][j];
    absphi=sqrt(absphi);
    for(j=0;j<Mp;j++)compbf[i][j]/=absphi;
    }
    //calculate chi^2
    sum1=0;
    for(j=0;j<Mp;j++){
      for(i=0;i<Nstar;i++){
	xij[i][j]=0;
	for(k=0;k<Nb;k++)xij[i][j]+=compbf[k][j]*coeff[k][i];
	sum1+=(xij[i][j]-images[i][j])*(xij[i][j]-images[i][j])*weight[i][j];
      }
    }
    detchi2=sum2/sum1-1;
    detchi2=fabs(detchi2);
    sum2=sum1;
    niter++;
    printf("niter=%d detchi2=%e   chi2=%e\n",niter,detchi2,sum1);
    
  }while (niter<nmax && detchi2>1.e-12);
  
  free_dmatrix(xij,0,Nstar,0,Mp);
  free_dmatrix(a,0,nbf,0,nbf);
  free_dvector(b,0,nbf);
  free_dvector(D,0,nbf);
  free_ivector(indx,0,nbf);
  free_d3tensor(sigmaka,0,Nstar,0,nbf,0,nbf); 
  free_dmatrix(sigmakb,0,Nstar,0,nbf);
  free_dmatrix(coff,0,Nstar,0,Nb);
  free_dmatrix(rbasisf,0,nbf,0,Mp);
  free_dmatrix(rcoeff,0,nbf,0,nbf);
  free_dmatrix(star,0,Nstar,0,Mp);
  free_dmatrix(basis,0,nbf,0,Mp);
}

void SPrepare_Dmatrix(double **weight,double **basisf,int Mp,int nimg,int nbf,
         double ***sigmaka)
{
  int i,j,k,l;
  double sum;
  for(j=0;j<nbf;j++){
    for(l=0;l<nbf;l++){
      for(i=0;i<nimg;i++){
	sigmaka[i][j][l]=0;
	for(k=0;k<Mp;k++)sigmaka[i][j][l]+=weight[i][k]*basisf[j][k]*basisf[l][k];
      }
    }
    /*for(i=0;i<nimg;i++){
      sigmakb[i][j]=0;
      for(k=0;k<Mp;k++)sigmakb[i][j]+=weight[i][k]*basisf[j][k]*images[i][k];
    }*/ 
  }

}

void Screat_Dmatrix(double ***sigmaka,int Mp,int nimg,int nbf,double *coeff,
       double **a,double *b,double **weight,double **images,double **basisf)
{
  int i,j,k,l;
  double sum,*tmp;
  tmp=dvector(0,Mp);
  for(k=0;k<Mp;k++){
    tmp[k]=0;
    for(i=0;i<nimg;i++)tmp[k]+=weight[i][k]*images[i][k]*coeff[i];
  }

  for(j=0;j<nbf;j++){
    for(l=j;l<nbf;l++){
      a[j][l]=0;
      for(i=0;i<nimg;i++){
	//sum=0;
	//for(k=0;k<Mp;k++)sum+=weight[i][k]*basisf[j][k]*basisf[l][k];
	a[j][l]+=sigmaka[i][j][l]*coeff[i]*coeff[i];
      }
      a[l][j]=a[j][l];
      //a[j][l]*=2;
    }

    b[j]=0;
    for(k=0;k<Mp;k++)b[j]+=basisf[j][k]*tmp[k];
    
    //b[j]*=2;
  }
  free_dvector(tmp,0,Mp);
}
