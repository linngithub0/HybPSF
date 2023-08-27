#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"
void krig(double **x,double *y,int npt,int ndim,double alph,double beta,
	  double nug,double **vo,double *voy, int *indx,double *err, int fit)
{
  void ludcmp(double **a, int n, int *indx, double *d);
  void lubksb(double **a, int n, int *indx, double *b);
  
  double nugsq, d,r,betat;
  int i,j,k;
  nugsq=nug*nug;
  betat=0.5*beta;
  for (i=0;i<npt;i++) {
    voy[i] = y[i];
    for (j=i;j<npt;j++) {
      r=0;
      for(k=0;k<ndim;k++)r += DSQR(x[i][k]-x[j][k]);
      vo[i][j] = vo[j][i] = nugsq+alph*pow(r,betat);
    }
    vo[i][npt] = vo[npt][i] = 1.;
  }
  vo[npt][npt] = voy[npt] = 0.;
  if (fit) for (i=0;i<npt;i++) vo[i][i] -= DSQR(err[i]);
  ludcmp(vo,npt+1,indx,&d);
  lubksb(vo,npt+1,indx,voy);
  
}

void krig_interp(double **x,double **vo,double *voy,int *indx,float *xt,
		 int npt,int ndim,double alph,double beta,double nug,
		 double *out,int esign)
{
  
  void lubksb(double **a, int n, int *indx, double *b);
  int i,k;
  double *vstar,r,nugsq,*dstar,lastval,lasterr,betat;
  nugsq=nug*nug;
  betat=0.5*beta;
  vstar=dvector(0,npt);
  for (i=0;i<npt;i++) {
    r=0;
    for(k=0;k<ndim;k++) r += DSQR(xt[k]-x[i][k]);
    vstar[i] = nugsq+alph*pow(r,betat);
  }
  vstar[npt] = 1.;
  lastval = 0.;
  for (i=0;i<=npt;i++) lastval += voy[i]*vstar[i];
  out[0]=lastval;
  if(esign){
    dstar=dvector(0,npt);
    for (i=0;i<=npt;i++)dstar[i]=vstar[i];    
    lubksb(vo,npt+1,indx,dstar);
    lasterr = 0;
    for (i=0;i<=npt;i++) lasterr += dstar[i]*vstar[i];
    if (lasterr<0)lasterr=0;
    out[1] = sqrt(lasterr);
    free_dvector(dstar,0,npt);
  }
  free_dvector(vstar,0,npt);
  
  
}

double Powvgram(double **x,double *y,int npt,int ndim,double beta,double nug)
{
  int i,j,k;
  double rb,num=0.,denom=0.,betat,nugsq,alph;
  betat=beta*0.5;
  nugsq=nug*nug;
  for (i=0;i<npt;i++) for (j=i+1;j<npt;j++) {
      rb=0;
      for(k=0;k<ndim;k++)rb += DSQR(x[i][k]-x[j][k]);
      rb = pow(rb,betat);
      num += rb*(0.5*DSQR(y[i]-y[j]) - nugsq);
      denom += DSQR(rb);
    }
  alph = num/denom;
  return alph;
}

