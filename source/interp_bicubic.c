#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"
#include "function.h"

void centriod_gaus2D(int nx,int ny,int nxt,int nyt,int nbound,
		double *para,float **yt)
{
  int i,j,ic,jc;
  double shift1,shift2,x,y;
  double kc,bgc,sigmac,xc,yc,sigma2v,det;
  kc=para[0];sigmac=para[1];bgc=para[2];xc=para[3];yc=para[4];
  sigma2v=-1./(2.*sigmac*sigmac);
  ic=nx*0.5;
  jc=ny*0.5;
  shift1=para[3]-ic;//-0.5;
  shift2=para[4]-jc;//-0.5;  
  for(i=0;i<nxt;i++){
    x=i+shift1+nbound;
    for(j=0;j<nyt;j++){
    y=j+shift2+nbound;
    det=DSQR(x-xc)+DSQR(y-yc);
    yt[i][j]=kc*exp(det*sigma2v)+bgc;
    }
  }
  
}
void Interp_bicubic(int nx,int ny,float **a0,int nbound,
		    int nxt,int nyt, float **at,double *xcen)
{

  void splie2(float **ya, int m, int n, float **y2a);
  void spline(float y[], int n, float yp1, float ypn, float y2[]);
  void splint(float ya[], float y2a[], int n, float x, float *y);
  int i,j,m,n;
  float *ytmp,*yytmp,**y2a,x1,x2,shift1,shift2,ic,jc;
  y2a=matrix(0,nx,0,ny);
  ytmp=vector(0,ny);
  yytmp=vector(0,ny);
  splie2(a0,nx,ny,y2a);
  ic=nx*0.5;
  jc=ny*0.5;
  shift1=xcen[0]-ic;
  shift2=xcen[1]-jc;
  if(fabs(shift1)>nbound || fabs(shift2)>nbound){
    printf("bicubic shifts too much %e %e %d\n",shift1,shift2,nbound);
    getchar();
  }
  for (n=0;n<nyt;n++){
    x2=n+nbound+shift2;
    for(i=0;i<nx;i++){
      splint(a0[i],y2a[i],ny,x2,&yytmp[i]);
      spline(yytmp,ny,1.0e30,1.0e30,ytmp);
    }
    for (m=0;m<nxt;m++){
      x1=m+nbound+shift1;
      splint(yytmp,ytmp,ny,x1,&at[m][n]);
    }
  }
  free_vector(yytmp,0,ny);
  free_vector(ytmp,0,ny);
  free_matrix(y2a,0,nx,0,ny);

}

void BInterp_bicubic(int nx,int ny,float **a0,int nbound,
        int nxt,int nyt, float **at,double *xcen)
{

  void splie2(float **ya, int m, int n, float **y2a);
  void spline(float y[], int n, float yp1, float ypn, float y2[]);
  void splint(float ya[], float y2a[], int n, float x, float *y);
  int i,j,m,n;
  float *ytmp,*yytmp,**y2a,x1,x2,shift1,shift2,ic,jc;
  y2a=matrix(0,nx,0,ny);
  ytmp=vector(0,ny);
  yytmp=vector(0,ny);
  splie2(a0,nx,ny,y2a);
  ic=nx*0.5;
  jc=ny*0.5;
  shift1=ic-(xcen[0]+nbound);
  shift2=jc-(xcen[1]+nbound);//printf("shift1=%f,shift1=%f\n",shift1,shift2 );
  if(fabs(shift1)>nbound || fabs(shift2)>nbound){
    printf("bicubic shifts too much %e %e %d\n",shift1,shift2,nbound);
    getchar();
  }
  for (n=0;n<nyt;n++){
    x2=n+nbound+shift2;//printf("x2=%f\n",x2 );
    for(i=0;i<nx;i++){
      splint(a0[i],y2a[i],ny,x2,&yytmp[i]);
      spline(yytmp,ny,1.0e30,1.0e30,ytmp);
    }
    for (m=0;m<nxt;m++){
      x1=m+nbound+shift1;
      splint(yytmp,ytmp,ny,x1,&at[m][n]);
    }
  }
  free_vector(yytmp,0,ny);
  free_vector(ytmp,0,ny);
  free_matrix(y2a,0,nx,0,ny);

}
void splie2(float **ya, int m, int n, float **y2a)
{
  void spline(float y[], int n, float yp1, float ypn, float y2[]);
  int j;
  
  for (j=0;j<m;j++)spline(ya[j],n,1.0e30,1.0e30,y2a[j]);
}

void spline(float y[], int n, float yp1, float ypn, float y2[])
{
  int i,k;
  float p,qn,sig,un,*u;
  
  u=vector(0,n-1);
  if (yp1 > 0.99e30)
    y2[0]=u[0]=0.0;
  else {
    y2[0] = -0.5;
    u[0]=3.0*(y[1]-y[0]-yp1);
  }
  sig=0.5;
  for (i=1;i<=n-2;i++) {
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-2.*y[i]+y[i-1]);
    u[i]=(3.0*u[i]-sig*u[i-1])/p;
  }
  if (ypn > 0.99e30)
    qn=un=0.0;
  else {
    qn=0.5;
    un=3.0*(ypn-y[n-1]+y[n-2]);
  }
  y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
  for (k=n-2;k>=0;k--)y2[k]=y2[k]*y2[k+1]+u[k];
  free_vector(u,0,n-1);
}
void splint(float ya[], float y2a[], int n, float x, float *y)
{
  void nrerror(char error_text[]);
  int klo,khi,k;
  float h,b,a;
  
  klo=x;
if( klo<0)klo=0;
if(klo==n-1)klo=n-2;
  khi=klo+1;//printf("klo=%d,khi=%d\n",klo,khi );
  if (klo<0 || khi>=n) nrerror("Bad xa input to routine splint");
  a=(khi-x);
  b=(x-klo);
  *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])/6.0;
}
