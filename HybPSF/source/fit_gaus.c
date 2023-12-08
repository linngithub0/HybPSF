#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"


void gasfit(double *x, double *y,int nbin,double *para)
{
  void search(double *x,double *y,int nbin,double kc,double kd,double sigmac,
              double sigmad,double meanc,double meand,double *para);
  int i,j,k,imax=0;
  double ymax,ymin,kc,kd,sigmac,sigmad,nuggetc,nuggetd;
  ymin=ymax=y[imax];
  for(i=1;i<nbin;i++){
    if(ymax<y[i]){imax=i;ymax=y[i];}
    if(ymin>y[i])ymin=y[i];
    
  }
  kc=ymax*2./3.;kd=ymax/6.;
  sigmac=0.4;sigmad=0.2;
  nuggetc=0.5*ymin;
  //nuggetc=0;
  nuggetd=0.5*nuggetc;
  
  for(i=0;i<5;i++){
    search(x,y,nbin,kc,kd,sigmac,sigmad,nuggetc,nuggetd,para);
    kd*=0.4;sigmad*=0.4;nuggetd*=0.4;
    kc=para[0];sigmac=para[1];nuggetc=para[2];
  }
}
void search(double *x,double *y,int nbin,double kc,double kd,double sigmac,
            double sigmad,double nuggetc,double nuggetd,double *para)
{
  double k,sigma,nugget,k2,k20,sigma2v,det;
  int i,j,l,m;
  sigma2v=-1./(2.*sigmac*sigmac);
  nugget=nuggetc;
  k20=0;
  for(m=0;m<nbin;m++){
    det=x[m];
    det=kc*(1.-exp(det*det*sigma2v))+nugget-y[m];
    k20+=det*det;
  }             
  for(i=-2;i<=2;i++){
    k=kc+i*kd;
    if(k>0){
      for(l=-2;l<=2;l++){
        nugget=nuggetc+l*nuggetd;
        if(nugget>0){
          for(j=-2;j<=2;j++){
            sigma=sigmac+j*sigmad;
            if(sigma>0){
              sigma2v=-1./(2.*sigma*sigma);
              k2=0;
              for(m=0;m<nbin;m++){
                det=x[m];
                det=k*(1-exp(det*det*sigma2v))+nugget-y[m];
                k2+=det*det;
              }
              //printf("k20=%e k2=%e\n",k20,k2);
              if(k2<=k20){k20=k2;para[0]=k;para[1]=sigma;para[2]=nugget;}
            }
          }
        }
      }
    }
  }
}
void star_gaus(float **y,int nx,int ny,float **yt,float **residu,double *para)
{
  void gasfit_2D(double **x, double *y,int np,double *para,double bg0, int fbg);
  double **xs,*ys,ymax;
  int i,j,np,npt,im,jm,**xi,k,ix,iy;
  double kc,bgc,sigmac,xc,yc,sigma2v,det;
  double bg0=0;
  double a,b,c,d,e,f,g,h;
  a=y[0][0];b=y[0][ny-1];c=y[nx-1][0];d=y[nx-1][ny-1];
  int fbg=0;
  np=nx*ny;
  xs=dmatrix(0,np-1,0,1);
  ys=dvector(0,np-1);
  xi=imatrix(0,np-1,0,1);
  a=(a<b)?a:b;c=(c<d)?c:d;ymax=(a<c)?a:c;
  for(i=-5;i<=5;i++)for(j=-5;j<=5;j++){
    ix=i+nx/2;iy=j+ny/2;
      if(ymax<y[ix][iy]){ymax=y[ix][iy];im=ix;jm=iy;}
    }
  //printf("maxpos:%d,%d\n",im,jm);
  npt=0;
  for(i=-3;i<=3;i++){
    for(j=-3;j<=3;j++){
      xs[npt][0]=xi[npt][0]=i+im;
      xs[npt][1]=xi[npt][1]=j+jm;
      ys[npt]=y[im+i][jm+j];
      npt++;
    }
  }
  gasfit_2D(xs, ys,npt,para,bg0,fbg);
  kc=para[0];sigmac=para[1];bgc=para[2];xc=para[3];yc=para[4];
  //printf("%e %e %e %e %e\n",kc,sigmac,bgc,xc,yc);
  sigma2v=-1./(2.*sigmac*sigmac);
  for(i=0;i<nx;i++){
    for(j=0;j<ny;j++){
    det=DSQR(i-xc)+DSQR(j-yc);
    yt[i][j]=kc*exp(det*sigma2v)+bgc;
    residu[i][j]=y[i][j]-yt[i][j];
    }
  }
  free_dmatrix(xs,0,np-1,0,1);
  free_imatrix(xi,0,np-1,0,1);
  free_dvector(ys,0,np-1);
}

void galaxy_gaus(float **y,int nx,int ny,float **yt,float **residu,double *para)
{
  void gasfit_2D(double **x, double *y,int np,double *para,double bg0, int fbg);
  double **xs,*ys,ymax;
  int i,j,np,npt,im,jm,**xi,k,ix,iy;
  
  double kc,bgc,sigmac,xc,yc,sigma2v,det;
  double bg0=0;
  double a,b,c,d,e,f,g,h;
  a=y[0][0];b=y[0][ny-1];c=y[nx-1][0];d=y[nx-1][ny-1];
  int fbg=0;
  np=nx*ny;
  xs=dmatrix(0,np-1,0,1);
  ys=dvector(0,np-1);
  xi=imatrix(0,np-1,0,1);
  a=(a<b)?a:b;c=(c<d)?c:d;ymax=(a<c)?a:c;
  //printf("ymax=%f\n",ymax,y[0][0] );
  int dpix;
  if(nx>50 && ny>50){
    dpix=20;
  }
  else{
    dpix=5;
  }
  for(i=-dpix;i<=dpix;i++)for(j=-dpix;j<=dpix;j++){
    ix=i+nx/2;iy=j+ny/2;
      if(ymax<y[ix][iy]){
        ymax=y[ix][iy];//printf("ymax=%f\n",ymax,y[0][0] );
        im=ix;jm=iy;
      }
    }//printf("im=%d,jm=%d\n",im,jm );
  npt=0;
  for(i=-5;i<=5;i++){
    for(j=-5;j<=5;j++){
      xs[npt][0]=xi[npt][0]=i+im;
      xs[npt][1]=xi[npt][1]=j+jm;
      ys[npt]=y[im+i][jm+j];
      npt++;
    }
  }
gasfit_2D(xs, ys,npt,para,bg0,fbg);
  kc=para[0];sigmac=para[1];bgc=para[2];xc=para[3];yc=para[4];
  //printf("%e %e %e %e %e\n",kc,sigmac,bgc,xc,yc);
  sigma2v=-1./(2.*sigmac*sigmac);
  for(i=0;i<nx;i++){
    for(j=0;j<ny;j++){
    det=DSQR(i-xc)+DSQR(j-yc);
    yt[i][j]=kc*exp(det*sigma2v)+bgc;
    residu[i][j]=y[i][j]-yt[i][j];
    }
  }
  free_dmatrix(xs,0,np-1,0,1);
  free_imatrix(xi,0,np-1,0,1);
  free_dvector(ys,0,np-1);
}


void burst_gaus(float **y,int nx,int ny,float **yt,float **residu,double *para)
{
  void gasfit_2D(double **x, double *y,int np,double *para,double bg0, int fbg);
  double **xs,*ys,ymax;
  int i,j,np,npt,im,jm,**xi,k,ix,iy;
  double kc,bgc,sigmac,xc,yc,sigma2v,det;
  double bg0=0;
  double a,b,c,d,e,f,g,h;
  a=y[0][0];b=y[0][ny-1];c=y[nx-1][0];d=y[nx-1][ny-1];
  int fbg=0;
  np=nx*ny;
  xs=dmatrix(0,np-1,0,1);
  ys=dvector(0,np-1);
  xi=imatrix(0,np-1,0,1);
  a=(a<b)?a:b;c=(c<d)?c:d;ymax=(a<c)?a:c;
  int dpix;
  if(nx>50 && ny>50){
    dpix=20;
  }
  else{
    dpix=5;
  }
  for(i=-dpix;i<=dpix;i++)for(j=-dpix;j<=dpix;j++){
    ix=i+nx/2;iy=j+ny/2;
      if(ymax<y[ix][iy]){
        ymax=y[ix][iy];//printf("ymax=%f\n",ymax,y[0][0] );
        im=ix;jm=iy;
      }
    }//printf("im=%d,jm=%d\n",im,jm );
  npt=0;
  for(i=-2;i<=2;i++){
    for(j=-2;j<=2;j++){
      xs[npt][0]=xi[npt][0]=i+im;
      xs[npt][1]=xi[npt][1]=j+jm;
      ys[npt]=y[im+i][jm+j];
      npt++;
    }
  }
gasfit_2D(xs, ys,npt,para,bg0,fbg);
  kc=para[0];sigmac=para[1];bgc=para[2];xc=para[3];yc=para[4];
  //printf("%e %e %e %e %e\n",kc,sigmac,bgc,xc,yc);
  sigma2v=-1./(2.*sigmac*sigmac);
  for(i=0;i<nx;i++){
    for(j=0;j<ny;j++){
    det=DSQR(i-xc)+DSQR(j-yc);
    yt[i][j]=kc*exp(det*sigma2v)+bgc;
    residu[i][j]=y[i][j]-yt[i][j];
    }
  }
  free_dmatrix(xs,0,np-1,0,1);
  free_imatrix(xi,0,np-1,0,1);
  free_dvector(ys,0,np-1);
}


void gasfit_2D(double **x, double *y,int np,double *para,double bg0,int fbg)
{
  void search_2D(double **x,double *y,int np,double kc,double kd,double sigmac,
              double sigmad,double bgc,double bgd,double xc, double xd, 
     double yc, double yd,double *para,double bg0,int fbg);
  int i,j,k,imax=0,isigma;
  double ymax,ymin,kc,kd,sigmac,sigmad,bgc,bgd,ysigma,xc,yc,xd,yd;
  double det,dett;
  ymin=ymax=y[imax];
  for(i=1;i<np;i++){
    if(ymax<y[i]){imax=i;ymax=y[i];}
    if(ymin>y[i])ymin=y[i];
    
  }
  //printf("ymax=%e,ymin=%e\n",ymax,ymin);
  kc=ymax;kd=ymax/12.;
  det=ysigma=kc*exp(-0.5);
  for(i=0;i<np;i++){
    dett=fabs(ysigma-y[i]);
    if(dett<det){det=dett;isigma=i;}
  }
  xc=x[imax][0];yc=x[imax][1];
  sigmac=sqrt(DSQR(xc-x[isigma][0])+DSQR(yc-x[isigma][1]));
  xd=yd=sigmac*0.25;
  sigmad=0.25*sigmac;
  bgc=0.;bgd=fabs(ymin);
  
  for(i=0;i<4;i++){
      
    //printf("%e %e %e %e %e k2=%e\n",kc,sigmac,bgc,xc,yc,para[5]);
    search_2D(x,y,np,kc,kd,sigmac,sigmad,bgc,bgd,xc,xd,yc,yd,para,bg0,fbg);
    kd*=0.33;sigmad*=0.33;bgd*=0.33;xd*=0.33;yd*=0.33;
    kc=para[0];sigmac=para[1];bgc=para[2];xc=para[3];yc=para[4]; 
    
  }
  if(fbg==0)para[2]=bg0;
}
void search_2D(double **x,double *y,int np,double kc,double kd,double sigmac,
         double sigmad,double bgc,double bgd,double xc, double xd, 
         double yc, double yd,double *para,double bg0,int fbg)
{
  double k,sigma,bg,k2,k20,sigma2v,det,xt,yt;
  int i,j,l,m,p,q;
  sigma2v=-1./(2.*sigmac*sigmac);
  bg=bgc;
  k20=0;
  for(m=0;m<np;m++){
    det=DSQR(x[m][0]-xc)+DSQR(x[m][1]-yc);
    det=kc*exp(det*sigma2v)+bgc-y[m];
    k20+=det*det;
  }             
  for(i=-4;i<=4;i++){
    k=kc+i*kd;
    if(k>0){
      for(j=-4;j<=4;j++){
  sigma=sigmac+j*sigmad;
  if(sigma>0){
    sigma2v=-1./(2.*sigma*sigma);
    for(p=-4;p<=4;p++){
      xt=xc+p*xd;
      for(q=-4;q<=4;q++){
        yt=yc+q*yd;
        k2=0;
        if(fbg==0){
    bg=bg0;
    for(m=0;m<np;m++){
      det=DSQR(x[m][0]-xt)+DSQR(x[m][1]-yt);
      det=k*exp(det*sigma2v)+bg-y[m];
      k2+=det*det;
    }
        }
        else{
    for(l=-4;l<=4;l++){
      bg=bgc+l*bgd;
      for(m=0;m<np;m++){
        det=DSQR(x[m][0]-xt)+DSQR(x[m][1]-yt);
        det=k*exp(det*sigma2v)+bg-y[m];
        k2+=det*det;
      }
    }
        }
        //printf("k20=%e k2=%e\n",k20,k2);
        if(k2<=k20){k20=k2;para[5]=k2;
    para[0]=k;para[1]=sigma;para[2]=bg;para[3]=xt;para[4]=yt;}
      }
    }
  }   
      }
    }
  }
}
