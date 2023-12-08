#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"
#include <string.h>
void moffat_fit(double **star,int m,int n,double *xcen,
		 double *para,double *gamma)
{
  void  star_ellip(double **I,int n1,int n2,double rg,double *xcen,
                 double *gamma,double *es);
  void moffat_moments(double **star,int m,int n,double *xcen,double *gamma,
		    double rd,double *moments);
  int i,j,k;
  double x,y,sum,scale,rt2,sigmaw=4,fwhm,rd,beta,es[2],moments[5];
  double weight,xc,yc,g1,g2,gscale,rd2v,r20,r42;
  g1=0;
  g2=0;
  gscale=1./(1.-g1*g1-g2*g2);
  sum=0;
  scale=0;
  weight=-0.5/(sigmaw*sigmaw);
  for(i=0;i<m;i++){
    x=(i-xcen[0]);
    for(j=0;j<n;j++){
      y=j-xcen[1];
      xc=(1.-g1)*x-g2*y;
      yc=(1.+g1)*y-g2*x;
      rt2=(xc*xc+yc*yc)*gscale;
      sum+=star[i][j]*rt2*exp(rt2*weight);
      scale+=star[i][j]*exp(rt2*weight);
    }
  }
  weight+=scale/sum;
  if(weight<=0){para[0]=-1;return;}
  else{sigmaw=1./sqrt(2.*weight);}
  star_ellip(star,m,n,sigmaw*1.25,xcen,gamma,es);
  g1=gamma[0];
  g2=gamma[1];
  gscale=1./(1.-g1*g1-g2*g2);
  //gauss estimate 1 done;
  sigmaw*=1.5;
  sum=0;
  scale=0;
  weight=-0.5/(sigmaw*sigmaw);
  for(i=0;i<m;i++){
    x=(i-xcen[0]);
    for(j=0;j<n;j++){
      y=j-xcen[1];
      xc=(1.-g1)*x-g2*y;
      yc=(1.+g1)*y-g2*x;
      rt2=(xc*xc+yc*yc)*gscale;
      sum+=star[i][j]*rt2*exp(rt2*weight);
      scale+=star[i][j]*exp(rt2*weight);
    }
  }
  weight+=scale/sum;
  if(weight<=0){para[0]=-1;return;}
  else{para[0]=sigmaw=1./sqrt(2.*weight);}
  //gauss estimate 2 done;
  fwhm=sigmaw*2.3548;
  //fwhm=5;
  beta=3.5;
  x=pow(2,1./beta);
  x=2*sqrt(x-1);
  rd=fwhm/x;
  for(i=0;i<3;i++){
    moffat_moments(star,m,n,xcen,gamma,rd,moments);
    r20=moments[1]/moments[0];
    r42=0.5*moments[2]/moments[1];
    rd=1./(1/r20-1/r42);
    if(rd<0){para[1]=-1;return;}
    rd=sqrt(fabs(rd));
  }
  beta=rd*rd/r42;
  if(beta<=1 || beta >100){para[2]=-1;return;}
  for(i=0;i<100;i++){
    moffat_moments(star,m,n,xcen,gamma,rd,moments);
    r20=moments[1]/moments[0];
    beta=1/((moments[2]*moments[0])/(2*moments[1]*moments[1])-1);
    rd=(beta+1)*r20; // beta'=beta+k;k=2;
    if(rd<0){para[1]=-1;return;}
    rd=sqrt(fabs(rd));
  }
  
  fwhm=2*rd*sqrt(pow(2,1/beta)-1);
  para[1]=rd;
  para[2]=beta;
  para[3]=fwhm;
  
}
void moffat_moments(double **star,int m,int n,double *xcen,double *gamma,
		    double rd,double *moments)
{ int i,j,k;
  double rd2v,r0,r2,r4,weight,x,y,xc,yc,rt2,g1,g2,g,gscale;
  
  g1=gamma[0];
  g2=gamma[1];
  gscale=1./(1.-g1*g1-g2*g2);
  rd2v=1/(rd*rd);
  r0=r2=r4=0;
  for(i=0;i<m;i++){
    x=(i-xcen[0]);
    for(j=0;j<n;j++){
      y=j-xcen[1];
      xc=(1.-g1)*x-g2*y;
      yc=(1.+g1)*y-g2*x;
      rt2=(xc*xc+yc*yc)*gscale;
      weight=(1+rt2*rd2v);
      weight=star[i][j]/(weight*weight*weight); //here k=2;
      r0+=weight;
      r2+=weight*rt2;
      r4+=weight*rt2*rt2;
    }
  }
  moments[0]=r0;
  moments[1]=r2;
  moments[2]=r4;
}  




void moffat_fit_help(double **star,int m,int n,double *xcen,
     double *para,double *gamma)
{
  void  star_ellip(double **I,int n1,int n2,double rg,double *xcen,
                 double *gamma,double *es);
  void moffat_moments(double **star,int m,int n,double *xcen,double *gamma,
        double rd,double *moments);
  void size_moment(double **I,int Ng,double *cen,double rg,double *se);
  int i,j,k;
  double x,y,sum,scale,rt2,sigmaw=4,fwhm,rd,beta,es[2],moments[5];
  double weight,xc,yc,g1,g2,gscale,rd2v,r20,r42;
  g1=0;
  g2=0;
  gscale=1./(1.-g1*g1-g2*g2);
  sum=0;
  scale=0;
  weight=-0.5/(sigmaw*sigmaw);
  for(i=0;i<m;i++){
    x=(i-xcen[0]);
    for(j=0;j<n;j++){
      y=j-xcen[1];
      xc=(1.-g1)*x-g2*y;
      yc=(1.+g1)*y-g2*x;
      rt2=(xc*xc+yc*yc)*gscale;
      sum+=star[i][j]*rt2*exp(rt2*weight);
      scale+=star[i][j]*exp(rt2*weight);
    }
  }
  weight=scale/sum;//printf("sigmaw=%f\n",weight );
  if(weight<=0){para[0]=-1;return;}
  else{sigmaw=sqrt(weight/2.);}
  //star_ellip(star,m,n,sigmaw*1.25,xcen,gamma,es);
  size_moment(star,m,xcen,sigmaw*1.25,gamma);
  g1=gamma[0];
  g2=gamma[1];//printf("g1=%f,g2=%f,fwhm=%f\n",g1,g2,sigmaw );
  g1=0;
  g2=0;
  gscale=1./(1.-g1*g1-g2*g2);
  //gauss estimate 1 done;
  sigmaw*=1.5;//printf("sigmaw=%f\n",sigmaw );
  sum=0;
  scale=0;
  weight=-0.5/(sigmaw*sigmaw);
  for(i=0;i<m;i++){
    x=(i-xcen[0]);
    for(j=0;j<n;j++){
      y=j-xcen[1];
      xc=(1.-g1)*x-g2*y;
      yc=(1.+g1)*y-g2*x;
      rt2=(xc*xc+yc*yc)*gscale;
      sum+=star[i][j]*rt2*exp(rt2*weight);
      scale+=star[i][j]*exp(rt2*weight);
    }
  }
  weight=scale/sum;
  if(weight<=0){para[0]=-1;return;}
  else{para[0]=sigmaw=sqrt(weight/2.);}
  //gauss estimate 2 done;
  fwhm=sigmaw*2.3548;//printf("sigmaw1=%f\n",sigmaw );
  //fwhm=5;
  beta=3.5;
  x=pow(2,1./beta);
  x=2*sqrt(x-1);
  rd=fwhm/x;
  for(i=0;i<3;i++){
    moffat_moments(star,m,n,xcen,gamma,rd,moments);
    r20=moments[1]/moments[0];
    r42=0.5*moments[2]/moments[1];
    rd=1./(1/r20-1/r42);
    if(rd<0){para[1]=-1;return;}
    rd=sqrt(fabs(rd));
  }
  beta=rd*rd/r42;
  if(beta<=1 || beta >100){para[2]=-1;return;}
  for(i=0;i<100;i++){
    moffat_moments(star,m,n,xcen,gamma,rd,moments);
    r20=moments[1]/moments[0];
    beta=1/((moments[2]*moments[0])/(2*moments[1]*moments[1])-1);
    if(beta<1){beta=1.;}if(beta>100){beta=100;}
    rd=fwhm/(2*sqrt(pow(2,1./beta)-1)); // beta'=beta+k;k=2;
    rd=sqrt(fabs(rd));
  }
  
  fwhm=2*rd*sqrt(pow(2,1/beta)-1);
  para[1]=rd;
  para[2]=beta;
  para[3]=fwhm;
  
}

void size_moment(double **I,int Ng,double *cen,double rg,double *se){

int i,j;
double scal,w,w11,w12,w22,r2,x,y;
w=0.;w11=0.;w12=0.;w22=0.;
scal=0.5/(rg*rg);
for(i=0;i<Ng;i++){
   for(j=0;j<Ng;j++){
      x=i-cen[0];y=j-cen[1];
      r2=x*x+y*y;
      w+=I[i][j]*exp(-r2*scal);
      w11+=x*x*I[i][j]*exp(-r2*scal);
      w12+=x*y*I[i][j]*exp(-r2*scal);
      w22+=y*y*I[i][j]*exp(-r2*scal);
     }
  }
w11/=w;w12/=w;w22/=w;
se[0]=(w11-w22)/(w11+w22);
se[1]=2.*w12/(w11+w22);

}