#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"
#include <string.h>
#define Pi 3.1415926535897932384626433
double Plr_moffat(double rd,double beta,double r,int l){
  double vr_moffat(double rd,double beta,double r);
  double x,f,L[100],cp1;
  int i;
  if(l>50){printf("moffatlets order is larger than 50,pause");getchar();}
  if(l<0){printf("moffatlets order is less than 0,pause");getchar();}
  x=vr_moffat(rd,beta,r);
  L[0]=1;
  L[1]=1-x;
  if(l==0)f=L[0];
  else if(l==1)f=L[1];
  else{
    cp1=-1-x;
    for(i=2;i<=l;i++)L[i]=(2.+cp1/i)*L[i-1]-(1-1./i)*L[i-2];
    f=L[l];
  }
  f=f*sqrt((2*beta-1)/acos(-1.))/rd*pow(1+r*r/(rd*rd),-beta);
  return f;
}
double vr_moffat(double rd,double beta,double r){
  double rs,beta2,f,tmp;
  rs=r/rd;
  beta2=2*beta;
  tmp=pow(1+rs*rs,-beta2);
  f=(1-beta2)/beta2*log(tmp);
  //getchar();
  return f;
}

double Plr_gauss(double sigma,double r,int l){
  double x,f,L[100],cp1,r2,sigma2;
  int i;
  if(l>50){printf("gaussianlets order is larger than 50,pause");getchar();}
  if(l<0){printf("gaussianlets order is less than 0,pause");getchar();}
  r2=r*r;sigma2=sigma*sigma;
  x=r2/sigma2;
  L[0]=1;
  L[1]=1-x;
  if(l==0)f=L[0];
  else if(l==1)f=L[1];
  else{
    cp1=-1-x;
    for(i=2;i<=l;i++)L[i]=(2.+cp1/i)*L[i-1]-(1-1./i)*L[i-2];
    f=L[l];
  }
  f=1./(sqrt(Pi)*sigma)*exp((-0.5)*x)*f;
  return f;
}




double Plr0_moffat(double rd,double beta,double r,int l){
  double vr_moffat(double rd,double beta,double r);
  double x,f,x2,x3,x4;
  if(l>8){printf("moffatlets order is larger than 8,pause");getchar();}
  if(l<0){printf("moffatlets order is less than 0,pause");getchar();}
  x=vr_moffat(rd,beta,r);
  if(l==0)f=1;
  else if(l==1)f=-x+1;
  else if(l==2)f=(x*x-4*x+2)*0.5;
  else if(l==3)f=(-x*x*x+9*x*x-18*x+6)/6.;
  else if(l==4)f=(x*x*x*x-16*x*x*x+72*x*x-96*x+24)/24.;
  else if(l==5)f=(-x*x*x*x*x+25*x*x*x*x-200*x*x*x+600*x*x-600*x+120)/120.;
  else if(l==6){x2=x*x;x3=x*x*x;
    f=(x3*x3-36*x2*x3+450*x2*x2-2400*x3+5400*x2-4320*x+720)/720.;}
  else if(l==7){x2=x*x;x3=x*x*x;x4=x2*x2;
    f=(-x3*x4+49*x3*x3-882*x2*x3+7350*x4-29400*x3+52920*x2-35280*x+5040)/5040;
  }
  else{x2=x*x;x3=x*x*x;x4=x2*x2;
    f=(x4*x4-64*x4*x3+1568*x3*x3-18816*x2*x3+117600*x4-376320*x3+564480*x2
       -322560*x+40320)/40320;
  }
  f=f*sqrt((2*beta-1)/acos(-1.))/rd*pow(1+r*r/(rd*rd),-beta);
  return f;
}
