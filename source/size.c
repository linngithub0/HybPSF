#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"

void size(float **I,int Ng,double *cen,double rg,double *se){

int i,j;
double scal,w,w11,w12,w22,r2,x,y,weight,nh=Ng/2.;
nh=nh*nh;
w=0.;w11=0.;w12=0.;w22=0.;
scal=0.5/(rg*rg);
for(i=0;i<Ng;i++){
   for(j=0;j<Ng;j++){
      x=i-cen[0];y=j-cen[1];
      r2=x*x+y*y;weight=exp(-r2*scal);
      if(r2<=nh){
        w+=I[i][j]*weight;
        w11+=x*x*I[i][j]*weight;
        w12+=x*y*I[i][j]*weight;
        w22+=y*y*I[i][j]*weight;
      }
     }
  }
w11/=w;w12/=w;w22/=w;
se[0]=w11+w22;
se[1]=(w11-w22)/(w11+w22);
se[2]=2.*w12/(w11+w22);

}

void size_rot(double **I,int Ng,double *cen,double rg,double *se,double *se2){

int i,j;
double scal,w,w11,w12,w22,r2,x,y,xr,yr,theta,pi=acos(-1.);
r2=sqrt(se[1]*se[1]+se[2]*se[2]);//printf("r2=%f\n",r2 );
if(r2!=0){
   theta=se[2]/r2;
   if(theta>1.)theta=1.;
   if(theta<-1.)theta=-1.;
   theta=acos(theta);
   if(se[2]<0)theta+=pi;
}else{
   theta=0;
}//printf("theta=%f,e1=%f,e2=%f\n",theta,se[1],se[2] );
theta/=2.;
w=0.;w11=0.;w12=0.;w22=0.;
scal=0.5/(rg*rg);
for(i=0;i<Ng;i++){
   for(j=0;j<Ng;j++){
      x=i-cen[0];y=j-cen[1];
      xr=cos(theta)*x+sin(theta)*y;yr=cos(theta)*y-sin(theta)*x;
      r2=xr*xr+yr*yr;
      w+=I[i][j]*exp(-r2*scal);
      w11+=xr*xr*I[i][j]*exp(-r2*scal);
      w12+=xr*yr*I[i][j]*exp(-r2*scal);
      w22+=yr*yr*I[i][j]*exp(-r2*scal);
     }
  }
w11/=w;w12/=w;w22/=w;
se2[0]=w11+w22;
se2[1]=(w11-w22)/(w11+w22);
se2[2]=2.*w12/(w11+w22);
//printf("e1=%f,e2=%f\n",se[1],se[2] );

}