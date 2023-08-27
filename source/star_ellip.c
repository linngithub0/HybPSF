#include<math.h>
#include<stdio.h>
#include"nrutil.h"
void  star_ellip(double **I,int n1,int n2,double rg,double *xcen,
		 double *gamma,double *gerr,double noise, int ierror )
{
 void  KSBc(double **I,int n1,int n2,double rg,double *xcen,double *gamma);
 void  KSBe(double **I, int n1, int n2,double rg,double *xcen,double *gamma);
 void estimate_gamma_err(double **I, int n1, int n2,double rg,double *xcen,
			 double noise,double *g,double *err);
 double g[2],detg,errf=0.0001;
 int i=0;
 KSBc(I,n1,n2,rg,xcen,gamma);
 do{
   g[0]=gamma[0];g[1]=gamma[1];
   KSBe(I,n1,n2,rg,xcen,gamma);
   detg=sqrt(DSQR(g[0]-gamma[0])+DSQR(g[1]-gamma[1]));
   i++;
 }while(i<10 && detg>errf);
 if(ierror!=0){
   estimate_gamma_err(I,n1,n2,rg,xcen,noise,gamma,gerr);
 }
}


void  KSBc(double **I,int n1,int n2,double rg,double *xcen,
		 double *gamma)
{
  void arr22_invers(double **a,double **at);
  void arr22_time_v(double **a,double *v,double *vt);
  double R2MAX(int n1, int n2,double *xcen,double *g);  
  double I11,I12,I22,esh[2],e[2],xsh11,xsh12,xsh21,xsh22;
  double w,w1,w2,scal,detx,dety,detx2,dety2,r2,r;
  double detxy,det2xy,T,esh1,esh2,R2max,gt[2],**Psh,**Pshv;
  int i,j,k;
  Psh=dmatrix(0,1,0,1);
  Pshv=dmatrix(0,1,0,1);
  
  gt[0]=0; gt[1]=0;
  R2max=R2MAX(n1,n2,xcen,gt);
  scal=0.5/(rg*rg);
  I11=0;I22=0;I12=0;xsh11=0;xsh12=0;xsh22=0;xsh21=0;esh1=0;esh2=0;
  for(j=0;j<n2;j++){
    dety=j-xcen[1];
    dety2=dety*dety;
    for(i=0;i<n1;i++){
      detx=i-xcen[0];
      detx2=detx*detx;
      detxy=2.*detx*dety;
      det2xy=detx2-dety2;
      r2=detx2+dety2;
      if(r2<R2max){
	w=I[i][j]*exp(-r2*scal);
	w1=-w*scal; //respect to r2;
	w2=-w1*scal;
	I11 += w*detx2;
	I22 += w*dety2; 
	I12 += w*detxy;
	xsh11 += w*r2+w1*det2xy*det2xy;
	xsh12 += w1*det2xy*detxy;
	xsh22 += w*r2+w1*detxy*detxy;
	esh1  += w1*r2*det2xy;
	esh2  += w1*r2*detxy;
      }
    }
  }
  T=1./(I11+I22);
  e[0]=(I11-I22)*T;
  e[1]=I12*T; // factor of 2 has been corrected
  xsh11=2.*xsh11*T;
  xsh12=2.*xsh12*T;
  xsh21=xsh12;
  xsh22=2.*xsh22*T;
  esh[0]=2.*(e[0]+esh1*T);
  esh[1]=2.*(e[1]+esh2*T);
    
  Psh[0][0]=xsh11-e[0]*esh[0];
  Psh[0][1]=xsh12-e[0]*esh[1];
  Psh[1][0]=xsh21-e[1]*esh[0];
  Psh[1][1]=xsh22-e[1]*esh[1];
  detx=Psh[0][0]+Psh[1][1];
  Psh[0][0]=Psh[1][1]=detx/2;
  Psh[0][1]=Psh[1][0]=0;
  arr22_invers(Psh,Pshv);
  arr22_time_v(Pshv,e,gamma);
  free_dmatrix(Psh,0,1,0,1);
  free_dmatrix(Pshv,0,1,0,1);

}

void  KSBe(double **I, int n1, int n2,double rg,double *xcen,double *gamma)
{ 
  double R2MAX(int n1, int n2,double *xcen,double *g);
  void arr22_invers(double **a,double **at);
  void arr22_time_v(double **a,double *v,double *vt);
  double I11,I12,I22,Iw11,Iw12,Iw22;
  double esh[2],esm[2],xsh11,xsh12,xsh21,xsh22;
  double w,w1,w2,scal,detx,dety,detx2,dety2,r2,r,R2,A,B,absg;
  double detxy,det2xy,T,esh1,esh2,ew1,ew2,ww,g1,g2,**Psh,**Pshv;
  double detm,detn,detmn,det2mn,R2max,gt[2],p1,p2,p3,p4,p5,p6,tmp,e[2];
  FILE   *fp;
  int i,j,k;
  
  Psh=dmatrix(0,1,0,1);
  Pshv=dmatrix(0,1,0,1);
  gt[0]=gamma[0]*0.8;
  gt[1]=gamma[1]*0.8;
  R2max=R2MAX(n1,n2,xcen,gt);
  g1=gamma[0];
  g2=gamma[1];
  A=(1.-g1)*(1.-g1)+g2*g2;
  B=(1.+g1)*(1.+g1)+g2*g2;
  scal=0.5/(rg*rg);
  I11=0;I22=0;I12=0;
  Iw11=0;Iw22=0;Iw12=0;
  xsh11=0;xsh12=0;xsh22=0;xsh21=0;
  p1=p2=p3=p4=p5=p6=0;
  esh1=0;esh2=0;
  for(j=0;j<n2;j++){
    dety=j-xcen[1];
    dety2=dety*dety;
    for(i=0;i<n1;i++){
      detx=i-xcen[0];
      detx2=detx*detx;
      detxy=2.*detx*dety;
      det2xy=detx2-dety2;
      r2=detx2+dety2;
      R2=A*detx2+B*dety2-2.*g2*detxy;
      if(R2<R2max){
	detm=A*detx-2.*g2*dety;
	detn=B*dety-2.*g2*detx;
	detmn=2.*detm*detn;
	det2mn=detm*detm-detn*detn;
	w=exp(-R2*scal);
	//w=expp(R2*scal);//the sign has been corrected in expp
	ww=w*w;//sqrt(w*I[i][j]);
	w *=I[i][j];
	w1=-w*scal; //respect to R2;
	w2=w*scal*scal;
	I11 += w*detx2;
	I22 += w*dety2; 
	I12 += w*detxy;
	Iw11 += ww*detx2;
	Iw22 += ww*dety2; 
	Iw12 += ww*detxy;
	xsh11 += w*r2+w1*(A*detx2*detx2+B*dety2*dety2-(A+B)*detx2*dety2);
	xsh12 += w1*det2xy*((A+B)*detx*dety-2.*g2*r2);
	xsh21 += w1*detxy*(A*detx2-B*dety2);
	xsh22 +=w*r2+2.*w1*(A+B)*detx2*dety2-2.*w1*detxy*g2*r2;
	esh1  +=w1*(A*detx2*detx2-B*dety2*dety2+(A-B)*detx2*dety2);
	esh2  +=w1*r2*((A+B)*detx*dety-2.*g2*r2);
	tmp=w1*((A-B)*detxy+4.*g2*det2xy);
	p2 += 2.*w*detxy+det2xy*tmp;
	p4 += detxy*tmp-2.*w*det2xy;
	p6  += tmp*r2;
	tmp=4*w + 2*w1*(A*detx2+B*dety2-2*g2*detxy);
	p1 += det2xy*tmp;
	p3 += detxy*tmp;
	p5 += r2*tmp;
      }
    }
  }
  T=1./(I11+I22);
  e[0]=(I11-I22)*T;
  e[1]=I12*T; // factor of 2 has been corrected

  xsh11=2.*xsh11*T+(-g1*p1-g2*p2)*T;
  xsh12=2.*xsh12*T+( g1*p2-g2*p1)*T;
  xsh21=2.*xsh21*T+(-g1*p3-g2*p4)*T;
  xsh22=2.*xsh22*T+( g1*p4-g2*p3)*T;
  esh[0]=2.*(e[0]+esh1*T)+(-g1*p5-g2*p6)*T;
  esh[1]=2.*(e[1]+esh2*T)+( g1*p6-g2*p5)*T;
    
  Psh[0][0]=xsh11-e[0]*esh[0];
  Psh[0][1]=xsh12-e[0]*esh[1];
  Psh[1][0]=xsh21-e[1]*esh[0];
  Psh[1][1]=xsh22-e[1]*esh[1];

  absg=g2*g2+g1*g1;
  absg=1./(1-absg);
  Psh[0][0] *=absg;Psh[0][1] *=absg;Psh[1][0] *=absg;Psh[1][1] *=absg;
  T=1./(Iw11+Iw22);
  ew1=(Iw11-Iw22)*T;
  ew2=Iw12*T; // factor of 2 has been corrected
  e[0] -=ew1;
  e[1] -=ew2;
  arr22_time_v(Psh,gamma,esm);
  e[0] +=esm[0];
  e[1] +=esm[1];
  arr22_invers(Psh,Pshv);
  arr22_time_v(Pshv,e,gamma);
  free_dmatrix(Psh,0,1,0,1);
  free_dmatrix(Pshv,0,1,0,1);
}

void estimate_gamma_err(double **I, int n1, int n2,double rg,double *xcen,
			double noise,double *g,double *err)
{ double expp(double x);
  double R2MAX(int n1, int n2,double *xcen,double *g);

  double I11,I12,I22,I1111,I2222,I1212,I1112,I2212;
  double w,w0,scal,detx,dety,detx2,dety2,r2,r,R2,A,B,absg;
  double detxy,det2xy,T,g1,g2,R2max,gt[2],tmp;
  double e1,e2,a11,a12,a22,a21,gh2,sigmae1,sigmae2,sigmae12;
  FILE   *fp;
  int i,j,k;
  gt[0]=g[0];
  gt[1]=g[1];
  R2max=R2MAX(n1,n2,xcen,gt);
  g1=g[0];
  g2=g[1];
  gh2=g1*g1+g2*g2;
  e1=2.*g1/(1.+gh2);
  e2=2.*g2/(1.+gh2);
  //////////////////////////
  scal=(1+gh2)/(1.-gh2)/2;
  a11=scal*(1+g1*g1-g2*g2);
  a22=scal*(1-g1*g1+g2*g2);
  a12=a21=scal*2.*g1*g2;
  
  //printf("%e %e %e\n",a11,a12,a22);
  A=(1.-g1)*(1.-g1)+g2*g2;
  B=(1.+g1)*(1.+g1)+g2*g2;
  I11=I12=I22=I1111=I2222=I1212=I1112=I2212=0;
  //  printf("g1=%5.2f, g2=%5.2f, A=%5.2f, B=%5.2f \n", g1, g2, A, B);
  scal=0.5/(rg*rg);
  for(j=0;j<n2;j++){
    dety=j-xcen[1];
    dety2=dety*dety;
    for(i=0;i<n1;i++){
      detx=i-xcen[0];
      detx2=detx*detx;
      detxy=2.*detx*dety;
      det2xy=detx2-dety2;
      r2=detx2+dety2;
      R2=A*detx2+B*dety2-2.*g2*detxy;
      if(R2<R2max){
	w=exp(-R2*scal);
	//w=expp(r2*scal); //the sign has been corrected in expp
	w0=w*I[i][j];
	I11 += w0*detx2;
	I22 += w0*dety2;
	I12 += w0*detxy;
	//w =w*w*(bsigma2+I[i][j]);
	w =w*w*(noise*noise);
	detxy *= 0.5;
	I1111+=detx2*detx2*w;
	I1112+=detx2*detxy*w;
	I1212+=detxy*detxy*w;
	I2212+=dety2*detxy*w;
	I2222+=dety2*dety2*w;
      }
    }
  }
  T=(I11+I22);
  e1=(I11-I22)/T;
  e2=I12/T;
  T=1./(T*T);
  sigmae1=T*((1.-e1)*(1.-e1)*I1111+(1+e1)*(1+e1)*I2222-2.*(1-e1*e1)*I1212);
  sigmae2=T*(e2*e2*(I1111+I2222)+(4.+2.*e2*e2)*I1212-4.*e2*(I1112+I2212));
  sigmae12=T*(e2*(e1-1.)*I1111+2*(1-e1)*I1112+
    2*e1*e2*I1212-2*(1+e1)*I2212+e2*(1+e1)*I2222);
  
  err[0]=sqrt(a11*a11*sigmae1+2*a11*a12*sigmae12+a12*a12*sigmae2);
  err[1]=sqrt(a21*a21*sigmae1+2*a22*a12*sigmae12+a22*a22*sigmae2);
  
  
}

double R2MAX(int n1, int n2,double *xcen,double *g)
{
  double xmin,xmax,ymax,ymin,x,y,r,r1,r2,r3,r4,A,B,g1,g2;
  g1=g[0];
  g2=g[1]; 
  xmin=-xcen[0];
  ymin=-xcen[1];
  xmax=n1-1-xcen[0];
  ymax=n2-1-xcen[1];
  A=(1.-g1)*(1.-g1)+g2*g2;
  B=(1.+g1)*(1.+g1)+g2*g2;
  x=2.*g2*ymin/A;
  r=A*x*x+B*ymin*ymin-4.*g2*x*ymin;
  r2=r;
  x=2.*g2*ymax/A;
  r=A*x*x+B*ymax*ymax-4.*g2*x*ymax;
  if(r<r2)r2=r;
  y=2.*g2*xmin/B;
  r=A*xmin*xmin+B*y*y-4.*g2*xmin*y;
  if(r<r2)r2=r;
  y=2.*g2*xmax/B;
  r=A*xmax*xmax+B*y*y-4.*g2*xmax*y;
  if(r<r2)r2=r;
  return r2;
  
}

void arr22_invers(double **a,double **at)
{
  double abs;
  abs=a[0][0]*a[1][1]-a[1][0]*a[0][1];
  abs=1./abs;
  at[0][0]=a[1][1]*abs;
  at[1][1]=a[0][0]*abs;
  at[1][0]=-a[1][0]*abs;
  at[0][1]=-a[0][1]*abs;
}



void arr22_time_v(double **a,double *v,double *vt)
{
  vt[0]=a[0][0]*v[0]+a[0][1]*v[1];
  vt[1]=a[1][0]*v[0]+a[1][1]*v[1];
}


