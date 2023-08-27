#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nrutil.h"
//l0    the max order
//msign if there is angular components
//mx    which moments we need   


void creat_moffatelets(double rd,double beta,int l0,int msign,int *mx,double **basisf0,int Mp,int ng,int *nbf)
{
  void basisfs_initial(int **lm,int l0,int msign,int *mx,int *nbf);
   double Plr_moffat(double rd,double beta,double r,int l);
   double Plr_gauss(double sigma,double r,int l);
   int **lm,nbft,k,i,j,l,m,kt,nbf0,ra,rb;
   double x,y,theta,tmpx,sum,r,pi2,dx,dy,xt,yt;
   pi2=2*acos(-1.);
   nbf0=(l0+1)*(l0+1);
   lm=imatrix(0,nbf0,0,1);
   basisfs_initial(lm,l0,msign,mx,&nbft);
   *nbf=nbft;
  
  for(k=0;k<nbft;k++){
    l=lm[k][0];m=lm[k][1];
    //printf("l=%d m=%d",l,m);
    //getchar();
    kt=0;
    sum=0;
    for(i=0;i<ng;i++){
      x=i-ng/2+0.5;
      for(j=0;j<ng;j++){
        y=j-ng/2+0.5;
        tmpx=0;
        for(ra=-5;ra<=5;ra++){
          xt=x+ra/11.;
          for(rb=-5;rb<=5;rb++){
            yt=y+rb/11.;
            r=sqrt(xt*xt+yt*yt);
            tmpx+=Plr_moffat(rd,beta,r,l);
          }
        }//tmpx/=121.;
        r=sqrt(x*x+y*y);
        //tmpx=Plr_moffat(rd,beta,r,l);
        if(m==0){
          basisf0[k][kt]=tmpx;
          sum+= basisf0[k][kt]*basisf0[k][kt];
          kt++;
        }
        else{
          if(r==0){basisf0[k][kt]=0;kt++;}
          else{
            theta=x/r;
            if(theta>1.)theta=1.;
            if(theta<-1.)theta=-1.;
            theta=acos(theta);
            if(y<0)theta=pi2-theta;
            if(m>0){
              basisf0[k][kt]=cos(m*theta)*tmpx;
              sum+= basisf0[k][kt]*basisf0[k][kt];
              kt++;
            }
            if(m<0){
              basisf0[k][kt]=sin(-m*theta)*tmpx;
              sum+= basisf0[k][kt]*basisf0[k][kt];
              kt++;
            }
          }
        }
      }
    }
    sum=sqrt(sum);
    for(i=0;i<kt;i++)basisf0[k][i]/=sum;
  }
  free_imatrix(lm,0,nbf0,0,1);
}

void creat_gaussianlets(double sigma,int l0,int msign,int *mx,double **basisf0,int Mp,int ng,int *nbf)
{
  void basisfs_initial(int **lm,int l0,int msign,int *mx,int *nbf);
   double Plr_gauss(double sigma,double r,int l);
   int **lm,nbft,k,i,j,l,m,kt,nbf0,ra,rb;
   double x,y,theta,tmpx,sum,r,pi2,dx,dy,xt,yt;
   pi2=2*acos(-1.);
   nbf0=(l0+1)*(l0+1);
   lm=imatrix(0,nbf0,0,1);
   basisfs_initial(lm,l0,msign,mx,&nbft);
   *nbf=nbft;
  
  for(k=0;k<nbft;k++){
    l=lm[k][0];m=lm[k][1];
    //printf("l=%d m=%d",l,m);
    //getchar();
    kt=0;
    sum=0;
    for(i=0;i<ng;i++){
      x=i-ng/2+0.5;
      for(j=0;j<ng;j++){
        y=j-ng/2+0.5;
        tmpx=0;
        for(ra=-5;ra<=5;ra++){
          xt=x+ra/11.;
          for(rb=-5;rb<=5;rb++){
            yt=y+rb/11.;
            r=sqrt(xt*xt+yt*yt);
            tmpx+=Plr_gauss(sigma,r,l);
          }
        }//tmpx/=169.;
        r=sqrt(x*x+y*y);
        if(m==0){
          basisf0[k][kt]=tmpx;
          sum+= basisf0[k][kt]*basisf0[k][kt];
          kt++;
        }
        else{
          if(r==0){basisf0[k][kt]=0;kt++;}
          else{
            theta=x/r;
            if(theta>1.)theta=1.;
            if(theta<-1.)theta=-1.;
            theta=acos(theta);
            if(y<0)theta=pi2-theta;
            if(m>0){
              basisf0[k][kt]=cos(m*theta)*tmpx;
              sum+= basisf0[k][kt]*basisf0[k][kt];
              kt++;
            }
            if(m<0){
              basisf0[k][kt]=sin(-m*theta)*tmpx;
              sum+= basisf0[k][kt]*basisf0[k][kt];
              kt++;
            }
          }
        }
      }
    }
    sum=sqrt(sum);
    for(i=0;i<kt;i++)basisf0[k][i]/=sum;
  }
  free_imatrix(lm,0,nbf0,0,1);
}




void basisfs_initial(int **lm,int l0,int msign,int *mx,int *nbf){
  int **ibig,i,j,l,m,k,mt[100],kt;
  double md,jd;
  if(l0>40){
    printf("the max order of basis function is too high, I won't di it in basisfs_initial\n");
    getchar();
  }
  for(i=0;i<=12;i++)mt[i]=0;
  if(mx[0])mt[0]=1;
  if(mx[1])mt[1]=2;
  if(mx[2])mt[2]=3;
  if(mx[3])mt[3]=5;
  if(mx[4])mt[4]=7;
  if(mx[5])mt[5]=11;
  if(mx[6])mt[6]=13;
  if(mx[7])mt[7]=17;
  if(mx[8])mt[8]=19;
  if(mx[9])mt[9]=23;
  if(mx[10])mt[10]=29;
  if(mx[11])mt[11]=31;
  if(mx[12])mt[12]=37;
  ibig=imatrix(0,l0,-l0,l0);
  for(i=0;i<=l0;i++)for(j=-l0;j<=l0;j++)ibig[i][j]=0;
  for(i=0;i<=l0;i++)ibig[i][0]=1;
  if(msign==1){
    if(mt[0])for(i=0;i<l0;i++)ibig[i][-1]=ibig[i][1]=1;
    for(k=1;k<=12;k++){
      if(mt[k]<l0 & mt[k]>0){
  for(i=0;i<=l0;i++)for(j=l0;j>=-l0;j--){
      kt=i+abs(j);
      if(kt<=l0){
        jd=j;
        md=jd/mt[k];
        m=md;
        if(m==md)ibig[i][j]=1;
      }
    }
      }
    }
  }
  /* for(j=l0;j>=-l0;j--){
    for(i=0;i<=l0;i++)printf("%d ",ibig[i][j]);
    printf("\n");
    }*/
  k=0;
  for(i=0;i<=l0;i++)for(j=i;j>=-i;j--){
      l=i-abs(j);
      if(ibig[l][j]){lm[k][0]=l;lm[k][1]=j;k++;}
    }
  *nbf=k;
  free_imatrix(ibig,0,l0,-l0,l0);

}

void creat_moffatelet(double rd,double beta,int l0,int msign,int *mx,double **basisf0,int Mp,int ng,int *nbf)
{
  void basisfs_initial(int **lm,int l0,int msign,int *mx,int *nbf);
   double Plr_moffat(double rd,double beta,double r,int l);
   int **lm,nbft,k,i,j,l,m,kt,nbf0,ra,rb;
   double x,y,theta,tmpx,sum,r,pi2,dx,dy,xt,yt;
   pi2=2*acos(-1.);
   nbf0=(l0+1)*(l0+1);
   lm=imatrix(0,nbf0,0,1);
   basisfs_initial(lm,l0,msign,mx,&nbft);
   *nbf=nbft;
  
  for(k=0;k<nbft;k++){
    l=lm[k][0];m=lm[k][1];
    //printf("l=%d m=%d",l,m);
    //getchar();
    kt=0;
    sum=0;
    for(i=0;i<ng;i++){
      x=i-ng/2;
      for(j=0;j<ng;j++){
        y=j-ng/2;tmpx=0;
        for(ra=-6;ra<=6;ra++){
          xt=x+ra/13.;
          for(rb=-6;rb<=6;rb++){
            yt=y+rb/13.;
            r=sqrt(xt*xt+yt*yt);
            tmpx+=Plr_moffat(rd,beta,r,l);
          }
        }tmpx/=169.;
        r=sqrt(x*x+y*y);
        //tmpx=Plr_moffat(rd,beta,r,l);
        if(m==0){
          basisf0[k][kt]=tmpx;
          sum+= basisf0[k][kt]*basisf0[k][kt];
          kt++;
        }
        else{
          if(r==0){basisf0[k][kt]=0;kt++;}
          else{
            theta=x/r;
            if(theta>1.)theta=1.;
            if(theta<-1.)theta=-1.;
            theta=acos(theta);
            if(y<0)theta=pi2-theta;
            if(m>0){
              basisf0[k][kt]=cos(m*theta)*tmpx;
              sum+= basisf0[k][kt]*basisf0[k][kt];
              kt++;
            }
            if(m<0){
              basisf0[k][kt]=sin(-m*theta)*tmpx;
              sum+= basisf0[k][kt]*basisf0[k][kt];
              kt++;
            }
          }
        }
      }
    }
    sum=sqrt(sum);
    for(i=0;i<kt;i++)basisf0[k][i]/=sum;
  }
  free_imatrix(lm,0,nbf0,0,1);
}

void creat_gaussianlet(double sigma,int l0,int msign,int *mx,double **basisf0,int Mp,int ng,int *nbf)
{
  void basisfs_initial(int **lm,int l0,int msign,int *mx,int *nbf);
   double Plr_gauss(double sigma,double r,int l);
   int **lm,nbft,k,i,j,l,m,kt,nbf0,ra,rb;
   double x,y,theta,tmpx,sum,r,pi2,dx,dy,xt,yt;
   pi2=2*acos(-1.);
   nbf0=(l0+1)*(l0+1);
   lm=imatrix(0,nbf0,0,1);
   basisfs_initial(lm,l0,msign,mx,&nbft);
   *nbf=nbft;
  
  for(k=0;k<nbft;k++){
    l=lm[k][0];m=lm[k][1];
    //printf("l=%d m=%d",l,m);
    //getchar();
    kt=0;
    sum=0;
    for(i=0;i<ng;i++){
      x=i-ng/2;
      for(j=0;j<ng;j++){
        y=j-ng/2;
        tmpx=0;
        for(ra=-5;ra<=5;ra++){
          xt=x+ra/11.;
          for(rb=-5;rb<=5;rb++){
            yt=y+rb/11.;
            r=sqrt(xt*xt+yt*yt);
            tmpx+=Plr_gauss(sigma,r,l);
          }
        }//tmpx/=169.;
        r=sqrt(x*x+y*y);
        if(m==0){
          basisf0[k][kt]=tmpx;
          sum+= basisf0[k][kt]*basisf0[k][kt];
          kt++;
        }
        else{
          if(r==0){basisf0[k][kt]=0;kt++;}
          else{
            theta=x/r;
            if(theta>1.)theta=1.;
            if(theta<-1.)theta=-1.;
            theta=acos(theta);
            if(y<0)theta=pi2-theta;
            if(m>0){
              basisf0[k][kt]=cos(m*theta)*tmpx;
              sum+= basisf0[k][kt]*basisf0[k][kt];
              kt++;
            }
            if(m<0){
              basisf0[k][kt]=sin(-m*theta)*tmpx;
              sum+= basisf0[k][kt]*basisf0[k][kt];
              kt++;
            }
          }
        }
      }
    }
    sum=sqrt(sum);
    for(i=0;i<kt;i++)basisf0[k][i]/=sum;
  }
  free_imatrix(lm,0,nbf0,0,1);
}
