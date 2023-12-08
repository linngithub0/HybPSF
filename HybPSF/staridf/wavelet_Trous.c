#include "nrutil.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>



void trous_filtering(double ***wij,double **cJ,double **Rn, int nx,int ny,
		     int J,int *sign)
{
  double k_sigma(double **data,int nx,int ny);
  double sigma0,w[7]={0.889,0.2,0.086,0.041,0.02,0.01,0.005};
  int i,j,k;
  double **wt,*wtt,significant;
  *sign=0;
  sigma0=k_sigma(wij[0],nx,ny)/0.889;
  significant=3*sigma0*w[J];
  for(i=0;i<nx;i++){
    wtt=cJ[i];
    for(j=0;j<ny;j++){
      if(fabs(wtt[j])<significant)wtt[j]=0;
      else *sign=1;
      Rn[i][j]=wtt[j];
    }
  }
  for(k=0;k<J;k++){
    wt=wij[k];
    significant=3.*sigma0*w[k];
    //printf("k=%d significant=%e\n",k,significant);
    for(i=0;i<nx;i++){
      wtt=wt[i];
      for(j=0;j<ny;j++){
	if(fabs(wtt[j])<significant)wtt[j]=0;
	else *sign=1;
	Rn[i][j] += wtt[j];
      }
    }   
  }
    
  
}

void trous_decomp(double **I0,double ***w,double **cJ,int nx,int ny,int J)
{
  double f[5][5],**ct,cJij;
  int i,j,k,it,jt,m,n,deti=1;
  f[0][0]=f[0][4]=f[4][0]=f[4][4]=1./256;
  f[0][1]=f[1][0]=f[0][3]=f[3][0]=f[1][4]=f[4][1]=f[3][4]=f[4][3]=1./64; 
  f[0][2]=f[2][0]=f[2][4]=f[4][2]=3./128;
  f[1][1]=f[1][3]=f[3][1]=f[3][3]=1./16; 
  f[1][2]=f[2][1]=f[2][3]=f[3][2]=3./32;
  f[2][2]=9./64;
  ct=dmatrix(0,nx-1,0,ny-1);
  for(i=0;i<nx;i++){for(j=0;j<ny;j++)ct[i][j]=I0[i][j];}
  for(k=0;k<J;k++){
    for(i=0;i<nx;i++){
      for(j=0;j<ny;j++){
	cJij=0;
	for(m=0;m<5;m++){
	  it=i+(m-2)*deti;
	  if(it<0) it += nx;
	  if(it>=nx)it -= nx;
	  for(n=0;n<5;n++){
	    jt=j+(n-2)*deti;
	    if(jt<0) jt += ny;
	    if(jt>=ny)jt -= ny;
	    cJij +=ct[it][jt]*f[m][n]; 
	  }
	}
	cJ[i][j]=cJij;
	w[k][i][j]=ct[i][j]-cJij;
      }
    }
    deti *=2;
    for(i=0;i<nx;i++){for(j=0;j<ny;j++)ct[i][j]=cJ[i][j];}
  }
  free_dmatrix(ct,0,nx-1,0,ny-1);
}


double k_sigma(double **data,int nx,int ny)
{
  void sigma_clean(double *data,double *r2,double *dms,int *np);
  int i,j,m,n,k,ktotal,np[2];
  double sum=0,dm,sigma,*r2,*d1,det,sigma2,dms[2];
  r2=dvector(0,nx*ny);
  d1=dvector(0,nx*ny);
  for(i=0;i<nx;i++){for(j=0;j<ny;j++)sum +=data[i][j];}
  dm=sum/(nx*ny);
  sum=0;
  k=0;
  for(i=0;i<nx;i++){
    for(j=0;j<nx;j++){
      det=data[i][j]-dm;
      r2[k]=det*det;
      d1[k]=data[i][j];
      sum +=r2[k];
      k++;
    }
  }
  ktotal=k;
  sigma=sqrt(sum/(k-1));
  //--------------------- m0
  dms[0]=dm;dms[1]=sigma;np[0]=np[1]=ktotal;
  k=0;
  do {
    np[0]=np[1];
    sigma_clean(d1,r2,dms,np);
    k++;
  }while((np[0]!=np[1]) && k<4);
 
  free_dvector(r2,0,nx*ny);
  free_dvector(d1,0,nx*ny);
 
  return dms[1];

}
void sigma_clean(double *data,double *r2,double *dms,int *np)
{
  int i,m,n,np0,ktotal,k;
  double det,sum,dm,sigma,sigma2;
  dm=dms[0];
  sigma=dms[1];
  np0=np[0];
  sigma2=9.*sigma*sigma;
  sum=0;
  k=0; 
  
  //printf("np0=%d\n",np[0]);
  for(i=0;i<np0;i++){
    if(r2[i]<sigma2){
      r2[k]=r2[i];
      data[k]=data[i];
      sum +=data[k];
      k++;
    }
  }
  np[1]=k;
  if(k==np0)return;
  dm=sum/k;
  sum=0;
  for(i=0;i<k;i++){
      det=data[i]-dm;
      r2[i]=det*det;
      sum +=r2[i];
  }
  sigma=sqrt(sum/(k-1.));
  dms[0]=dm;
  dms[1]=sigma;

}
