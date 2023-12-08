#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"
#include "string.h"

double k_sigma_1D(double *data,int np0,double *dms,int *indx)
{
  void sigma_clean_1D(double *data,double *r2,double nsigma,double *dms,int *np,
		   int *indx);
  int i,j,m,n,k,ktotal,np[2];
  double sum=0,dm,sigma,*r2,*d1,det,sigma2;
  r2=dvector(0,np0);
  d1=dvector(0,np0);
  for(i=0;i<np0;i++)sum +=data[i];
  dm=sum/np0;
  sum=0;
  k=0;
  for(i=0;i<np0;i++){
      det=data[i]-dm;
      r2[i]=det*det;
      d1[i]=data[i];
      sum +=r2[i];
  }
  ktotal=np0;
  sigma=sqrt(sum/(np0-1));
  //--------------------- m0
  dms[0]=dm;dms[1]=sigma;np[0]=np[1]=ktotal;
  k=0;
  do {
    np[0]=np[1];
    //printf("mean=%e,sigma=%e\n",dms[0],dms[1]);
    sigma_clean_1D(d1,r2,2.,dms,np,indx);
    k++;
  }while((np[0]!=np[1]) && k<2);
  indx[np0]=np[1];
  free_dvector(r2,0,np0);
  free_dvector(d1,0,np0);
 
  return dms[1];

}
void sigma_clean_1D(double *data,double *r2,double nsigma,double *dms,int *np,int *indx)
{
  int i,m,n,np0,ktotal,k;
  double det,sum,dm,sigma,sigma2;
  dm=dms[0];
  sigma=dms[1];
  np0=np[0];
  sigma2=nsigma*nsigma*sigma*sigma;
  sum=0;
  k=0; 
  
  for(i=0;i<np0;i++){
    if(r2[i]<sigma2){
      r2[k]=r2[i];
      data[k]=data[i];
      sum +=data[k];
      indx[k]=indx[i];
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
