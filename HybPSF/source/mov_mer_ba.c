#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"
#include <string.h>
void mov_mer_pc(double **erbasisf,double ***basis,int Nke,int Ng,
                int num,double **cent,int Nstar,int nbf){

  int Nge;
  int i,j,l,f,ra,rb,ic,pl,r,k,t;
  double x,y,dx,dy;
  int write_fits_3D(char *argv,float ***stars,int *dim);
  void fft_int(float **image,int n1_image,int n2_image,float **result,
    double dx,double dy);
  int dim[3];char ftmp[200];
  double ***basisl,***recive0,***basis0,pp,w1,w2,w3,w4,sum;
  float ***recive1,***recive;
  Nge=Ng*num;
  r=(Nke-Nge)/2;
  basisl=d3tensor(0,nbf,0,Ng-1,0,Ng-1);
  recive1=f3tensor(0,nbf,0,Nke-1,0,Nke-1);
  recive=f3tensor(0,nbf,0,Nke-1,0,Nke-1);

  //printf("starting slips!\n");

  for(l=0;l<nbf;l++){
      f=0;
      for(i=0;i<Nke;i++){
        for(j=0;j<Nke;j++){
          recive1[l][i][j]=erbasisf[l][f];
          f++;
        }
      } 
    }//printf("slips done&start merging!\n");
  for(ic=0;ic<Nstar;ic++){
    ra=pp=(Nke/2)-cent[ic][0];//t=pp*10;x=t%10;
    x=(pp-ra)*10.;
    rb=pp=(Nke/2)-cent[ic][1];//t=pp*10;y=t%10;
    y=(pp-rb)*10.;
    w1=(10-x)*(10-y);w1/=100.;w2=x*(10-y);w2/=100.;w3=(10-x)*y;w3/=100.;w4=x*y;w4/=100.;
      for(l=0;l<nbf;l++){
        k=0;f=0;
        for(i=0;i<Nge;i++){
          for(j=0;j<Nge;j++){
            recive[l][i][j]=recive1[l][ra+i][rb+j]*w1+recive1[l][ra+i+1][rb+j]*w2+
            recive1[l][ra+i][rb+j+1]*w3+recive1[l][ra+i+1][rb+j+1]*w4;
            k++;
          }
        }
      }
     for(l=0;l<nbf;l++){
      sum=0;
        for(i=0;i<Ng;i++){
          for(j=0;j<Ng;j++){
            basisl[l][i][j]=0;
            for(pl=0;pl<num;pl++){
              ra=i*num+pl;
              for(f=0;f<num;f++){
                rb=j*num+f;
                basisl[l][i][j]+=recive[l][ra][rb];//merging divided pixels
               } 
              } 
              sum+=basisl[l][i][j]*basisl[l][i][j];
            }
          }
        for(i=0;i<Ng;i++){
          for(j=0;j<Ng;j++){
            //basisl[l][i][j]/=sqrt(sum);
          }
        }
     }//printf("merging is done&staring wite!Nstar=%d\n",ic);
    for(l=0;l<nbf;l++){
      f=0;
      for(i=0;i<Ng;i++){
        for(j=0;j<Ng;j++){
          basis[ic][l][f]=basisl[l][i][j];f++;
        }
      }
    }
  }


  free_d3tensor(basisl,0,nbf,0,Ng-1,0,Ng-1);
  free_f3tensor(recive1,0,nbf,0,Nke-1,0,Nke-1);
  free_f3tensor(recive,0,nbf,0,Nke-1,0,Nke-1);

}

void mop(double *erbasisf,double **basis,int Nke,int Ng,
                int num,double **cent,int Nstar){

  int Nge;
  int i,j,l,f,ra,rb,ic,pl,r,k,t;
  double x,y,dx,dy;
  int write_fits_3D(char *argv,float ***stars,int *dim);
  void fft_int(float **image,int n1_image,int n2_image,float **result,
    double dx,double dy);
  int dim[3];char ftmp[200];
  double **basisl,**recive0,**basis0,pp,w1,w2,w3,w4,sum;
  float **recive1,**recive;
  Nge=Ng*num;
  r=(Nke-Nge)/2;
  basisl=dmatrix(0,Ng-1,0,Ng-1);
  recive1=matrix(0,Nke-1,0,Nke-1);
  recive=matrix(0,Nke-1,0,Nke-1);

  float ***write;
  write=f3tensor(0,Nstar,0,Nke,0,Nke);
  char *tmp;


  //printf("starting slips!\n");

      f=0;
      for(i=0;i<Nke;i++){
        for(j=0;j<Nke;j++){
          recive1[i][j]=erbasisf[f];
          write[0][i][j]=recive1[i][j];
          f++;
        }
      }
      tmp="fits/recive.fits";
      dim[0]=dim[1]=Nke;dim[2]=1;
      //write_fits_3D(tmp,write,dim); 
    //printf("slips done&start merging!\n");
  for(ic=0;ic<Nstar;ic++){
    ra=pp=(Nke/2)-cent[ic][0];//t=pp*10;x=t%10;
    x=(pp-ra)*10.;
    rb=pp=(Nke/2)-cent[ic][1];//t=pp*10;y=t%10;
    y=(pp-rb)*10.;
    w1=(10.-x)*(10.-y);w1/=100.;w2=x*(10.-y);w2/=100.;
    w3=(10.-x)*y;w3/=100.;w4=x*y;w4/=100.;
        k=0;
        for(i=0;i<Nge;i++){
          for(j=0;j<Nge;j++){
            recive[i][j]=recive1[ra+i][rb+j]*w1+recive1[ra+i+1][rb+j]*w2+
                         recive1[ra+i][rb+j+1]*w3+recive1[ra+i+1][rb+j+1]*w4;
            k++;
          }
        }
      
      //printf("moving is done&start merging!Nstar=%d\n",ic);
      sum=0;
        for(i=0;i<Ng;i++){
          for(j=0;j<Ng;j++){
            basisl[i][j]=0;
            for(pl=0;pl<num;pl++){
              ra=i*num+pl;
              for(f=0;f<num;f++){
                rb=j*num+f;
                basisl[i][j]+=recive[ra][rb];//merging divided pixels
               } 
              } 
              sum+=basisl[i][j]*basisl[i][j];
            }
          }
        for(i=0;i<Ng;i++){
          for(j=0;j<Ng;j++){
            //basisl[i][j]/=sqrt(sum);
          }
        }
     //printf("merging is done&staring wite!Nstar=%d\n",ic);
      for(i=0;i<Ng;i++){
        for(j=0;j<Ng;j++){
          f=i*Ng+j;
          basis[ic][f]=basisl[i][j];
        }
      }
    
  }

  free_dmatrix(basisl,0,Ng-1,0,Ng-1);
  free_matrix(recive1,0,Nke-1,0,Nke-1);
  free_matrix(recive,0,Nke-1,0,Nke-1);

}
