#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"
#include <string.h>
#include "function.h"
void iSPCA(double **images,double **weight,double **basisf,int Mp,int Nstar,int nbf,double **cent,
      double ***compbf,double **coeff,int Nb,int Nk,int num,double **compbfr,int Ng){

  void Prepare_Dmatrix(double **weight,double ***basis,int Mp,int nimg,int nbf,
           double ***sigmaka);
  void creat_Dmatrix(double ***sigmaka,int Mp,int nimg,int nbf,double *coeff,
         double **a,double *b,double **weight,double ***basis,double **images);
  void lubksb(double **a, int n, int *indx, double b[]);
  void ludcmp(double **a, int n, int *indx, double *d);
  void Rebasis(double **basisf0,int nbf0,int Mp,double **basisf,double **coeff,int *nbf);
  void mov_mer_pc(double **erbasisf,double ***basis,int Nke,int Ng,
                int num,double **cent,int Nstar,int nbf);
  void EC(double *xj,double **phi,double *weight,int Mp, int nbf,double *coeff);
  int write_fits_3D(char *argv,float ***stars,int *dim);
  void mop(double *erbasisf,double **basis,int Nke,int Ng,
                int num,double **cent,int Nstar);
  float ran1(long *idum);

    int i,j,k,l,ipc,nbfs,ic,k0,Nke=Nk*num,Mpke=Nke*Nke;
    double **star,**xij,**coff,absphi;
    double **a,*b,d,*compbf0,*D,sum,sum1,sum2=0,detbf,sumt,***basis,**dbasisf;
    double ***sigmaka,**rbasisf, **rcoeff,detchi2,**coefferr,**com;
    int *indx,niter=0,nmax=700,rnbf,dim[3];
    float ***ebasisf;
    char *ftmp,fimage[200],dat[200];
    long iseed=-3;

    xij=dmatrix(0,Nstar,0,Mp);
    a=dmatrix(0,nbf,0,nbf);
    b=dvector(0,nbf);
    D=dvector(0,nbf);
    indx=ivector(0,nbf);
    sigmaka=d3tensor(0,Nstar,0,nbf,0,nbf); 
    
    coff=dmatrix(0,Nstar,0,Nb);
    rbasisf=dmatrix(0,nbf,0,Mpke);
    rcoeff=dmatrix(0,nbf,0,nbf);
    star=dmatrix(0,Nstar,0,Mp);
    basis=d3tensor(0,Nstar,0,nbf,0,Mp);
    coefferr=dmatrix(0,Nstar,0,Nstar);
    dbasisf=dmatrix(0,nbf,0,Mpke);
    ebasisf=f3tensor(0,nbf,0,Nke,0,Nke);
    com=dmatrix(0,Nstar,0,Mp);
    
    //printf("Ng=%d,Nk=%d,num=%d,Nb=%d,Mp=%d,Nstar=%d,nbf=%d\n",Ng,Nk,num,Nb,Mp,Nstar,nbf);

    for(ic=0;ic<Nstar;ic++){
      for(ipc=0;ipc<Nb;ipc++){
        for(k=0;k<Mp;k++)compbf[ic][ipc][k]=1+ran1(&iseed);//images[ic][k];
      }
    }
    for(i=0;i<nbf;i++){
      for(k=0;k<Mpke;k++)dbasisf[i][k]=basisf[i][k];
    }

    /*float ***basisfs;
    basisfs=f3tensor(0,nbf,0,Nke,0,Nke);
    for(ic=0;ic<nbf;ic++){
      for(i=0;i<Nke;i++){
        for(j=0;j<Nke;j++){
          basisfs[ic][i][j]=basisf[ic][i*Nke+j];
        }
      }
    }
    dim[0]=dim[1]=Nke;dim[2]=nbf;
    write_fits_3D("fits/basisf.fits",basisfs,dim);*/

    Rebasis(dbasisf,nbf,Mpke,rbasisf,rcoeff,&rnbf); // Given nbf basis function basisf[][],this subroutine
    // creats rnbf orthogonal basis functions, rbasisf[][].
    printf("erbasis is done!\n");

    mov_mer_pc(rbasisf,basis,Nke,Ng,num,cent,Nstar,rnbf);
    Prepare_Dmatrix(weight,basis,Mp,Nstar,rnbf,sigmaka);//provide sigmaka only
    printf("basisf is done\n");
    
    niter=0;
    do{
      for(i=0;i<Nstar;i++){
      for(j=0;j<Mp;j++){
        star[i][j]=xij[i][j]=images[i][j];
      }
    }
      for(i=0;i<Nstar;i++){
        EC(xij[i],compbf[i],weight[i],Mp,Nb,coff[i]);
      }
      for(i=0;i<Nstar;i++){
        for(j=0;j<Nb;j++)coeff[j][i]=coff[i][j];
      }

      for(ipc=0;ipc<Nb;ipc++){
        creat_Dmatrix(sigmaka,Mp,Nstar,rnbf,coeff[ipc],a,b,weight,basis,star);

        ludcmp(a, rnbf, indx, &d);
        lubksb(a, rnbf, indx, b);


        sumt=0;
        for(k=0;k<Mpke;k++){
            compbfr[ipc][k]=0;
            for(i=0;i<rnbf;i++)compbfr[ipc][k]+=rbasisf[i][k]*b[i];
            //sumt+=compbfr[ipc][k]*compbfr[ipc][k];
         }
          //sumt=sqrt(sumt);
          //for(k=0;k<Mpke;k++)compbfr[ipc][k]/=sumt;
        //use the calculated PCs
        //mop(compbfr[ipc],com,Nke,Ng,num,cent,Nstar);
       for(ic=0;ic<Nstar;ic++){
          for(k=0;k<Mp;k++){
            com[ic][k]=0;
            for(l=0;l<rnbf;l++){
              com[ic][k]+=basis[ic][l][k]*b[l];
            }
          }
        }


         //substract PCs from stars
         
        for(i=0;i<Nstar;i++)for(k=0;k<Mp;k++)star[i][k]-=com[i][k]*coeff[ipc][i];
         
      }
      //force orthogonality PCs
      for(k0=0;k0<Nb-1;k0++){
        absphi=0;
        for(j=0;j<Mpke;j++)absphi+=compbfr[k0][j]*compbfr[k0][j]; //normalization
        absphi=sqrt(absphi);
        for(j=0;j<Mpke;j++)compbfr[k0][j]/=absphi;
        for(ic=0;ic<Nstar;ic++)coeff[k0][ic]/=absphi;
        for(k=k0+1;k<Nb;k++){
          sum=0;
          for(j=0;j<Mpke;j++)sum+=compbfr[k0][j]*compbfr[k][j];
          for(j=0;j<Mpke;j++)compbfr[k][j]-=compbfr[k0][j]*sum;
        }
      }
        
      //renormalize the last PC
      absphi=0;
      for(j=0;j<Mpke;j++)absphi+=compbfr[Nb-1][j]*compbfr[Nb-1][j];
      absphi=sqrt(absphi);
      for(j=0;j<Mpke;j++)compbfr[Nb-1][j]/=absphi;
      for(ic=0;ic<Nstar;ic++)coeff[Nb-1][ic]/=absphi;


      mov_mer_pc(compbfr,compbf,Nke,Ng,num,cent,Nstar,Nb);
      //calculate chi^2
      sum1=0;
      for(j=0;j<Mp;j++){
        for(i=0;i<Nstar;i++){
          xij[i][j]=0;
          for(k=0;k<Nb;k++)xij[i][j]+=compbf[i][k][j]*coeff[k][i];
          sum1+=(xij[i][j]-images[i][j])*(xij[i][j]-images[i][j])*weight[i][j];
        }
      }
      detchi2=sum2/sum1-1;
      detchi2=fabs(detchi2);
      sum2=sum1;
      niter++;
      if(niter%50==0){
        printf("niter=%d detchi2=%e   chi2=%e\n",niter,detchi2,sum1);
      }

      
    }while(niter<nmax && detchi2>1.e-12);

    free_dmatrix(xij,0,Nstar,0,Mp);
    free_dmatrix(a,0,nbf,0,nbf);
    free_dvector(b,0,nbf);
    free_dvector(D,0,nbf);
    free_ivector(indx,0,nbf);
    free_d3tensor(sigmaka,0,Nstar,0,nbf,0,nbf); 
    
    free_dmatrix(coff,0,Nstar,0,Nb);
    free_dmatrix(rbasisf,0,nbf,0,Mpke);
    free_dmatrix(rcoeff,0,nbf,0,nbf);
    free_dmatrix(star,0,Nstar,0,Mp);
    free_d3tensor(basis,0,Nstar,0,nbf,0,Mp);
    free_dmatrix(coefferr,0,Nstar,0,Nstar);
    free_dmatrix(dbasisf,0,nbf,0,Mpke);
    free_f3tensor(ebasisf,0,nbf,0,Nke,0,Nke);
    free_dmatrix(com,0,Nstar,0,Mp);



}


void Prepare_Dmatrix(double **weight,double ***basis,int Mp,int nimg,int nbf,
         double ***sigmaka)
{
  int i,j,k,l;
  int ic;
  double sum;
  int write_fits_3D(char *argv,float ***stars,int *dim);
  int dim[3];char ftmp[200];
    for(j=0;j<nbf;j++){
      for(l=0;l<nbf;l++){
        for(i=0;i<nimg;i++){
    sigmaka[i][j][l]=0;
    for(k=0;k<Mp;k++)sigmaka[i][j][l]+=weight[i][k]*basis[i][j][k]*basis[i][l][k];
        }
      }
    } 

}

void creat_Dmatrix(double ***sigmaka,int Mp,int nimg,int nbf,double *coeff,
       double **a,double *b,double **weight,double ***basis,double **images)
{
  int i,j,k,l;
  double sum,**sigmakb;
  sigmakb=dmatrix(0,nbf,0,nimg);
  /*for(j=0;j<nbf;j++){
      for(i=0;i<nimg;i++){
        sigmakb[j][i]=0;
        for(k=0;k<Mp;k++)sigmakb[j][i]+=weight[i][k]*basis[i][j][k]*images[i][k];
      }
  }*/


  for(j=0;j<nbf;j++){
    for(l=j;l<nbf;l++){
      a[j][l]=0;
      for(i=0;i<nimg;i++){
  //sum=0;
  //for(k=0;k<Mp;k++)sum+=weight[i][k]*basisf[j][k]*basisf[l][k];
  a[j][l]+=sigmaka[i][j][l]*coeff[i]*coeff[i];
      }
      a[l][j]=a[j][l];
      //a[j][l]*=2;
    }

    b[j]=0;
    for(i=0;i<nimg;i++){
      sum=0;
      for(k=0;k<Mp;k++)sum+=weight[i][k]*basis[i][j][k]*images[i][k];
      b[j]+=sum*coeff[i];
    } 
    //b[j]*=2;
  }

  free_dmatrix(sigmakb,0,nbf,0,nimg);

}


void EC(double *xj,double **phi,double *weight,int Mp, int nbf,double *coeff)
{
  void matrix_matrix(double **a,int m, int n,double **b,int k, int l,double **c);
  void matrix_vector(double **a,int m, int n,double *b,int k,double *c);
  void matrix_invers(double **a,int N,double **b);
  double **PtV,**P3,**P3in,*tmp;
  int i,j,k;
  PtV=dmatrix(0,Mp,0,nbf);
  P3=dmatrix(0,nbf,0,nbf);
  P3in=dmatrix(0,nbf,0,nbf);
  tmp=dvector(0,nbf);
  
  for(i=0;i<Mp;i++)for(j=0;j<nbf;j++)PtV[i][j]=phi[j][i]*weight[i];
  for(j=0;j<nbf;j++){
    tmp[j]=0;
    for(i=0;i<Mp;i++)tmp[j]+=PtV[i][j]*xj[i];
  }
  //printf("s1\n");
  //matrix_matrix(PtV,Mp,nbf,phi,nbf,Mp,P3);
  for(i=0;i<nbf;i++){
    for(j=0;j<nbf;j++){
      P3[i][j]=0;
      for(k=0;k<Mp;k++)P3[i][j]+=PtV[k][i]*phi[j][k];
    }
  }
  matrix_invers(P3,nbf,P3in);
  for(j=0;j<nbf;j++){
    coeff[j]=0;
    for(i=0;i<nbf;i++)coeff[j]+=P3in[i][j]*tmp[i];
  }
  //matrix_vector(P3in,nbf, nbf,tmp,nbf,coeff);
  
  free_dmatrix(PtV,0,Mp,0,nbf);
  free_dmatrix(P3,0,nbf,0,nbf);
  free_dmatrix(P3in,0,nbf,0,nbf);
  free_dvector(tmp,0,nbf);
}







void hiSPCA(double **images,double **weight,double **basisf,int Mp,int Nstar,int nbf,double **cent,
      double ***compbf,double **coeff,int Nb,int Nk,int num,double **compbfr,int Ng,
      double **wbpsf,double *wbcoeff){


  void Prepare_Dmatrix(double **weight,double ***basis,int Mp,int nimg,int nbf,
           double ***sigmaka);
  void creat_Dmatrix(double ***sigmaka,int Mp,int nimg,int nbf,double *coeff,
         double **a,double *b,double **weight,double ***basis,double **images);
  void lubksb(double **a, int n, int *indx, double b[]);
  void ludcmp(double **a, int n, int *indx, double *d);
  void Rebasis(double **basisf0,int nbf0,int Mp,double **basisf,double **coeff,int *nbf);
  void mov_mer_pc(double **erbasisf,double ***basis,int Nke,int Ng,
                int num,double **cent,int Nstar,int nbf);
  void EC(double *xj,double **phi,double *weight,int Mp, int nbf,double *coeff);
  int write_fits_3D(char *argv,float ***stars,int *dim);
  void mop(double *erbasisf,double **basis,int Nke,int Ng,
                int num,double **cent,int Nstar);
  float ran1(long *idum);
  void mcmc_coeff(double **stars,double ***PCs,double **weight,int Nstar,
        int Mp, int nbf,double **coeff);

    int i,j,k,l,ipc,nbfs,ic,k0,Nke=Nk*num,Mpke=Nke*Nke;
    double **star,**xij,**coff,absphi;
    double **a,*b,d,*compbf0,*D,sum,sum1,sum2=0,detbf,sumt,***basis,**dbasisf;
    double ***sigmaka,**rbasisf, **rcoeff,detchi2,**coefferr,**com;
    int *indx,niter=0,nmax=700,rnbf,dim[3];
    float ***ebasisf,**wa;
    char *ftmp,fimage[200],dat[200];
    long iseed=-3;

    double ***hcompbf,**hcoff;
    hcompbf=d3tensor(0,Nstar,0,Nb+1,0,Mp);
    hcoff=dmatrix(0,Nstar,0,Nb+1);
    printf("in hiSPCA\n");
    xij=dmatrix(0,Nstar,0,Mp);
    a=dmatrix(0,nbf,0,nbf);
    b=dvector(0,nbf);
    wa=matrix(0,nbf+1,0,nbf);
    D=dvector(0,nbf);
    indx=ivector(0,nbf);
    sigmaka=d3tensor(0,Nstar,0,nbf,0,nbf); 
    
    coff=dmatrix(0,Nstar,0,Nb);
    rbasisf=dmatrix(0,nbf,0,Mpke);
    rcoeff=dmatrix(0,nbf,0,nbf);
    star=dmatrix(0,Nstar,0,Mp);
    basis=d3tensor(0,Nstar,0,nbf,0,Mp);
    coefferr=dmatrix(0,Nstar,0,Nstar);
    dbasisf=dmatrix(0,nbf,0,Mpke);
    ebasisf=f3tensor(0,nbf,0,Nke,0,Nke);
    com=dmatrix(0,Nstar,0,Mp);
    
    //printf("Ng=%d,Nk=%d,num=%d,Nb=%d,Mp=%d,Nstar=%d,nbf=%d\n",Ng,Nk,num,Nb,Mp,Nstar,nbf);

    for(ic=0;ic<Nstar;ic++){
      /*for(ipc=0;ipc<Nb;ipc++){
        for(k=0;k<Mp;k++)compbf[ic][ipc][k]=1.+ran1(&iseed);//images[ic][k];
      }*/
      for(ipc=0;ipc<Nb;ipc++){
        for(k=0;k<Mp;k++)hcompbf[ic][ipc+1][k]=1.+ran1(&iseed);//images[ic][k];
        /*for(k=0;k<Mp;k++)hcompbf[ic][ipc+1][k]=0;
        for(i=0;i<Ng;i++){
          hcompbf[ic][ipc+1][i*Ng+i]=1.;
        }*/
      }
      for(k=0;k<Mp;k++)hcompbf[ic][0][k]=wbpsf[ic][k];
    }
    for(i=0;i<nbf;i++){
      for(k=0;k<Mpke;k++)dbasisf[i][k]=basisf[i][k];
    }


    Rebasis(dbasisf,nbf,Mpke,rbasisf,rcoeff,&rnbf); // Given nbf basis function basisf[][],this subroutine
    // creats rnbf orthogonal basis functions, rbasisf[][].
    printf("erbasis is done!\n");

    mov_mer_pc(rbasisf,basis,Nke,Ng,num,cent,Nstar,rnbf);
    Prepare_Dmatrix(weight,basis,Mp,Nstar,rnbf,sigmaka);//provide sigmaka only
    printf("basisf is done\n");
    
    niter=0;
    do{
      for(i=0;i<Nstar;i++){
        for(j=0;j<Mp;j++){
          star[i][j]=xij[i][j]=images[i][j];
        }
        EC(xij[i],hcompbf[i],weight[i],Mp,Nb+1,hcoff[i]);
        wbcoeff[i]=hcoff[i][0];//printf("wbcoeff=%f\n", wbcoeff[i]);
        for(j=0;j<Nb;j++){
          coeff[j][i]=hcoff[i][j+1];
        }
        for(k=0;k<Mp;k++){
          star[i][k]-=wbcoeff[i]*wbpsf[i][k];
        }
    }


      for(ipc=0;ipc<Nb;ipc++){
        creat_Dmatrix(sigmaka,Mp,Nstar,rnbf,coeff[ipc],a,b,weight,basis,star);

        ludcmp(a, rnbf, indx, &d);
        lubksb(a, rnbf, indx, b);


        sumt=0;
        for(k=0;k<Mpke;k++){
            compbfr[ipc][k]=0;
            for(i=0;i<rnbf;i++)compbfr[ipc][k]+=rbasisf[i][k]*b[i];
            //sumt+=compbfr[ipc][k]*compbfr[ipc][k];
         }
          //sumt=sqrt(sumt);
          //for(k=0;k<Mpke;k++)compbfr[ipc][k]/=sumt;
        //use the calculated PCs
        //mop(compbfr[ipc],com,Nke,Ng,num,cent,Nstar);
       for(ic=0;ic<Nstar;ic++){
          for(k=0;k<Mp;k++){
            com[ic][k]=0;
            for(l=0;l<rnbf;l++){
              com[ic][k]+=basis[ic][l][k]*b[l];
            }
          }
        }
         //substract PCs from stars
         
        for(i=0;i<Nstar;i++)for(k=0;k<Mp;k++)star[i][k]-=com[i][k]*coeff[ipc][i];
         
      }
      //force orthogonality PCs
      for(k0=0;k0<Nb-1;k0++){
        absphi=0;
        for(j=0;j<Mpke;j++)absphi+=compbfr[k0][j]*compbfr[k0][j]; //normalization
        absphi=sqrt(absphi);
        for(j=0;j<Mpke;j++)compbfr[k0][j]/=absphi;
        for(ic=0;ic<Nstar;ic++)coeff[k0][ic]/=absphi;
        for(k=k0+1;k<Nb;k++){
          sum=0;
          for(j=0;j<Mpke;j++)sum+=compbfr[k0][j]*compbfr[k][j];
          for(j=0;j<Mpke;j++)compbfr[k][j]-=compbfr[k0][j]*sum;
        }
      }
        
      //renormalize the last PC
      absphi=0;
      for(j=0;j<Mpke;j++)absphi+=compbfr[Nb-1][j]*compbfr[Nb-1][j];
      absphi=sqrt(absphi);
      for(j=0;j<Mpke;j++)compbfr[Nb-1][j]/=absphi;
      for(ic=0;ic<Nstar;ic++)coeff[Nb-1][ic]/=absphi;


      mov_mer_pc(compbfr,compbf,Nke,Ng,num,cent,Nstar,Nb);
      for(ic=0;ic<Nstar;ic++){
        for(ipc=0;ipc<Nb;ipc++){
          for(k=0;k<Mp;k++)hcompbf[ic][ipc+1][k]=compbf[ic][ipc][k];//images[ic][k];
        }
      }
      //calculate chi^2
      sum1=0;
      for(j=0;j<Mp;j++){
        for(i=0;i<Nstar;i++){
          xij[i][j]=0;
          xij[i][j]=wbpsf[i][j]*wbcoeff[i];
          for(k=0;k<Nb;k++)xij[i][j]+=compbf[i][k][j]*coeff[k][i];
          sum1+=(xij[i][j]-images[i][j])*(xij[i][j]-images[i][j])*weight[i][j];
        }
      }
      detchi2=sum2/sum1-1;
      detchi2=fabs(detchi2);
      sum2=sum1;
      niter++;
      if(niter%50==0){
        printf("niter=%d detchi2=%e   chi2=%e\n",niter,detchi2,sum1);
      }
      
      
    }while(niter<nmax && detchi2>1.e-12);
    
    /*niter=0;
    do{
      for(i=0;i<Nstar;i++){
      for(j=0;j<Mp;j++){
        star[i][j]=xij[i][j]=images[i][j];
      }
    }

      mcmc_coeff(xij,hcompbf,weight,Nstar,Mp,Nb+1,hcoff);

      for(i=0;i<Nstar;i++){
        wbcoeff[i]=hcoff[i][0];//printf("wbcoeff=%f\n", wbcoeff[i]);
        for(j=0;j<Nb;j++)coeff[j][i]=hcoff[i][j+1];
        for(k=0;k<Mp;k++){
          star[i][k]-=wbcoeff[i]*wbpsf[i][k];
        }
      }

      for(ipc=0;ipc<Nb;ipc++){
        creat_Dmatrix(sigmaka,Mp,Nstar,rnbf,coeff[ipc],a,b,weight,basis,star);

        ludcmp(a, rnbf, indx, &d);
        lubksb(a, rnbf, indx, b);


        sumt=0;
        for(k=0;k<Mpke;k++){
            compbfr[ipc][k]=0;
            for(i=0;i<rnbf;i++)compbfr[ipc][k]+=rbasisf[i][k]*b[i];
            //sumt+=compbfr[ipc][k]*compbfr[ipc][k];
         }
          //sumt=sqrt(sumt);
          //for(k=0;k<Mpke;k++)compbfr[ipc][k]/=sumt;
        //use the calculated PCs
        //mop(compbfr[ipc],com,Nke,Ng,num,cent,Nstar);
       for(ic=0;ic<Nstar;ic++){
          for(k=0;k<Mp;k++){
            com[ic][k]=0;
            for(l=0;l<rnbf;l++){
              com[ic][k]+=basis[ic][l][k]*b[l];
            }
          }
        }
         //substract PCs from stars
         
        for(i=0;i<Nstar;i++)for(k=0;k<Mp;k++)star[i][k]-=com[i][k]*coeff[ipc][i];
         
      }
      //force orthogonality PCs
      for(k0=0;k0<Nb-1;k0++){
        absphi=0;
        for(j=0;j<Mpke;j++)absphi+=compbfr[k0][j]*compbfr[k0][j]; //normalization
        absphi=sqrt(absphi);
        for(j=0;j<Mpke;j++)compbfr[k0][j]/=absphi;
        for(ic=0;ic<Nstar;ic++)coeff[k0][ic]/=absphi;
        for(k=k0+1;k<Nb;k++){
          sum=0;
          for(j=0;j<Mpke;j++)sum+=compbfr[k0][j]*compbfr[k][j];
          for(j=0;j<Mpke;j++)compbfr[k][j]-=compbfr[k0][j]*sum;
        }
      }
        
      //renormalize the last PC
      absphi=0;
      for(j=0;j<Mpke;j++)absphi+=compbfr[Nb-1][j]*compbfr[Nb-1][j];
      absphi=sqrt(absphi);
      for(j=0;j<Mpke;j++)compbfr[Nb-1][j]/=absphi;
      for(ic=0;ic<Nstar;ic++)coeff[Nb-1][ic]/=absphi;


      mov_mer_pc(compbfr,compbf,Nke,Ng,num,cent,Nstar,Nb);
      for(ic=0;ic<Nstar;ic++){
        for(ipc=0;ipc<Nb;ipc++){
          for(k=0;k<Mp;k++)hcompbf[ic][ipc+1][k]=compbf[ic][ipc][k];//images[ic][k];
        }
      }
      //calculate chi^2
      sum1=0;
      for(j=0;j<Mp;j++){
        for(i=0;i<Nstar;i++){
          xij[i][j]=0;
          xij[i][j]=wbpsf[i][j]*wbcoeff[i];
          for(k=0;k<Nb;k++)xij[i][j]+=compbf[i][k][j]*coeff[k][i];
          sum1+=(xij[i][j]-images[i][j])*(xij[i][j]-images[i][j])*weight[i][j];
        }
      }
      detchi2=sum2/sum1-1;
      detchi2=fabs(detchi2);
      sum2=sum1;
      niter++;
      printf("niter=%d detchi2=%e   chi2=%e\n",niter,detchi2,sum1);

      
    }while(niter<20);*/

    free_dmatrix(xij,0,Nstar,0,Mp);
    free_dmatrix(a,0,nbf,0,nbf);
    free_dvector(b,0,nbf);
    free_dvector(D,0,nbf);
    free_ivector(indx,0,nbf);
    free_d3tensor(sigmaka,0,Nstar,0,nbf,0,nbf); 
    
    free_dmatrix(coff,0,Nstar,0,Nb);
    free_dmatrix(rbasisf,0,nbf,0,Mpke);
    free_dmatrix(rcoeff,0,nbf,0,nbf);
    free_dmatrix(star,0,Nstar,0,Mp);
    free_d3tensor(basis,0,Nstar,0,nbf,0,Mp);
    free_dmatrix(coefferr,0,Nstar,0,Nstar);
    free_dmatrix(dbasisf,0,nbf,0,Mpke);
    free_f3tensor(ebasisf,0,nbf,0,Nke,0,Nke);
    free_dmatrix(com,0,Nstar,0,Mp);



}
