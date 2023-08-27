#include "function.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"
#include <string.h>

void method_branches(float ***stars,float **spos,int npc,int Nstar,int Ng0,int Ng,
  float ***PCs,float ***coeff,float gain,int method,int *Nobj,int osam,float snrs,float ***wgt){

	switch(method){
    case(1):iSPCA_entr(stars,spos,npc,Nstar,Ng0,Ng,PCs,coeff,gain,Nobj,snrs,osam,wgt);break;
    case(2):SPCA_entr(stars,spos,npc,Nstar,Ng0,Ng,PCs,coeff,gain,Nobj,snrs,wgt);break;
    case(3):EMPCA_entr(stars,spos,npc,Nstar,Ng0,Ng,PCs,coeff,gain,Nobj,snrs,wgt);break;
    //case(4):PCA(stars,spos,npc,Nstar,Ng0,PCs,coeff,coefferr);break;
    default:{
    	printf("please choose the method you want to use to construct PCs for PSF modeling: 1->iSPCA, 2->SPCA, 3->EMPCA, 4->PCA");
    	exit(EXIT_FAILURE);
    }
  } 
}


void psf_fit(float *get_stars,float *get_spos,int npc,int Nstar,int Ng0,int Ng,float *get_PCs,
    float *get_coeff,float gain,int method,int *get_Nobj,int osam,float snrs,float *wgt){
    float ***stars,**spos,***PCs,***coeff,***mask;
    int k,ic,i,j,Nobj;
    printf("Nstar=%d\n",Nstar );
    FILE *fp,*wp;
    stars=f3tensor(0,Nstar,0,Ng0,0,Ng0);
    mask=f3tensor(0,Nstar,0,Ng0,0,Ng0);
    spos=matrix(0,Nstar,0,2);
    PCs=f3tensor(0,npc,0,Ng*osam,0,Ng*osam);
    coeff=f3tensor(0,Nstar,0,npc,0,2);
    for(ic=0;ic<Nstar;ic++){
      for(i=0;i<Ng0;i++)for(j=0;j<Ng0;j++){
        k=ic*Ng0*Ng0+i*Ng0+j;
        stars[ic][i][j]=get_stars[k];
        mask[ic][i][j]=wgt[k];
      }
    }

    int dim[5];char ftmp[100];
    strcpy(ftmp,"fits/branchstars.fits");
    dim[0]=dim[1]=Ng0;dim[2]=Nstar;
    //write_fits_3D(ftmp,mask,dim);
    for(ic=0;ic<Nstar;ic++){
      k=ic*2;
      spos[ic][0]=get_spos[k];spos[ic][1]=get_spos[k+1];
    }

    method_branches(stars,spos,npc,Nstar,Ng0,Ng,PCs,coeff,gain,method,&Nobj,osam,snrs,mask);

    get_Nobj[0]=Nobj;
    printf("get_Nobj=%d\n",get_Nobj[0] );
    //fp=fopen("results/ocoeff.dat","w");
    //wp=fopen("results/oweight.dat","w");
    for(ic=0;ic<get_Nobj[0];ic++){
      k=ic*2;
      get_spos[k]=spos[ic][0];
      get_spos[k+1]=spos[ic][1];
      for(i=0;i<npc;i++){
        j=ic*npc*2+i*2;
        get_coeff[j]=coeff[ic][i][0];
        get_coeff[j+1]=coeff[ic][i][1];
        //fprintf(fp, "%f\t", get_coeff[j]);
        //fprintf(wp, "%f\t", get_coeff[j+1]);
      }//fprintf(fp, "\n" );fprintf(wp, "\n" );
    }//fclose(fp);fclose(wp);
    for(ic=0;ic<npc;ic++){
      for(i=0;i<Ng*osam;i++)for(j=0;j<Ng*osam;j++){
        k=ic*Ng*Ng*osam*osam+i*Ng*osam+j;
        get_PCs[k]=PCs[ic][i][j];
      }
    }
    for(ic=0;ic<Nobj;ic++){
      for(i=0;i<Ng0;i++)for(j=0;j<Ng0;j++){
        k=ic*Ng0*Ng0+i*Ng0+j;
        get_stars[k]=stars[ic][i][j];
      }
    }
    strcpy(ftmp,"fits/consPCs.fits");
    dim[0]=dim[1]=Ng;dim[2]=npc;
    //write_fits_3D(ftmp,PCs,dim);
    printf("out method_branches.c: psf_fit\n");
    free_f3tensor(stars,0,Nstar,0,Ng0,0,Ng0);
    free_f3tensor(mask,0,Nstar,0,Ng0,0,Ng0);
    free_matrix(spos,0,Nstar,0,2);
    free_f3tensor(PCs,0,npc,0,Ng*osam,0,Ng*osam);
    free_f3tensor(coeff,0,Nstar,0,npc,0,2);
}


void web_method_branches(float ***stars,float **spos,int npc,int Nstar,int Ng0,int Ng,
  float ***PCs,float ***coeff,float gain,int method,int *Nobj,int osam,
  float snrs,float ***imdoel,float ***wgt,float **wcoff){

  switch(method){
    case(1):webiSPCA_entr(stars,spos,npc,Nstar,Ng0,Ng,PCs,coeff,gain,Nobj,snrs,osam,imdoel,wgt,wcoff);break;
    case(2):webSPCA_entr(stars,spos,npc,Nstar,Ng0,Ng,PCs,coeff,gain,Nobj,snrs,imdoel,wgt);break;
    case(3):webEMPCA_entr(stars,spos,npc,Nstar,Ng0,Ng,PCs,coeff,gain,Nobj,snrs,imdoel,wgt);break;
    //case(4):PCA(stars,spos,npc,Nstar,Ng0,PCs,coeff,coefferr);break;
    default:{
      printf("please choose the method you want to use to construct PCs for PSF modeling: 1->iSPCA, 2->SPCA, 3->EMPCA, 4->PCA");
      exit(EXIT_FAILURE);
    }
  } 
}


void web_psf_fit(float *get_stars,float *get_spos,int npc,int Nstar,int Ng0,int Ng,
  float *get_PCs,float *get_coeff,float gain,int method,int *get_Nobj,int osam,
  float snrs,float *imodel,float *wgt,float *get_wcoff){
    printf("Nstar=%d\n",Nstar );
    float ***stars,**spos,***PCs,***coeff,***models,***mask,**wcoff;
    int k,ic,i,j,Nobj,Ngm=Ng0+5;
    printf("Nstar=%d\n",Nstar );
    FILE *fp,*wp;
    stars=f3tensor(0,Nstar,0,Ng0,0,Ng0);
    mask=f3tensor(0,Nstar,0,Ng0,0,Ng0);
    spos=matrix(0,Nstar,0,2);
    PCs=f3tensor(0,npc,0,Ng*osam,0,Ng*osam);
    coeff=f3tensor(0,Nstar,0,npc,0,2);
    wcoff=matrix(0,Nstar,0,2);
    models=f3tensor(0,Nstar,0,Ngm*osam,0,Ngm*osam);
    for(ic=0;ic<Nstar;ic++){
      for(i=0;i<Ng0;i++)for(j=0;j<Ng0;j++){
        k=ic*Ng0*Ng0+i*Ng0+j;
        stars[ic][i][j]=get_stars[k];
        mask[ic][i][j]=wgt[k];
      }
    }printf("mask\n");
    for(ic=0;ic<Nstar;ic++){
      for(i=0;i<Ngm*osam;i++)for(j=0;j<Ngm*osam;j++){
        k=ic*Ngm*osam*Ngm*osam+i*Ngm*osam+j;
        models[ic][i][j]=imodel[k];
      }
    }printf("imodel\n");

    int dim[5];char ftmp[100];
    strcpy(ftmp,"fits/branchstars.fits");
    dim[0]=dim[1]=Ngm*osam;dim[2]=Nstar;
    //write_fits_3D(ftmp,models,dim);
    for(ic=0;ic<Nstar;ic++){
      k=ic*2;
      spos[ic][0]=get_spos[k];spos[ic][1]=get_spos[k+1];
    }

    web_method_branches(stars,spos,npc,Nstar,Ng0,Ng,PCs,coeff,gain,method,&Nobj,osam,snrs,models,mask,wcoff);

    get_Nobj[0]=Nobj;
    printf("get_Nobj=%d\n",get_Nobj[0] );
    //fp=fopen("results/ocoeff.dat","w");
    //wp=fopen("results/oweight.dat","w");
    for(ic=0;ic<get_Nobj[0];ic++){
      k=ic*2;
      get_spos[k]=spos[ic][0];get_spos[k+1]=spos[ic][1];
      for(i=0;i<npc;i++){
        j=ic*npc*2+i*2;
        get_coeff[j]=coeff[ic][i][0];
        get_coeff[j+1]=coeff[ic][i][1];
        //fprintf(fp, "%f\t", get_coeff[j]);
        //fprintf(wp, "%f\t", get_coeff[j+1]);
      }//fprintf(fp, "\n" );fprintf(wp, "\n" );
      get_wcoff[k]=wcoff[ic][0];get_wcoff[k+1]=wcoff[ic][1];
    }//fclose(fp);fclose(wp);
    for(ic=0;ic<npc;ic++){
      for(i=0;i<Ng*osam;i++)for(j=0;j<Ng*osam;j++){
        k=ic*Ng*osam*Ng*osam+i*Ng*osam+j;
        get_PCs[k]=PCs[ic][i][j];
      }
    }
    for(ic=0;ic<Nobj;ic++){
      for(i=0;i<Ng0;i++)for(j=0;j<Ng0;j++){
        k=ic*Ng0*Ng0+i*Ng0+j;
        get_stars[k]=stars[ic][i][j];
        wgt[k]=mask[ic][i][j];
      }
    }
    for(ic=0;ic<Nobj;ic++){
      for(i=0;i<Ngm*osam;i++)for(j=0;j<Ngm*osam;j++){
        k=ic*Ngm*osam*Ngm*osam+i*Ngm*osam+j;
        imodel[k]=models[ic][i][j];
      }
    }printf("imodel\n");
    strcpy(ftmp,"fits/consPCs.fits");
    dim[0]=dim[1]=Ng*osam;dim[2]=npc;
    //write_fits_3D(ftmp,PCs,dim);
    printf("out method_branches.c: psf_fit\n");

    free_f3tensor(stars,0,Nstar,0,Ng0,0,Ng0);
    free_f3tensor(mask,0,Nstar,0,Ng0,0,Ng0);
    free_matrix(spos,0,Nstar,0,2);
    free_f3tensor(PCs,0,npc,0,Ng*osam,0,Ng*osam);
    free_f3tensor(coeff,0,Nstar,0,npc,0,2);
    free_f3tensor(models,0,Nstar,0,Ngm*osam,0,Ngm*osam);
    free_matrix(wcoff,0,Nstar,0,2);

}





