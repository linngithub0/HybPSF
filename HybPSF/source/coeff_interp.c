#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nrutil.h"
#include "function.h"
void coeff_interp(int order_PC,int order_poly,int Nstar,float **spos,float **Coeff,
	  float **weight,float ***PCs,int ngp,float ***rPSFs,int Nobj,float **psfpos){
	void npsf(int order_PC,int order_poly,int Nstar,float **pos,float **Coeff,
    float **weight,double **dcoeff,float ***PCs,int ngp);
	void psf_interp(int Nobj,int order_poly,int order_PC,int ngp,float ***rPSFs,
    float **gpos,double **dcoeff,float ***PCs);
	int i,l0=0;
	for(i=order_poly+1;i>0;i--)l0+=i;
  double **dcoeff;
  dcoeff=dmatrix(0,order_PC,0,l0);

    printf("start npsf\n");
	npsf(order_PC,order_poly,Nstar,spos,Coeff,weight,dcoeff,PCs,ngp);
	printf("npsf down start psf_interp\n");
	psf_interp(Nobj,order_poly,order_PC,ngp,rPSFs,psfpos,dcoeff,PCs);
    printf("psf_interp down\n");
  free_dmatrix(dcoeff,0,order_PC,0,l0);
}

void web_coeff_interp(int order_PC,int order_poly,int Nstar,float **spos,float **Coeff,
	  float **weight,float ***PCs,int ngp,float ***rPSFs,int Nobj,float **psfpos){
	void web_npsf(int order_PC,int order_poly,int Nstar,float **pos,float **Coeff,
    float **weight,double **dcoeff,int ngp);
	void web_psf_interp(int Nobj,int order_poly,int order_PC,int ngp,float ***rPSFs,
    float **gpos,double **dcoeff,float ***PCs);
	int i,l0=0;
	for(i=order_poly+1;i>0;i--)l0+=i;
  double **dcoeff;
  dcoeff=dmatrix(0,order_PC,0,l0);

  printf("start npsf\n");
	web_npsf(order_PC,order_poly,Nstar,spos,Coeff,weight,dcoeff,ngp);
	printf("npsf down start psf_interp\n");
	web_psf_interp(Nobj,order_poly,order_PC,ngp,rPSFs,psfpos,dcoeff,PCs);
    printf("psf_interp down\n");
  free_dmatrix(dcoeff,0,order_PC,0,l0);
}


void pix_interp(int Nstar,float **spos,float ***stars,int ngp,int order_poly,float ***rPSFs,int Nobj,float **psfpos){
	void pix_polypre(int Nstar,float **spos,float ***stars,int ngp,int order_poly,double **dcoeff);
	void pix_polyint(int ngp,float ***rPSFs,int Nobj,float **psfpos,int order_poly,double **dcoeff);
	int i,l0=0;
	for(i=order_poly+1;i>0;i--)l0+=i;
  double **dcoeff;
  dcoeff=dmatrix(0,ngp*ngp,0,l0);
  printf("pix_polypre\n");printf("Nstar=%d,Nobj=%d,ngp=%d,order_poly=%d\n",Nstar,Nobj,ngp,order_poly );
  pix_polypre(Nstar,spos,stars,ngp,order_poly,dcoeff);
  printf("pix_polypre down start pix_polyint\n");
  pix_polyint(ngp,rPSFs,Nobj,psfpos,order_poly,dcoeff);
  printf("pix_polyint down\n");

  free_dmatrix(dcoeff,0,ngp*ngp,0,l0);

}


void Krig_interp(int order_PC,int Nstar,float **spos,float **Coeff,float **Coefferr,float ***PCs,
	int ngp,float ***rPSFs,int Nobj,float **psfpos){
	void krig_prepare(float **x,float **y,float **error,double ***vo,double **voy,
      int **indx,double *alpha,int npt,int nbf);
  void creat_star(float **basisf,float **x,double ***vo,double **voy,int **indx,
    double *alpha,int npt, int nbf,int ng,int Mp,float *xt,float **out_star);
	double ***vo,**voy,*alpha;
	float **basisf;
	int **indx,Mp=ngp*ngp,i,j,ic;
  vo=d3tensor(0,order_PC,0,Nstar,0,Nstar);
  voy=dmatrix(0,order_PC,0,Nstar);
  indx=imatrix(0,order_PC,0,Nstar);      
  alpha=dvector(0,order_PC);
  basisf=matrix(0,order_PC,0,Mp);
  for(ic=0;ic<order_PC;ic++){
  	for(i=0;i<ngp;i++)for(j=0;j<ngp;j++){
  		basisf[ic][i*ngp+j]=PCs[ic][i][j];
  	}
  }
	krig_prepare(spos,Coeff,Coefferr,vo,voy,indx,alpha,Nstar,order_PC);
  for(ic=0;ic<Nobj;ic++){
    creat_star(basisf,spos,vo,voy,indx,alpha,Nstar,order_PC,ngp,Mp,psfpos[ic],rPSFs[ic]);
  }
	

	free_d3tensor(vo,0,order_PC,0,Nstar,0,Nstar);
  free_dmatrix(voy,0,order_PC,0,Nstar);
  free_imatrix(indx,0,order_PC,0,Nstar);      
  free_dvector(alpha,0,order_PC);
  free_matrix(basisf,0,order_PC,0,Mp);
}

void psf_rec(int npc,int order_poly,int Nstar,float *get_spos,float *get_Coeff,
	  float *get_PCs,int ngp,float *get_rPSFs,int Nobj,float *get_psfpos,int method){
	void coeff_interp(int order_PC,int order_poly,int Nstar,float **spos,float **Coeff,
	  float **weight,float ***PCs,int ngp,float ***rPSFs,int Nobj,float **psfpos);
	void pix_interp(int Nstar,float **spos,float ***stars,int ngp,int order_poly,
		float ***rPSFs,int Nobj,float **psfpos);
	void Krig_interp(int order_PC,int Nstar,float **spos,float **Coeff,float **Coefferr,float ***PCs,
	int ngp,float ***rPSFs,int Nobj,float **psfpos);
	int ic,i,j,k,l;
	float **spos,***PCs,***rPSFs,**psfpos;
	float **Coeff,**weight,**Coefferr,***stars;
	FILE *fp,*wp;
	spos=matrix(0,Nstar,0,2);
	Coeff=matrix(0,Nstar,0,npc);
	weight=matrix(0,Nstar,0,npc);
	Coefferr=matrix(0,Nstar,0,npc);
	PCs=f3tensor(0,npc,0,ngp,0,ngp);
	rPSFs=f3tensor(0,Nobj,0,ngp,0,ngp);
	psfpos=matrix(0,Nobj,0,2);
	stars=f3tensor(0,Nstar,0,ngp,0,ngp);
	//fp=fopen("results/gcoeff.dat","w");
	//wp=fopen("results/gweight.dat","w");
	for(ic=0;ic<Nstar;ic++){
	  for(i=0;i<npc;i++){
	  	k=ic*npc*2+i*2;
	  	Coeff[ic][i]=get_Coeff[k];
	  	weight[ic][i]=1./get_Coeff[k+1];
	  	Coefferr[ic][i]=get_Coeff[k+1];
	  	//fprintf(fp, "%f\t",Coeff[ic][i] );
	  	//fprintf(wp, "%f\t",weight[ic][i] );
	  }//fprintf(fp, "\n");fprintf(wp, "\n" );
	  j=ic*2;
	  spos[ic][0]=get_spos[j];
	  spos[ic][1]=get_spos[j+1];
	}//fclose(fp);fclose(wp);
	for(ic=0;ic<npc;ic++){
	  for(i=0;i<ngp;i++)for(j=0;j<ngp;j++){
	  	k=ic*ngp*ngp+i*ngp+j;
	  	PCs[ic][i][j]=get_PCs[k];
	  }
	}

	int dim[5];char ftmp[100];
    strcpy(ftmp,"fits/get_PCs.fits");
    dim[0]=dim[1]=ngp;dim[2]=npc;
    //write_fits_3D(ftmp,PCs,dim);
	for(ic=0;ic<Nobj;ic++){
	  k=ic*2;
	  psfpos[ic][0]=get_psfpos[k];
	  psfpos[ic][1]=get_psfpos[k+1];
	}
	for(ic=0;ic<Nstar;ic++){
		for(i=0;i<ngp;i++)for(j=0;j<ngp;j++){
			stars[ic][i][j]=0;
			for(l=0;l<npc;l++){
				stars[ic][i][j]+=PCs[l][i][j]*Coeff[ic][l];
			}
		}
	}

    if(method==1){
    	coeff_interp(npc,order_poly,Nstar,spos,Coeff,weight,PCs,ngp,rPSFs,Nobj,psfpos);
    }
    if(method==2){
    	pix_interp(Nstar,spos,stars,ngp,order_poly,rPSFs,Nobj,psfpos);
    }
    if(method==3){
    	Krig_interp(npc,Nstar,spos,Coeff,Coefferr,PCs,ngp,rPSFs,Nobj,psfpos);
    }

	for(ic=0;ic<Nobj;ic++){
	  for(i=0;i<ngp;i++)for(j=0;j<ngp;j++){
	  	k=ic*ngp*ngp+i*ngp+j;
	  	get_rPSFs[k]=rPSFs[ic][i][j];
	  }
	}

	strcpy(ftmp,"fits/rPSFs.fits");
    dim[0]=dim[1]=ngp;dim[2]=Nobj;
    //write_fits_3D(ftmp,rPSFs,dim);
  free_matrix(spos,0,Nstar,0,2);
	free_matrix(Coeff,0,Nstar,0,npc);
	free_matrix(weight,0,Nstar,0,npc);
	free_matrix(Coefferr,0,Nstar,0,npc);
	free_f3tensor(PCs,0,npc,0,ngp,0,ngp);
	free_f3tensor(rPSFs,0,Nobj,0,ngp,0,ngp);
	free_matrix(psfpos,0,Nobj,0,2);
	free_f3tensor(stars,0,Nstar,0,ngp,0,ngp);
}


void web_psf_rec(int npc,int order_poly,int Nstar,float *get_spos,float *get_Coeff,
	  float *get_PCs,int ngp,float *get_rPSFs,int Nobj,float *get_psfpos,int method){
	void web_coeff_interp(int order_PC,int order_poly,int Nstar,float **spos,float **Coeff,
	  float **weight,float ***PCs,int ngp,float ***rPSFs,int Nobj,float **psfpos);
	void pix_interp(int Nstar,float **spos,float ***stars,int ngp,int order_poly,
		float ***rPSFs,int Nobj,float **psfpos);
	void Krig_interp(int order_PC,int Nstar,float **spos,float **Coeff,float **Coefferr,float ***PCs,
	int ngp,float ***rPSFs,int Nobj,float **psfpos);
	int ic,i,j,k,l;
	float **spos,***PCs,***rPSFs,**psfpos;
	float **Coeff,**weight,**Coefferr,***stars;
	FILE *fp,*wp;
	spos=matrix(0,Nstar,0,2);
	Coeff=matrix(0,Nstar,0,npc);
	weight=matrix(0,Nstar,0,npc);
	Coefferr=matrix(0,Nstar,0,npc);
	PCs=f3tensor(0,npc,0,ngp,0,ngp);
	rPSFs=f3tensor(0,Nobj,0,ngp,0,ngp);
	psfpos=matrix(0,Nobj,0,2);
	stars=f3tensor(0,Nstar,0,ngp,0,ngp);
	//fp=fopen("results/gcoeff.dat","w");
	//wp=fopen("results/gweight.dat","w");
	for(ic=0;ic<Nstar;ic++){
	  for(i=0;i<npc;i++){
	  	k=ic*npc*2+i*2;
	  	Coeff[ic][i]=get_Coeff[k];
	  	weight[ic][i]=1./get_Coeff[k+1];
	  	Coefferr[ic][i]=get_Coeff[k+1];
	  	//fprintf(fp, "%f\t",Coeff[ic][i] );
	  	//fprintf(wp, "%f\t",weight[ic][i] );
	  }//fprintf(fp, "\n");fprintf(wp, "\n" );
	  j=ic*2;
	  spos[ic][0]=get_spos[j];
	  spos[ic][1]=get_spos[j+1];
	}//fclose(fp);fclose(wp);
	for(ic=0;ic<npc;ic++){
	  for(i=0;i<ngp;i++)for(j=0;j<ngp;j++){
	  	k=ic*ngp*ngp+i*ngp+j;
	  	PCs[ic][i][j]=get_PCs[k];
	  }
	}

	int dim[5];char ftmp[100];
    strcpy(ftmp,"fits/get_PCs.fits");
    dim[0]=dim[1]=ngp;dim[2]=npc;
    //write_fits_3D(ftmp,PCs,dim);
	for(ic=0;ic<Nobj;ic++){
	  k=ic*2;
	  psfpos[ic][0]=get_psfpos[k];
	  psfpos[ic][1]=get_psfpos[k+1];
	}
	for(ic=0;ic<Nstar;ic++){
		for(i=0;i<ngp;i++)for(j=0;j<ngp;j++){
			stars[ic][i][j]=0;
			for(l=0;l<npc;l++){
				stars[ic][i][j]+=PCs[l][i][j]*Coeff[ic][l];
			}
		}
	}

    if(method==1){
    	web_coeff_interp(npc,order_poly,Nstar,spos,Coeff,weight,PCs,ngp,rPSFs,Nobj,psfpos);
    }
    if(method==2){
    	pix_interp(Nstar,spos,stars,ngp,order_poly,rPSFs,Nobj,psfpos);
    }
    if(method==3){
    	Krig_interp(npc,Nstar,spos,Coeff,Coefferr,PCs,ngp,rPSFs,Nobj,psfpos);
    }

	for(ic=0;ic<Nobj;ic++){
	  for(i=0;i<ngp;i++)for(j=0;j<ngp;j++){
	  	k=ic*ngp*ngp+i*ngp+j;
	  	get_rPSFs[k]=rPSFs[ic][i][j];
	  }
	}

	strcpy(ftmp,"fits/rPSFs.fits");
    dim[0]=dim[1]=ngp;dim[2]=Nobj;
    //write_fits_3D(ftmp,rPSFs,dim);
  free_matrix(spos,0,Nstar,0,2);
	free_matrix(Coeff,0,Nstar,0,npc);
	free_matrix(weight,0,Nstar,0,npc);
	free_matrix(Coefferr,0,Nstar,0,npc);
	free_f3tensor(PCs,0,npc,0,ngp,0,ngp);
	free_f3tensor(rPSFs,0,Nobj,0,ngp,0,ngp);
	free_matrix(psfpos,0,Nobj,0,2);
	free_f3tensor(stars,0,Nstar,0,ngp,0,ngp);
}

void estimate_parapy(double *stars,double *cent,int Nstar,int Ng, int nbound,
       double *paras){
	void estimate_para(double **star,double **center,int Nstar,int Ng, int nbound,
       double *para);
	//printf("Nstar=%d,Ng=%d,nbound=%d\n",Nstar,Ng,nbound );
  double **star,**center,*para;
  float ***wstar;
  int ic,i,j,k;
  star=dmatrix(0,Nstar,0,Ng*Ng);
  center=dmatrix(0,Nstar,0,2);
  para=dvector(0,10);
  //wstar=f3tensor(0,Nstar,0,Ng,0,Ng);
  for(ic=0;ic<Nstar;ic++){
  	for(i=0;i<Ng;i++)for(j=0;j<Ng;j++){
  		star[ic][i*Ng+j]=stars[ic*Ng*Ng+i*Ng+j];
  		//wstar[ic][i][j]=(float)stars[ic*Ng*Ng+i*Ng+j];
  	}
  	center[ic][0]=cent[ic*2];
  	center[ic][1]=cent[ic*2+1];
  }
    int dim[5];
    dim[0]=dim[1]=Ng;dim[2]=Nstar;
    //write_fits_3D("fits/instar.fits",wstar,dim);


  estimate_para(star,center,Nstar,Ng,nbound,para);
  paras[0]=para[0];
  paras[1]=para[1];
  paras[2]=para[2];
  
  free_dmatrix(star,0,Nstar,0,Ng*Ng);
  free_dmatrix(center,0,Nstar,0,2);
  free_dvector(para,0,10);
}






