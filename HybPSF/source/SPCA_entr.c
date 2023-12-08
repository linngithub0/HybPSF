#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"
#include <string.h>
#include "function.h"
void SPCA_entr(float ***stars,float **spos,int npc,int Nstar0,int Ng0,int Ng,
  float ***PCs,float ***coeff,float gain,int *Nobj,float snrs,float ***wgt){

	int i,j,k,ic,Mp,Mp0,r,count,pl,ra,rb,f,ipc;
	int nbond=(Ng0-Ng)/2+1,mx[13],msign=1,l0=20,nbf,NPCmax,rnbf0,Nstar;
	float ***star0,**start,**starf,**residu,***star,***residual,***ebasisfs;
	double mean,sigma,detx,dety,sum,sumt,rd,beta,para[10],xcen[2],sigmasum;
	double **error0,**weight,SNR[Nstar0],**coefferr,**a;
	double **cent,**center,**ebasisf,***compbf,**compcoeff,**compbfr,***basis;
  FILE *fp,*wp;
	Mp0=Ng0*Ng0;Mp=Ng*Ng;
  NPCmax=Nstar0;
  rnbf0=(l0+1)*(l0+1);
  if(NPCmax<rnbf0)NPCmax=rnbf0;
	star0=f3tensor(0,Nstar0,0,Ng0,0,Ng0);
	star=f3tensor(0,Nstar0,0,Ng0,0,Ng0);
  residual=f3tensor(0,Nstar0,0,Ng0,0,Ng0);
	a=dmatrix(0,Nstar0,0,Mp0);
	error0=dmatrix(0,Nstar0,0,Mp0);
	weight=dmatrix(0,Nstar0,0,Mp0);
	start=matrix(0,Ng0,0,Ng0);
	starf=matrix(0,Ng0,0,Ng0);
	residu=matrix(0,Ng0,0,Ng0);
	cent=dmatrix(0,Nstar0,0,Ng0);
	center=dmatrix(0,Nstar0,0,Ng0);
	compbf=d3tensor(0,Nstar0,0,npc,0,Mp0);
	compcoeff=dmatrix(0,npc,0,Nstar0);
	compbfr=dmatrix(0,npc,0,Mp0);
	basis=d3tensor(0,npc,0,Ng0,0,Ng0);
  coefferr=dmatrix(0,Nstar0,0,npc);
  ebasisf=dmatrix(0,NPCmax,0,Mp);
  printf("method SPCA\n");
	

	for(ic=0;ic<Nstar0;ic++){
	  for(i=0;i<Ng0;i++)for(j=0;j<Ng0;j++){
		  star0[ic][i][j]=stars[ic][i][j];
	  }
	}

  int dim[5];char ftmp[100];
  
  
	for(ic=0;ic<Nstar0;ic++){
	  gaus_estimate(star0[ic],Ng0,Ng0,&mean,&sigma); //estimate the background noise as gasussian randoms
    printf("ic=%d,mean=%f,sigma=%f\n",ic+1,mean,sigma );
    k=0;
	  for(i=0;i<Ng0;i++){
	  	for(j=0;j<Ng0;j++){
	  	  a[ic][k]=star0[ic][i][j]-mean;//subtract the bacground sky intensity
        start[i][j]=a[ic][k];
        k++;
      }
	  }
	  star_gaus(start,Ng0,Ng0,starf,residu,para);
	  xcen[0]=para[3];xcen[1]=para[4];printf("x=%f,y=%f\n",xcen[0],xcen[1] );
    //center[ic][0]=xcen[0];center[ic][1]=xcen[1];
    k=0;sum=0;
    for(i=0;i<Ng0;i++){
      for(j=0;j<Ng0;j++){
        start[i][j]=star0[ic][i][j]-mean;
        detx=sigma*sigma*gain*gain;
        dety=gain*fabs(start[i][j]);
        if(dety<detx)error0[ic][k]=sigma;
        else{error0[ic][k]=sqrt(dety+detx)/gain;}//estimate poisson noise
        sum +=start[i][j];
        k++;
      }
    }
    sigmasum=0;
    for(k=0;k<Mp0;k++){
      sigmasum+=error0[ic][k]*error0[ic][k];
    }SNR[ic]=sum/sqrt(sigmasum);//calculate the SNR of each image
    //printf("snr=%f\n",SNR[ic] );
	}

  count=0;
  for(ic=0;ic<Nstar0;ic++){
    if(SNR[ic]>=snrs){   //select the high snr images
      for(i=0;i<Ng0;i++)for(j=0;j<Ng0;j++){
        star[count][i][j]=a[ic][i*Ng0+j];
        wgt[count][i][j]=wgt[ic][i][j];
      }
      //center[count][0]=center[ic][0];
      //center[count][1]=center[ic][1];
      spos[count][0]=spos[ic][0];
      spos[count][1]=spos[ic][1];
      count++;
    }
  }Nstar=count;*Nobj=Nstar;
  printf("selected star number=%d\n",Nstar);
  strcpy(ftmp,"fits/instars.fits");
  dim[0]=dim[1]=Ng0;dim[2]=Nstar;
  //write_fits_3D(ftmp,star,dim);


  nbond=(Ng0-Ng)/2+1;
  for(ic=0;ic<Nstar;ic++){
    for(i=0;i<Ng0;i++)for(j=0;j<Ng0;j++){
      error0[ic][i*Ng0+j]=0;
      start[i][j]=star[ic][i][j];
    }
    star_gaus(start,Ng0,Ng0,starf,residu,para);
    //center[ic][0]=xcen[0];center[ic][1]=xcen[1];
    xcen[0]=para[3];xcen[1]=para[4];//printf("x=%f,y=%f\n",xcen[0],xcen[1] );
    //printf("x=%f,y=%f\n",xcen[0],xcen[1] );
    Interp_bicubic(Ng0,Ng0,start,nbond,Ng,Ng,starf,xcen);
    gaus_estimate(starf,Ng,Ng,&mean,&sigma);
    k=0;sum=0;
    for(i=0;i<Ng;i++){
      for(j=0;j<Ng;j++){
        a[ic][k]=starf[i][j];
        detx=sigma*sigma*gain*gain;
        dety=gain*fabs(a[ic][k]);
        if(dety<detx)error0[ic][k]=sigma;
        else{error0[ic][k]=sqrt(dety+detx)/gain;}//estimate poisson noise
        sum +=a[ic][k];
        k++;
      }
    }
    for(i=0;i<Ng0;i++)for(j=0;j<Ng0;j++){
      start[i][j]=wgt[ic][i][j];
    }
    Interp_bicubic(Ng0,Ng0,start,nbond,Ng,Ng,starf,xcen);
    for(k=0;k<Mp;k++){
      error0[ic][k] /= sum;a[ic][k] /= sum;
      weight[ic][k]=1./(error0[ic][k]*error0[ic][k]);
    }
    for(i=0;i<Ng;i++)for(j=0;j<Ng;j++){
      if(starf[i][j]<0.5)starf[i][j]=0;
      weight[ic][i*Ng+j]*=starf[i][j]*starf[i][j];
      star0[ic][i][j]=a[ic][i*Ng+j];
      stars[ic][i][j]=a[ic][i*Ng+j];
    }

  }
  strcpy(ftmp,"fits/star_ctr.fits");
  dim[0]=dim[1]=Ng;dim[2]=Nstar;
  //write_fits_3D(ftmp,star0,dim);

  sigmasum=0;
  for(ic=0;ic<Nstar;ic++){
    star_gaus(star0[ic],Ng,Ng,starf,residu,para); 
    center[ic][0]=para[3];center[ic][1]=para[4];sigmasum+=para[1]*1;
  }sigmasum/=Nstar;
    estimate_para(a,center,Nstar,Ng,0,para);
    if(para[1]<=0 || para[2]<=0){
      para[1]=3.0;para[2]=3.0;
    }
    rd=para[1];beta=para[2];printf("rd=%e\tbeta=%e\n",para[1],beta );

  //EMPCA(a,Nstar,Mp,compbfr,compcoeff,coefferr,weight,npc);

  for(k=0;k<13;k++)mx[k]=0;
  //mx[0]=mx[1]=1; // Here we only chose the basis function with theta, 2n*theta and 3n*theta.
  mx[0]=mx[1]=mx[2]=mx[3]=mx[4]=1;
  creat_moffatelet(rd,beta,l0,msign,mx,ebasisf,Mp,Ng,&nbf);printf("creat is done!\n");
  //creat_gaussianlet(sigmasum,l0,msign,mx,ebasisf,Mp,Ng,&nbf);printf("nbf=%d\n", nbf);
  

  SPCA(a,weight,ebasisf,Mp,Nstar,nbf,compbfr,compcoeff,npc);
  printf("SPCA down\n");
    for(k=0;k<npc;k++){
      for(i=0;i<Ng0;i++){
        for(j=0;j<Ng0;j++){
          PCs[k][i][j]=0;
        }
      }
    }
    for(k=0;k<npc;k++){
      for(i=0;i<Ng;i++){
        for(j=0;j<Ng;j++){
          PCs[k][i][j]=compbfr[k][i*Ng+j];
        }
      }
    }
    strcpy(ftmp,"fits/PCs.fits");
    dim[0]=dim[1]=Ng,dim[2]=npc;
    //write_fits_3D(ftmp,PCs,dim);

    for(ic=0;ic<Nstar;ic++){
      for(ipc=0;ipc<npc;ipc++){
        sum=0;
        for(i=0;i<Ng;i++)for(j=0;j<Ng;j++){
          sum+=compbfr[ipc][i*Ng+j]*compbfr[ipc][i*Ng+j]*weight[ic][i*Ng+j];
        }sum=sqrt(sum);sum=1./sum;coefferr[ic][ipc]=sum;
        sumt=fabs(compcoeff[ipc][ic]/sum);
        if(sumt<1.5)compcoeff[ipc][ic]=0.;
      }
    }
    fp=fopen("results/coeff.dat","w");
    wp=fopen("results/weight.dat","w");
    for(ic=0;ic<Nstar;ic++){
      for(k=0;k<npc;k++){
    	  coeff[ic][k][0]=compcoeff[k][ic];
    	  coeff[ic][k][1]=coefferr[ic][k];
        fprintf(fp, "%f\t", coeff[ic][k][0]);
        fprintf(wp, "%f\t", coeff[ic][k][1]);
      }fprintf(fp, "\n" );fprintf(wp, "\n" );
    }fclose(fp);fclose(wp);
    
    printf("recostruction\n");
    strcpy(ftmp,"results/");
    strcat(ftmp,"chi.dat");
    fp=fopen(ftmp,"w");
    for(ic=0;ic<Nstar;ic++){
      sum=0;//printf("ic=%d\n",ic );
      for(i=0;i<Ng;i++)for(j=0;j<Ng;j++){
        start[i][j]=0;
        for(ipc=0;ipc<npc;ipc++){
          start[i][j]+=compbfr[ipc][i*Ng+j]*compcoeff[ipc][ic];
        }
        star[ic][i][j]=start[i][j];
        residual[ic][i][j]=star0[ic][i][j]-start[i][j];
        sum+=residual[ic][i][j]*residual[ic][i][j]*weight[ic][i*Ng+j];
      }
      sum/=Mp;
      fprintf(fp, "%f\n", sum);
    }fclose(fp);

    strcpy(ftmp,"fits/residual.fits");
    dim[0]=dim[1]=Ng,dim[2]=Nstar;
    //write_fits_3D(ftmp,residual,dim);
    strcpy(ftmp,"fits/restar.fits");
    dim[0]=dim[1]=Ng,dim[2]=Nstar;
    //write_fits_3D(ftmp,star,dim);printf("SPCA down\n");

    double oparas[5],rparas[5],cmpara[5];
    for(ic=0;ic<Nstar;ic++){
      star_gaus(star0[ic],Ng,Ng,starf,residu,para);
      center[ic][0]=para[3];center[ic][1]=para[4];
    }
    estimate_para(a,center,Nstar,Ng,0,para);
    double rg=para[0]*3.5;printf("rg=%f\n",rg);
    //rg=1.422417;printf("Ng=%d\n",Ng );
    count=0;
    for(ic=0;ic<Nstar;ic++){
      size(star0[ic],Ng,center[ic],rg,oparas);
      star_gaus(star[ic],Ng,Ng,starf,residu,para);
      center[ic][0]=para[3];center[ic][1]=para[4];
      size(star[ic],Ng,center[ic],rg,rparas);
      cmpara[1]=fabs(oparas[1]-rparas[1]);
      cmpara[2]=fabs(oparas[2]-rparas[2]);
      cmpara[0]=fabs((oparas[0]-rparas[0])/oparas[0]);
      if(cmpara[1]<=0.1 && cmpara[2]<=0.1 && cmpara[0]<=0.1){
        //printf("%f,%f,%f\n",rparas[0],rparas[1],rparas[2] );
        for(k=0;k<npc;k++){//copy coeff
          coeff[count][k][0]=coeff[ic][k][0];
          coeff[count][k][1]=coeff[ic][k][1];
        }
        spos[count][0]=spos[ic][0];
        spos[count][1]=spos[ic][1];
        for(i=0;i<Ng;i++)for(j=0;j<Ng;j++){
          stars[count][i][j]=stars[ic][i][j];
          star[count][i][j]=star[ic][i][j];
        }
        count++;
      }
    }*Nobj=count;
    printf("selected Nstar=%d\n",count );

    strcpy(ftmp,"fits/slecstar.fits");
    dim[0]=dim[1]=Ng,dim[2]=count;
    //write_fits_3D(ftmp,star,dim);


  


  free_f3tensor(star0,0,Nstar0,0,Ng0,0,Ng0);
  free_f3tensor(star,0,Nstar0,0,Ng0,0,Ng0);
  free_f3tensor(residual,0,Nstar0,0,Ng0,0,Ng0);
  free_dmatrix(a,0,Nstar0,0,Mp0);
  free_dmatrix(error0,0,Nstar0,0,Mp0);
  free_dmatrix(weight,0,Nstar0,0,Mp0);
  free_matrix(start,0,Ng0,0,Ng0);
  free_matrix(starf,0,Ng0,0,Ng0);
  free_matrix(residu,0,Ng0,0,Ng0);
  free_dmatrix(cent,0,Nstar0,0,Ng0);
  free_dmatrix(center,0,Nstar0,0,Ng0);
  free_d3tensor(compbf,0,Nstar0,0,npc,0,Mp0);
  free_dmatrix(compcoeff,0,Nstar0,0,npc);
  free_dmatrix(compbfr,0,npc,0,Mp0);
  free_d3tensor(basis,0,npc,0,Ng0,0,Ng0);
  free_dmatrix(coefferr,0,Nstar0,0,npc);
  free_dmatrix(ebasisf,0,NPCmax,0,Mp);

}



void webSPCA_entr(float ***stars,float **spos,int npc,int Nstar0,int Ng0,int Ng,
  float ***PCs,float ***coeff,float gain,int *Nobj,float snrs,float ***models,
  float ***wgt){

  int i,j,k,ic,Mp,Mp0,r,count,pl,ra,rb,f,ipc;
  int nbond=(Ng0-Ng)/2+1,mx[13],msign=1,l0=20,nbf,NPCmax,rnbf0,Nstar;
  float ***star0,**start,**starf,**residu,***star,***residual,***ebasisfs;
  double mean,sigma,detx,dety,sum,sumt,rd,beta,para[10],xcen[2],sigmasum;
  double **error0,**weight,SNR[Nstar0],**coefferr,**a;
  double **cent,**center,**ebasisf,***compbf,**compcoeff,**compbfr,***basis;
  FILE *fp,*wp;
  Mp0=Ng0*Ng0;Mp=Ng*Ng;
  NPCmax=Nstar0;
  rnbf0=(l0+1)*(l0+1);
  if(NPCmax<rnbf0)NPCmax=rnbf0;
  star0=f3tensor(0,Nstar0,0,Ng0,0,Ng0);
  star=f3tensor(0,Nstar0,0,Ng0,0,Ng0);
  residual=f3tensor(0,Nstar0,0,Ng0,0,Ng0);
  a=dmatrix(0,Nstar0,0,Mp0);
  error0=dmatrix(0,Nstar0,0,Mp0);
  weight=dmatrix(0,Nstar0,0,Mp0);
  start=matrix(0,Ng0,0,Ng0);
  starf=matrix(0,Ng0,0,Ng0);
  residu=matrix(0,Ng0,0,Ng0);
  cent=dmatrix(0,Nstar0,0,Ng0);
  center=dmatrix(0,Nstar0,0,Ng0);
  compbf=d3tensor(0,Nstar0,0,npc,0,Mp0);
  compcoeff=dmatrix(0,npc,0,Nstar0);
  compbfr=dmatrix(0,npc,0,Mp0);
  basis=d3tensor(0,npc,0,Ng0,0,Ng0);
  coefferr=dmatrix(0,Nstar0,0,npc);
  ebasisf=dmatrix(0,NPCmax,0,Mp);
  ebasisfs=f3tensor(0,rnbf0,0,Ng0,0,Ng0);
  printf("method SPCA\n");
  

  for(ic=0;ic<Nstar0;ic++){
    for(i=0;i<Ng0;i++)for(j=0;j<Ng0;j++){
      star0[ic][i][j]=stars[ic][i][j];
    }
  }

  int dim[5];char ftmp[100];
  
  
  for(ic=0;ic<Nstar0;ic++){
    gaus_estimate(star0[ic],Ng0,Ng0,&mean,&sigma); //estimate the background noise as gasussian randoms
    printf("ic=%d,mean=%f,sigma=%f\n",ic+1,mean,sigma );
    k=0;
    for(i=0;i<Ng0;i++){
      for(j=0;j<Ng0;j++){
        a[ic][k]=star0[ic][i][j]-mean;//subtract the bacground sky intensity
        start[i][j]=a[ic][k];
        k++;
      }
    }
    star_gaus(start,Ng0,Ng0,starf,residu,para);
    xcen[0]=para[3];xcen[1]=para[4];printf("x=%f,y=%f\n",xcen[0],xcen[1] );
    //center[ic][0]=xcen[0];center[ic][1]=xcen[1];
    k=0;sum=0;
    for(i=0;i<Ng0;i++){
      for(j=0;j<Ng0;j++){
        start[i][j]=star0[ic][i][j]-mean;
        detx=sigma*sigma*gain*gain;
        dety=gain*fabs(start[i][j]);
        if(dety<detx)error0[ic][k]=sigma;
        else{error0[ic][k]=sqrt(dety+detx)/gain;}//estimate poisson noise
        sum +=start[i][j];
        k++;
      }
    }
    sigmasum=0;
    for(k=0;k<Mp0;k++){
      sigmasum+=error0[ic][k]*error0[ic][k];
    }SNR[ic]=sum/sqrt(sigmasum);//calculate the SNR of each image
    //printf("snr=%f\n",SNR[ic] );
  }

  count=0;
  for(ic=0;ic<Nstar0;ic++){
    if(SNR[ic]>=snrs){   //select the high snr images
      for(i=0;i<Ng0;i++)for(j=0;j<Ng0;j++){
        star[count][i][j]=a[ic][i*Ng0+j];
        wgt[count][i][j]=wgt[ic][i][j];
      }
      for(i=0;i<Ng0-1;i++)for(j=0;j<Ng0-1;j++){
        models[count][i][j]=models[ic][i][j];
      }
      //center[count][0]=center[ic][0];
      //center[count][1]=center[ic][1];
      spos[count][0]=spos[ic][0];
      spos[count][1]=spos[ic][1];
      count++;
    }
  }Nstar=count;*Nobj=Nstar;
  printf("selected star number=%d\n",Nstar);
  strcpy(ftmp,"fits/instars.fits");
  dim[0]=dim[1]=Ng0;dim[2]=Nstar;
  //write_fits_3D(ftmp,star,dim);


  nbond=(Ng0-Ng)/2+1;
  for(ic=0;ic<Nstar;ic++){
    for(i=0;i<Ng0;i++)for(j=0;j<Ng0;j++){
      error0[ic][i*Ng0+j]=0;
      start[i][j]=star[ic][i][j];
    }
    star_gaus(start,Ng0,Ng0,starf,residu,para);
    //center[ic][0]=xcen[0];center[ic][1]=xcen[1];
    xcen[0]=para[3];xcen[1]=para[4];//printf("x=%f,y=%f\n",xcen[0],xcen[1] );
    //printf("x=%f,y=%f\n",xcen[0],xcen[1] );
    Interp_bicubic(Ng0,Ng0,start,nbond,Ng,Ng,starf,xcen);
    gaus_estimate(starf,Ng,Ng,&mean,&sigma);
    k=0;sum=0;
    for(i=0;i<Ng;i++){
      for(j=0;j<Ng;j++){
        a[ic][k]=starf[i][j];
        detx=sigma*sigma*gain*gain;
        dety=gain*fabs(a[ic][k]);
        if(dety<detx)error0[ic][k]=sigma;
        else{error0[ic][k]=sqrt(dety+detx)/gain;}//estimate poisson noise
        sum +=a[ic][k];
        k++;
        residual[ic][i][j]=models[ic][i+nbond-1][j+nbond-1];//resize models 
      }
    }
    for(i=0;i<Ng0;i++)for(j=0;j<Ng0;j++){
      start[i][j]=wgt[ic][i][j];
    }
    Interp_bicubic(Ng0,Ng0,start,nbond,Ng,Ng,starf,xcen);
    for(k=0;k<Mp;k++){
      error0[ic][k] /= sum;a[ic][k] /= sum;
      weight[ic][k]=1./(error0[ic][k]*error0[ic][k]);
    }
    sum=0;
    for(i=0;i<Ng;i++)for(j=0;j<Ng;j++)sum+=residual[ic][i][j];
    if(sum!=0)for(i=0;i<Ng;i++)for(j=0;j<Ng;j++)residual[ic][i][j]/=sum;//normalise models
    for(i=0;i<Ng;i++)for(j=0;j<Ng;j++){
      weight[ic][i*Ng+j]*=starf[i][j]*starf[i][j];
      stars[ic][i][j]=a[ic][i*Ng+j];
      star0[ic][i][j]=a[ic][i*Ng+j];
    }

  }
  strcpy(ftmp,"fits/star_ctr.fits");
  dim[0]=dim[1]=Ng;dim[2]=Nstar;
  //write_fits_3D(ftmp,stars,dim);
  

  sigmasum=0;
  for(ic=0;ic<Nstar;ic++){
    star_gaus(star0[ic],Ng,Ng,starf,residu,para); 
    center[ic][0]=para[3];center[ic][1]=para[4];sigmasum+=para[1]*1;
  }sigmasum/=Nstar;
    estimate_para(a,center,Nstar,Ng,0,para);
    if(para[1]<=0 || para[2]<=0){
      para[1]=3.0;para[2]=3.0;
    }
    rd=para[1];beta=para[2];printf("rd=%e\tbeta=%e\n",para[1],beta );

  for(ic=0;ic<Nstar;ic++){
    for(i=0;i<Ng;i++)for(j=0;j<Ng;j++){
      a[ic][i*Ng+j]-=residual[ic][i][j];//subtract model from observed data
      star0[ic][i][j]=a[ic][i*Ng+j];
    }
  }
  strcpy(ftmp,"fits/substar.fits");
  dim[0]=dim[1]=Ng;dim[2]=Nstar;
  //write_fits_3D(ftmp,star0,dim);

  for(k=0;k<13;k++)mx[k]=0;
  //mx[0]=mx[1]=1; // Here we only chose the basis function with theta, 2n*theta and 3n*theta.
  mx[0]=mx[1]=mx[2]=mx[3]=mx[4]=1;
  creat_moffatelet(rd,beta,l0,msign,mx,ebasisf,Mp,Ng,&nbf);printf("creat is done!\n");
  //creat_gaussianlet(sigmasum,l0,msign,mx,ebasisf,Mp,Ng,&nbf);printf("nbf=%d\n", nbf);


  SPCA(a,weight,ebasisf,Mp,Nstar,nbf,compbfr,compcoeff,npc);
  printf("SPCA down\n");
    for(k=0;k<npc;k++){
      for(i=0;i<Ng0;i++){
        for(j=0;j<Ng0;j++){
          PCs[k][i][j]=0;
        }
      }
    }
    for(k=0;k<npc;k++){
      for(i=0;i<Ng;i++){
        for(j=0;j<Ng;j++){
          PCs[k][i][j]=compbfr[k][i*Ng+j];
        }
      }
    }
    strcpy(ftmp,"fits/PCs.fits");
    dim[0]=dim[1]=Ng,dim[2]=npc;
    //write_fits_3D(ftmp,PCs,dim);

    for(ic=0;ic<Nstar;ic++){
      for(ipc=0;ipc<npc;ipc++){
        sum=0;
        for(i=0;i<Ng;i++)for(j=0;j<Ng;j++){
          sum+=compbfr[ipc][i*Ng+j]*compbfr[ipc][i*Ng+j]*weight[ic][i*Ng+j];
        }sum=sqrt(sum);sum=1./sum;coefferr[ic][ipc]=sum;
        sumt=fabs(compcoeff[ipc][ic]/sum);
        if(sumt<1.5)compcoeff[ipc][ic]=0.;
      }
    }
    fp=fopen("results/coeff.dat","w");
    wp=fopen("results/weight.dat","w");
    for(ic=0;ic<Nstar;ic++){
      for(k=0;k<npc;k++){
        coeff[ic][k][0]=compcoeff[k][ic];
        coeff[ic][k][1]=coefferr[ic][k];
        fprintf(fp, "%f\t", coeff[ic][k][0]);
        fprintf(wp, "%f\t", coeff[ic][k][1]);
      }fprintf(fp, "\n" );fprintf(wp, "\n" );
    }fclose(fp);fclose(wp);
    
    printf("recostruction\n");
    strcpy(ftmp,"results/");
    strcat(ftmp,"chi.dat");
    fp=fopen(ftmp,"w");
    for(ic=0;ic<Nstar;ic++){
      sum=0;//printf("ic=%d\n",ic );
      for(i=0;i<Ng;i++)for(j=0;j<Ng;j++){
        start[i][j]=0;
        for(ipc=0;ipc<npc;ipc++){
          start[i][j]+=compbfr[ipc][i*Ng+j]*compcoeff[ipc][ic];
        }
        star[ic][i][j]=start[i][j];
        residual[ic][i][j]=star0[ic][i][j]-start[i][j];
        sum+=residual[ic][i][j]*residual[ic][i][j]*weight[ic][i*Ng+j];
      }
      sum/=Mp;
      fprintf(fp, "%f\n", sum);
    }fclose(fp);

    strcpy(ftmp,"fits/residual.fits");
    dim[0]=dim[1]=Ng,dim[2]=Nstar;
    //write_fits_3D(ftmp,residual,dim);
    strcpy(ftmp,"fits/restar.fits");
    dim[0]=dim[1]=Ng,dim[2]=Nstar;
    //write_fits_3D(ftmp,star,dim);printf("SPCA down\n");

    /*double oparas[5],rparas[5],cmpara[5];
    for(ic=0;ic<Nstar;ic++){
      star_gaus(star0[ic],Ng,Ng,starf,residu,para);
      center[ic][0]=para[3];center[ic][1]=para[4];
    }
    estimate_para(a,center,Nstar,Ng,0,para);
    double rg=para[0]*3.5;printf("rg=%f\n",rg);
    //rg=1.422417;printf("Ng=%d\n",Ng );
    count=0;
    for(ic=0;ic<Nstar;ic++){
      size(star0[ic],Ng,center[ic],rg,oparas);
      star_gaus(star[ic],Ng,Ng,starf,residu,para);
      center[ic][0]=para[3];center[ic][1]=para[4];
      size(star[ic],Ng,center[ic],rg,rparas);
      cmpara[1]=fabs(oparas[1]-rparas[1]);
      cmpara[2]=fabs(oparas[2]-rparas[2]);
      cmpara[0]=fabs((oparas[0]-rparas[0])/oparas[0]);
      if(cmpara[1]<=0.1 && cmpara[2]<=0.1 && cmpara[0]<=0.1){
        //printf("%f,%f,%f\n",rparas[0],rparas[1],rparas[2] );
        for(k=0;k<npc;k++){//copy coeff
          coeff[count][k][0]=coeff[ic][k][0];
          coeff[count][k][1]=coeff[ic][k][1];
        }
        spos[count][0]=spos[ic][0];
        spos[count][1]=spos[ic][1];
        for(i=0;i<Ng;i++)for(j=0;j<Ng;j++){
          stars[count][i][j]=stars[ic][i][j];
          star[count][i][j]=star[ic][i][j];
        }
        count++;
      }
    }*Nobj=count;
    printf("selected Nstar=%d\n",count );*/
    *Nobj=Nstar;

    strcpy(ftmp,"fits/slecstar.fits");
    dim[0]=dim[1]=Ng,dim[2]=count;
    //write_fits_3D(ftmp,star,dim);

}
