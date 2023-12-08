#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"
#include <string.h>
#include "function.h"
void iSPCA_entr(float ***stars,float **spos,int npc,int Nstar0,int Ng0,int Ng,
  float ***PCs,float ***coeff,float gain,int *Nobj,float snrs,int osam,float ***wgt){

	int i,j,k,ic,Mp,scale=2,num=osam*scale,Nk,Nke,Mpke,Nge,r,count,pl,ra,rb,f,ipc;
	int mx[13],msign=1,l0=9,nbf,NPCmax,rnbf0,Nstar,nbond=2;
	float ***star0,**start,**starf,**residu,***star,***residual,***ebasisfs;
	double mean,sigma,detx,dety,sum,sumt,rd,beta,para[10],xcen[2],sigmasum;
	double **error0,**weight,SNR[Nstar0],**coefferr,**a;
	double **cent,**center,**ebasisf,***compbf,**compcoeff,**compbfr,***recept,***basis;
  float **mask;
  FILE *fp;
	Mp=Ng*Ng;Nk=Ng+4;Nke=Nk*num;Mpke=Nke*Nke;Nge=Ng*num;
	NPCmax=Nstar0;
	rnbf0=(l0+1)*(l0+1);
  if(NPCmax<rnbf0)NPCmax=rnbf0;
	star0=f3tensor(0,Nstar0,0,Ng0,0,Ng0);
	star=f3tensor(0,Nstar0,0,Ng0,0,Ng0);
  residual=f3tensor(0,Nstar0,0,Ng0,0,Ng0);
	a=dmatrix(0,Nstar0,0,Mp);
	error0=dmatrix(0,Nstar0,0,Mp);
	weight=dmatrix(0,Nstar0,0,Mp);
	start=matrix(0,Ng0,0,Ng0);
	starf=matrix(0,Ng0,0,Ng0);
	residu=matrix(0,Ng0,0,Ng0);
	cent=dmatrix(0,Nstar0,0,Ng0);
	center=dmatrix(0,Nstar0,0,Ng0);
	ebasisf=dmatrix(0,NPCmax,0,Mpke);
	compbf=d3tensor(0,Nstar0,0,npc,0,Mp);
	compcoeff=dmatrix(0,npc,0,Nstar0);
	compbfr=dmatrix(0,npc,0,Mpke);
	recept=d3tensor(0,npc,0,Nke,0,Nke);
	basis=d3tensor(0,npc,0,Nge,0,Nge);
  coefferr=dmatrix(0,Nstar0,0,npc);
  ebasisfs=f3tensor(0,rnbf0,0,Nke,0,Nke);
  mask=matrix(0,Nstar0,0,Mp);
  printf("method iSPCA\n");
	

	for(ic=0;ic<Nstar0;ic++){
	  for(i=0;i<Ng0;i++)for(j=0;j<Ng0;j++){
		  star0[ic][i][j]=stars[ic][i][j];
	  }
	}

  int dim[5];char ftmp[100];
  strcpy(ftmp,"fits/instars.fits");
  dim[0]=dim[1]=Ng0;dim[2]=Nstar0;
  //write_fits_3D(ftmp,star0,dim);

	for(ic=0;ic<Nstar0;ic++){
	  gaus_estimate(star0[ic],Ng0,Ng0,&mean,&sigma); //estimate the background noise as gasussian randoms
    //printf("ic=%d,mean=%f,sigma=%f\n",ic+1,mean,sigma );
	  for(i=0;i<Ng;i++){
	  	for(j=0;j<Ng;j++){
	  	  start[i][j]=star0[ic][i+nbond][j+nbond]-mean;//subtract the bacground sky intensity
          mask[ic][i*Ng+j]=wgt[ic][i+nbond][j+nbond];
	  	}
	  }
	  star_gaus(start,Ng,Ng,starf,residu,para);
	  xcen[0]=para[3];xcen[1]=para[4];//printf("x=%f,y=%f\n",xcen[0],xcen[1] );
	  cent[ic][0]=xcen[0]*num+num/2;
    cent[ic][1]=xcen[1]*num+num/2;
    k=0;sum=0;
    for(i=0;i<Ng;i++){
      for(j=0;j<Ng;j++){
      	a[ic][k]=start[i][j];
      	detx=sigma*sigma*gain*gain;
      	dety=gain*fabs(a[ic][k]);
      	if(dety<detx)error0[ic][k]=sigma;
      	else{error0[ic][k]=sqrt(dety+detx)/gain;}//estimate poisson noise
      	sum +=a[ic][k];
      	k++;
      }
    }
    sigmasum=0;
    for(k=0;k<Mp;k++){
      sigmasum+=error0[ic][k]*error0[ic][k];
    }SNR[ic]=sum/sqrt(sigmasum);//calculate the SNR of each image
    //printf("snr=%f\n",SNR[ic] );
    for(k=0;k<Mp;k++){
    	error0[ic][k] /= sum;a[ic][k] /= sum;
    	weight[ic][k]=1./(error0[ic][k]*error0[ic][k])*mask[ic][k];
    }
    for(i=0;i<Ng;i++)for(j=0;j<Ng;j++){
    	star[ic][i][j]=a[ic][i*Ng+j];
    }
	}
  
  count=0;
  for(ic=0;ic<Nstar0;ic++){
    if(SNR[ic]>=snrs){   //select the high snr images
      for(k=0;k<Mp;k++){
        error0[count][k]=error0[ic][k];
        a[count][k]=a[ic][k];
        weight[count][k]=weight[ic][k];
      }

      for(i=0;i<Ng;i++)for(j=0;j<Ng;j++){
        star[count][i][j]=a[count][i*Ng+j];
        stars[count][i][j]=a[count][i*Ng+j];
      }
      cent[count][0]=cent[ic][0];
      cent[count][1]=cent[ic][1];
      spos[count][0]=spos[ic][0];
      spos[count][1]=spos[ic][1];
      count++;
    }
  }Nstar=count;*Nobj=Nstar;
  printf("selected star number=%d\n",Nstar);

  strcpy(ftmp,"fits/star_ctr.fits");
  dim[0]=dim[1]=Ng;dim[2]=Nstar;
  //write_fits_3D(ftmp,star,dim);


  sigmasum=0;
	for(ic=0;ic<Nstar;ic++){
    star_gaus(star[ic],Ng,Ng,starf,residu,para); 
    center[ic][0]=para[3];center[ic][1]=para[4];sigmasum+=para[1]*2;
  }sigmasum/=Nstar;sigmasum*=num;
  printf("sigma=%f\n",sigmasum );printf("a mask here\n");
  estimate_para(a,center,Nstar,Ng,0,para);
  printf("para1=%e\t para2=%e\n",para[1],para[2] );
    if(para[1]<=0 || para[2]<=0){
      printf("para1=%e\tpara2=%e\n",para[1],para[2] );
      para[1]=3.0;para[2]=3.0;
    }
    rd=para[1]*num*2.5;beta=para[2];printf("rd=%e\tbeta=%e\n",para[1],beta );


	for(k=0;k<13;k++)mx[k]=0;
  //mx[0]=mx[1]=1; // Here we only chose the basis function with theta, 2n*theta and 3n*theta.
  mx[0]=mx[1]=mx[2]=mx[3]=mx[4]=1;
  creat_moffatelets(rd,beta,l0,msign,mx,ebasisf,Mpke,Nke,&nbf);printf("creat is done!\n");
  //creat_gaussianlets(sigmasum,l0,msign,mx,ebasisf,Mpke,Nke,&nbf);printf("nbf=%d\n",nbf );


  iSPCA(a,weight,ebasisf,Mp,Nstar,nbf,cent,compbf,compcoeff,npc,Nk,num,compbfr,Ng);printf("iSPCA down\n");

    for(ic=0;ic<npc;ic++){
      for(i=0;i<Nke;i++){
        for(j=0;j<Nke;j++){
        	k=i*Nke+j;
        	recept[ic][i][j]=compbfr[ic][k];
        }
      }
    }
    r=(Nke-Nge)/2.;
    
    for(k=0;k<npc;k++){
      sum=0;
      for(i=0;i<Ng*osam;i++){
        for(j=0;j<Ng*osam;j++){
          basis[k][i][j]=0;
          for(pl=0;pl<scale;pl++){
            ra=i*scale+pl+r;
            for(f=0;f<scale;f++){
              rb=j*scale+f+r;
              sum+=basis[k][i][j]+=recept[k][ra][rb];
            } 
          } 
        }
      }
      for(i=0;i<Ng*osam;i++){
        for(j=0;j<Ng*osam;j++){
          PCs[k][i][j]=basis[k][i][j];
        }
      }
    }

    for(ic=0;ic<Nstar;ic++){
      for(ipc=0;ipc<npc;ipc++){
        sum=0;
        for(i=0;i<Ng;i++)for(j=0;j<Ng;j++){
          sum+=compbf[ic][ipc][i*Ng+j]*compbf[ic][ipc][i*Ng+j]*weight[ic][i*Ng+j];
        }sum=sqrt(sum);sum=1./sum;coefferr[ic][ipc]=sum;
        sumt=fabs(compcoeff[ipc][ic]/sum);
        if(sumt<1.5)compcoeff[ipc][ic]=0.;
      }
    }

    for(ic=0;ic<Nstar;ic++)for(k=0;k<npc;k++){
    	coeff[ic][k][0]=compcoeff[k][ic];
    	coeff[ic][k][1]=coefferr[ic][k];
    }
    
    printf("recostruction\n");
    strcpy(ftmp,"results/chi.dat");
    fp=fopen(ftmp,"w");
    for(ic=0;ic<Nstar;ic++){
      sum=0;//printf("ic=%d\n",ic );
      for(i=0;i<Ng;i++)for(j=0;j<Ng;j++){
        start[i][j]=0;
        for(ipc=0;ipc<npc;ipc++){
          start[i][j]+=compbf[ic][ipc][i*Ng+j]*compcoeff[ipc][ic];
        }
        star0[ic][i][j]=start[i][j];
        residual[ic][i][j]=star[ic][i][j]-start[i][j];
        sum+=residual[ic][i][j]*residual[ic][i][j]*weight[ic][i*Ng+j];
      }
      sum/=Ng*Ng;
      fprintf(fp, "%f\n", sum);
    }fclose(fp);

    strcpy(ftmp,"fits/residual.fits");
    dim[0]=dim[1]=Ng,dim[2]=Nstar;
    //write_fits_3D(ftmp,residual,dim);
    strcpy(ftmp,"fits/restar.fits");
    dim[0]=dim[1]=Ng,dim[2]=Nstar;
    //write_fits_3D(ftmp,star0,dim);

    double oparas[5],rparas[5],cmpara[5];
    for(ic=0;ic<Nstar;ic++){
      star_gaus(stars[ic],Ng,Ng,starf,residu,para);
      center[ic][0]=para[3];center[ic][1]=para[4];
    }
    estimate_para(a,center,Nstar,Ng,0,para);
    double rg=para[0]*3.5;printf("rg=%f\n",rg);
    count=0;
    for(ic=0;ic<Nstar;ic++){
      size(stars[ic],Ng,center[ic],rg,oparas);
      star_gaus(star0[ic],Ng,Ng,starf,residu,para);
      center[ic][0]=para[3];center[ic][1]=para[4];
      size(star0[ic],Ng,center[ic],rg,rparas);
      cmpara[1]=fabs(oparas[1]-rparas[1]);
      cmpara[2]=fabs(oparas[2]-rparas[2]);
      cmpara[0]=fabs((oparas[0]-rparas[0])/oparas[0]);
      if(cmpara[1]<=0.1 && cmpara[2]<=0.1 && cmpara[0]<=0.1){
        for(k=0;k<npc;k++){//copy coeff
          coeff[count][k][0]=coeff[ic][k][0];
          coeff[count][k][1]=coeff[ic][k][1];
        }
        spos[count][0]=spos[ic][0];
        spos[count][1]=spos[ic][1];
        for(i=0;i<Ng;i++)for(j=0;j<Ng;j++){
          stars[count][i][j]=stars[ic][i][j];
          star0[count][i][j]=star0[ic][i][j];
        }
        count++;
      }
    }*Nobj=count;
    printf("selected Nstar=%d\n",count );

    strcpy(ftmp,"fits/slecstar.fits");
    dim[0]=dim[1]=Ng,dim[2]=count;
    //write_fits_3D(ftmp,star0,dim);

  free_f3tensor(star0,0,Nstar0,0,Ng0,0,Ng0);
  free_f3tensor(star,0,Nstar0,0,Ng0,0,Ng0);
  free_f3tensor(residual,0,Nstar0,0,Ng0,0,Ng0);
  free_dmatrix(a,0,Nstar0,0,Mp);
  free_dmatrix(error0,0,Nstar0,0,Mp);
  free_dmatrix(weight,0,Nstar0,0,Mp);
  free_matrix(start,0,Ng0,0,Ng0);
  free_matrix(starf,0,Ng0,0,Ng0);
  free_matrix(residu,0,Ng0,0,Ng0);
  free_dmatrix(cent,0,Nstar0,0,Ng0);
  free_dmatrix(center,0,Nstar0,0,Ng0);
  free_dmatrix(ebasisf,0,NPCmax,0,Mpke);
  free_d3tensor(compbf,0,Nstar0,0,npc,0,Mp);
  free_dmatrix(compcoeff,0,npc,0,Nstar0);
  free_dmatrix(compbfr,0,npc,0,Mpke);
  free_d3tensor(recept,0,npc,0,Nke,0,Nke);
  free_d3tensor(basis,0,npc,0,Nge,0,Nge);
  free_dmatrix(coefferr,0,Nstar0,0,npc);
  free_f3tensor(ebasisfs,0,rnbf0,0,Nke,0,Nke);
  free_matrix(mask,0,Nstar0,0,Mp);

}



void webiSPCA_entr(float ***stars,float **spos,int npc,int Nstar0,int Ng0,int Ng,
  float ***PCs,float ***coeff,float gain,int *Nobj,float snrs,int osam,
  float ***models,float ***wgt, float **wcoff,int big){
  void galaxy_centra(float **y,int nx,int ny,double *pa);

  //printf("(Ng0=%d,Ng=%d) and no sigma output\n",Ng0,Ng );
  void align_model(float **models,double *center,int Ng0,int Ng,int osam,float **aligned);
  int i,j,k,ic,Mp,scale=2,num,Nk,Nke,Mpke,Nge,r,count,pl,ra,rb,f,ipc;
  int mx[13],msign=1,l0=9,nbf,NPCmax,rnbf0,Nstar,nbond=(Ng0-Ng)/2,cutNg=10,cutbond=Ng/2-cutNg/2;
  float ***star0,**start,**starf,**residu,***star,***residual,***ebasisfs,**cutimg;
  double mean,sigma,detx,dety,sum,sumt,rd,beta,para[10],xcen[2],sigmasum;
  double **error0,**weight,SNR[Nstar0],**coefferr,**a,***recept1;;
  double **cent,**center,**ebasisf,***compbf,**compcoeff,**compbfr,***recept,***basis;
  FILE *fp;
  if(osam<=2){
    scale=2;
    
  }
  else{
    scale=1;
  }
  num=osam*scale;
  Mp=Ng*Ng;Nk=Ng+4;Nke=Nk*num;Mpke=Nke*Nke;Nge=Ng*num;
  NPCmax=Nstar0;
  rnbf0=(l0+1)*(l0+1);
  if(NPCmax<rnbf0)NPCmax=rnbf0;
  star0=f3tensor(0,Nstar0,0,Ng0,0,Ng0);
  star=f3tensor(0,Nstar0,0,Ng0,0,Ng0);
  residual=f3tensor(0,Nstar0,0,Ng0,0,Ng0);
  a=dmatrix(0,Nstar0,0,Mp);
  error0=dmatrix(0,Nstar0,0,Mp);
  weight=dmatrix(0,Nstar0,0,Mp);
  start=matrix(0,Ng0,0,Ng0);
  starf=matrix(0,Ng0,0,Ng0);
  residu=matrix(0,Ng0,0,Ng0);
  cent=dmatrix(0,Nstar0,0,Ng0);
  center=dmatrix(0,Nstar0,0,Ng0);
  ebasisf=dmatrix(0,NPCmax,0,Mpke);
  compbf=d3tensor(0,Nstar0,0,npc,0,Mp);
  compcoeff=dmatrix(0,npc,0,Nstar0);
  compbfr=dmatrix(0,npc,0,Mpke);
  recept=d3tensor(0,NPCmax,0,Nke,0,Nke);
  recept1=d3tensor(0,NPCmax,0,Ng0,0,Ng0);
  basis=d3tensor(0,npc,0,Nge,0,Nge);
  coefferr=dmatrix(0,Nstar0,0,npc);
  ebasisfs=f3tensor(0,rnbf0,0,Nke,0,Nke);
  cutimg=matrix(0,cutNg,0,cutNg);
  printf("method iSPCA\n");
  

  for(ic=0;ic<Nstar0;ic++){
    for(i=0;i<Ng0;i++)for(j=0;j<Ng0;j++){
      star0[ic][i][j]=stars[ic][i][j]*wgt[ic][i][j];
      //stars[ic][i][j]=star0[ic][i][j];
    }
  }

  int dim[5];char ftmp[100];
  strcpy(ftmp,"fits/instars.fits");
  dim[0]=dim[1]=(Ng0+5)*osam;dim[2]=Nstar0;
  //write_fits_3D(ftmp,models,dim);
  printf("gain=%f\n",gain );

  for(ic=0;ic<Nstar0;ic++){
    gaus_estimate(star0[ic],Ng0,Ng0,&mean,&sigma); //estimate the background noise as gasussian randoms
    //printf("ic=%d,mean=%f,sigma=%f\n",ic+1,mean,sigma );
    for(i=0;i<Ng;i++){
      for(j=0;j<Ng;j++){
        start[i][j]=star0[ic][i+nbond][j+nbond]-mean;//subtract the bacground sky intensity
        residual[ic][i][j]=wgt[ic][i+nbond][j+nbond];
      }
    }
    star_gaus(start,Ng,Ng,starf,residu,para);
    xcen[0]=para[3];xcen[1]=para[4];
    //printf("x=%f,y=%f\n",xcen[0],xcen[1] );
    cent[ic][0]=xcen[0]*num+num/2;
    cent[ic][1]=xcen[1]*num+num/2;
    k=0;sum=0;
    count=0;
    for(i=0;i<Ng;i++){
      for(j=0;j<Ng;j++){
        a[ic][k]=start[i][j];
        detx=sigma*sigma*gain*gain;
        dety=gain*fabs(a[ic][k]);
        if(dety<detx){
        //if(dety<2.*sigma*gain){
          error0[ic][k]=sigma;//count++;
        }
        else{
          error0[ic][k]=sqrt(dety+detx)/gain;
        }//estimate poisson noise
        sum +=a[ic][k]*residual[ic][i][j];
        k++;
      }
    }//printf("sum=%f\n",sum );
    sigmasum=0;
    for(k=0;k<Mp;k++){
      sigmasum+=error0[ic][k]*error0[ic][k];
    }SNR[ic]=sum/sqrt(sigmasum);//calculate the SNR of each image
    //printf("snr=%f\n",SNR[ic] );
    for(k=0;k<Mp;k++){
      error0[ic][k] /= sum;a[ic][k] /= sum;
      weight[ic][k]=1./(error0[ic][k]*error0[ic][k]);
    }
    for(i=0;i<Ng;i++)for(j=0;j<Ng;j++){
      weight[ic][i*Ng+j]*=residual[ic][i][j]*residual[ic][i][j];
      //star[ic][i][j]=a[ic][i*Ng+j];
    }
    
    if(sum!=0){
      //printf("align_model\n");
      align_model(models[ic],xcen,Ng0+5,Ng,osam,star0[ic]);
      //printf("align_model down\n");
      //shift the model to the observed positions
      sum=0;
      for(i=0;i<Ng;i++)for(j=0;j<Ng;j++)sum+=star0[ic][i][j];//*star0[ic][i][j];//residual[ic][i][j];
      if(sum!=0){
        for(i=0;i<Ng;i++)for(j=0;j<Ng;j++){
          star0[ic][i][j]/=sum;//normalise models
          //stars[ic][i][j]=error0[ic][i*Ng+j];
          //star0[ic][i][j]*=residual[ic][i][j];
        }
      }
    }
  }
  strcpy(ftmp,"fits/star_ctr.fits");
  dim[0]=dim[1]=Ng;dim[2]=Nstar0;
  //write_fits_3D(ftmp,star0,dim);
  strcpy(ftmp,"fits/weight.fits");
  dim[0]=dim[1]=Ng;dim[2]=Nstar0;
  //write_fits_3D(ftmp,stars,dim);

  
  
  count=0;
  for(ic=0;ic<Nstar0;ic++){
    if(SNR[ic]>=snrs){   //select the high snr images
      for(k=0;k<Mp;k++){
        a[count][k]=a[ic][k];
        weight[count][k]=weight[ic][k];
      }
      //printf("ic=%d\n",ic );
      for(i=0;i<Ng;i++)for(j=0;j<Ng;j++){
        a[count][i*Ng+j]*=residual[ic][i][j];
        star[count][i][j]=a[count][i*Ng+j];
        recept[count][i][j]=a[count][i*Ng+j];
        recept1[count][i][j]=weight[count][i*Ng+j];
        star0[count][i][j]=star0[ic][i][j];
      }
      for(i=0;i<Ng0;i++)for(j=0;j<Ng0;j++){
        stars[count][i][j]=stars[ic][i][j];
        wgt[count][i][j]=wgt[ic][i][j];
      }
      for(i=0;i<(Ng0+5)*osam;i++)for(j=0;j<(Ng0+5)*osam;j++){
        models[count][i][j]=models[ic][i][j];
      }
      cent[count][0]=cent[ic][0];
      cent[count][1]=cent[ic][1];
      spos[count][0]=spos[ic][0];
      spos[count][1]=spos[ic][1];
      count++;
    }
  }Nstar=count;*Nobj=Nstar;
  printf("selected star number=%d\n",Nstar);
  /*if(Nstar<7){
    printf("star number<7, programm will quit!\n");
    exit(EXIT_FAILURE);
  }*/

  sigmasum=0;
  for(ic=0;ic<Nstar;ic++){
    star_gaus(star[ic],Ng,Ng,starf,residu,para); 
    center[ic][0]=para[3];center[ic][1]=para[4];sigmasum+=para[1]*2.;
  }sigmasum/=Nstar;sigmasum*=num;
  printf("sigma=%f\n",sigmasum );
  estimate_para(a,center,Nstar,Ng,0,para);
    if(para[1]<=0 || para[2]<=0){
      para[1]=3.0;para[2]=3.0;
    }
    rd=para[1]*num*2.;beta=para[2];printf("rd=%e\tbeta=%e\n",para[1],beta );


  /*********************************************/
  //only cut the center region to calculate???
  /*********************************************/ 
  double **new_a,**new_weight,**new_ebasisf,**new_cent,**new_a1;
  double ***new_compbf,**new_compbfr,**wbpsf, *wbcoeff, *wbcoefferr;
  int new_Mp,new_Nk,new_Ng,new_Mpke,new_Nke,dpix=Ng/4,new_Nge;
  //new_Ng=Ng/2;//if((new_Ng%2)==0)new_Ng+=1;
  new_Ng=Ng*0.6;if((new_Ng%2)!=0)new_Ng-=1;printf("new_Ng=%d\n",new_Ng );
  dpix=(Ng-new_Ng)/2;
  new_Nk=new_Ng+4;new_Mp=new_Ng*new_Ng;
  new_Nke=new_Nk*num;new_Mpke=new_Nke*new_Nke;
  new_Nge=new_Ng*num;

  new_a=dmatrix(0,Nstar,0,new_Mp);
  new_a1=dmatrix(0,Nstar,0,new_Mp);
  wbpsf=dmatrix(0,Nstar,0,new_Mp);
  wbcoeff=dvector(0,Nstar);
  wbcoefferr=dvector(0,Nstar);
  new_weight=dmatrix(0,Nstar,0,new_Mp);
  new_ebasisf=dmatrix(0,NPCmax,0,new_Mpke);
  new_cent=dmatrix(0,Nstar,0,2);
  new_compbf=d3tensor(0,Nstar,0,npc,0,new_Mp);
  new_compbfr=dmatrix(0,npc,0,new_Mpke);
  for(ic=0;ic<Nstar;ic++){
    for(i=0;i<new_Ng;i++){
      for(j=0;j<new_Ng;j++){
        ra=i+dpix;rb=j+dpix;
        new_a[ic][i*new_Ng+j]=recept[ic][ra][rb];
        star[ic][i][j]=new_a[ic][i*new_Ng+j];
        new_weight[ic][i*new_Ng+j]=recept1[ic][ra][rb];
        //star[ic][i][j]=new_weight[ic][i*new_Ng+j];
      }
    }
  }printf("down\n");
  strcpy(ftmp,"fits/new_star_ctr.fits");
  dim[0]=dim[1]=new_Ng;dim[2]=Nstar;
  //write_fits_3D(ftmp,star,dim);
  
  sigmasum=0;
  for(ic=0;ic<Nstar;ic++){
    star_gaus(star[ic],new_Ng,new_Ng,starf,residu,para); 
    center[ic][0]=para[3];center[ic][1]=para[4];sigmasum+=para[1]*2.;
    //printf("x=%f,y=%f\n",para[3],para[4] );
    new_cent[ic][0]=para[3]*num+num/2;
    new_cent[ic][1]=para[4]*num+num/2;
  }sigmasum/=Nstar;sigmasum*=num;
 
  for(ic=0;ic<Nstar;ic++){
    sum=0;sumt=0;
    for(i=0;i<new_Ng;i++)for(j=0;j<new_Ng;j++){
      //new_a[ic][i*new_Ng+j]-=star0[ic][i+dpix][j+dpix];
      star[ic][i][j]=new_a[ic][i*new_Ng+j];
      wbpsf[ic][i*new_Ng+j]=star0[ic][i+dpix][j+dpix];
      star0[ic][i][j]=star[ic][i][j]-wbpsf[ic][i*new_Ng+j];
      //star[ic][i][j]=star0[ic][i+dpix][j+dpix];
      //sum+=star[ic][i][j];//sumt+=wbpsf[ic][i*new_Ng+j]*wbpsf[ic][i*new_Ng+j];
    }
    //for(i=0;i<new_Ng;i++)for(j=0;j<new_Ng;j++){
      //star[ic][i][j]/=sum;
      //wbpsf[ic][i*new_Ng+j]/=sqrt(sumt);
    //}
  }
  strcpy(ftmp,"fits/substar.fits");
  dim[0]=dim[1]=new_Ng;dim[2]=Nstar;
  //write_fits_3D(ftmp,star0,dim);

  for(ic=0;ic<Nstar;ic++){
    /*sum=0;sumt=0;
    for(k=0;k<new_Ng*new_Ng;k++){
      sum+=wbpsf[ic][k]*new_a[ic][k]*new_weight[ic][k];
      sumt+=wbpsf[ic][k]*wbpsf[ic][k]*new_weight[ic][k];
    }
    wbcoeff[ic]=sum/sumt;//printf("new coeff:%f\t",wbcoeff[ic]);
    */
    wbcoeff[ic]=1.;
    for(i=0;i<new_Ng;i++){
      for(j=0;j<new_Ng;j++){
        //star0[ic][i][j]=star[ic][i][j]-wbcoeff[ic]*wbpsf[ic][i*new_Ng+j];
        new_a1[ic][i*new_Ng+j]=star[ic][i][j]-wbcoeff[ic]*wbpsf[ic][i*new_Ng+j];
      }
    }
  }

  strcpy(ftmp,"fits/substar1.fits");
  dim[0]=dim[1]=new_Ng;dim[2]=Nstar;
  //write_fits_3D(ftmp,star0,dim);
  //exit(EXIT_FAILURE);
  for(k=0;k<13;k++)mx[k]=0;
  //mx[0]=mx[1]=1; // Here we only chose the basis function with theta, 2n*theta and 3n*theta.
  mx[0]=mx[1]=mx[2]=mx[3]=mx[4]=1;
  //creat_gaussianlets(sigmasum,l0,msign,mx,new_ebasisf,new_Mpke,new_Nke,&nbf);printf("nbf=%d\n",nbf );
  creat_moffatelets(rd,beta,l0,msign,mx,new_ebasisf,new_Mpke,new_Nke,&nbf);

  float ***basisf;
  basisf=f3tensor(0,nbf,0,new_Nke,0,new_Nke);
  for(ic=0;ic<nbf;ic++){
    for(i=0;i<new_Nke;i++){
      for(j=0;j<new_Nke;j++){
        basisf[ic][i][j]=new_ebasisf[ic][i*new_Nke+j];
      }
    }
  }
  strcpy(ftmp,"fits/basis.fits");
  dim[0]=dim[1]=new_Nke;dim[2]=nbf;
  //write_fits_3D(ftmp,basisf,dim);

  void hiSPCA(double **images,double **weight,double **basisf,int Mp,int Nstar,int nbf,double **cent,
      double ***compbf,double **coeff,int Nb,int Nk,int num,double **compbfr,int Ng,
      double **wbpsf,double *wbcoeff);
  if(big==0){
    printf("fitting for normal star images\n");
    hiSPCA(new_a,new_weight,new_ebasisf,new_Mp,Nstar,nbf,new_cent,
    new_compbf,compcoeff,npc,new_Nk,num,new_compbfr,new_Ng,wbpsf,wbcoeff);
  }
  else{
    printf("fitting for very bright star images\n");
    iSPCA(new_a1,new_weight,new_ebasisf,new_Mp,Nstar,nbf,new_cent,
    new_compbf,compcoeff,npc,new_Nk,num,new_compbfr,new_Ng);
  }
  
  

  //printf("out SPCA\n");
  for(ic=0;ic<npc;ic++){
      //printf("npc=%d\n",npc);
      for(i=0;i<new_Nke;i++){
        for(j=0;j<new_Nke;j++){
          k=i*new_Nke+j;
          recept[ic][i][j]=new_compbfr[ic][k];//transform 2d PCs to cubic 
        }
      }
    }printf("get basis\n");
    r=(new_Nke-new_Nge)/2.;printf("r=%d\n",r );
    
    for(k=0;k<npc;k++){
      for(i=0;i<new_Ng*osam;i++){
        for(j=0;j<new_Ng*osam;j++){
          basis[k][i][j]=0;//printf("i=%d,j=%d\n",i,j );
          for(pl=0;pl<scale;pl++){
            ra=i*scale+pl+r;
            for(f=0;f<scale;f++){
              rb=j*scale+f+r;
              basis[k][i][j]+=recept[k][ra][rb];
              //mergering oversampled PCs to required resolutions
            } 
          }
          ebasisfs[k][i][j]=basis[k][i][j]; 
        }
      }//printf("merger basis\n");
      for(i=0;i<Ng*osam;i++)for(j=0;j<Ng*osam;j++)PCs[k][i][j]=0;
      for(i=0;i<new_Ng*osam;i++){
        for(j=0;j<new_Ng*osam;j++){
          PCs[k][i+dpix*osam][j+dpix*osam]=basis[k][i][j];//copy in PCs
          
        }
      }
    }printf("get down\n");
    strcpy(ftmp,"fits/oPCs.fits");
    dim[0]=dim[1]=new_Ng*osam,dim[2]=npc;
    //write_fits_3D(ftmp,ebasisfs,dim);

    for(ic=0;ic<Nstar;ic++){
      for(ipc=0;ipc<npc;ipc++){
        sum=0;
        for(i=0;i<new_Ng;i++)for(j=0;j<new_Ng;j++){
          sum+=new_compbf[ic][ipc][i*new_Ng+j]*new_compbf[ic][ipc][i*new_Ng+j]*new_weight[ic][i*new_Ng+j];
        }sum=sqrt(sum);sum=1./sum;coefferr[ic][ipc]=sum;
        sumt=fabs(compcoeff[ipc][ic]/sum);
        //if(sumt<1.5)compcoeff[ipc][ic]=0.;
        //printf("%f\t",compcoeff[ipc][ic] );
      }
      sum=0;
      for(i=0;i<new_Ng;i++)for(j=0;j<new_Ng;j++){
          sum+=wbpsf[ic][i*new_Ng+j]*wbpsf[ic][i*new_Ng+j]*new_weight[ic][i*new_Ng+j];
        }sum=sqrt(sum);sum=1./sum;wbcoefferr[ic]=sum;
        sumt=fabs(wbcoeff[ic]/sum);
        //if(sumt<1.2)wbcoeff[ic]=0.;
        //printf("%f\t",wbcoeff[ic] );
        //printf("\n");
    }

    for(ic=0;ic<Nstar;ic++){
      wcoff[ic][0]=wbcoeff[ic];wcoff[ic][1]=wbcoefferr[ic];
      for(k=0;k<npc;k++){
        coeff[ic][k][0]=compcoeff[k][ic];
        coeff[ic][k][1]=coefferr[ic][k];
      }
    }

    printf("recostruction\n");
    strcpy(ftmp,"results/chi.dat");
    //fp=fopen(ftmp,"w");
    for(ic=0;ic<Nstar;ic++){
      sum=0;//printf("ic=%d\n",ic );
      for(i=0;i<new_Ng;i++)for(j=0;j<new_Ng;j++){
        start[i][j]=0;
        start[i][j]=wbpsf[ic][i*new_Ng+j]*wbcoeff[ic];
        for(ipc=0;ipc<npc;ipc++){
          start[i][j]+=new_compbf[ic][ipc][i*new_Ng+j]*compcoeff[ipc][ic];
        }
        star0[ic][i][j]=start[i][j];
        residual[ic][i][j]=star[ic][i][j]-start[i][j];
        sum+=residual[ic][i][j]*residual[ic][i][j]*new_weight[ic][i*new_Ng+j];
      }
      sum/=new_Ng*new_Ng;
      //fprintf(fp, "%f\n", sum);
    }//fclose(fp);

    strcpy(ftmp,"fits/residual.fits");
    dim[0]=dim[1]=new_Ng,dim[2]=Nstar;
    //write_fits_3D(ftmp,residual,dim);
    strcpy(ftmp,"fits/restar.fits");
    dim[0]=dim[1]=new_Ng,dim[2]=Nstar;
    //write_fits_3D(ftmp,star0,dim);
    strcpy(ftmp,"fits/ostar.fits");
    dim[0]=dim[1]=new_Ng,dim[2]=Nstar;
    //write_fits_3D(ftmp,star,dim);

    
    
    double oparas[5],rparas[5],cmpara[5],cents[2];
    for(ic=0;ic<Nstar;ic++){
      star_gaus(star[ic],new_Ng,new_Ng,starf,residu,para);
      center[ic][0]=para[3];center[ic][1]=para[4];
      for(i=0;i<new_Ng;i++){
        for(j=0;j<new_Ng;j++){
          new_a[ic][i*new_Ng+j]=star[ic][i][j];
        }
      }
    }
    estimate_para(new_a,center,Nstar,new_Ng,0,para);
    double rg=para[0]*3.5;printf("rg=%f\n",rg);rg=3.;
    count=0;
    for(ic=0;ic<Nstar;ic++){
      size(star[ic],Ng,center[ic],rg,oparas);
      star_gaus(star0[ic],new_Ng,new_Ng,starf,residu,para);
      cents[0]=para[3];cents[1]=para[4];
      size(star0[ic],new_Ng,cents,rg,rparas);
      cmpara[1]=fabs(oparas[1]-rparas[1]);
      cmpara[2]=fabs(oparas[2]-rparas[2]);
      cmpara[0]=fabs((oparas[0]-rparas[0])/oparas[0]);
      //printf("%f\t%f\t%f\n",cmpara[1],cmpara[2],cmpara[0]);
      if(cmpara[1]<=0.05 && cmpara[2]<=0.05 && cmpara[0]<=0.1){
        for(k=0;k<npc;k++){//copy coeff
          coeff[count][k][0]=coeff[ic][k][0];
          coeff[count][k][1]=coeff[ic][k][1];
        }
        wcoff[count][0]=wcoff[ic][0];
        wcoff[count][1]=wcoff[ic][1];
        spos[count][0]=spos[ic][0];
        spos[count][1]=spos[ic][1];
        for(i=0;i<Ng;i++)for(j=0;j<Ng;j++){
          star[count][i][j]=star[ic][i][j];
          star0[count][i][j]=star0[ic][i][j];
        }
        for(i=0;i<Ng0;i++)for(j=0;j<Ng0;j++){
          stars[count][i][j]=stars[ic][i][j];
          wgt[count][i][j]=wgt[ic][i][j];
        }
        for(i=0;i<(Ng0+5)*osam;i++)for(j=0;j<(Ng0+5)*osam;j++){
          models[count][i][j]=models[ic][i][j];
        }
        count++;
      }
    }*Nobj=count;Nstar=count;//printf("Nobj=%d\n",count);
  //////////////////////////////////////////////// 
    for(ic=0;ic<Nstar;ic++){
      sum=0;sumt=0;
      for(i=0;i<new_Ng;i++)for(j=0;j<new_Ng;j++){
        sum+=star[ic][i][j];sumt+=star0[ic][i][j];
      }
      //printf("fraction sum=%f\n",sumt/sum );
    }


    *Nobj=Nstar;//printf("Nstar=%d\n",Nstar);

    strcpy(ftmp,"fits/slecstar.fits");
    dim[0]=dim[1]=Ng,dim[2]=count;
    //write_fits_3D(ftmp,star0,dim);

  free_f3tensor(star0,0,Nstar0,0,Ng0,0,Ng0);
  free_f3tensor(star,0,Nstar0,0,Ng0,0,Ng0);
  free_f3tensor(residual,0,Nstar0,0,Ng0,0,Ng0);
  free_dmatrix(a,0,Nstar0,0,Mp);
  free_dmatrix(error0,0,Nstar0,0,Mp);
  free_dmatrix(weight,0,Nstar0,0,Mp);
  free_matrix(start,0,Ng0,0,Ng0);
  free_matrix(starf,0,Ng0,0,Ng0);
  free_matrix(residu,0,Ng0,0,Ng0);
  free_dmatrix(cent,0,Nstar0,0,Ng0);
  free_dmatrix(center,0,Nstar0,0,Ng0);
  free_dmatrix(ebasisf,0,NPCmax,0,Mpke);
  free_d3tensor(compbf,0,Nstar0,0,npc,0,Mp);
  free_dmatrix(compcoeff,0,npc,0,Nstar0);
  free_dmatrix(compbfr,0,npc,0,Mpke);
  free_d3tensor(recept,0,NPCmax,0,Nke,0,Nke);
  free_d3tensor(basis,0,npc,0,Nge,0,Nge);
  free_dmatrix(coefferr,0,Nstar0,0,npc);
  free_f3tensor(ebasisfs,0,rnbf0,0,Nke,0,Nke);

}


void align_model(float **models,double *center,int Ng0,int Ng,int osam,float **aligned){
  int i,j,it,jt,k,ra,rb,eNg0=Ng0*osam,eNg=Ng*osam;
  double cent[2],nbond=(eNg0-eNg)/2.+1,sum;
  float **ealigned;
  ealigned=matrix(0,eNg,0,eNg);

  cent[0]=center[0]*osam+osam/2.;
  cent[1]=center[1]*osam+osam/2.;
  BInterp_bicubic(eNg0,eNg0,models,nbond,eNg,eNg,ealigned,cent);
  sum=0;
  for(i=0;i<Ng;i++){
    for(j=0;j<Ng;j++){
      aligned[i][j]=0;
      for(it=0;it<osam;it++){
        ra=i*osam+it;
        for(jt=0;jt<osam;jt++){
          rb=j*osam+jt;
          aligned[i][j]+=ealigned[ra][rb];
        }
      }
    }
  }

  free_matrix(ealigned,0,eNg,0,eNg);
}

void galaxy_centra(float **y,int nx,int ny,double *pa){
  int i,j,k;
  double sum,sumx,sumy,ix,iy,sum1,**weight,e1,e2;
  double sum2,ra,rb,r,sigma,mean,para[10],center[2],it,jt;
  float **I;
  void size(float **I,int Ng,double *cen,double rg,double *se);
  weight=dmatrix(0,nx,0,ny);
  I=matrix(0,nx,0,ny);
  sum=sumx=sumy=0;
  for(i=0;i<nx;i++){
    for(j=0;j<ny;j++){
      k=i*nx+j;
      sum+=y[i][j];
    }
  }mean=sum/nx/ny;
  for(i=0;i<nx;i++){
    for(j=0;j<ny;j++){
      k=i*nx+j;
      sum+=DSQR(y[i][j]-mean);
    }
  }sigma=sum/(nx*ny-1);

  sum=sumx=sumy=0;
  for(i=0;i<nx;i++){
    for(j=0;j<ny;j++){
      sum+=y[i][j];
      sumx+=y[i][j]*i;
      sumy+=y[i][j]*j;
    }
  }
  ix=sumx/sum;iy=sumy/sum;
  center[0]=ix;center[1]=iy;
  sigma=DSQR(3);
  size(I,nx,center,sigma,para);
  e1=para[1];e2=para[2];
  for(i=0;i<nx;i++)for(j=0;j<ny;j++){
    it=i-ix;jt=j-iy;
    ra=(1.-e1)*it-e2*jt;
    rb=(1.+e1)*jt-e2*it;
    r=ra*ra+rb*rb;
    weight[i][j]=exp(-r/sigma/2.);
  }
  do{
    sum=sumx=sumy=0;
    for(i=0;i<nx;i++){
      for(j=0;j<ny;j++){
        ra=i-ix;rb=j-iy;
        r=ra*ra+rb*rb;
        //sum+=y[i][j]*weight[i][j];
        //sumx+=y[i][j]*i*weight[i][j];
        //sumy+=y[i][j]*j*weight[i][j];
        sum+=y[i][j]*exp(-r/sigma/2);
        sumx+=y[i][j]*i*exp(-r/sigma/2);
        sumy+=y[i][j]*j*exp(-r/sigma/2);
      }
    }
    sum1=sumx/sum;sum2=sumy/sum;
    mean=(sum1-ix)*(sum1-ix)+(sum2-iy)*(sum2-iy);
    ix=sum1;iy=sum2;
    center[0]=ix;center[1]=iy;
    size(I,nx,center,sigma,para);
    e1=para[1];e2=para[2];
    for(i=0;i<nx;i++)for(j=0;j<ny;j++){
      it=i-ix;jt=j-iy;
      ra=(1.-e1)*it-e2*jt;
      rb=(1.+e1)*jt-e2*it;
      r=ra*ra+rb*rb;
      weight[i][j]=exp(-r/sigma/2.);
    }
  }while(mean<1.0e-14);
  pa[3]=ix;pa[4]=iy;

  free_dmatrix(weight,0,nx,0,ny);
  free_matrix(I,0,nx,0,ny);
}

