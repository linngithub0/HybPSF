#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"
#include <string.h>
#include <Python.h>
#include "function.h"
#define max( a, b) ( (a) > (b) ? (a) : (b) )
#define min( a, b) ( (a) < (b) ? (a) : (b) )
struct filte{
  int ndetecor;
  char filte_name[100];
  char slice_name[100];
  char detector[8][100];
  float gain[8];
};
int main(int argc,char *argv[]){
  if(argc!=2){
    printf("the input should be arranged as: ./KSB+ num\n");
    exit(EXIT_FAILURE);
  }printf("input num=%s\n",argv[1]);
  void get_data(char *stamp_name,char *ofname,char *dirs,char *PCdirs,char *det,char *filte);
  void pre_HybPSF(char *stamp_name,char *ofname,char *dirs,char *PCdirs,
    char *det,char *filte,float ***stars,float ***imodel,
    float **spos,float ***wgt,int *Nstar,int *Ng0,int Nim);
  void dignose(char *stamp_name,char *ofname,char* dirs,char *PCdirs,char *det,
    char *filte,float ***PCs,float **pos,float ***coeff,int getNobj,
    float **wcoeff,float ***imodel,int Ng0,int Ng,int osam,int Nb,double **center,int Nim);
  int send_name(char *stamp_name);
  Py_Initialize();
  if (!Py_IsInitialized()){
    printf("nitialising failed!\n");
    Py_Finalize();
    exit(EXIT_FAILURE);
  }
  PyRun_SimpleString("import sys");
  PyRun_SimpleString("sys.path.append('./')");

  long iseed=-1;
  int t=atof(argv[1]),l,f,m,n,l0=0,Nhead,dim[5],i,j,k;//preparing for the ploynomials;
  int f_num=1,ic,colum=9,Nstar0,Nobj,Nstari;
  int Ngp,hgp,Ngi,Mp,chip_num=36,NPC,Ex,*Ngal,Ncom=1,charnum=200;
  double *NAXIS1,*NAXIS2,*CRPIX1,*CRPIX2,keyva,*date_obs;
  double *CD1_1,*CD1_2,*CD2_1,*CD2_2,*CRVAL1,*CRVAL2,*GAIN;
  char fname_image[charnum],psf_catalogue[charnum],kywod[charnum],time_obs[charnum],fname_image_ex[charnum],ftmp[charnum];
  char delim[charnum];
  float ***stars,**spos,rx,ry,**starf,**residu,***wgt;
  float ***imodel,***coeff,***PCs,**wcoff,snrs=0.0001;
  int npc=4,Ng0,Nstar,Ng,Nim=5,osam=2;
  double **paras,**a0,rg,sum,**oparas,para[10];
  double mean,sigma,**center;
  FILE *fp,*fp1,*fp2;
  int detec_long_num=2,detec_short_num=8,nbands=3;
  float get_gain;
  char *ldetec[]={"nrcalong","nrcblong"};
  char *sdetec[]={"nrca1","nrca2","nrca3","nrca4","nrcb1","nrcb2","nrcb3","nrcb4"};
  float sgain[]={2.08,2.02,2.17,2.02,2.01,2.14,1.94,2.03};
  float lgain[]={1.84,1.80};
  char *filt_sname[]={"jw02736001001_02101_","jw02736001001_02103_","jw02736001001_02105_"};
  char *filters[]={"F090W","F150W","F200W","F277W","F356W","F444W"};
  char odirs[100];
  char ofname[500],stamp_name[500],cat_name[500],mask_name[500];
  char dirs[100];
  char PCdirs[100];
  char stamp_dirs[500],PCdir[500];
  char indx[100],indx_tmp[10],tmp[100];
  strcpy(indx,"0000");
  stars=f3tensor(0,100,0,100,0,100);
  wgt=f3tensor(0,100,0,100,0,100);
  imodel=f3tensor(0,100,0,200,0,200);
  spos=matrix(0,100,0,2);
  coeff=f3tensor(0,100,0,10,0,2);
  PCs=f3tensor(0,10,0,200,0,200);
  wcoff=matrix(0,100,0,2);
  residu=matrix(0,100,0,100);
  starf=matrix(0,100,0,100);
  center=dmatrix(0,100,0,2);
  strcpy(odirs,"/data/SMACS0723_NIRCam_pipeline1.8.1/");
  strcpy(dirs,"/home/lnie/data/JWST_rev/");
  strcpy(PCdirs,"/home/lnie/data/PC/JWST/");

  struct filte sbands[3],lbands[3];
  //set the parameter for short filter
  for(ic=0;ic<nbands;ic++){
    strcpy(sbands[ic].filte_name, filters[ic]);
    strcpy(sbands[ic].slice_name, filt_sname[ic]);
    sbands[ic].ndetecor=detec_short_num;
    for(i=0;i<sbands[ic].ndetecor;i++){
      strcpy(sbands[ic].detector[i],sdetec[i]);
      sbands[ic].gain[i]=sgain[i];
    }
    /*printf("%s.slice_name: %s\t",sbands[ic].filte_name,sbands[ic].slice_name );
    for(i=0;i<sbands[ic].ndetecor;i++){
      printf("%s.detector: %s\t",sbands[ic].filte_name,sbands[ic].detector[i] );
      printf("%s.gain: %f\t",sbands[ic].filte_name,sbands[ic].gain[i] );
    }printf("\n");*/
  }
  //set the parameter for long filter
  for(ic=0;ic<nbands;ic++){
    strcpy(lbands[ic].filte_name, filters[ic+3]);
    strcpy(lbands[ic].slice_name, filt_sname[ic]);
    lbands[ic].ndetecor=detec_long_num;
    for(i=0;i<lbands[ic].ndetecor;i++){
      strcpy(lbands[ic].detector[i],ldetec[i]);
      lbands[ic].gain[i]=lgain[i];
    }
    printf("%s.slice_name: %s\t",lbands[ic].filte_name,lbands[ic].slice_name );
    for(i=0;i<lbands[ic].ndetecor;i++){
      printf("%s.detector: %s\t",lbands[ic].filte_name,lbands[ic].detector[i] );
      printf("%s.gain: %f\t",lbands[ic].filte_name,lbands[ic].gain[i] );
    }printf("\n");
  }
  
  
  /*for(Ex=0;Ex<nbands;Ex++){
  for(ic=0;ic<sbands[Ex].ndetecor;ic++){
  //for(ic=0;ic<1;ic++){
    get_gain=sbands[Ex].gain[ic];
    strcpy(stamp_dirs,dirs);strcat(stamp_dirs,"NIRCam/");
    strcat(stamp_dirs,sbands[Ex].filte_name);strcat(stamp_dirs,"/");
    strcat(stamp_dirs,sbands[Ex].detector[ic]);strcat(stamp_dirs,"/");
    strcat(stamp_dirs,"catalogue/star_stamps/");
    //printf("%s\n",stamp_dirs );
    for(int Tim=0;Tim<9;Tim++){
      int expos=Tim+1;
      itoa(expos,tmp,10);
      strcpy(ofname,odirs);strcat(ofname,sbands[Ex].filte_name);strcat(ofname,"/");
      strcat(ofname,sbands[Ex].slice_name);
      strcat(ofname,indx);strcat(ofname,tmp);strcat(ofname,"_");
      strcat(ofname,sbands[Ex].detector[ic]);
      strcat(ofname,"_cal.fits");
      strcpy(stamp_name,stamp_dirs);strcat(stamp_name,sbands[Ex].slice_name);
      strcat(stamp_name,indx);strcat(stamp_name,tmp);strcat(stamp_name,"_");
      strcat(stamp_name,sbands[Ex].detector[ic]);
      strcat(stamp_name,"_cal_pv181_rev4_star.fits");
      printf("stamp_name: %s\n", stamp_name);
      printf("ofname:%s\n",ofname);
      
      k=send_name(stamp_name);
      
      printf("ofname:%d\n",k);
      if(k==1){
        printf("do_HybPSF\n");fflush(stdout);
        pre_HybPSF(stamp_name,ofname,dirs,PCdirs,
          sbands[Ex].detector[ic],sbands[Ex].filte_name,stars,imodel,
          spos,wgt,&Nstar,&Ng0,Nim);
        printf("Nstar=%d,Ng0=%d\n",Nstar,Ng0 );
        printf("%s\n",sbands[Ex].detector[ic]);
        Ng=Ng0-4;
        webiSPCA_entr(stars,spos,npc,Nstar,Ng0,Ng,PCs,coeff,
          sbands[Ex].gain[ic],&Nobj,snrs,osam,imodel,wgt,wcoff);
        //estimate center
        for(f=0;f<Nstar;f++){
          gaus_estimate(stars[f],Ng0,Ng0,&mean,&sigma); //estimate the background noise as gasussian randoms
          //printf("ic=%d,mean=%f,sigma=%f\n",f+1,mean,sigma );
          for(i=0;i<Ng0;i++){
            for(j=0;j<Ng0;j++){
              stars[f][i][j]-=mean;//subtract the bacground sky intensity
            }
          }
          star_gaus(stars[f],Ng0,Ng0,starf,residu,para);
          center[f][0]=para[3];center[f][1]=para[4];
          //printf("x=%f,y=%f\n",para[3],para[4]);
        }
          //printf("start dignose\n");printf("%s\n",sbands[Ex].detector[ic]);
        dignose(stamp_name,ofname,dirs,PCdirs,sbands[Ex].detector[ic],
          sbands[Ex].filte_name,PCs,spos,coeff,Nobj,wcoff,imodel,Ng0,Ng,osam,npc,center,Nim);
        
      }
      
    }
  }printf("bands\n");
  }*/
    
    
  for(Ex=0;Ex<1;Ex++){
  //for(ic=0;ic<lbands[Ex].ndetecor;ic++){
  for(ic=0;ic<1;ic++){
    get_gain=lbands[Ex].gain[ic];
    strcpy(stamp_dirs,dirs);strcat(stamp_dirs,"NIRCam/");
    strcat(stamp_dirs,lbands[Ex].filte_name);strcat(stamp_dirs,"/");
    strcat(stamp_dirs,lbands[Ex].detector[ic]);strcat(stamp_dirs,"/");
    strcat(stamp_dirs,"catalogue/star_stamps/");
    //printf("%s\n",stamp_dirs );
    for(int Tim=0;Tim<1;Tim++){
      int expos=Tim+1;
      itoa(expos,tmp,10);
      strcpy(ofname,odirs);strcat(ofname,lbands[Ex].filte_name);strcat(ofname,"/");
      strcat(ofname,lbands[Ex].slice_name);
      strcat(ofname,indx);strcat(ofname,tmp);strcat(ofname,"_");
      strcat(ofname,lbands[Ex].detector[ic]);
      strcat(ofname,"_cal.fits");
      strcpy(stamp_name,stamp_dirs);strcat(stamp_name,lbands[Ex].slice_name);
      strcat(stamp_name,indx);strcat(stamp_name,tmp);strcat(stamp_name,"_");
      strcat(stamp_name,lbands[Ex].detector[ic]);
      strcat(stamp_name,"_cal_pv181_rev4_star.fits");
      printf("stamp_name: %s\n", stamp_name);
      printf("ofname:%s\n",ofname);
      
      k=send_name(stamp_name);
      
      printf("ofname:%d\n",k);
      if(k==1){
        printf("do_HybPSF\n");fflush(stdout);
        pre_HybPSF(stamp_name,ofname,dirs,PCdirs,
          lbands[Ex].detector[ic],lbands[Ex].filte_name,stars,imodel,
          spos,wgt,&Nstar,&Ng0,Nim);
        printf("Nstar=%d,Ng0=%d\n",Nstar,Ng0 );
        printf("%s\n",lbands[Ex].detector[ic]);
        Ng=Ng0-4;
        webiSPCA_entr(stars,spos,npc,Nstar,Ng0,Ng,PCs,coeff,
          lbands[Ex].gain[ic],&Nobj,snrs,osam,imodel,wgt,wcoff);
        //estimate center
        for(f=0;f<Nstar;f++){
          gaus_estimate(stars[f],Ng0,Ng0,&mean,&sigma); //estimate the background noise as gasussian randoms
          //printf("ic=%d,mean=%f,sigma=%f\n",f+1,mean,sigma );
          for(i=0;i<Ng0;i++){
            for(j=0;j<Ng0;j++){
              stars[f][i][j]-=mean;//subtract the bacground sky intensity
            }
          }
          star_gaus(stars[f],Ng0,Ng0,starf,residu,para);
          center[f][0]=para[3];center[f][1]=para[4];
          //printf("x=%f,y=%f\n",para[3],para[4]);
        }
          //printf("start dignose\n");printf("%s\n",lbands[0].detector[ic]);
        dignose(stamp_name,ofname,dirs,PCdirs,lbands[Ex].detector[ic],
          lbands[Ex].filte_name,PCs,spos,coeff,Nobj,wcoff,imodel,Ng0,Ng,osam,npc,center,Nim);
        
      }
      
    }
  }printf("bands\n");
  }

  Py_Finalize();

}

