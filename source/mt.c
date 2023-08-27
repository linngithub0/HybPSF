#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"
#include <string.h>
#include "function.h"
#include "dirent.h" 
#define max( a, b) ( (a) > (b) ? (a) : (b) )
#define min( a, b) ( (a) < (b) ? (a) : (b) )
int main(int argc,char *argv[]){
  if(argc!=2){
    printf("the input should be arranged as: ./cut num\n");
    exit(EXIT_FAILURE);
  }printf("input num=%s\n",argv[1]);

  long iseed=-1;
  int t=atof(argv[1]),l,f,m,n,l0=0,Nhead,dim[5],i,j,k,it,jt;//preparing for the ploynomials;
  int f_num=1,ic,colum=9,Nstar0,Nstar,Nobj,Nstari,Nobj_deli,file_num;
  int Ngp,hgp,Ngi,Ng=23,Mp,chip_num=36,NPC,Ex,*Ngal,Ncom=1,charnum=500,Ng0=Ng+5;
  double *NAXIS1,*NAXIS2,*CRPIX1,*CRPIX2,keyva,*date_obs;
  double *CD1_1,*CD1_2,*CD2_1,*CD2_2,*CRVAL1,*CRVAL2,*GAIN;
  char **fname_image,**psf_catalogue,kywod[charnum],time_obs[charnum];
  char fname_image_ex[charnum],ftmp[charnum],num[charnum];
  char delim[charnum],tmp[charnum],fname[charnum*1000],dirs[charnum];
  char chip[charnum],fsuffix[charnum],file[charnum];
  float ***stars,**mpara,**spos,rx,ry,**starf,**residu,***start;
  double paras[10],**a0,**center,rg,sum,**oparas,para[10],cmpara[5];
  FILE *fp;
  printf("double=%lu\n",sizeof(double) );


  GAIN=dvector(0,f_num);
  starf=matrix(0,Ng0,0,Ng0);
  residu=matrix(0,Ng0,0,Ng0);
  file_num=20;
  fname_image=(char **)calloc(file_num,sizeof(char *));
  psf_catalogue=(char **)calloc(file_num,sizeof(char *));
  for(i=0;i<file_num;i++){
    fname_image[i]=(char *)calloc(charnum,sizeof(char));
    psf_catalogue[i]=(char *)calloc(charnum,sizeof(char));
  }


  char fcata[5][charnum];
  char keyna[100]={"GAIN1"};
  int cNobj=0;
  chip_num=10;
  if(t<10){
    strcpy(chip,"0");itoa(t,num,10);strcat(chip,num);
  }
  else{
    itoa(t,num,10);strcpy(chip,num);
  }
  //strcpy(chip,"09");
  strcpy(fsuffix,chip);strcat(fsuffix,"_raw.fits");printf("%s\n",fsuffix );
      
      for(i=1;i<=chip_num;i++){
        strcpy(dirs,"/Volumes/SHIN/Nie/psfTest/TEST_150/MSC_00000");
        if(i<10){
          strcpy(ftmp,"0");
          itoa(i,tmp,10);strcat(ftmp,tmp);strcpy(tmp,ftmp);
        }
        else{
          itoa(i,tmp,10);
        }
        strcat(dirs,tmp);strcat(dirs,"/");
        printf("%s\n",dirs);
        read_files(dirs,fname,&n,fsuffix,charnum);
        //printf("return:%s\n",fname );
        for(ic=0;ic<n;ic++){
          for(j=0;j<charnum;j++){
            fcata[ic][j]=fname[ic*charnum+j];
          }//printf("get fname:%s\n",fcata[ic] );
          if(fcata[ic][0]=='r' && fcata[ic][1]=='e'){
            strcpy(ftmp,fcata[ic]);
          }
        }printf("get fname1:%s\n",ftmp );
        strcpy(psf_catalogue[i],dirs);strcat(psf_catalogue[i],"catalogue/star/");
        strcat(psf_catalogue[i],ftmp);
        strcat(psf_catalogue[i],".cat");
        strcpy(fname_image[i],dirs);
        strcat(fname_image[i],ftmp);
        strcat(fname_image[i],"[raw]");
        printf("%s\n",psf_catalogue[i] );
        printf("%s\n",fname_image[i] );
        funcpsfs_size(psf_catalogue[i],&Nstar0,&Nhead);colum=Nhead;
        cNobj+=Nstar0;
        //append the suffix finally
      }printf("cNobj=%d\n",cNobj );
      if(cNobj<40){
        printf("The star number is less than 40, exit!\n");
        exit(EXIT_FAILURE);
      }
      Nstar=120;//only select part of the stars to do reconstruction
      //Nstar=cNobj;

      a0=dmatrix(0,cNobj,0,Ng0*Ng0);
      center=dmatrix(0,cNobj,0,2);
      mpara=matrix(0,cNobj,0,colum);
      stars=f3tensor(0,cNobj,0,Ng0,0,Ng0);
      start=f3tensor(0,cNobj,0,Ng0,0,Ng0);
      spos=matrix(0,500,0,2);
      int nimgx,nimgy,imgsize,Nh=Ng0/2;
      float *image0,**image;
      get_fits_size(fname_image[1],dim);
      get_header_key(fname_image[1],&keyva,keyna);//printf("gain=%f\n",keyva );
      nimgx=dim[0];nimgy=dim[1];imgsize=nimgx*nimgy;
      image0=vector(0,imgsize-1);
      image=matrix(0,nimgx-1,0,nimgy-1);
      printf("nx=%d,ny=%d\n",nimgx,nimgy );
      strcpy(delim,"\t");
      m=0;
      for(ic=1;ic<=chip_num;ic++){
        read_fits_2D(fname_image[ic],image0,imgsize);
        for(j=0;j<nimgy;j++){
          for(i=0;i<nimgx;i++){k=j*nimgx+i;image[i][j]=image0[k];}
        }//printf("read_fits_2D\n");
        funcpsfs_size(psf_catalogue[ic],&Nstar0,&Nhead);
        //printf("Nstar0=%d\n",colum );
        funcpsfs_read(psf_catalogue[ic],Nstar0,Nhead,mpara,colum,delim);
        printf("funcpsfs_read\n");
        for(i=0;i<Nstar0;i++){
          //printf("i=%d\n",i );
          spos[m][0]=mpara[i][1];
          spos[m][1]=mpara[i][2];
          rx=spos[m][0];ry=spos[m][1];
          for(l=0;l<Ng0;l++)for(f=0;f<Ng0;f++){
            it=rx+l-Nh;jt=ry+f-Nh;
            if(it>=0 && it<nimgx && jt>=0 && jt<nimgy){
              stars[m][l][f]=image[it][jt];
            }
          }m++;
        }printf("m=%d,ic=%d\n",m,ic );
      }printf("m=%d\n",m );
      strcpy(ftmp,"fits/instars0.fits");
      dim[0]=dim[1]=Ng0,dim[2]=cNobj;
      write_fits_3D(ftmp,stars,dim);
      //exit(EXIT_FAILURE);




      void method_branches(float ***stars,float **spos,int npc,int Nstar,int Ng0,float ***PCs,
        float ***coeff,float gain,int method,int *Nobj,float *means);
      int npc=15;
      float ***PCs,***coeff,gain=keyva,means[cNobj];
      PCs=f3tensor(0,npc,0,Ng0,0,Ng0);
      coeff=f3tensor(0,cNobj,0,npc,0,2);printf("here\n");
      method_branches(stars,spos,npc,Nstar,Ng,PCs,coeff,gain,2,&Nobj,means);
      strcpy(ftmp,"fits/PCs.fits");
      dim[0]=dim[1]=Ng,dim[2]=npc;
      write_fits_3D(ftmp,PCs,dim);printf("Nobj=%d\n",Nobj );
      
      printf("method_branches out\n");
      ///////////////////////////////////////////////////////////////////////
      //read clean images
      for(ic=0;ic<Nobj;ic++){
        for(i=0;i<Ng;i++)for(j=0;j<Ng;j++){
          start[ic][i][j]=stars[ic][i][j];
          a0[ic][i*Ng+j]=stars[ic][i][j];
        }
      }
      strcpy(ftmp,"fits/instars.fits");
      dim[0]=dim[1]=Ng,dim[2]=Nobj;
      write_fits_3D(ftmp,start,dim);
      oparas=dmatrix(0,Nobj,0,5);
      for(ic=0;ic<Nobj;ic++){
        star_gaus(start[ic],Ng,Ng,starf,residu,paras);
        center[ic][0]=paras[3];center[ic][1]=paras[4];
      }
      estimate_para(a0,center,Nobj,Ng,0,para);
      rg=para[0]*3.5;printf("rg=%f\n", rg);
      //rg=3.5;
      strcpy(ftmp,"data/TEST150/");
      strcat(ftmp,chip);strcat(ftmp,"oe.dat");printf("%s\n",ftmp );//original stars
      fp=fopen(ftmp,"w");
      for(ic=0;ic<Nobj;ic++){
        size(start[ic],Ng,center[ic],rg,oparas[ic]);
        fprintf(fp, "%e\t%e\t%e\t%e\t%e\n",oparas[ic][0],oparas[ic][1],oparas[ic][2],spos[ic][0],spos[ic][1] );
      }fclose(fp);

      printf("recostruction\n");
      //reconstruct PSFs
      for(ic=0;ic<Nobj;ic++){
        for(i=0;i<Ng;i++)for(j=0;j<Ng;j++){
          start[ic][i][j]=0;
          for(l=0;l<npc;l++){
            start[ic][i][j]+=PCs[l][i][j]*coeff[ic][l][0];
          }
        }
      }
      for(ic=0;ic<Nobj;ic++){
        star_gaus(start[ic],Ng,Ng,starf,residu,paras);
        center[ic][0]=paras[3];center[ic][1]=paras[4];
      }
      strcpy(ftmp,"data/TEST150/");
      strcat(ftmp,chip);strcat(ftmp,"re.dat");printf("%s\n",ftmp );//reconstructed stars
      fp=fopen(ftmp,"w");printf("Ng=%d\n",Ng );
      for(ic=0;ic<Nobj;ic++){
        size(start[ic],Ng,center[ic],rg,paras);
        fprintf(fp, "%e\t%e\t%e\t%e\t%e\n",paras[0],paras[1],paras[2],spos[ic][0],spos[ic][1] );
        cmpara[0]=fabs(oparas[ic][0]-paras[0])/oparas[ic][0];
        cmpara[1]=fabs(oparas[ic][1]-paras[1]);
        cmpara[2]=fabs(oparas[ic][2]-paras[2]);
        //printf("%f,%f,%f\n",paras[0],paras[1],paras[2] );
      }fclose(fp);
      for(ic=0;ic<Nobj;ic++){
        for(i=0;i<Ng;i++)for(j=0;j<Ng;j++){
          start[ic][i][j]=stars[ic][i][j];
        }
      }
      strcpy(ftmp,"fits/restar.fits");
      dim[0]=dim[1]=Ng;dim[2]=Nobj;
      write_fits_3D(ftmp,start,dim); 
      

      void coeff_interp(int order_PC,int order_poly,int Nstar,float **spos,float **Coeff,
        float **weight,float ***PCs,int ngp,float ***rPSFs,int Nobj,float **psfpos);
      void Krig_interp(int order_PC,int Nstar,float **spos,float **Coeff,float **Coefferr,
        float ***PCs,int ngp,float ***rPSFs,int Nobj,float **psfpos);
      void pix_interp(int Nstar,float **spos,float ***stars,int ngp,int order_poly,
        float ***rPSFs,int Nobj,float **psfpos);
      float **Coeff,**weight,**Coefferr;
      float **psfpos,***rPSFs;
      int order_poly=6;
      double imean,isigma;
      l0=0;
      for(i=order_poly+1;i>0;i--)l0+=i;
      Coeff=matrix(0,cNobj,0,npc);
      Coefferr=matrix(0,cNobj,0,npc);
      weight=matrix(0,cNobj,0,npc);
      psfpos=matrix(0,cNobj,0,2);
      rPSFs=f3tensor(0,cNobj,0,Ng,0,Ng);

      strcpy(ftmp,"results/coeffs.dat");
      fp=fopen(ftmp,"w");
      for(ic=0;ic<Nobj;ic++){
        for(l=0;l<npc;l++){
          Coeff[ic][l]=coeff[ic][l][0];
          weight[ic][l]=1;//./coeff[ic][l][1];
          Coefferr[ic][l]=coeff[ic][l][1];
          fprintf(fp, "%f\t",Coeff[ic][l] );
        }fprintf(fp, "\n" );
      }fclose(fp);
      
      coeff_interp(npc,order_poly,Nobj,spos,Coeff,weight,PCs,Ng,rPSFs,Nobj,spos);
      //Using constructed polynomial to interpolate PSFs at stars,\
       which are used to reconstructions
      //Krig_interp(npc,Nobj,spos,Coeff,Coefferr,PCs,Ng,rPSFs,Nobj,spos);
      //pix_interp(Nobj,spos,stars,Ng,order_poly,rPSFs,Nobj,spos);
      printf("coeff_interp down\n");

      for(ic=0;ic<Nobj;ic++){
        star_gaus(rPSFs[ic],Ng,Ng,starf,residu,paras);
        center[ic][0]=paras[3];center[ic][1]=paras[4];
      }
      strcpy(ftmp,"data/TEST150/");
      strcat(ftmp,chip);strcat(ftmp,"ioe.dat");printf("%s\n",ftmp );
      fp=fopen(ftmp,"w");
      for(ic=0;ic<Nobj;ic++){
        size(rPSFs[ic],Ng,center[ic],rg,paras);
        fprintf(fp, "%e\t%e\t%e\t%e\t%e\n",paras[0],paras[1],paras[2],spos[ic][0],spos[ic][1] );
      }fclose(fp);
      printf("calculated\n");
      for(ic=0;ic<Nobj;ic++){
        for(i=0;i<Ng;i++)for(j=0;j<Ng;j++){
          rPSFs[ic][i][j]-=start[ic][i][j];
        }
      }

      strcpy(ftmp,"fits/restars.fits");
      dim[0]=dim[1]=Ng,dim[2]=Nobj;
      write_fits_3D(ftmp,rPSFs,dim); 

      free_matrix(mpara,0,Nstar0,0,colum);
      free_matrix(psfpos,0,Nstar0,0,2);
      free_f3tensor(stars,0,Nstar0,0,Ng0,0,Ng0);
      //free_f3tensor(start,0,Nstar0,0,Ng0,0,Ng0);
      free_dmatrix(a0,0,Nstar0,0,Ng0*Ng0);
      free_dmatrix(center,0,Nstar0,0,2);
      free_f3tensor(rPSFs,0,Nstar0,0,Ng,0,Ng);


      //testing for interpolated PSFs
      ///////////////////////////////////////////////////////////////////////
      cNobj=0;
      int ichip_num=file_num-chip_num;//ichip_num=file_num;
      printf("ichip_num=%d,%d,%d\n",ichip_num,file_num,chip_num );
      strcpy(fsuffix,chip);strcat(fsuffix,"_raw.fits");printf("%s\n",fsuffix );
      
      for(i=0;i<ichip_num;i++){
        strcpy(dirs,"/Volumes/SHIN/Nie/psfTest/TEST_300/MSC_00000");
        //strcpy(dirs,"/Volumes/SHIN/Nie/psfTest/TEST_150/MSC_00000");
        it=i+chip_num;//printf("it=%d\n",it );
        //it=i+1;
        if(it<10){
          strcpy(ftmp,"0");
          itoa(it,tmp,10);strcat(ftmp,tmp);strcpy(tmp,ftmp);
        }
        else{
          itoa(it,tmp,10);
        }
        strcat(dirs,tmp);strcat(dirs,"/");
        printf("dirs: %s\n",dirs);
        read_files(dirs,fname,&n,fsuffix,charnum);
        //printf("return:%s\n",fname );
        for(ic=0;ic<n;ic++){
          for(j=0;j<charnum;j++){
            fcata[ic][j]=fname[ic*charnum+j];
          }//printf("get fname:%s\n",fcata[ic] );
          //if(fcata[ic][0]=='r' && fcata[ic][1]=='e'){
            
          //}
        }
        strcpy(ftmp,fcata[0]);
        //printf("get fname: %s\n",ftmp );
        strcpy(psf_catalogue[i],dirs);
        strcat(psf_catalogue[i],"300/catalogue/star/");
        //strcat(psf_catalogue[i],"catalogue/star/re");
        strcat(psf_catalogue[i],ftmp);
        strcat(psf_catalogue[i],".cat");//printf("%s\n",psf_catalogue[i] );
        strcpy(fname_image[i],dirs);
        strcat(fname_image[i],ftmp);
        strcat(fname_image[i],"[raw]");
        //printf("%s\n",psf_catalogue[i] );
        //printf("%s\n",fname_image[i] );
        funcpsfs_size(psf_catalogue[i],&Nstar0,&Nhead);colum=Nhead;
        printf("Nstar0=%d\n",Nstar0 );
        cNobj+=Nstar0;
        //append the suffix finally
      }printf("cNobj=%d\n",cNobj );
      Nstari=cNobj;
      printf("read catalogue down\n");

      a0=dmatrix(0,Nstari,0,Ng*Ng);
      center=dmatrix(0,Nstari,0,2);
      mpara=matrix(0,Nstari,0,colum);
      stars=f3tensor(0,Nstari,0,Ng,0,Ng);
      //start=f3tensor(0,Nstari,0,Ng,0,Ng);
      psfpos=matrix(0,Nstari,0,2);
      rPSFs=f3tensor(0,Nstari,0,Ng,0,Ng);
      Nh=Ng/2;
      strcpy(delim,"\t");
      m=0;
      for(ic=0;ic<ichip_num;ic++){
        read_fits_2D(fname_image[ic],image0,imgsize);
        for(j=0;j<nimgy;j++){
          for(i=0;i<nimgx;i++){k=j*nimgx+i;image[i][j]=image0[k];}
        }//printf("read_fits_2D\n");
        funcpsfs_size(psf_catalogue[ic],&Nstar0,&Nhead);
        //printf("Nstar0=%d\n",colum );
        funcpsfs_read(psf_catalogue[ic],Nstar0,Nhead,mpara,colum,delim);
        //printf("funcpsfs_read\n");
        for(i=0;i<Nstar0;i++){
          //printf("i=%d\n",i );
          psfpos[m][0]=mpara[i][1]-0.5;
          psfpos[m][1]=mpara[i][2]-0.5;
          rx=psfpos[m][0];ry=psfpos[m][1];
          for(l=0;l<Ng;l++)for(f=0;f<Ng;f++){
            it=rx+l-Nh;jt=ry+f-Nh;
            if(it>=0 && it<nimgx && jt>=0 && jt<nimgy){
              stars[m][l][f]=image[it][jt];
            }
          }m++;
        }//printf("m=%d,ic=%d\n",m,ic );
      }printf("m=%d\n",m );
      strcpy(ftmp,"fits/stars.fits");
      dim[0]=dim[1]=Ng,dim[2]=Nstari;
      write_fits_3D(ftmp,stars,dim);

      for(ic=0;ic<Nstari;ic++){
        //gaus_estimate(stars[ic],Ng,Ng,&imean,&isigma);
        //printf("ic=%d,means=%f\n",ic,imean );
        //for(i=0;i<Ng;i++)for(j=0;j<Ng;j++)stars[ic][i][j]-=imean;
        star_gaus(stars[ic],Ng,Ng,starf,residu,paras);
        center[ic][0]=paras[3];center[ic][1]=paras[4];
      }
      strcpy(ftmp,"data/TEST150/");
      strcat(ftmp,chip);strcat(ftmp,"oie.dat");printf("%s\n",ftmp );
      fp=fopen(ftmp,"w");
      for(ic=0;ic<Nstari;ic++){
        size(stars[ic],Ng,center[ic],rg,paras);
        fprintf(fp, "%e\t%e\t%e\t%e\t%e\n",paras[0],paras[1],paras[2],psfpos[ic][0],psfpos[ic][1] );
      }fclose(fp);
      printf("read down\n");


      printf("start coeff_interp\n");
      coeff_interp(npc,order_poly,Nobj,spos,Coeff,weight,PCs,Ng,rPSFs,Nstari,psfpos);
      //Krig_interp(npc,Nobj,spos,Coeff,Coefferr,PCs,Ng,rPSFs,Nstari,psfpos);
      //pix_interp(Nobj,spos,start,Ng,order_poly,rPSFs,Nstari,psfpos);
      printf("coeff_interp down\n");
      strcpy(ftmp,"fits/rPSFs.fits");
      dim[0]=dim[1]=Ng,dim[2]=Nstari;
      write_fits_3D(ftmp,rPSFs,dim);

      for(ic=0;ic<Nstari;ic++){
        star_gaus(rPSFs[ic],Ng,Ng,starf,residu,paras);
        center[ic][0]=paras[3];center[ic][1]=paras[4];
        //printf("x=%f,y=%f\n",paras[3],paras[4] );
      }
      strcpy(ftmp,"data/TEST150/");
      strcat(ftmp,chip);strcat(ftmp,"ie.dat");printf("%s\n",ftmp );
      fp=fopen(ftmp,"w");
      for(ic=0;ic<Nstari;ic++){
        size(rPSFs[ic],Ng,center[ic],rg,paras);
        fprintf(fp, "%e\t%e\t%e\t%e\t%e\n",paras[0],paras[1],paras[2],psfpos[ic][0],psfpos[ic][1] );
      }fclose(fp);
      printf("calculated\n");

      exit(EXIT_FAILURE);


      

}
