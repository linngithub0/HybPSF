#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"
#include <string.h>
#include "function.h"
#define max( a, b) ( (a) > (b) ? (a) : (b) )
#define min( a, b) ( (a) < (b) ? (a) : (b) )
int main(){
  

  long iseed=-1;
  int t,l,f,m,n,l0=0,Nhead,dim[5],i,j,k;//preparing for the ploynomials;
  int f_num=1,ic,colum=9,Nstar0,Nstar,Nobj,Nstari,Nobj_deli;
  int Ngp,hgp,Ngi,Ng=23,Mp,chip_num=36,NPC,Ex,*Ngal,Ncom=1,charnum=200,Ng0=Ng+5;
  double *NAXIS1,*NAXIS2,*CRPIX1,*CRPIX2,**keyva,*date_obs;
  double *CD1_1,*CD1_2,*CD2_1,*CD2_2,*CRVAL1,*CRVAL2,*GAIN;
  char fname_image[charnum],psf_catalogue[charnum],kywod[charnum],time_obs[charnum],fname_image_ex[charnum],ftmp[charnum];
  char delim[charnum],tmp[charnum];
  float ***stars,**mpara,**spos,rx,ry,**starf,**residu,***start;
  double paras[10],**a0,**center,rg,sum,oparas[10],para[10];
  FILE *fp;
  printf("double=%lu\n",sizeof(double) );

  keyva=dmatrix(0,10,0,4);
  date_obs=dvector(0,f_num);
  NAXIS1=dvector(0,f_num);
  NAXIS2=dvector(0,f_num);
  CRPIX1=dvector(0,f_num);
  CRPIX2=dvector(0,f_num);
  CD1_1=dvector(0,f_num);
  CD1_2=dvector(0,f_num);
  CD2_1=dvector(0,f_num);
  CD2_2=dvector(0,f_num);
  CRVAL1=dvector(0,f_num);
  CRVAL2=dvector(0,f_num);
  GAIN=dvector(0,f_num);
  starf=matrix(0,Ng0,0,Ng0);
  residu=matrix(0,Ng0,0,Ng0);

  
      t=6;
      if(t<10){
        strcpy(ftmp,"0");
        itoa(t,tmp,10);strcat(ftmp,tmp);strcpy(tmp,ftmp);
      }
      else{
        itoa(t,tmp,10);
      }
      strcpy(fname_image,"data/TEST150/MSC_0000009/MSC_MS_210525145000_100000009_");
      strcat(fname_image,tmp);strcat(fname_image,"_raw.fits");
      strcat(fname_image,"[raw]");
      
      //strcpy(fname_image,"data/MSC_210304093000_0000010_07_raw_img.fits");
      //printf("%s\n",fname_image );
      //strcpy(fname_image_ex,fname_image);strcat(fname_image_ex,"[1]");
      get_fits_header(fname_image,keyva,time_obs);

      //get_fits_img_data(fname_image);


      //printf("fname_image\n");
      strcpy(psf_catalogue,"data/TEST150/50/catalogue/star/MSC_MS_210525145000_100000009_");
      strcat(psf_catalogue,tmp);strcat(psf_catalogue,"_raw.fits.cat");
      printf("%s\n",psf_catalogue);
      //strcpy(psf_catalogue,"data/MSC_210304093000_0000010_07.cat");
      funcpsfs_size(psf_catalogue,&Nstar0,&Nhead);//the catalogue should be arranged as "id x y SNR"
      printf("Nstar0=%d,Nhead=%d\n",Nstar0,Nhead );colum=Nhead;
      mpara=matrix(0,Nstar0,0,colum);
      spos=matrix(0,Nstar0,0,2);
      printf("colum=%d\n",colum );
      strcpy(delim,"\t");
      funcpsfs_read(psf_catalogue,Nstar0,Nhead,mpara,colum,delim);//Nstar0=250;
      for(ic=0;ic<Nstar0;ic++){
        spos[ic][0]=mpara[ic][1]+1;
        spos[ic][1]=mpara[ic][2]+1;
        //printf("%e\t%e\tic=%d\n",mpara[ic][1],mpara[ic][2],ic );
      }
      if(Nstar0<40){
        printf("The star number is less than 40, exit!\n");
        exit(EXIT_FAILURE);
      } Nstar=20;//Nobj_deli=Nstar0-Nstar;
      /*if(Nstar0>40){
        Nstar=Nstar0/2;Nstari=Nstar0-Nstar;
      }
      if(Nstar0>20 && Nstar0<=40){
        Nstar=20;Nstari=Nstar0-Nstar;
      }*/

      //Nstar=150;Nstari=Nstar0-Nstar;
      //only select half stars to do reconstruction

      a0=dmatrix(0,Nstar0,0,Ng0*Ng0);
      center=dmatrix(0,Nstar0,0,2);

     
      int nimgx,nimgy,imgsize,Nh,it,jt;
      float *image0,**image;
      get_fits_size(fname_image,dim);
      nimgx=dim[0];nimgy=dim[1];imgsize=nimgx*nimgy;
      image0=vector(0,imgsize-1);
      image=matrix(0,nimgx-1,0,nimgy-1);
      stars=f3tensor(0,Nstar0,0,Ng0,0,Ng0);
      start=f3tensor(0,Nstar0,0,Ng0,0,Ng0);
      read_fits_2D(fname_image,image0,imgsize);
      for(j=0;j<nimgy;j++){
        for(i=0;i<nimgx;i++){k=j*nimgx+i;image[i][j]=image0[k];}
      }
      
      Nh=Ng0/2;

      
      
      for(ic=0;ic<Nstar;ic++){
        rx=spos[ic][0]-1;ry=spos[ic][1]-1;
        for(i=0;i<Ng0;i++)for(j=0;j<Ng0;j++){
          it=rx+i-Nh;jt=ry+j-Nh;
          if(it>=0 && it<nimgx && jt>=0 && jt<nimgy){
            stars[ic][i][j]=image[it][jt];
            
          }
        }
      }

      strcpy(ftmp,"fits/stars.fits");
      dim[0]=dim[1]=Ng0,dim[2]=Nstar;
      write_fits_3D(ftmp,stars,dim);

      

      /*double mean,sigma;
      for(ic=0;ic<Nstar;ic++){
        //gaus_estimate(stars[ic],Ng0,Ng0,&mean,&sigma);//printf("sigma=%f,mean=%f,ic=%d\n",sigma,mean,ic+1 );
        mean=0;
        for(i=0;i<Ng0;i++)for(j=0;j<Ng0;j++){
          stars[ic][i][j]-=mean;
          a0[ic][i*Ng0+j]=stars[ic][i][j];
        }
        //printf("mean=%f,sigma=%f\n", mean,sigma);
        star_gaus(stars[ic],Ng0,Ng0,starf,residu,paras);
        center[ic][0]=paras[3];center[ic][1]=paras[4];
        //printf("x=%f\ty=%f\n",center[ic][0],center[ic][1] );
      }
      estimate_para(a0,center,Nstar,Ng0,0,paras);
      rg=paras[0]*2;printf("rg=%f,rd=%f,beta=%f\n",rg,paras[1],paras[2] );

      for(ic=0;ic<Nstar;ic++){
        star_gaus(stars[ic],Ng0,Ng0,starf,residu,paras);
        center[ic][0]=paras[3];center[ic][1]=paras[4];
      }
      strcpy(ftmp,"data/cdata/");
      strcat(ftmp,tmp);strcat(ftmp,"e.dat");printf("%s\n",ftmp );
      fp=fopen(ftmp,"w");
      for(ic=0;ic<Nstar;ic++){
        size(stars[ic],Ng0,center[ic],rg,paras);
        fprintf(fp, "%e\t%e\t%e\t%e\t%e\n",paras[0],paras[1],paras[2],spos[ic][0],spos[ic][1] );
      }fclose(fp);*/
      

      void method_branches(float ***stars,float **spos,int npc,int Nstar,int Ng0,float ***PCs,
        float ***coeff,float gain,int method,int *Nobj,float *means);
      int npc=10;
      float ***PCs,***coeff,gain=keyva[4][0],means[Nstar0];
      PCs=f3tensor(0,npc,0,Ng0,0,Ng0);
      coeff=f3tensor(0,Nstar0,0,npc,0,2);
      method_branches(stars,spos,npc,Nstar,Ng,PCs,coeff,gain,3,&Nobj,means);
      strcpy(ftmp,"fits/PCs.fits");
      dim[0]=dim[1]=Ng,dim[2]=npc;
      write_fits_3D(ftmp,PCs,dim);printf("Nobj=%d\n",Nobj );

      ///////////////////////////////////////////////////////////////////////
      //read clean images
      strcpy(fname_image,"data/TEST150/MSC_0000009/MSC_MS_210525145000_100000009_");
      strcat(fname_image,tmp);strcat(fname_image,"_raw.fits");
      strcat(fname_image,"[raw]");
      read_fits_2D(fname_image,image0,imgsize);
      for(j=0;j<nimgy;j++){
        for(i=0;i<nimgx;i++){k=j*nimgx+i;image[i][j]=image0[k];}
      }
      
      Nh=Ng/2;k=0;
      for(ic=0;ic<Nobj;ic++){
        rx=spos[ic][0]-1.5;ry=spos[ic][1]-1.5;sum=0;
        for(i=0;i<Ng;i++)for(j=0;j<Ng;j++){
          it=rx+i-Nh;jt=ry+j-Nh;
          if(it>=0 && it<nimgx && jt>=0 && jt<nimgy){
            stars[ic][i][j]=image[it][jt]-means[ic];
            //sum+=stars[ic][i][j];
          }
        }
        for(i=0;i<Ng;i++)for(j=0;j<Ng;j++){
          //stars[ic][i][j]/=sum;
          a0[ic][i*Ng+j]=stars[ic][i][j];
        }
      }
      strcpy(ftmp,"fits/ostar.fits");
      dim[0]=dim[1]=Ng,dim[2]=Nobj;
      write_fits_3D(ftmp,stars,dim);
      for(ic=0;ic<Nobj;ic++){
        star_gaus(stars[ic],Ng,Ng,starf,residu,oparas);
        center[ic][0]=oparas[3];center[ic][1]=oparas[4];
      }
      estimate_para(a0,center,Nstar,Ng0,0,para);rg=para[0]*1.2;printf("rg=%f\n",rg );
      strcpy(ftmp,"data/TEST150/");
      strcat(ftmp,tmp);strcat(ftmp,"oe.dat");printf("%s\n",ftmp );//original stars
      fp=fopen(ftmp,"w");
      for(ic=0;ic<Nobj;ic++){
        size(stars[ic],Ng,center[ic],rg,oparas);
        fprintf(fp, "%e\t%e\t%e\t%e\t%e\n",oparas[0],oparas[1],oparas[2],spos[ic][0],spos[ic][1] );
      }fclose(fp);
      ///////////////////////////////////////////////////////////////////////

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
        //printf("x=%f\ty=%f\n",center[ic][0],center[ic][1] );
      }
      strcpy(ftmp,"data/TEST150/");
      strcat(ftmp,tmp);strcat(ftmp,"re.dat");printf("%s\n",ftmp );//reconstructed stars
      fp=fopen(ftmp,"w");
      for(ic=0;ic<Nobj;ic++){
        size(start[ic],Ng,center[ic],rg,paras);
        fprintf(fp, "%e\t%e\t%e\t%e\t%e\n",paras[0],paras[1],paras[2],spos[ic][0],spos[ic][1] );
      }fclose(fp);
      strcpy(ftmp,"fits/restar.fits");
      dim[0]=dim[1]=Ng,dim[2]=Nobj;
      write_fits_3D(ftmp,start,dim);  
 
      void coeff_interp(int order_PC,int order_poly,int Nstar,float **spos,float **Coeff,
        float **weight,float ***PCs,int ngp,float ***rPSFs,int Nobj,float **psfpos);
      float **Coeff,**weight,npsf;
      float **psfpos,***rPSFs;
      int order_poly=3;
      double imean,isigma;
      for(i=order_poly+1;i>0;i--)l0+=i;
      Coeff=matrix(0,Nstar0,0,npc);
      weight=matrix(0,Nstar0,0,npc);
      psfpos=matrix(0,Nstar0,0,2);
      rPSFs=f3tensor(0,Nstar0,0,Ng,0,Ng);


      strcpy(ftmp,"results/coeffs.dat");
      fp=fopen(ftmp,"w");
      for(ic=0;ic<Nobj;ic++){
        for(l=0;l<npc;l++){
          Coeff[ic][l]=coeff[ic][l][0];
          weight[ic][l]=1./coeff[ic][l][1];
          fprintf(fp, "%f\t",Coeff[ic][l] );
        }fprintf(fp, "\n" );
      }fclose(fp);

      coeff_interp(npc,order_poly,Nobj,spos,Coeff,weight,PCs,Ng,rPSFs,Nobj,spos);
      //Using constructed polynomial to interpolate PSFs at stars,\
       which are used to reconstructions
      printf("coeff_interp down\n");

      for(ic=0;ic<Nobj;ic++){
        star_gaus(rPSFs[ic],Ng,Ng,starf,residu,paras);
        center[ic][0]=paras[3];center[ic][1]=paras[4];
      }
      strcpy(ftmp,"data/TEST150/");
      strcat(ftmp,tmp);strcat(ftmp,"ioe.dat");printf("%s\n",ftmp );
      fp=fopen(ftmp,"w");
      for(ic=0;ic<Nobj;ic++){
        size(rPSFs[ic],Ng,center[ic],rg,paras);
        fprintf(fp, "%e\t%e\t%e\t%e\t%e\n",paras[0],paras[1],paras[2],spos[ic][0],spos[ic][1] );
      }fclose(fp);
      printf("calculated\n");
      free_matrix(mpara,0,Nstar0,0,colum);
      free_matrix(psfpos,0,Nstar0,0,2);
      free_f3tensor(stars,0,Nstar0,0,Ng0,0,Ng0);
      free_f3tensor(start,0,Nstar0,0,Ng0,0,Ng0);
      free_dmatrix(a0,0,Nstar0,0,Ng0*Ng0);
      free_dmatrix(center,0,Nstar0,0,2);
      free_f3tensor(rPSFs,0,Nstar0,0,Ng,0,Ng);






      printf("start interp,Nstari=%d\n",Nstari);rg=rg/1.2*1.5;
      ///////////////////////////////////////////////////////////////////////
      strcpy(psf_catalogue,"data/TEST300/200/catalogue/star/MSC_MS_210525150000_100000010_");
      strcat(psf_catalogue,tmp);strcat(psf_catalogue,"_raw.fits.cat");
      printf("%s\n",psf_catalogue);
      //strcpy(psf_catalogue,"data/MSC_210304093000_0000010_07.cat");
      funcpsfs_size(psf_catalogue,&Nstari,&Nhead);//the catalogue should be arranged as "id x y SNR"
      printf("Nstar0=%d,Nhead=%d\n",Nstari,Nhead );colum=Nhead;
      mpara=matrix(0,Nstari,0,colum);
      psfpos=matrix(0,Nstari,0,2);
      stars=f3tensor(0,Nstari,0,Ng0,0,Ng0);
      start=f3tensor(0,Nstari,0,Ng0,0,Ng0);
      center=dmatrix(0,Nstari,0,2);
      rPSFs=f3tensor(0,Nstari,0,Ng,0,Ng);
      printf("colum=%d\n",colum );
      strcpy(delim,"\t");
      funcpsfs_read(psf_catalogue,Nstari,Nhead,mpara,colum,delim);//Nstar0=250;
      Nstari=Nobj_deli;
      for(ic=0;ic<Nstari;ic++){
        //psfpos[ic][0]=mpara[ic][1]+1;
        //psfpos[ic][1]=mpara[ic][2]+1;
        psfpos[ic][0]=spos[ic+Nstar][0];
        psfpos[ic][1]=spos[ic+Nstar][1];
        //printf("%e\t%e\tic=%d\n",mpara[ic][1],mpara[ic][2],ic );
      }
      strcpy(fname_image,"data/TEST300/MSC_0000010/MSC_MS_210525150000_100000010_");
      strcat(fname_image,tmp);strcat(fname_image,"_raw.fits");
      strcat(fname_image,"[raw]");
      read_fits_2D(fname_image,image0,imgsize);
      for(j=0;j<nimgy;j++){
        for(i=0;i<nimgx;i++){k=j*nimgx+i;image[i][j]=image0[k];}
      }
      printf("start read\n");
      Nh=Ng/2;k=0;
      for(ic=0;ic<Nstari;ic++){
        rx=psfpos[ic][0]-1.5;ry=psfpos[ic][1]-1.5;sum=0;
        for(i=0;i<Ng;i++)for(j=0;j<Ng;j++){
          it=rx+i-Nh;jt=ry+j-Nh;
          if(it>=0 && it<nimgx && jt>=0 && jt<nimgy){
            stars[ic][i][j]=image[it][jt];//-means[ic];
            //sum+=stars[ic][i][j];
          }
        }
        //gaus_estimate(stars[ic],Ng,Ng,&imean,&isigma);printf("ic=%d,means=%f\n",ic,imean );
      }
      strcpy(ftmp,"fits/oistar.fits");
      dim[0]=dim[1]=Ng,dim[2]=Nstari;
      write_fits_3D(ftmp,stars,dim);

      for(ic=0;ic<Nstari;ic++){
        gaus_estimate(stars[ic],Ng,Ng,&imean,&isigma);
        //printf("ic=%d,means=%f\n",ic,imean );
        for(i=0;i<Ng;i++)for(j=0;j<Ng;j++)stars[ic][i][j]-=imean;
        star_gaus(stars[ic],Ng,Ng,starf,residu,oparas);
        center[ic][0]=oparas[3];center[ic][1]=oparas[4];
      }
      strcpy(ftmp,"data/TEST150/");
      strcat(ftmp,tmp);strcat(ftmp,"oie.dat");printf("%s\n",ftmp );
      fp=fopen(ftmp,"w");
      for(ic=0;ic<Nstari;ic++){
        size(stars[ic],Ng,center[ic],rg,oparas);
        fprintf(fp, "%e\t%e\t%e\t%e\t%e\n",oparas[0],oparas[1],oparas[2],psfpos[ic][0],psfpos[ic][1] );
      }fclose(fp);
      printf("read down\n");


      printf("start coeff_interp\n");
      coeff_interp(npc,order_poly,Nobj,spos,Coeff,weight,PCs,Ng,rPSFs,Nstari,psfpos);
      printf("coeff_interp down\n");
      strcpy(ftmp,"fits/rPSFs.fits");
      dim[0]=dim[1]=Ng,dim[2]=Nstari;
      write_fits_3D(ftmp,rPSFs,dim);

      for(ic=0;ic<Nstari;ic++){
        star_gaus(rPSFs[ic],Ng,Ng,starf,residu,paras);
        center[ic][0]=paras[3];center[ic][1]=paras[4];
      }
      strcpy(ftmp,"data/TEST150/");
      strcat(ftmp,tmp);strcat(ftmp,"ie.dat");printf("%s\n",ftmp );
      fp=fopen(ftmp,"w");
      for(ic=0;ic<Nstari;ic++){
        size(rPSFs[ic],Ng,center[ic],rg,paras);
        fprintf(fp, "%e\t%e\t%e\t%e\t%e\n",paras[0],paras[1],paras[2],psfpos[ic][0],psfpos[ic][1] );
      }fclose(fp);
      printf("calculated\n");
      ///////////////////////////////////////////////////////////////////////



      

}
