#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"
int main()
{
  int get_fits_size(char *argv,int *dim);
  int write_fits_2D(char *argv,float **stars,int *dim);
  int write_fits_3D(char *argv,float ***stars,int *dim);
  int read_fits_2D(char *argv,float *galaxy,int imagesize);
  int Ntotal,Nstar,i,j,ix,iy,m,n,k,Npix=11,nt,nimgx,nimgy,imgsize,Nh,dim[3];
  float **pos,**postar,*FWHM,*R50,*RMS,*RMSstar,*flux,class;
  float minx=0,miny=0,maxx=100,maxy=100;
  int flag;
  float *image0,**image,***Stars;
  char *fname0="check.fits ";
  char *fnames="star/stars.fit";
  FILE *fp;
  get_fits_size(fname0,dim);
  nimgx=dim[0];nimgy=dim[1];imgsize=nimgx*nimgy;
  image0=vector(0,imgsize-1);
  image=matrix(0,nimgx-1,0,nimgy-1);
  read_fits_2D(fname0,image0,imgsize);
  for(j=0;j<nimgy;j++){
    for(i=0;i<nimgx;i++){k=j*nimgx+i;image[i][j]=image0[k];}
  }
  free_vector(image0,0,imgsize-1);
  fp=fopen("listtotal.d","r");
  fscanf(fp,"%d\n",&Ntotal);
  pos=matrix(0,Ntotal-1,0,1);
  postar=matrix(0,Ntotal-1,0,1);
  FWHM=vector(0,Ntotal-1);
  RMS=vector(0,Ntotal-1);
  RMSstar=vector(0,Ntotal-1);
  R50=vector(0,Ntotal-1);
  flux=vector(0,Ntotal-1);
  for(i=0;i<Ntotal;i++){
     fscanf(fp,"%d %e %e %e %e %e %e %d %e\n",&j,&pos[i][0],&pos[i][1],
	    &FWHM[i],&R50[i],&flux[i],&RMS[i],&flag,&class);
     //printf("%d %e %e %e %e %e %e %d %e\n",j,pos[i][0],pos[i][1],
     //	    FWHM[i],R50[i],flux[i],RMS[i],flag,class);
  }
  fclose(fp);
  for(nt=0;nt<Ntotal;nt++){
    Npix=R50[nt]*3+1;
    Npix=Npix*2+1;    
    Stars=f3tensor(0,1,0,Npix-1,0,Npix-1);
    Nh=Npix*0.5;
    ix=pos[nt][0]+0.5;iy=pos[nt][1]+0.5;
    ix--;iy--;
    printf("%d %e %e\n",nt,FWHM[nt],R50[nt]);
    for(j=0;j<Npix;j++){
      for(i=0;i<Npix;i++){
        m=ix+i-Nh;
        n=iy+j-Nh;
	if(m>=0 & m<nimgx & n>=0 & n<nimgx)
        Stars[0][i][j]=image[m][n];
	else{Stars[0][i][j]=0;}
	
      }
    }
    if(Npix>100){
    dim[0]=dim[1]=Npix;
    write_fits_2D(fnames,Stars[0],dim);   
    //free_f3tensor(Stars,0,Nstar-1,0,Npix-1,0,Npix-1);
    getchar();
    } 
    //free_f3tensor(Stars,0,1,0,Npix-1,0,Npix-1);
  }


  free_matrix(pos,0,Ntotal-1,0,1);
  free_matrix(postar,0,Ntotal-1,0,1);
  free_vector(FWHM,0,Ntotal-1);
  free_vector(RMS,0,Ntotal-1);
  free_vector(RMSstar,0,Ntotal-1);  
  free_matrix(image,0,nimgx-1,0,nimgy-1);

}
