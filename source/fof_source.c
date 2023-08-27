#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"
#include <string.h>
#include "dirent.h" 
#define max( a, b) ( (a) > (b) ? (a) : (b) )
#define min( a, b) ( (a) < (b) ? (a) : (b) )
#define Pi 3.14159256
void fof_image(float **image,double sigma,int nimgx,int nimgy,int *groupnum,
  int *tails,int **pgroup){
  int write_fits_2D(char *argv,float **stars,int *dim);
  int write_fits_3D(char *argv,float ***stars,int *dim);
  int Nstar,Nhead,ic,i,j,k,Nb,m,n,Nobj,f;
  int ra,rb,ix,iy;
  int *isop,Ncount,*tmps,*tgroup,*tcount;
  int pin,tail,ngroup;
  double r,sum,sumx,sumy,centx,centy,r0,r2,peak;
  float **tmpimage;
  isop=ivector(0,nimgx*nimgy);
  tcount=ivector(0,nimgx*nimgy);
  tgroup=ivector(0,nimgx*nimgy);
  tmps=ivector(0,nimgx*nimgy);
  tmpimage=matrix(0,nimgx,0,nimgy);
  

  printf("sigma=%f\n",sigma*3. );
  for(i=0;i<nimgx;i++)for(j=0;j<nimgy;j++){
    if(image[i][j]<(2.*sigma))tmpimage[i][j]=-1;
    else{tmpimage[i][j]=image[i][j];}
  }
  f=0;
  for(i=0;i<nimgx;i++)for(j=0;j<nimgy;j++){
    k=0;
    if(tmpimage[i][j]>0)for(m=-1;m<=1;m++)for(n=-1;n<=1;n++){
      if(i+m>=0 && i+m<nimgx && j+n>=0 &&j+n<nimgy){
        if(tmpimage[i+m][j+n]>0){k++;}
      }
    }if(k<2){
      if(tmpimage[i][j]<15*sigma)tmpimage[i][j]=-1;//  filte isolated pixels
      else{image[i][j]=0;}
    }else{
      isop[f]=i*nimgy+j;f++;
    }
  }Ncount=f;printf("f=%d\n",f );

  char fname[200];
  int dim[4];
  strcpy(fname,"fits/test.fits");
  dim[0]=nimgx;dim[1]=nimgy;
  write_fits_2D(fname,tmpimage,dim);

  for(i=0;i<Ncount;i++)tcount[i]=isop[i];
  ngroup=0;
  if(Ncount>=5){
    for(i=0;i<Ncount;i++)tcount[i]=isop[i];
    //ngroup=0;
    peak=0;

    //printf("Ncount=%d\n",Ncount );
    do{
      pin=tail=0;
      tgroup[pin]=tcount[0];//for the initial base pixel
      ra=tgroup[pin]/nimgy;rb=tgroup[pin]%nimgy;//construct the pixel position
      tail++;
      f=0;
      for(i=1;i<Ncount;i++){
        ix=tcount[i]/nimgy;iy=tcount[i]%nimgy;
        r=(ra-ix)*(ra-ix)+(rb-iy)*(rb-iy);
        if(r<2.5){//find the connected pixel
          tgroup[tail]=tcount[i];tail++;//add this pixel in the potential fofgroup list
        }else{
          tmps[f]=tcount[i];f++;//the rest pixels
        }
      }//printf("\n");
      Ncount=f;pin++;//shift the base pixel
      do{
        for(i=0;i<Ncount;i++)tcount[i]=tmps[i];//find the rest fof pixels
        ra=tgroup[pin]/nimgy;rb=tgroup[pin]%nimgy;//construct the pixel position
        f=0;
        for(i=0;i<Ncount;i++){
          ix=tcount[i]/nimgy;iy=tcount[i]%nimgy;
          r=(ra-ix)*(ra-ix)+(rb-iy)*(rb-iy);
          if(r<2.5){
            tgroup[tail]=tcount[i];tail++;//printf("ix=%d,iy=%d\t",ix,iy );
          }else{
            tmps[f]=tcount[i];f++;
          }
        }pin++;
        Ncount=f;//printf("Ncount=%d\n",Ncount[ic]);
      }while(pin!=tail);//untill the base pixel go to the last one

      if(tail>=3 && tail <6){ //hot pixel detection
        for(i=0;i<tail;i++){
          ra=tgroup[i]/nimgy;rb=tgroup[i]%nimgy;
          peak=max(peak,tmpimage[ra][rb]);
        }
        if(peak>15*sigma){
          for(i=0;i<tail;i++)pgroup[ngroup][i]=tgroup[i];
          tails[ngroup]=tail;
          ngroup++;
        }
      }
    
      //printf("tail=%d,Ncount left=%d\n",tail,Ncount );
      f=0;
      if(tail>=5){//the minimum pixel for a group
        for(i=0;i<tail;i++){
          ra=tgroup[i]/nimgy;rb=tgroup[i]%nimgy;
          if(tmpimage[ra][rb]>3.*sigma){ //threashhold selection
            f++;
          }
        }
        if(f>=2){
          for(i=0;i<tail;i++)pgroup[ngroup][i]=tgroup[i];
          tails[ngroup]=tail;//printf("group=%d,tail=%d\n",ngroup,tails[ngroup] );
          ngroup++;
        }
      }
    }while(Ncount>=5);
    *groupnum=ngroup;


    float ***tstar;
    tstar=f3tensor(0,ngroup,0,nimgx,0,nimgy);
    for(ic=0;ic<ngroup;ic++){
      for(i=0;i<nimgx;i++)for(j=0;j<nimgy;j++)tstar[ic][i][j]=0;
      for(i=0;i<tails[ic];i++){
        ra=pgroup[ic][i]/nimgy;rb=pgroup[ic][i]%nimgy;//printf("x=%d,y=%d\n",ra,rb );
        tstar[ic][ra][rb]=image[ra][rb];
      }
    }
    strcpy(fname,"fits/test1.fits");
    dim[0]=nimgx;dim[1]=nimgy;dim[2]=ngroup;
    write_fits_3D(fname,tstar,dim);
    }
    else{
      *groupnum=0;
      printf("the effective pixel number is small:<10, Quit!\n");
      //exit(EXIT_FAILURE);
    }
  /*int tip=12;//printf("ngroup=\n");
  for(ic=0;ic<ngroups[tip];ic++){
    printf("tails=%d\n",tails[tip][ic] );
    for(i=0;i<Npix;i++)for(j=0;j<Npix;j++)cutstar[ic][i][j]=0;
    for(k=0;k<tails[tip][ic];k++){
      ra=pgroup[tip][ic][k]/Npix;rb=pgroup[tip][ic][k]%Npix;
      cutstar[ic][ra][rb]=wstar[tip][ra][rb];
    }
  }
  strcpy(fnames,"test/groups.fits");
  dim[0]=dim[1]=Npix;dim[2]=ngroups[tip];
  write_fits_3D(fnames,cutstar,dim);*/



  free_ivector(isop,0,nimgx*nimgy);
  free_ivector(tcount,0,nimgx*nimgy);
  free_ivector(tgroup,0,nimgx*nimgy);
  free_ivector(tmps,0,nimgx*nimgy);
  free_matrix(tmpimage,0,nimgx,0,nimgy);

}

void subfind(float **image,float *lthresh,int nthresh,int nimgx,int nimgy,
  int tails,int *pgroup,int **psubs,int *nsub,int *subtail,int Npix,float sigma){
  int Nstar,Nhead,ic,i,j,k,Nb,m,n,Nobj,l,f,dim[5];
  int ra,rb,ix,iy;
  int *isop,*Ncount,*tmps,*tgroup,*tcount,**tpsub,*ttail,**tsubs;
  int pin,tail,ngroup,ntmp,thresh;
  double r,sum,sumx,sumy,centx,centy,r0,r2;
  isop=ivector(0,Npix*Npix);
  tcount=ivector(0,Npix*Npix);
  tgroup=ivector(0,Npix*Npix);
  tmps=ivector(0,Npix*Npix);
  Ncount=ivector(0,10);
  
  tpsub=imatrix(0,10,0,Npix*Npix);
  ttail=ivector(0,10);
  tsubs=imatrix(0,10,0,Npix*Npix);
  
  ntmp=1;thresh=1;Ncount[0]=tails;
  for(k=0;k<tails;k++)tpsub[0][k]=pgroup[k];

  n=1;//number of all source
  do{
  l=0;//number for potential group
  for(ic=0;ic<ntmp;ic++){
    f=0;
    for(k=0;k<Ncount[ic];k++){
      ra=tpsub[ic][k]/nimgy;rb=tpsub[ic][k]%nimgy;//printf("thresh=%f\n",lthresh[thresh] );
      if(image[ra][rb]>lthresh[thresh]){
      	isop[f]=tpsub[ic][k];f++;
      }
    }//printf("f=%d\n",f );
    for(i=0;i<f;i++)tcount[i]=isop[i];
    ngroup=0;
    Ncount[ic]=f;

    //printf("Ncount=%d\n",f );
    do{
      pin=tail=0;
      tgroup[pin]=tcount[0];//for the initial base pixel
      ra=tgroup[pin]/nimgy;rb=tgroup[pin]%nimgy;
      tail++;
      f=0;
      for(i=1;i<Ncount[ic];i++){
        ix=tcount[i]/nimgy;iy=tcount[i]%nimgy;
        r=(ra-ix)*(ra-ix)+(rb-iy)*(rb-iy);
        if(r<2.01){//find the connected pixel
          tgroup[tail]=tcount[i];tail++;//add this pixel in the potential fofgroup list
          //printf("1ra=%d,rb=%d,ix=%d,iy=%d,r=%f,tail=%d\n",ra,rb,ix,iy,r,tail );
        }else{
          tmps[f]=tcount[i];f++;//the rest pixels

        }
      }//printf("else 1ra=%d,rb=%d,ix=%d,iy=%d,r=%f,tail=%d\n",ra,rb,ix,iy,r,tail);
      Ncount[ic]=f;pin++;//shift the base pixel
      for(i=0;i<Ncount[ic];i++)tcount[i]=tmps[i];
      //printf("pin=%d,tail=%d\n",pin,tail );
      while(pin!=tail){
        for(i=0;i<Ncount[ic];i++)tcount[i]=tmps[i];//find the rest fof pixels
        ra=tgroup[pin]/nimgy;rb=tgroup[pin]%nimgy;
        f=0;
        for(i=0;i<Ncount[ic];i++){
          ix=tcount[i]/nimgy;iy=tcount[i]%nimgy;
          r=(ra-ix)*(ra-ix)+(rb-iy)*(rb-iy);
          if(r<2.01){
            tgroup[tail]=tcount[i];tail++;
            //printf("2ra=%d,rb=%d,ix=%d,iy=%d,r=%f,tail=%d\n",ra,rb,ix,iy,r,tail );
          }else{
            tmps[f]=tcount[i];f++;
          }
        }pin++;
        Ncount[ic]=f;
      }//untill the base pixel go to the last one
      for(i=0;i<Ncount[ic];i++)tcount[i]=tmps[i];

      if(tail>=7){//the minimum pixel for a group
        for(i=0;i<tail;i++)tsubs[ngroup][i]=tgroup[i];
        ttail[ngroup]=tail;
        ngroup++;
      }
      else{
        if(tail>=4 && lthresh[thresh]>10.*sigma){
          for(i=0;i<tail;i++)tsubs[ngroup][i]=tgroup[i];
          ttail[ngroup]=tail;
          ngroup++;
        }
      }
      //printf("Ncount=%d\n",Ncount[ic] );
    }while(Ncount[ic]>=4);
    //printf("n=%d\n",n );
    if(ngroup>0)n+=ngroup-1;//if(ngroup-1!=0)printf("lthresh=%f,subtail=%d\n",lthresh[thresh],tail );
    //printf("n=%d\n",n );
    for(i=0;i<ngroup;i++){
      if(ttail[i]>=7){
      	for(k=0;k<ttail[i];k++)tpsub[l][k]=tsubs[i][k];
      	Ncount[l]=ttail[i];
      	l++;
      }
    }

  }
  ntmp=l;thresh++;//printf("l=%d,thresh=%d\n",l,thresh );
  }while(l!=0);
  *nsub=n;
  //printf("*n=%d\n",n );



  free_ivector(isop,0,Npix*Npix);
  free_ivector(tcount,0,Npix*Npix);
  free_ivector(tgroup,0,Npix*Npix);
  free_ivector(tmps,0,Npix*Npix);
  free_ivector(Ncount,0,10);
  //free_imatrix(psubs,0,10,0,Npix*Npix);
  free_imatrix(tpsub,0,10,0,Npix*Npix);
  free_ivector(ttail,0,10);
  free_imatrix(tsubs,0,10,0,Npix*Npix);

}

void deblend(float **image,int nimgx,int nimgy,int ngroups,int *tails,
	int **pgroup,int *nsub,float sigma){
  void subfind(float **image,float *lthresh,int nthresh,int nimgx,int nimgy,
    int tails,int *pgroup,int **psubs,int *nsub,int *subtail,int Npix,float sigma);
  void sort(int n, float arr[]);
  void sort2(unsigned long n, float arr[], float brr[]);
  int write_fits_2D(char *argv,float **stars,int *dim);
  int write_fits_3D(char *argv,float ***stars,int *dim);

  int i,j,k,l,f,m,n,ra,rb,subcount,subs;
  double sigs,sigs1,peak,theta0,theta,dh,rx,ry,rleast,di,dc,rmin;
  int u,v,pleast[2],nbigfof;
  int subpix,**psubs,*subtail,Npix=nimgy;
  float **tmpimage,**lthresh,**rlthresh,*rlist,fr;
  int pin,tail,ngroup,s,*nthresh,**pedge,*nedge,**pthresh;
  double r,sum,sumx,sumy,centx,centy,r0,r2,mean;
  char tmp[200];
  int dim[4];
  long int countf;

  tmpimage=matrix(0,Npix,0,Npix);
  nthresh=ivector(0,10);
  lthresh=matrix(0,10,1,Npix);
  rlthresh=matrix(0,10,1,Npix);
  pedge=imatrix(0,10,0,Npix*2);
  nedge=ivector(0,10);
  pthresh=imatrix(0,10,0,Npix*2);//assum only 10 groups
  rlist=vector(1,Npix*Npix);
  subtail=ivector(0,10);
  psubs=imatrix(0,10,0,Npix*Npix);
    l=0;subcount=0;
    for(i=0;i<ngroups;i++){
      peak=0;sum=0;sumx=0;sumy=0;
      for(m=0;m<Npix;m++)for(n=0;n<Npix;n++)tmpimage[m][n]=-1;
      //printf("ngroups=%d,tails=%d,ic=%d\n",ngroups[ic],tails[ic][i],ic );
      if(tails[i]>25){
        for(k=0;k<tails[i];k++){
          ra=pgroup[i][k]/Npix;rb=pgroup[i][k]%Npix;
          tmpimage[ra][rb]=image[ra][rb];
          peak=max(peak,image[ra][rb]);
          if(image[ra][rb]==peak){//peak detection
            centx=ra;centy=rb;
          }
        }//printf("x=%f,y=%f\n",centx,centy);
        //center pick
        for(k=0;k<tails[i];k++){
          ra=pgroup[i][k]/Npix;rb=pgroup[i][k]%Npix;
          r=(ra-centx)*(ra-centx)+(rb-centy)*(rb-centy);
          sum+=image[ra][rb]*exp(-r/9/2);
          sumx+=image[ra][rb]*ra*exp(-r/9/2);
          sumy+=image[ra][rb]*rb*exp(-r/9/2);
        }
        centx=sumx/sum;centy=sumy/sum;//printf("x=%f,y=%f,ic=%d\n",centx,centy,ic+1 );

        pin=0;
        for(k=0;k<tails[i];k++){
          m=pgroup[i][k]/Npix;n=pgroup[i][k]%Npix;
          s=0;
          for(u=-1;u<=1;u++)for(v=-1;v<=1;v++){
            if(m+u>=0 && m+u<nimgx && n+v>=0 &&n+v<nimgy){
              if(tmpimage[m+u][n+v]>0){s++;}
            }
          }if(s<9){//edge detection
            pedge[l][pin]=pgroup[i][k];pin++;
          }
        }nedge[l]=pin;//printf("pin=%d\n",pin );

        sum=2.*Npix*Npix;sumx=0;
        for(k=0;k<nedge[l];k++){
          ra=pedge[l][k]/Npix;rb=pedge[l][k]%Npix;
          rx=ra-centx;ry=rb-centy;
          r=rx*rx+ry*ry;
          rlist[k+1]=r;
        }sort(nedge[l],rlist);sum=rlist[k-k/4];//printf("sum=%f\n",sum );
        //for(m=0;m<k;m++){printf("thresh=%f\n",rlist[m] );}
        //find a proper point that not the longest,which could cross the subpeak(more efficiency at most)
        //and not the shortest,which may miss the subpeak
        for(k=0;k<nedge[l];k++){
          ra=pedge[l][k]/Npix;rb=pedge[l][k]%Npix;
          rx=ra-centx;ry=rb-centy;
          fr=rx*rx+ry*ry;//printf("r1=%f\n",r );
          //sum=min(sum,r);
          sumx=max(sumx,r);
          if(fr==sum){//pick the minimum pixel to the center
            pleast[0]=ra;pleast[1]=rb;theta=atan(rx/ry);
            if(ry<0){theta+=Pi;}
            //printf("ra=%d,rb=%d\n",ra,rb );
          }
        }
        rmin=sum;subpix=(int)sumx+1;//find the longest distance, which could contain all the pixels
        //printf("ra1=%d,rb1=%d\n",pleast[0],pleast[1] );
        f=0;
        for(k=0;k<tails[i];k++){
          ra=pgroup[i][k]/Npix;rb=pgroup[i][k]%Npix;
          rx=ra-centx;ry=rb-centy;
          r=rx*rx+ry*ry;
          di=(ra-pleast[0])*(ra-pleast[0])+(rb-pleast[1])*(rb-pleast[1]);
          theta0=atan(rx/ry)-theta;if(ry<0)theta0+=Pi;
          dh=r*sin(theta0)*sin(theta0);
          if(dh<.5 && di<rmin && r>0 && r<=(rmin)){//dh<.5: distance to the line; di<rmin:distance to the edge point
            //r>0: avoid the peak point; r<=(rmin): constrain the areas
            pthresh[l][f]=pgroup[i][k];f++;
            lthresh[l][f]=image[ra][rb];
          }
        }nthresh[l]=f;countf=f;//printf("f=%d\n",f );
        sort(countf,lthresh[l]);//for(k=0;k<f;k++)printf("thresh=%f\n",lthresh[l][k] );
        subfind(image,lthresh[l],nthresh[l],Npix,Npix,tails[i],
          pgroup[i],psubs,&subs,subtail,subpix,sigma);
        //if(nsub>1)printf("nsub=%d\n",nsub );
        l++;subcount+=subs;
      }
    }nbigfof=l;//printf("nbigfof=%d\n",nbigfof[ic] );
   
   *nsub=subcount;

   for(i=0;i<Npix;i++)for(j=0;j<Npix;j++)tmpimage[i][j]=0;
   for(k=0;k<nthresh[0];k++){
     ra=pthresh[0][k]/Npix;rb=pthresh[0][k]%Npix;
     tmpimage[ra][rb]=image[ra][rb];
   }
   strcpy(tmp,"fits/thresh.fits");
   dim[0]=dim[1]=Npix;
   //write_fits_2D(tmp,tmpimage,dim);
   //for(i=0;i<Npix;i++)for(j=0;j<Npix;j++)tmpimage[i][j]=-1;
   	//printf("subtail=%d\n",subs);
   for(k=0;k<nedge[0];k++){
   	 ra=pedge[0][k]/Npix;rb=pedge[0][k]%Npix;
   	 tmpimage[ra][rb]=image[ra][rb];
   }
   strcpy(tmp,"fits/subs.fits");
   dim[0]=dim[1]=Npix;
   //write_fits_2D(tmp,tmpimage,dim);

  free_matrix(tmpimage,0,Npix,0,Npix);
  free_ivector(nthresh,0,10);
  free_matrix(lthresh,0,10,1,Npix);
  free_matrix(rlthresh,0,10,1,Npix);
  free_imatrix(pedge,0,10,0,Npix*2);
  free_ivector(nedge,0,10);
  free_imatrix(pthresh,0,10,0,Npix*2);
  free_vector(rlist,1,Npix*Npix);
  free_ivector(subtail,0,10);
  free_imatrix(psubs,0,10,0,Npix*Npix);

}



void deblend1(float **image,int nimgx,int nimgy,int ngroups,int *tails,
  int **pgroup,int *nsub,float sigma){
  void subfind(float **image,float *lthresh,int nthresh,int nimgx,int nimgy,
    int tails,int *pgroup,int **psubs,int *nsub,int *subtail,int Npix,float sigma);
  void sort(int n, float arr[]);
  void sort2(unsigned long n, float arr[], float brr[]);
  int write_fits_2D(char *argv,float **stars,int *dim);
  int write_fits_3D(char *argv,float ***stars,int *dim);

  int i,j,k,l,f,m,n,ra,rb,subcount,subs;
  double sigs,sigs1,peak,theta0,theta,dh,rx,ry,rleast,di,dc,rmin;
  int u,v,pleast[2],nbigfof;
  int subpix,**psubs,*subtail,Npix=nimgy;
  float **tmpimage,**lthresh,**rlthresh,*rlist,fr;
  int pin,tail,ngroup,s,*nthresh,**pedge,*nedge,**pthresh;
  double r,sum,sumx,sumy,centx,centy,r0,r2,mean;
  char tmp[200];
  int dim[4];
  long int countf;

  tmpimage=matrix(0,Npix,0,Npix);
  nthresh=ivector(0,10);
  lthresh=matrix(0,10,1,Npix);
  rlthresh=matrix(0,10,1,Npix);
  pedge=imatrix(0,10,0,Npix*2);
  nedge=ivector(0,10);
  pthresh=imatrix(0,10,0,Npix*2);//assum only 10 groups
  rlist=vector(1,Npix*Npix);
  subtail=ivector(0,10);
  psubs=imatrix(0,10,0,Npix*Npix);
    l=0;subcount=0;
    for(i=0;i<ngroups;i++){
      peak=0;sum=0;sumx=0;sumy=0;
      for(m=0;m<Npix;m++)for(n=0;n<Npix;n++)tmpimage[m][n]=-1;
      //printf("ngroups=%d,tails=%d,ic=%d\n",ngroups[ic],tails[ic][i],ic );
      if(tails[i]>25){
        for(k=0;k<tails[i];k++){
          ra=pgroup[i][k]/Npix;rb=pgroup[i][k]%Npix;
          tmpimage[ra][rb]=image[ra][rb];
          peak=max(peak,image[ra][rb]);
          if(image[ra][rb]==peak){//peak detection
            centx=ra;centy=rb;
          }
        }//printf("x=%f,y=%f\n",centx,centy);
        //center pick
        for(k=0;k<tails[i];k++){
          ra=pgroup[i][k]/Npix;rb=pgroup[i][k]%Npix;
          r=(ra-centx)*(ra-centx)+(rb-centy)*(rb-centy);
          sum+=image[ra][rb]*exp(-r/9/2);
          sumx+=image[ra][rb]*ra*exp(-r/9/2);
          sumy+=image[ra][rb]*rb*exp(-r/9/2);
        }
        centx=sumx/sum;centy=sumy/sum;//printf("x=%f,y=%f,ic=%d\n",centx,centy,ic+1 );

        pin=0;
        for(k=0;k<tails[i];k++){
          m=pgroup[i][k]/Npix;n=pgroup[i][k]%Npix;
          s=0;
          for(u=-1;u<=1;u++)for(v=-1;v<=1;v++){
            if(m+u>=0 && m+u<nimgx && n+v>=0 &&n+v<nimgy){
              if(tmpimage[m+u][n+v]>0){s++;}
            }
          }if(s<9){//edge detection
            pedge[l][pin]=pgroup[i][k];pin++;
          }
        }nedge[l]=pin;//printf("pin=%d\n",pin );

        sum=2.*Npix*Npix;sumx=0;
        for(k=0;k<nedge[l];k++){
          ra=pedge[l][k]/Npix;rb=pedge[l][k]%Npix;
          rx=ra-centx;ry=rb-centy;
          r=rx*rx+ry*ry;
          rlist[k+1]=r;
        }sort(nedge[l],rlist);sum=rlist[k-1];//printf("sum=%f\n",sum );
        //for(m=0;m<k;m++){printf("thresh=%f\n",rlist[m] );}
        //find a proper point that not the longest,which could cross the subpeak(more efficiency at most)
        //and not the shortest,which may miss the subpeak
        for(k=0;k<nedge[l];k++){
          ra=pedge[l][k]/Npix;rb=pedge[l][k]%Npix;
          rx=ra-centx;ry=rb-centy;
          fr=rx*rx+ry*ry;//printf("r1=%f\n",r );
          //sum=min(sum,r);
          sumx=max(sumx,r);
          if(fr==sum){//pick the minimum pixel to the center
            pleast[0]=ra;pleast[1]=rb;theta=atan(rx/ry);
            if(ry<0){theta+=Pi;}
            //printf("ra=%d,rb=%d\n",ra,rb );
          }
        }
        rmin=sum;subpix=(int)sumx+1;//find the longest distance, which could contain all the pixels
        //printf("ra1=%d,rb1=%d\n",pleast[0],pleast[1] );
        f=0;
        for(k=0;k<tails[i];k++){
          ra=pgroup[i][k]/Npix;rb=pgroup[i][k]%Npix;
          rx=ra-centx;ry=rb-centy;
          r=rx*rx+ry*ry;
          di=(ra-pleast[0])*(ra-pleast[0])+(rb-pleast[1])*(rb-pleast[1]);
          theta0=atan(rx/ry)-theta;if(ry<0)theta0+=Pi;
          dh=r*sin(theta0)*sin(theta0);
          if(dh<.5 && di<rmin && r>0 && r<=(rmin)){//dh<.5: distance to the line; di<rmin:distance to the edge point
            //r>0: avoid the peak point; r<=(rmin): constrain the areas
            pthresh[l][f]=pgroup[i][k];f++;
            lthresh[l][f]=image[ra][rb];
          }
        }nthresh[l]=f;countf=f;//printf("f=%d\n",f );
        sort(countf,lthresh[l]);//for(k=0;k<f;k++)printf("thresh=%f\n",lthresh[l][k] );
        subfind(image,lthresh[l],nthresh[l],Npix,Npix,tails[i],
          pgroup[i],psubs,&subs,subtail,subpix,sigma);
        //if(nsub>1)printf("nsub=%d\n",nsub );
        l++;subcount+=subs;
      }
    }nbigfof=l;//printf("nbigfof=%d\n",nbigfof[ic] );
   
   *nsub=subcount;

   for(i=0;i<Npix;i++)for(j=0;j<Npix;j++)tmpimage[i][j]=0;
   for(k=0;k<nthresh[0];k++){
     ra=pthresh[0][k]/Npix;rb=pthresh[0][k]%Npix;
     tmpimage[ra][rb]=image[ra][rb];
   }
   strcpy(tmp,"test/thresh.fits");
   dim[0]=dim[1]=Npix;
   //write_fits_2D(tmp,tmpimage,dim);
   //for(i=0;i<Npix;i++)for(j=0;j<Npix;j++)tmpimage[i][j]=-1;
    //printf("subtail=%d\n",subs);
   /*for(k=0;k<nedge[0];k++){
     ra=pedge[0][k]/Npix;rb=pedge[0][k]%Npix;
     tmpimage[ra][rb]=image[ra][rb];
   }
   strcpy(tmp,"test/subs.fits");
   dim[0]=dim[1]=Npix;
   write_fits_2D(tmp,tmpimage,dim);*/

  free_matrix(tmpimage,0,Npix,0,Npix);
  free_ivector(nthresh,0,10);
  free_matrix(lthresh,0,10,1,Npix);
  free_matrix(rlthresh,0,10,1,Npix);
  free_imatrix(pedge,0,10,0,Npix*2);
  free_ivector(nedge,0,10);
  free_imatrix(pthresh,0,10,0,Npix*2);
  free_vector(rlist,1,Npix*Npix);
  free_ivector(subtail,0,10);
  free_imatrix(psubs,0,10,0,Npix*Npix);

}





void sigma_neibor(float **image,int Nstar,int Npix,double **var,float **pos,
  int nimgx,int nimgy){
  int sign[4],ic,i,j,k,ra,rb,nh;
  double mean[4],sigma[4],peak[4],media,sum;
  nh=Npix/2;
  for(ic=0;ic<Nstar;ic++){
    //printf("ic=%d,",ic );
    ra=pos[ic][0];rb=pos[ic][1];
    //printf("x=%d,y=%d\n",ra,rb );if(ra-Npix>0)printf("hhhhh\n");
    if(ra-Npix>0){//printf("here\n");
      sign[0]=1;mean[0]=0;peak[0]=0;sum=0;k=0;
      for(i=0;i<nh;i++)for(j=0;j<Npix;j++){
        if(rb-nh+j>=0 && rb-nh+j<nimgy){
          mean[0]+=image[ra-Npix+i][rb-nh+j];
          peak[0]=max(peak[0],image[ra-Npix+i][rb-nh+j]);k++;
        }
      }mean[0]/=k;
      for(i=0;i<nh;i++)for(j=0;j<Npix;j++){
      	if(rb-nh+j>=0 && rb-nh+j<nimgy){
          sum+=DSQR(image[ra-Npix+i][rb-nh+j]-mean[0]);
        }
      }sum/=(k-1);
      sigma[0]=sqrt(sum);//printf("peak=%f,sigma=%f,mean=%f\n", peak[0],sqrt(sum),mean);
    }else{sign[0]=0;peak[0]=1.0e+10;}
    if(rb+Npix<nimgy){//printf("here1,%d\n",rb+Npix);
      sign[1]=1;mean[1]=0;peak[1]=0;sum=0;k=0;
      for(i=0;i<Npix;i++)for(j=0;j<nh;j++){
      	if(ra-nh+i>=0 && ra-nh+i<nimgx){
          mean[1]+=image[ra-nh+i][rb+nh+j];
          peak[1]=max(peak[1],image[ra-nh+i][rb+nh+j]);k++;
        }
      }mean[1]/=k;
      for(i=0;i<Npix;i++)for(j=0;j<nh;j++){
      	if(ra-nh+i>=0 && ra-nh+i<nimgx){
          sum+=DSQR(image[ra-nh+i][rb+nh+j]-mean[1]);
        }
      }sum/=(k-1);
      sigma[1]=sqrt(sum);//printf("peak=%f,sigma=%f,mean=%f\n", peak[1],sqrt(sum),mean);
    }else{sign[1]=0;peak[1]=1.0e+10;}
    if(ra+Npix<nimgx){//printf("here2\n");
      sign[2]=1;mean[2]=0;peak[2]=0;sum=0;k=0;
      for(i=0;i<nh;i++)for(j=0;j<Npix;j++){
      	if(rb-nh+j>=0 && rb-nh+j<nimgy){
          mean[2]+=image[ra+nh+i][rb-nh+j];
          peak[2]=max(peak[2],image[ra+nh+i][rb-nh+j]);k++;
        }
      }mean[2]/=k;
      for(i=0;i<nh;i++)for(j=0;j<Npix;j++){
      	if(rb-nh+j>=0 && rb-nh+j<nimgy){
          sum+=DSQR(image[ra+nh+i][rb-nh+j]-mean[2]);
        }
      }sum/=(k-1);
      sigma[2]=sqrt(sum);//printf("peak=%f,sigma=%f,mean=%f\n", peak[2],sqrt(sum),mean);
    }else{sign[2]=0;peak[2]=1.0e+10;}
    if(rb-Npix>0){//printf("here3\n");
      sign[3]=1;mean[3]=0;peak[3]=0;sum=0;k=0;
      for(i=0;i<Npix;i++)for(j=0;j<nh;j++){
      	if(ra-nh+i>=0 && ra-nh+i<nimgx){
          mean[3]+=image[ra-nh+i][rb-Npix+j];
          peak[3]=max(peak[3],image[ra-nh+i][rb-Npix+j]);k++;
        }
      }mean[3]/=k;
      for(i=0;i<Npix;i++)for(j=0;j<nh;j++){
      	if(ra-nh+i>=0 && ra-nh+i<nimgx){
          sum+=DSQR(image[ra-nh+i][rb-Npix+j]-mean[3]);
        }
      }sum/=(k-1);
      sigma[3]=sqrt(sum);//printf("peak=%f,sigma=%f,mean=%f\n", peak[3],sqrt(sum),mean);
    }else{sign[3]=0;peak[3]=1.0e+10;}

    media=min(min(peak[0],peak[1]),min(peak[2],peak[3]));
    
    //printf("media=%f\n",media );
    for(k=0;k<4;k++){
      if(peak[k]==media){
        var[ic][0]=sigma[k];
        var[ic][1]=mean[k];
        //printf("peak=%f,ic=%d\n", mean[k],ic);
      }
    }
  }
}


void mask_groups(int count,int ngroup,int *tails,int **pgroup,float **image,
  double sigma,int nimgx,int nimgy,float **mask){
  int ic,i,j,k,ra,rb;
  long iseed;
  float gasdev(long *idum);
  float ***tstar;
  double sumx,sum,cx,cy,sumy,cenx,ceny,csumx,csumy,csum;
  int num=0,dim[5];
  iseed=(long)(-count);//printf("nimgx=%d\n",nimgx );
  tstar=f3tensor(0,ngroup,0,nimgx,0,nimgy);
  for(ic=0;ic<ngroup;ic++){
    sumx=sumy=sum=csumx=csumy=csum=0;
    for(i=0;i<nimgx;i++)for(j=0;j<nimgy;j++)tstar[ic][i][j]=0;
    for(i=0;i<tails[ic];i++){
      ra=pgroup[ic][i]/nimgy;rb=pgroup[ic][i]%nimgy;
      tstar[ic][ra][rb]=image[ra][rb];
      sumx+=1*ra;sumy+=1*rb;sum+=1;
      csumx+=image[ra][rb]*ra;csumy+=image[ra][rb]*rb;csum+=image[ra][rb];

    }
    //printf("csumx=%f,csumy=%f\n",csumx,csumy );
    cx=sumx/sum+1.;cy=sumy/sum+1.;//geometry center
    cenx=csumx/csum+1.;ceny=csumy/csum+1.;//weight center
    printf("sx=%f,sy=%f,cx=%f,cy=%f\n",cx,cy,cenx,ceny );
    if(fabs(cx-nimgx/2.)>=3 || fabs(cy-nimgy/2.)>=3 || 
       fabs(cenx-nimgx/2.)>=3 || fabs(ceny-nimgy/2.)>=3){
      num+=1;//printf("bias\n");
      for(i=0;i<tails[ic];i++){
        ra=pgroup[ic][i]/nimgy;rb=pgroup[ic][i]%nimgy;
        image[ra][rb]=0;//gasdev(&iseed)*sigma;
        mask[ra][rb]=0;
        sumx+=1*ra;sumy+=1*rb;
        sum+=1;
      }
    }
    /*else{
      if(tails[ic]<30){ //exclude the stars that contain too less pixels
        for(i=0;i<tails[ic];i++){
          ra=pgroup[ic][i]/nimgy;rb=pgroup[ic][i]%nimgy;
          image[ra][rb]=0;//gasdev(&iseed)*sigma;
          mask[ra][rb]=0;
          sumx+=1*ra;sumy+=1*rb;
          sum+=1;
        }
      }
    }*/

  }
  /*if(num==ngroup){
    // *flag=0;
  }*/
  int write_fits_3D(char *argv,float ***stars,int *dim);
  char tmp[100];
   strcpy(tmp,"fits/mask.fits");
   dim[0]=dim[1]=nimgx;dim[2]=ngroup;
   //write_fits_3D(tmp,tstar,dim);

}

/*double sigs,sigs1,peak,theta0,theta,dh,rx,ry,rleast,di,dc;
  int ***pthresh,*nbigfof,**nthresh,***pedge,u,v,**nedge,pleast[2];
  int subpix,**psubs,nsub,*subtail;
  float ***lthresh,**tmpimage,*tlist;
  pthresh=i3tensor(0,Nstar,0,10,0,Npix*2);//assum only 10 groups
  pedge=i3tensor(0,Nstar,0,10,0,Npix*2);
  nedge=imatrix(0,Nstar,0,10);
  nthresh=imatrix(0,Nstar,0,10);
  nbigfof=ivector(0,Nstar);
  lthresh=f3tensor(0,Nstar,0,10,1,Npix);
  tlist=vector(1,Npix);
  tmpimage=matrix(0,Npix,0,Npix);

  for(ic=0;ic<Nstar;ic++){
    l=0;
    for(i=0;i<ngroups[ic];i++){
      peak=0;sum=0;sumx=0;sumy=0;
      for(m=0;m<Npix;m++)for(n=0;n<Npix;n++)tmpimage[m][n]=-1;
      //printf("ngroups=%d,tails=%d,ic=%d\n",ngroups[ic],tails[ic][i],ic );
      if(tails[ic][i]>20){
        for(k=0;k<tails[ic][i];k++){
          ra=pgroup[ic][i][k]/Npix;rb=pgroup[ic][i][k]%Npix;
          tmpimage[ra][rb]=cutstar[ic][ra][rb];
          peak=max(peak,cutstar[ic][ra][rb]);
          if(cutstar[ic][ra][rb]==peak){//peak detection
            centx=ra;centy=rb;
          }
        }//printf("x=%f,y=%f\n",centx,centy);
        //center pick
        for(k=0;k<tails[ic][i];k++){
          ra=pgroup[ic][i][k]/Npix;rb=pgroup[ic][i][k]%Npix;
          r=(ra-centx)*(ra-centx)+(rb-centy)*(rb-centy);
          sum+=cutstar[ic][ra][rb]*exp(-r/9/2);
          sumx+=cutstar[ic][ra][rb]*ra*exp(-r/9/2);
          sumy+=cutstar[ic][ra][rb]*rb*exp(-r/9/2);
        }
        centx=sumx/sum;centy=sumy/sum;//printf("x=%f,y=%f,ic=%d\n",centx,centy,ic+1 );

        pin=0;
        for(k=0;k<tails[ic][i];k++){
          m=pgroup[ic][i][k]/Npix;n=pgroup[ic][i][k]%Npix;
          s=0;
          for(u=-1;u<=1;u++)for(v=-1;v<=1;v++){
            if(m+u>=0 && m+u<nimgx && n+v>=0 &&n+v<nimgy){
              if(tmpimage[m+u][n+v]>0){s++;}
            }
          }if(s<9){//edge detection
            pedge[ic][l][pin]=pgroup[ic][i][k];pin++;
          }
        }nedge[ic][l]=pin;//printf("pin=%d\n",pin );

        sum=2.*Npix*Npix;sumx=0;
        for(k=0;k<nedge[ic][l];k++){
          ra=pedge[ic][l][k]/Npix;rb=pedge[ic][l][k]%Npix;
          rx=ra-centx;ry=rb-centy;
          r=rx*rx+ry*ry;
          sum=min(sum,r);sumx=max(sumx,r);
          if(r==sum){//pick the least pixel to the center
            pleast[0]=ra;pleast[1]=rb;theta=atan(rx/ry);
          }
        }
        rleast=sum;subpix=(int)sumx+1;

        f=0;
        for(k=0;k<tails[ic][i];k++){
          ra=pgroup[ic][i][k]/Npix;rb=pgroup[ic][i][k]%Npix;
          rx=ra-centx;ry=rb-centy;
          r=rx*rx+ry*ry;
          di=(ra-pleast[0])*(ra-pleast[0])+(rb-pleast[1])*(rb-pleast[1]);
          theta0=atan(rx/ry)-theta;
          dh=r*sin(theta0)*sin(theta0);
          if(dh<0.5 && di<rleast && r>0 && r<=rleast){
            pthresh[ic][l][f]=pgroup[ic][i][k];f++;
            lthresh[ic][l][f]=cutstar[ic][ra][rb];
          }
        }nthresh[ic][l]=f;
        sort(f,lthresh[ic][l]);
        subfind(cutstar[ic],lthresh[ic][l],nthresh[ic][l],Npix,Npix,tails[ic][i],
          pgroup[ic][i],psubs,&nsub,subtail,subpix);
        if(nsub>1)printf("nsub=%d,ic=%d\n",nsub,ic+1 );
        l++;
      }
    }nbigfof[ic]=l;//printf("nbigfof=%d\n",nbigfof[ic] );

  }
  //for(ic=0;ic<Nstar;ic++)printf("nbigfof[ic]=%d,ic=%d\n", nbigfof[ic],ic);

  */
