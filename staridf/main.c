#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"
#include "fitsio.h"
int main(int argc,char *argv[])
{
  int get_fits_size(char *argv,int *dim);
  int write_fits_2D(char *argv,float **stars,int *dim);
  int write_fits_3D(char *argv,float ***stars,int *dim);
  int read_fits_2D(char *argv,float *galaxy,int imagesize);
  void starbox(float *flux,float *size1,float **box1,
	       float *size2,float **box,int Ntotal,char *boxname);
  void trous_decomp(double **I0,double ***w,double **cJ,int nx,int ny,int J);
  double k_sigma(double **data,int nx,int ny);
  char* itoa(int num,char *str,int radix);
  int pathexist(char *path);
  void sigma_neibor(float **image,int Nstar,int Npix,double **var,float **pos,
    int nimgx,int nimgy);
  void fof_image(float **image,double sigma,int nimgx,int nimgy,int *groupnum,
    int *tails,int **pgroup);
  void deblend(float **image,int nimgx,int nimgy,int *ngroups,int *tails,
    int **pgroup,int *nsub,float sigma);
  void deblend1(float **image,int nimgx,int nimgy,int *ngroups,int *tails,
    int **pgroup,int *nsub,float sigma);
  double poidev(double xm,long *idum);
  void funcpsfs_size(char *fname,int *Nf,int *khead);
  void funcpsfs_read(char *fname,int Nf,int khead,float **para,int npara,char *delim);
  void mask_groups(int count,int ngroup,int *tails,int **pgroup,float **image,
    double sigma,int nimgx,int nimgy,float **mask);
  void  gaus_estimate(float **star0,int Ng1,int Ng2,double *mean,
         double *sigma);
  float **pos,*FWHM,*R50,*RMS,*flux,*class;
  float minx,maxx,min50,max50,minf,maxf,x,y,x1,y1,yf,y50,**box50,**boxf;
  float *image0,**image,***Stars,gain,**mpara;
  long iseed=-1;
  double **stamp,***wij,**cJ;
  int *flag,J=3,ir,is,ic,slen=1000;
  int Ntotal,Nstar,i,j,ix,iy,m,n,k,Npix=80,nt,nimgx,nimgy,imgsize,Nh,dim[3];
  int Nhead,colum;
  char *fmid0="R43/",*cend="output",fmid[20],fmid1[9][20],delim[20];
  char head0[slen],fname0[slen],fnames[slen],fnamel[slen],fnameln[slen],fbox[slen],ctmp[slen],
    flist0[slen],fgalaxy[slen],ffix[slen],cata[slen],header_app[slen],fweight[slen],slicename[slen]; 
  FILE *fp,*xp;
  box50=matrix(0,1,0,1);boxf=matrix(0,1,0,1);//R50 & FWHM

  if(argc!=8 && argc!=6 && argc!=9 && argc!=10){
    printf("parameter input error\n");
    printf("should be arranged as:\n");
    printf("or :%s <raw cata file path> <input cata name> <fits path> <input fitsname> <output path> [<left> <right>] [<gain>] [<weight fits>]\n",
      argv[0]);
    exit(EXIT_FAILURE);
  }
  strcpy(header_app,"[SCI]");

  strcpy(flist0,argv[1]);strcat(flist0,argv[2]);  //input catalogue
  strcpy(fname0,argv[3]);strcat(fname0,argv[4]);//input fits file
  strcpy(cata,argv[5]);strcat(cata,"/catalogue");
  printf("catalogue:%s\n",flist0);printf("fits:%s\n",fname0);printf("cata path:%s\n",cata);
  pathexist(cata);
  is=strlen(argv[4]);
  for(i=0;i<is-5;i++)slicename[i]=argv[4][i];
  slicename[i]='\0';printf("slicename=%s\n",slicename );
  /*pathexist("catalogue/star_stamps");
  pathexist("catalogue/star");
  pathexist("catalogue/source");
  pathexist("catalogue/box");*/
  strcpy(fnames,cata);strcat(fnames,"/star_stamps/");//output star stamps
  pathexist(fnames);
  strcpy(fnamel,cata);strcat(fnamel,"/star/");//output star list
  pathexist(fnamel);
  strcpy(fnameln,cata);strcat(fnameln,"/source/");//output source catalogue
  pathexist(fnameln);
  strcpy(fbox,cata);strcat(fbox,"/box/");//output starbox
  pathexist(fbox);

  
  

  strcat(fnames,slicename); strcat(fnames,"_star"); //star stamps; append the file name at the path
  strcat(fnamel,slicename); strcat(fnamel,"_star"); //star list
  strcat(fnameln,slicename); strcat(fnameln,"_source"); //source catalogue
  strcat(fbox,slicename);  strcat(fbox,"_box"); //starbox 

  strcpy(ffix,".fits");
  strcat(fnames,ffix);
  strcpy(ffix,".cat");
  strcat(fnamel,ffix);
  strcat(fnameln,ffix);
  strcat(fbox,ffix);

  printf("%s\t%s\t%s\t%s\t%s\t%s\n",flist0,fname0,fnames,fnameln,fnameln,fbox );

  
  
    //Transform flux to adu value
    char comment[100],comment_DATE[100],DATE[100]={"DATE_BEG"},value_DATE[100],stmp[100];    
    float xposure,photmjsr;
    fitsfile *afptr;
    int status = 0;
    int  anaxis,keysexist,morekeys;
    char XPOSURE[100]={"XPOSURE"},PHOTMJSR[100]={"PHOTMJSR"};
    strcpy(ctmp,fname0);strcat(ctmp,header_app);printf("%s\n",ctmp );
    double value;
    fits_open_file(&afptr, ctmp, READONLY, &status);
    if (status) {
      fits_report_error(stderr, status); /* print error message */
      return(status);exit(EXIT_FAILURE);
    }
    fits_read_key(afptr,TDOUBLE,XPOSURE,&value,comment,&status);//printf("%d\n",status );
    xposure=value;printf("xposure=%f\n", xposure);
    fits_read_key(afptr,TDOUBLE,PHOTMJSR,&value,comment,&status);//printf("%d\n",status );
    photmjsr=value;printf("photmjsr=%f\n", photmjsr);
    fits_close_file (afptr,&status);
    
    
    
    
  //read gain from command line or from fits file
  if(argc==9){
    gain=atof(argv[8]);
  }
  else{
    fitsfile *afptr;
    int status = 0;
    int  anaxis,keysexist,morekeys;
    char comment[100],date[100],GAIN[100]={"GAIN"};
    strcpy(ctmp,fname0);strcat(ctmp,header_app);printf("%s\n",ctmp );
    double value;
    fits_open_file(&afptr, ctmp, READONLY, &status);
    if (status) {
      fits_report_error(stderr, status); /* print error message */
      return(status);exit(EXIT_FAILURE);
    }
    //fits_get_hdrspace(afptr,&keysexist,&morekeys,&status);
    fits_read_key(afptr,TDOUBLE,GAIN,&value,comment,&status);printf("%d\n",status );
    gain=value;printf("gain=%f\n", gain);
    fits_close_file (afptr,&status);
  }printf("gain=%f\n",gain);

  //read the OPD information
    
    /*fits_open_file(&afptr, ctmp, READONLY, &status);
    if (status) {
      fits_report_error(stderr, status); // print error message 
      return(status);exit(EXIT_FAILURE);
    }
    //fits_read_key(afptr,TSTRING,DATE,&value,comment,&status);printf("%d\n",status );
    keysexist=fits_read_keyword(afptr,DATE, stmp,comment_DATE,&status);
    //printf("%d\n",keysexist);
    if(keysexist!=KEY_NO_EXIST){
    printf("len:%d\n",ic=strlen(stmp) );
    for(i=1;i<ic-1;i++){
      value_DATE[i-1]=stmp[i];
    }
    value_DATE[i-1]='\0';
    printf("DATE:%s,comment:%s\n",value_DATE,comment_DATE );
    //fits_close_file (afptr,&status);
    }
    else{
      strcpy(value_DATE,"2022-07-12T17:32:18.400");
      strcpy(comment_DATE,"/ UTC date file created");
    }
    fits_close_file(afptr,&status);
    printf("OPD:%s\n",value_DATE);*/

  //get the magnitude section
  if(argc==6){
    minx=boxf[0][0]=box50[0][0]=-13;maxx=boxf[1][0]=box50[1][0]=-8;
  }
  else{
    minx=boxf[0][0]=box50[0][0]=atof(argv[6]);maxx=boxf[1][0]=box50[1][0]=atof(argv[7]);//left & right
  }
  boxf[0][1]=2;box50[0][1]=1;//bottom
  stamp=dmatrix(0,Npix-1,0,Npix-1);
  cJ   =dmatrix(0,Npix-1,0,Npix-1);
  wij  =d3tensor(0,J-1,0,Npix-1,0,Npix-1);

  /*fp=fopen(flist0,"r");printf("%s\n",flist0 );
  if(fp==NULL)printf("NULL, line 81\n");
  k=1;
  while(k!=0){
    fscanf(fp,"%s\n",ctmp); 
    k=strcmp(ctmp,cend);  }
  do{fscanf(fp,"%d %e %e %e %e %e %e %d %e\n",&i,&x,&x,&x,&x,&x,&x,&j,&x); 
  }while(!feof(fp));
  fclose(fp);printf("i=%d\n",i);Ntotal=i;*/
  funcpsfs_size(flist0,&Ntotal,&Nhead);colum=Nhead;
  printf("Ntotal=%d\n",Ntotal );
  pos   =matrix(0,Ntotal-1,0,1);
  FWHM  =vector(0,Ntotal-1);
  R50   =vector(0,Ntotal-1);
  RMS   =vector(0,Ntotal-1);
  flux  =vector(0,Ntotal-1);
  class =vector(0,Ntotal-1);
  flag  =ivector(0,Ntotal-1);
  mpara=matrix(0,Ntotal,0,colum);
  strcpy(delim,"\t");printf("column=%d\n",colum);
  funcpsfs_read(flist0,Ntotal,Nhead,mpara,colum,delim);
  printf("read func, total source:%d\n",Ntotal);


  //k=0;
  //fp=fopen(flist0,"r");
  k=0;
  //while(k!=0){
  //  fscanf(fp,"%s\n",ctmp); 
  //  k=strcmp(ctmp,cend);
  //}
  for(i=0;i<Ntotal;i++){
    //fscanf(fp,"%d %e %e %e %e %e %e %d %e\n",&j,&pos[i][0],&pos[i][1],
    //        &FWHM[i],&R50[i],&flux[i],&RMS[i],&flag[i],&class[i]);  
    //printf("%d %e %e %e %e %e %e %d %e\n",j,pos[i][0],pos[i][1],FWHM[i],R50[i],flux[i],RMS[i],flag[i],class[i]); 
    pos[i][0]=mpara[i][1];pos[i][1]=mpara[i][2];
    FWHM[i]=mpara[i][3];R50[i]=mpara[i][4];flux[i]=mpara[i][5];
    RMS[i]=mpara[i][6];flag[i]=mpara[i][7];class[i]=mpara[i][8];
    if(flux[i]>0 && FWHM[i]>1.2 && R50[i]>0.7){
      pos[k][0]=pos[i][0];
      pos[k][1]=pos[i][1];
      R50[k]=R50[i];
      FWHM[k]=FWHM[i];
      flux[k]=flux[i];
      RMS[k]=RMS[i];
      flag[k]=flag[i];
      class[k]=class[i];
      k++;
    }
  }
  //fclose(fp);
  
  Ntotal=k;
  
  starbox(flux,FWHM,boxf,R50,box50,Ntotal,fbox);printf("finish starbox\n");
  minf=boxf[0][1];maxf=boxf[1][1];
  min50=box50[0][1];max50=box50[1][1];
  
  strcpy(ctmp,fname0);strcat(ctmp,header_app);
  get_fits_size(ctmp,dim);//printf("x=%d,y=%d\n",dim[0],dim[1] );
  nimgx=dim[0];nimgy=dim[1];
  k=0;
  for(i=0;i<Ntotal;i++){
    x=-2.5*log10(flux[i]);yf=FWHM[i],y50=R50[i];
    if(x>minx & x<maxx & yf>minf & yf<maxf & y50>min50 & y50<max50 &&
      pos[i][0]>Npix && pos[i][1]>Npix && pos[i][0]<nimgx-Npix && 
      pos[i][1]<nimgx-Npix){
      pos[k][0]=pos[i][0];
      pos[k][1]=pos[i][1];
      R50[k]=R50[i];
      FWHM[k]=FWHM[i];
      flux[k]=flux[i];
      RMS[k]=RMS[i];
      flag[k]=flag[i];
      class[k]=class[i];
      k++;
    }
  }
  Nstar=k;
  
  if(Nstar==0){
    printf("No star is found,change the window\n");exit(EXIT_FAILURE);
  }

  //////////////////////////////////////////////////////////////////
  

  fp=fopen(fnameln,"w");
  //fprintf(fp,"%d\n",Ntotal);
  fprintf(fp, "#   1 NUMBER                 Running object number\n" );
  fprintf(fp, "#   2 X_IMAGE                Object position along x      [pixel]\n" );
  fprintf(fp, "#   3 Y_IMAGE                Object position along y      [pixel]\n" );
  fprintf(fp, "#   4 FWHM_IMAGE             FWHM assuming a gaussian core[pixel]\n" );
  fprintf(fp, "#   5 FLUX_RADIUS            Fraction-of-light radii      [pixel]\n" );
  fprintf(fp, "#   6 FLUX_AUTO              Flux within a Kron-like elliptical aperture [counts]\n" );
  fprintf(fp, "#   7 FLUXERR_AUTO           RMS error for AUTO flux      [counts]\n" );
  fprintf(fp, "#   8 FLAGS                  Extraction flags\n" );
  fprintf(fp, "#   9 CLASS_STAR             S/G classifier output\n" );
  m=0;
  for(i=0;i<Nstar;i++){
    k=0;//printf("i=%d\n",i );
    x=pos[i][0];y=pos[i][1];
    x--;y--;
    for(j=0;j<Nstar;j++){
      x1=fabs(pos[j][0]-x)*2.;y1=fabs(pos[j][1]-y)*2.;
      if(x1<Npix/4. && y1<Npix/4.){
        //if(i!=j)printf("i=%d,j=%d\n",i,j );
        k++;
      }
    }
    if(k==1){//exclude the closing objects
      pos[m][0]=pos[i][0];pos[m][1]=pos[i][1];FWHM[m]=FWHM[i];
      R50[m]=R50[i];flux[m]=flux[i];RMS[m]=RMS[i];flag[m]=flag[i];class[m]=class[i];
      fprintf(fp,"%d\t%e\t%e\t%e\t%e\t%e\t%e\t%d\t%e\n",m+1,pos[m][0],pos[m][1],FWHM[m],
        R50[m],flux[m],RMS[m],flag[m],class[m]);
      m++;
    }
  }
  fclose(fp);
  printf("m=%d\n", m);Nstar=Ntotal=m;
  printf("start starbox\n");

  //////////////////////////////////////////////////////////////////
  // printf("test......\n");
  strcpy(ctmp,fname0);strcat(ctmp,header_app);
  printf("%s\n",ctmp );
  get_fits_size(ctmp,dim);//printf("x=%d,y=%d\n",dim[0],dim[1] );
  nimgx=dim[0];nimgy=dim[1];imgsize=nimgx*nimgy;
  image0=vector(0,imgsize-1);
  image=matrix(0,nimgx-1,0,nimgy-1);
  Stars=f3tensor(0,Nstar,0,Npix-1,0,Npix-1);
  read_fits_2D(ctmp,image0,imgsize);
  for(j=0;j<nimgy;j++){
    for(i=0;i<nimgx;i++){k=j*nimgx+i;image[i][j]=image0[k]*xposure*photmjsr;}
  }
  //free_vector(image0,0,imgsize-1);
  Nh=Npix/2.;
  for(nt=0;nt<Nstar;nt++){
    ix=pos[nt][0]+0.5;iy=pos[nt][1]+0.5; 
    //ix--;iy--;
    for(j=0;j<Npix;j++){
      for(i=0;i<Npix;i++){
        m=ix+i-Nh;
        n=iy+j-Nh;
    if(m>=0 & m<nimgx & n>=0 & n<nimgy)
        stamp[i][j]=Stars[nt][i][j]=image[m][n];
    else stamp[i][j]=Stars[nt][i][j]=0;
      }
    }
    trous_decomp(stamp,wij,cJ,Npix,Npix,J);
    RMS[nt]=k_sigma(wij[0],Npix,Npix)/0.899;
  }
  dim[0]=dim[1]=Npix;dim[2]=Nstar;
  //write_fits_3D("fits/test.fits",Stars,dim);

  /////////////////////////////////////////////////////////////////////
  if(argc==10){
    printf("Using weight to select the non-contaminated stamps\n");
    float **imagew,***weight;
    imagew=(float **)calloc(nimgx,sizeof(float *));
    weight=(float ***)calloc(Nstar,sizeof(float **));
    for(i=0;i<nimgx;i++){
      imagew[i]=(float *)calloc(nimgy,sizeof(float));
    }
    for(i=0;i<Nstar;i++){
      weight[i]=(float **)calloc(Npix,sizeof(float *));
      for(j=0;j<Npix;j++){
        weight[i][j]=(float *)calloc(Npix,sizeof(float));
      }
    }//printf("x=%d,y=%d\n",nimgx,nimgy );
    strcpy(fweight,argv[3]);strcat(fweight,argv[8]);printf("%s\n",fweight );
    get_fits_size(fweight,dim);//printf("x=%d,y=%d\n",dim[0],dim[1] );
    nimgx=dim[0];nimgy=dim[1];imgsize=nimgx*nimgy;
    read_fits_2D(fweight,image0,imgsize);
    for(j=0;j<nimgy;j++){
      for(i=0;i<nimgx;i++){k=j*nimgx+i;imagew[i][j]=image0[k];}
    }//printf("read down\n");
    //dim[0]=nimgx;dim[1]=nimgy;
    //write_fits_2D("fits/test.fits",imagew,dim);
    ic=0;is=0;
    for(nt=0;nt<Nstar;nt++){
      k=0;
      ix=pos[nt][0]+0.5;iy=pos[nt][1]+0.5; 
      //ix--;iy--;//printf("x=%d,y=%d\n",ix,iy );
      for(i=0;i<Npix/2;i++){
        for(j=0;j<Npix/2;j++){
          m=ix+i-Nh/2;
          n=iy+j-Nh/2;
          weight[nt][i][j]=imagew[m][n];
          if(imagew[m][n]<0){
            k++;
          }
        }   
      }//printf("nt=%d\n",nt );
      if(k<2){
        for(ix=0;ix<Npix;ix++)for(iy=0;iy<Npix;iy++){
          Stars[ic][ix][iy]=Stars[nt][ix][iy];    
        }
        ic++;
      }
      else{
        for(ix=0;ix<Npix/2;ix++)for(iy=0;iy<Npix/2;iy++){
          weight[is][ix][iy]=weight[nt][ix][iy];
        }is++;
      }
    }
    Nstar=ic;printf("Nstar=%d\n",Nstar );
    strcpy(ctmp,cata);strcat(ctmp,"/star_stamps/");strcat(ctmp,slicename);strcat(ctmp,"_wgt.fits");
    dim[0]=dim[1]=Npix/2;dim[2]=is;printf("start write\n");
    write_fits_3D(ctmp,weight,dim);
  }
  
  
  printf("star snrs\n");
  float SNRs_estimate(float **image,int nx,int ny,float gain);
  float ***wstar;
  wstar=f3tensor(0,Nstar,0,Npix,0,Npix);
  float snr[Nstar];
  for(i=0;i<Nstar;i++){
    snr[i]=SNRs_estimate(Stars[i],Npix,Npix,gain);
    //printf("SNRs=%f\n",snr[i] );
  }
  m=0;
  for(i=0;i<Nstar;i++){
    if(snr[i]>=0.){
      for(ix=0;ix<Npix;ix++)for(iy=0;iy<Npix;iy++){
        wstar[m][ix][iy]=Stars[i][ix][iy];
      }
      pos[m][0]=pos[i][0];pos[m][1]=pos[i][1];FWHM[m]=FWHM[i];
      R50[m]=R50[i];flux[m]=flux[i];RMS[m]=RMS[i];flag[m]=flag[i];class[m]=class[i];
      m++;
      //wstar[m-1]=Stars[i];//exclude the low SNRs images
    }
  }
  

  Nstar=m;
  if(Nstar==0){
    printf("No star's snr is lager than 10, press Enter and exit");
    getchar();exit(EXIT_FAILURE);
  }printf("I found %d prestars\n",Nstar);
  
  double **sigm,*sigma;
  sigm=dmatrix(0,Nstar,0,1);
  sigma=dvector(0,Nstar);
  sigma_neibor(image,Nstar,Npix,sigm,pos,nimgx,nimgy);
  for(ic=0;ic<Nstar;ic++){
    sigma[ic]=sigm[ic][0];
    for(i=0;i<Npix;i++)for(j=0;j<Npix;j++){
      wstar[ic][i][j]-=sigm[ic][1];//subtract the background
    }
  }
  dim[0]=dim[1]=Npix;dim[2]=m;
  //write_fits_3D("fits/wstar.fits",wstar,dim);

  int **isop,*Ncount,*tmps,*tgroup,*tcount,**tails,***pgroup;
  int pin,tail,ngroup,*ngroups,s,*nsub,*nsub1,ra,rb,flags;
  float ***cutstar,***mask;
  double r,sum,sumx,sumy,centx,centy,r0,r2,mean;
  isop=imatrix(0,Nstar,0,Npix*Npix);
  pgroup=i3tensor(0,Nstar,0,Npix*Npix,0,Npix*Npix);
  cutstar=f3tensor(0,Nstar,0,Npix,0,Npix);
  mask=f3tensor(0,Nstar,0,Npix,0,Npix);
  tcount=ivector(0,Npix*Npix);
  tgroup=ivector(0,Npix*Npix);
  ngroups=ivector(0,Nstar);
  nsub=ivector(0,Nstar);
  nsub1=ivector(0,Nstar);
  tails=imatrix(0,Nstar,0,Npix*Npix);
  Ncount=ivector(0,Nstar);
  tmps=ivector(0,Npix*Npix);//printf("Nstar=%d\n",Nstar );
  int group_count=0;
  for(ic=0;ic<Nstar;ic++){
    //printf("in fof_image\n");
    for(i=0;i<Npix;i++)for(j=0;j<Npix;j++)mask[ic][i][j]=1.;
    fof_image(wstar[ic],sigma[ic],Npix,Npix,&ngroup,tails[ic],pgroup[ic]);
    //printf("ngroups=%d,sigma=%f,ic=%d\n",ngroups[ic],sigma[ic],ic);
    if(ngroup>1)group_count+=1;
    mask_groups(ic+1,ngroup,tails[ic],pgroup[ic],wstar[ic],sigma[ic],Npix,Npix,mask[ic]);
    //printf("in subfind\n");
    gaus_estimate(wstar[ic],Npix,Npix,&mean,&sigma[ic]);
    fof_image(wstar[ic],sigma[ic],Npix,Npix,&ngroup,tails[ic],pgroup[ic]);
    ngroups[ic]=ngroup;
    //printf("ngroups=%d,sigma=%f,ic=%d\n",ngroups[ic],sigma[ic],ic);
    mask_groups(ic+1,ngroups[ic],tails[ic],pgroup[ic],wstar[ic],sigma[ic],Npix,Npix,mask[ic]);
    snr[ic]=SNRs_estimate(wstar[ic],Npix,Npix,gain);
  }//printf("Nstar=%d\n",Nstar );
  fp=fopen(fnamel,"w");//xp=fopen("test/tmp.dat","w");
  fprintf(fp, "#   1 NUMBER                 Running object number\n" );
  fprintf(fp, "#   2 X_IMAGE                Object position along x      [pixel]\n" );
  fprintf(fp, "#   3 Y_IMAGE                Object position along y      [pixel]\n" );
  fprintf(fp, "#   4 FWHM_IMAGE             FWHM assuming a gaussian core[pixel]\n" );
  fprintf(fp, "#   5 FLUX_RADIUS            Fraction-of-light radii      [counts]\n" );
  fprintf(fp, "#   6 FLUX_AUTO              Flux within a Kron-like elliptical aperture [counts]\n" );
  fprintf(fp, "#   7 FLUXERR_AUTO           RMS error for AUTO flux      [pixel]\n" );
  fprintf(fp, "#   8 SNR                    signal to noise ratio\n" );
  m=0;n=0;
  for(i=0;i<Nstar;i++){
    if(snr[i]>10. && ngroups[i]!=0){
      
      //printf("ngroups=%d,nsub=%d\n",ngroups[ic],nsub[ic] );
      for(ix=0;ix<Npix;ix++)for(iy=0;iy<Npix;iy++){
        cutstar[m][ix][iy]=wstar[i][ix][iy];
        Stars[m][ix][iy]=mask[i][ix][iy];
      }
      m++;
      fprintf(fp,"%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",m,
            pos[i][0],pos[i][1],FWHM[i],R50[i],flux[i],RMS[i],snr[i]);
      //fprintf(xp, "%d\t%d\n",m,i+1 );
      }
      /*else{
        for(ix=0;ix<Npix;ix++)for(iy=0;iy<Npix;iy++)Stars[n][ix][iy]=wstar[i][ix][iy];
        n++;
      }*/
  }
  fclose(fp);
  printf("I found %d stars, %f stars are blended\n",m, ((group_count+0.)/Nstar));
  strcpy(ctmp,cata);strcat(ctmp,"/star_stamps/");strcat(ctmp,slicename);strcat(ctmp,"_mask.fits");
  //strcpy(fnames,"test/catalogue/star_stamps/selected.fits");
  dim[0]=dim[1]=Npix;dim[2]=m;
  if(m>0){
    write_fits_3D(fnames,cutstar,dim);
    write_fits_3D(ctmp,Stars,dim);
    if(value_DATE){printf("start write OPD\n");
    status=0;
    fits_open_file(&afptr, fnames, READWRITE, &status); 
    //printf("status:%d\n",status);
    if (status) {
      //printf("read succes\n");
      fits_report_error(stderr, status); /* print error message */
      return(status);
      exit(EXIT_FAILURE);
    }
    //printf("DATE:%s,value_DATE:%s,commment_DATE:%s\n",DATE,value_DATE,comment_DATE);
    //fits_write_key(afptr, TSTRING, DATE, value_DATE,comment_DATE,&status);
    printf("write done!\n");
    fits_close_file (afptr,&status);
    }
  }



  free_dmatrix(stamp,0,Npix-1,0,Npix-1);
  free_dmatrix(cJ,0,Npix-1,0,Npix-1);
  free_d3tensor(wij,0,J-1,0,Npix-1,0,Npix-1);
  free_matrix(pos,0,Ntotal-1,0,1);
  free_vector(FWHM,0,Ntotal-1);
  free_vector(RMS,0,Ntotal-1); 
  free_vector(R50,0,Ntotal-1);
  free_vector(flux,0,Ntotal-1);
  free_vector(class,0,Ntotal-1);
  free_ivector(flag,0,Ntotal-1);
  free_matrix(image,0,nimgx-1,0,nimgy-1);
  free_f3tensor(Stars,0,Nstar,0,Npix-1,0,Npix-1);
  free_f3tensor(mask,0,Nstar,0,Npix,0,Npix);
  }


void starbox(float *flux,float *size1,float **box1,
	     float *size2,float **box2,int Ntotal,char *boxname)
{
  void k_sigma_1D(double *data,int np0,double *dms);
  void sort(int n, float arr[]);
  int i,j,k;
  double minx,miny,maxx,maxy,x,y,w=3.,ktmp;
  double *sizet,dms[2];
  float *arr1,*arr2;
  FILE *fp;

  sizet=dvector(0,Ntotal);
  arr1=vector(1,Ntotal);
  arr2=vector(1,Ntotal);
  
  minx=box1[0][0];  
  maxx=box1[1][0]-2.5;
  k=0;
  for(i=0;i<Ntotal;i++){
    x=-2.5*log10(flux[i]);
    if(x>minx & x<maxx ){k++;arr1[k]=size1[i],arr2[k]=size2[i];}
  }
  ktmp=k;
  sort(ktmp,arr1);
  sort(ktmp,arr2);
  y=arr1[5];
  miny=y-box1[0][1]*0.5;
  maxy=y+box1[0][1]*0.5;
  k=0;
  for(i=0;i<ktmp;i++){
    y=arr1[i+1];
    if(y>miny & y<maxy ){sizet[k]=y;k++;}
  }
  k_sigma_1D(sizet,k,dms);printf("k_sigma_1D\n");
  box1[0][1]=dms[0]-w*dms[1];
  box1[1][1]=dms[0]+w*dms[1];
 
  y=arr2[5];
  miny=y-box2[0][1]*0.5;
  maxy=y+box2[0][1]*0.5;
  k=0;
  for(i=0;i<ktmp;i++){
    y=arr2[i+1];
    if(y>miny & y<maxy ){sizet[k]=y;k++;}
  }
  k_sigma_1D(sizet,k,dms);printf("k_sigma\n");
  box2[0][1]=dms[0]-w*dms[1];
  box2[1][1]=dms[0]+w*dms[1];
 
  //free_dvector(sizet,0,Ntotal);
  //free_vector(arr1,1,Ntotal);
  //free_vector(arr2,1,Ntotal);printf("%s\n",boxname);
  fp=fopen(boxname,"w");
  fprintf(fp,"%e %e\n",box1[0][0],box1[0][1]);
  fprintf(fp,"%e %e\n",box1[0][0],box1[1][1]);
  fprintf(fp,"%e %e\n",box1[1][0],box1[1][1]);
  fprintf(fp,"%e %e\n",box1[1][0],box1[0][1]);
  fprintf(fp,"%e %e\n",box1[0][0],box1[0][1]);
  fprintf(fp,"%e %e\n",box2[0][0],box2[0][1]);
  fprintf(fp,"%e %e\n",box2[0][0],box2[1][1]);
  fprintf(fp,"%e %e\n",box2[1][0],box2[1][1]);
  fprintf(fp,"%e %e\n",box2[1][0],box2[0][1]);
  fprintf(fp,"%e %e\n",box2[0][0],box2[0][1]);
  fclose(fp);printf("write down\n");

}


void k_sigma_1D(double *data,int np0,double *dms)
{
  void sigma_clean(double *data,double *r2,double *dms,int *np);
  int i,j,m,n,k,ktotal,np[2];
  double sum=0,dm,sigma,*r2,*d1,det,sigma2;
  r2=dvector(0,np0);
  d1=dvector(0,np0);
  for(i=0;i<np0;i++)sum +=data[i];
  dm=sum/(np0);
  sum=0;
  k=0;
  for(i=0;i<np0;i++){
      det=data[i]-dm;
      r2[i]=det*det;
      d1[i]=data[i];
      sum +=r2[i];
  }
  ktotal=np0;
  sigma=sqrt(sum/(np0-1));
  //--------------------- m0
  dms[0]=dm;dms[1]=sigma;np[0]=np[1]=ktotal;
  k=0;
  do {
    np[0]=np[1];
    sigma_clean(d1,r2,dms,np);
    k++;
  }while((np[0]!=np[1]) && k<5);
 
  free_dvector(r2,0,np0);
  free_dvector(d1,0,np0);
 
}

void CloseP(float **pos,float *flux,float *R50,int Ntotal,float **posgal,
	    float *fluxgal,float *R50gal,float *SNgal,int *Ngal)
{
  int i,j,k,n0,*sign;
  double x,y,detx,dety,r2,r,r50g,fluxg,nr=1.5,thres;
  n0=*Ngal;
  sign=ivector(0,n0);
  for(i=0;i<n0;i++)sign[i]=1;
  for(i=0;i<n0;i++){
    x=posgal[i][0];
    y=posgal[i][1];
    r50g=R50gal[i]*nr;
    fluxg=fluxgal[i]*0.1;
    for(j=0;j<Ntotal;j++){
      detx=x-pos[j][0];
      dety=y-pos[j][1];
      r2=detx*detx+dety*dety;
      thres=r50g+R50[j]*nr;
      thres=thres*thres;
      if(r2>0.0001){
	if(r2<thres)sign[i]=0; 
      }
    }
  }
  k=0;
  for(i=0;i<n0;i++){
    if(sign[i]==1){
      posgal[k][0]=posgal[i][0];
      posgal[k][1]=posgal[i][1];
      fluxgal[k]  =fluxgal[i];
      R50gal[k]   =R50gal[i];
      SNgal[k]    =SNgal[i];
      k++;
    }
  }
  printf("k=%d\n",k);
  *Ngal=k;
}
