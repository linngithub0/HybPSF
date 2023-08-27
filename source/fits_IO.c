/*-------------------------------------------------------------
 Program to read a GREAT08 fits file

 Uses cfitsio which can be downloaded from
 http://heasarc.gsfc.nasa.gov/docs/software/fitsio/fitsio.html

 Compile with a C compiler (cc for example) (makefile attached)

 To compile:
 make -f Makefile_Example
 need to change the paths to point to the directory containing 
 cfitsio libraries

 To run:
 ./read_GREAT08_fits $GREAT08DIR/set0001.fit postagestamps.fit

 Author: Thomas Kitching on behalf of GREAT08
 C fits code questions?  e-mail tdk@astro.ox.ac.uk
 GREAT08 questions?   e-mail sarah@sarahbridle.net

 http://www.great08challenge.info/
 
-------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "fitsio.h"
#include "longnam.h"

int get_fits_size(char *argv,int *dim)
{

  /* fits image*/
  fitsfile *afptr; /*fitsfile is a C variable defined by cfitsio*/
  int status = 0;  /*  CFITSIO status value must be initialized to zero  */
  int  anaxis;
  long anaxes[2] = {1,1}, fpixel[2]={1,1}; //dimensions of the input fitsfile
  /* open input image */
  fits_open_file(&afptr, argv, READONLY, &status); 

  /* read input image dimensions */
  fits_get_img_dim(afptr, &anaxis, &status);  
  fits_get_img_size(afptr, 2, anaxes, &status);

  if (status) {
    fits_report_error(stderr, status); /* print error message */
    return(status);
  }
  
  if (anaxis != 2) {
    printf("Error: images with other than 2 dimensions are not supported\n");
    exit (2);
  }

  /*define input image dimensions*/
  dim[0] = (int)anaxes[0]; 
  dim[1] = (int)anaxes[1];
  /* close main image file */
  fits_close_file (afptr,&status);
  if (status) {
    fits_report_error(stderr, status); /* print error message */
    return(status);
  }
  return 1;
}

int read_fits_2D(char *argv,float *galaxy,int imagesize)
{
  /* fits image*/
  fitsfile *afptr; /*fitsfile is a C variable defined by cfitsio*/
  int status = 0;  /*  CFITSIO status value must be initialized to zero  */

  /*fits image properties*/
  int  anaxis, size3d;
  long anaxes[2] = {1,1}, fpixel[2]={1,1}; //dimensions of the input fitsfile
  /* open input image */
  fits_open_file(&afptr, argv, READONLY, &status); 
  /* read input image dimensions */
  fits_get_img_dim(afptr, &anaxis, &status);  
  fits_get_img_size(afptr, 2, anaxes, &status);

  if (status) {
    fits_report_error(stderr, status); /* print error message */
    return(status);
  }
  
  if (anaxis != 2) {
    printf("Error: images with other than 2 dimensions are not supported\n");
    exit (2);
  }

 /*define input image dimensions*/
  if (galaxy==NULL){
      fprintf(stderr," error allocating memory for image\n");
      exit (1);
    }
  
  /* read input data into image array */
  if (fits_read_pix(afptr, TFLOAT, fpixel, imagesize, NULL, galaxy,
		    NULL, &status) )
    {
      printf(" error reading pixel data \n");
      exit (2);
    }

  /* close main image file */
  fits_close_file (afptr,&status);
  if (status) {
    fits_report_error(stderr, status); /* print error message */
    return(status);
  }
  return 1;
}

int write_fits_2D( char *argv,float **stars,int *dim)
{

  int i,j,n,nobj,pwidth,pheight,k;
  /* fits image*/
  fitsfile *afptr; /*fitsfile is a C variable defined by cfitsio*/
  int status = 0;  /*  CFITSIO status value must be initialized to zero  */
  /*fits image properties*/
  int  anaxis, size2d;
  long anaxes_out[2], fpixel_out[2]={1,1}; //dimensions of the output fitfile
  int  bitpix=-32;  //IEEE single precision floating point (LT) 
  /* postage stamp data*/
  float *array2d;
  
  remove(argv); //remove file if it already exists
  fits_create_file(&afptr, argv, &status); 
  if (status) {
    fprintf (stderr," error creating output FITS table \n");
    fits_report_error(stderr, status); /*  print error message  */
    exit (2);
  }
  /*define the size of the 2d data cube*/
  pwidth=dim[0];pheight=dim[1];
  anaxis = 2;
  anaxes_out[0] = pwidth;
  anaxes_out[1] = pheight;
  size2d = anaxes_out[0]*anaxes_out[1];
  /* Need to fill the 3D data cube correctly*/ 
  array2d = (float*)calloc(size2d,sizeof(float));
  k=0;
  for (j=0; j<pheight; j++){
    for(i=0;i<pwidth;i++){
      array2d[k] = stars[i][j];
      k++;
    }
  }
  /*first create the image, bits, and dimensionsm, to be filled with data*/
  fits_create_img(afptr, bitpix, anaxis, anaxes_out, &status);
  if (status) {
    fits_report_error(stderr, status); /* print error message */
    return(status);
  }
  /*then fill the define image with data from array3d*/
  if (fits_write_pix(afptr, TFLOAT, fpixel_out, size2d, array2d, &status))
    {
      printf(" error reading pixel data \n");
      exit (2);
    }
  /* close the fits file */
  fits_close_file(afptr, &status);
  free(array2d);
  return 1;
}



int write_fits_3D(char *argv,float ***stars,int *dim)
{

  int i,j,n,nobj,pwidth,pheight,k;
  /* fits image*/
  fitsfile *afptr; /*fitsfile is a C variable defined by cfitsio*/
  int status = 0;  /*  CFITSIO status value must be initialized to zero  */
  /*fits image properties*/
  int  anaxis, size3d;
  long anaxes_out[3], fpixel_out[3]={1,1,1}; //dimensions of the output fitfile
  int  bitpix=-32;  //IEEE single precision floating point (LT) 
  /* postage stamp data*/
  float *array3d;
  
  remove(argv); //remove file if it already exists
  fits_create_file(&afptr, argv, &status); 
  if (status) {
    fprintf (stderr," error creating output FITS table \n");
    fits_report_error(stderr, status); /*  print error message  */
    exit (2);
  }
  /*define the size of the 3d data cube*/
  nobj=dim[2];pwidth=dim[0];pheight=dim[1];
  anaxis = 3;
  anaxes_out[2] = nobj;
  anaxes_out[0] = pwidth;
  anaxes_out[1] = pheight;
  size3d = anaxes_out[0]*anaxes_out[1]*anaxes_out[2];
  /* Need to fill the 3D data cube correctly*/ 
  array3d = (float*)calloc(size3d,sizeof(float));
  k=0;
  for (n=0; n<nobj; n++){
      /* write postage stamps into 3D array */
    for (j=0; j<pheight; j++){
      for(i=0;i<pwidth;i++){
	array3d[k] = stars[n][i][j];
	k++;
      }
    }
  }
  /*first create the image, bits, and dimensionsm, to be filled with data*/
  fits_create_img(afptr, bitpix, anaxis, anaxes_out, &status);
  if (status) {
    fits_report_error(stderr, status); /* print error message */
    return(status);
  }
  /*then fill the define image with data from array3d*/
  if (fits_write_pix(afptr, TFLOAT, fpixel_out, size3d, array3d, &status))
    {
      printf(" error reading pixel data \n");
      exit (2);
    }
  /* close the fits file */
  fits_close_file(afptr, &status);
  free(array3d);
  return 1;
}

int get_fits_header(char *argv,double **keyva,char *time_obs)
//edit by Linn 2018/12/12 22:27:40,this function only fit for the CFHTLens data
{

  /* fits image*/
  fitsfile *afptr; /*fitsfile is a C variable defined by cfitsio*/
  int status = 0;  /*  CFITSIO status value must be initialized to zero  */
  int  anaxis,keysexist,morekeys,hdutype=ASCII_TBL,hdupos,hdunum,keynum,extver;
  long anaxes[2] = {1,1}, fpixel[2]={1,1}; //dimensions of the input fitsfile
  char CRPIX1[100]={"CRPIX1"},CRPIX2[100]={"CRPIX2"};//pixel position of the reference
  char NAXIS1[100]={"NAXIS1"},NAXIS2[100]={"NAXIS2"};
  char CD1_1[100]={"CD1_1"},CD1_2[100]={"CD1_2"},CD2_1[100]={"CD2_1"},CD2_2[100]={"CD2_2"};
  char CRVAL1[100]={"CRVAL1"},CRVAL2[100]={"CRVAL2"},GAIN[100]={"GAIN"};//coordinate of reference
  char CTYPE1[100]={"CTYPE1"},CTYPE2[100]={"CTYPE2"};
  char MJD_OBS[100]={"MJD-OBS"};
  /* open input image */
  printf("%s\n",argv );
  fits_open_image(&afptr, argv, READONLY, &status); 
  if (status) {
    fits_report_error(stderr, status); // print error message
    return(status);exit(EXIT_FAILURE);
  }

  /*fits_open_data(&afptr, argv, READONLY, &status);
  if (status) {
    fits_report_error(stderr, status); // print error message
    return(status);exit(EXIT_FAILURE);
  }*/

  //fits_get_num_hdus(afptr, &hdunum, &status);printf("hdunum=%d\n",hdunum );
  //fits_get_hdu_num(afptr, &hdupos);printf("hdupos=%d\n",hdupos );
  //fits_movrel_hdu(afptr, 0, &hdutype, &status);
  //fits_movrel_hdu(afptr, 1, NULL, &status);

  //fits_get_hdrspace(afptr,&keysexist,&morekeys,&status);
  //printf("morekeys=%d,keysexist=%d\n",morekeys,keysexist );
  //fits_get_hdrpos(afptr,&keysexist, &keynum, &status);
  //printf("keysexist=%d,keynum=%d\n",keysexist,keynum );
  /*int extend=TRUE,naxis,simple=TRUE,maxdim,bitpix;
  long gcount=0,pcount=0,naxes;
  fits_read_imghdr(afptr,maxdim,&simple,&bitpix,&naxis,&naxes,&pcount,&gcount,&extend,&status);
  printf("maxdim=%d,simple=%d,bitpix=%d,naxis=%d,naxes=%d,pcount=%d,gcount=%d,extend=%d\n",
    maxdim,simple,bitpix,naxis,naxes,pcount,gcount,extend );*/
  //char extname[100]={"RAW,07  "};extver=0;
  //fits_movnam_hdu(afptr, ASCII_TBL, extname, extver,&status);
  //fits_get_hdrspace(afptr,&keysexist,&morekeys,&status);
  //printf("morekeys=%d,keysexist=%d\n",morekeys,keysexist );

  //fits_get_hdrpos(afptr,&keysexist, &keynum, &status);printf("keysexist=%d,keynum=%d\n",keysexist,keynum );



  
  char comment[100],date[100],card[100];double value;
  fits_read_key(afptr,TINT,NAXIS1,&anaxis,comment,&status);
  keyva[0][0]=anaxis;//printf("NAXIS1=%d\n",anaxis );
  fits_read_key(afptr,TINT,NAXIS2,&anaxis,comment,&status);
  keyva[0][1]=anaxis;//printf("NAXIS2=%d\n",anaxis );
  fits_read_key(afptr,TDOUBLE,CRPIX1,&value,comment,&status);
  keyva[1][0]=value;//printf("CRPIX1=%e\n",value );
  fits_read_key(afptr,TDOUBLE,CRPIX2,&value,comment,&status);
  keyva[1][1]=value;//printf("CRPIX2=%e\n",value );
  fits_read_key(afptr,TDOUBLE,CD1_1,&value,comment,&status);
  keyva[2][0]=value;//printf("CD1_1=%e\n",value );
  fits_read_key(afptr,TDOUBLE,CD1_2,&value,comment,&status);
  keyva[2][1]=value;//printf("CD1_2=%e\n",value );
  fits_read_key(afptr,TDOUBLE,CD2_1,&value,comment,&status);
  keyva[2][2]=value;//printf("CD2_1=%e\n",value );
  fits_read_key(afptr,TDOUBLE,CD2_2,&value,comment,&status);
  keyva[2][3]=value;//printf("CD2_2=%e\n",value );
  fits_read_key(afptr,TDOUBLE,CRVAL1,&value,comment,&status);
  keyva[3][0]=value;//printf("CRVAL1=%e\n",value );
  fits_read_key(afptr,TDOUBLE,CRVAL2,&value,comment,&status);
  keyva[3][1]=value;//printf("CRVAL2=%e\n",value );
  fits_read_key(afptr,TDOUBLE,GAIN,&value,comment,&status);
  keyva[4][0]=value;printf("GAIN=%e\n",value );
  //fits_read_key(afptr,TDOUBLE,MJD_OBS,&value,comment,&status);
  //keyva[5][0]=value;printf("MJD_OBS=%e\n",value );
  //strcpy(time_obs,date);

  
  /* close main image file */
  fits_close_file (afptr,&status);
  if (status) {
    fits_report_error(stderr, status); /* print error message */
    return(status);
  }
  return 1;
}

int get_header_key(char *argv,double *keyva,char *keyna)
//edit by Linn 2018/12/12 22:27:40,this function only fit for the CFHTLens data
{

  /* fits image*/
  fitsfile *afptr; /*fitsfile is a C variable defined by cfitsio*/
  int status = 0;  /*  CFITSIO status value must be initialized to zero  */
  int  anaxis,keysexist,morekeys,hdutype=ASCII_TBL,hdupos,hdunum,keynum,extver;
  long anaxes[2] = {1,1}, fpixel[2]={1,1}; //dimensions of the input fitsfile
  /* open input image */
  printf("%s\n",argv );
  fits_open_image(&afptr, argv, READONLY, &status); 
  if (status) {
    fits_report_error(stderr, status); // print error message
    return(status);exit(EXIT_FAILURE);
  }

  
  char comment[100],date[100],card[100];double value;
  fits_read_key(afptr,TDOUBLE,keyna,&value,comment,&status);
  *keyva=value;printf("GAIN=%e\n",value );

  /* close main image file */
  fits_close_file (afptr,&status);
  if (status) {
    fits_report_error(stderr, status); /* print error message */
    return(status);
  }
  return 1;
}

int get_fits_img_data(char *argv){
  /* fits image*/
  fitsfile *afptr; /*fitsfile is a C variable defined by cfitsio*/
  int status = 0,naxis;  /*  CFITSIO status value must be initialized to zero  */
  int  anaxis,keysexist,morekeys,hdutype=ASCII_TBL,hdupos,hdunum,keynum,extver=0;
  fits_open_file(&afptr, argv, READONLY, &status); 
  if (status) {
    fits_report_error(stderr, status); // print error message
    return(status);exit(EXIT_FAILURE);
  }
  fits_get_img_dim(afptr,&naxis,&status);
  printf("naxis=%d\n",naxis );
  return 1;
}

int read_header_strkey(char *fname,int datatype,char *keyname,char *sval,char *comment){

  /* fits image*/
  fitsfile *afptr; /*fitsfile is a C variable defined by cfitsio*/
  int status = 0;  /*  CFITSIO status value must be initialized to zero  */
  int  anaxis,keysexist,morekeys,hdutype=ASCII_TBL,hdupos,hdunum,keynum,extver;
  long anaxes[2] = {1,1}, fpixel[2]={1,1}; //dimensions of the input fitsfile
  /* open input image */
  //printf("%s\n",argv );
  fits_open_image(&afptr, fname, READONLY, &status); 
  if (status) {
    fits_report_error(stderr, status); // print error message
    return(status);exit(EXIT_FAILURE);
  }

  
  char date[100],card[100];
  fits_read_keyword(afptr, keyname, sval, comment,&status);
  /* close main image file */
  fits_close_file (afptr,&status);
  if (status) {
    fits_report_error(stderr, status); /* print error message */
    return(status);
  }
  return 1;

} 


int write_header_strkey(char *fname,int datatype,char *keyname,char *sval,char *comment){

  /* fits image*/
  fitsfile *afptr; /*fitsfile is a C variable defined by cfitsio*/
  int status = 0;  /*  CFITSIO status value must be initialized to zero  */
  int  anaxis,keysexist,morekeys,hdutype=ASCII_TBL,hdupos,hdunum,keynum,extver;
  long anaxes[2] = {1,1}, fpixel[2]={1,1}; //dimensions of the input fitsfile
  /* open input image */
  //printf("%s\n",argv );
  fits_open_image(&afptr, fname, READONLY, &status); 
  if (status) {
    fits_report_error(stderr, status); // print error message
    return(status);exit(EXIT_FAILURE);
  }

  
  char date[100],card[100];
  fits_write_key(afptr, TSTRING, keyname,sval,comment,&status);

  /* close main image file */
  fits_close_file (afptr,&status);
  if (status) {
    fits_report_error(stderr, status); /* print error message */
    return(status);
  }
  return 1;

} 




/******************************************************************/
void printerror( int status)
{
    /*****************************************************/
    /* Print out cfitsio error messages and exit program */
    /*****************************************************/

    if (status)
    {
       fits_report_error(stderr, status); /* print error report */
       exit( status );    /* terminate the program, returning error status */
    }
    return;
}

