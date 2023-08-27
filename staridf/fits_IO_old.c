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

int get_fits_size(char *argv,int *dim)
{printf("get_fits_size,%s\n",argv );

  /* fits image*/
  fitsfile *afptr; /*fitsfile is a C variable defined by cfitsio*/
  int status = 0;  /*  CFITSIO status value must be initialized to zero  */
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
  dim[0] = (int)anaxes[0]; 
  dim[1] = (int)anaxes[1];
  /* close main image file */
  fits_close_file (afptr,&status);
  if (status) {
    fits_report_error(stderr, status); /* print error message */
    return(status);
  }

}

int read_fits_2D(char *argv,float *galaxy,int imagesize)
{printf("read_fits_2D,%s\n",argv );
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
}

int write_fits_2D(char *argv,float **stars,int *dim)
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
