#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nrutil.h"
#include "function.h"
void npsf(int order_PC,int order_poly,int Nstar,float **pos,float **Coeff,
  float **weight,double **dcoeff,float ***PCs,int ngp){
  void creat_matrix_coeff(double **a,double *b,int order_PC,int l0,int Nstar,
    double **poly,float **Coeff,float **weight,float ***PCs,int ngp);
  void lubksb(double **a, int n, int *indx, double b[]);
  void ludcmp(double **a, int n, int *indx, double *d);

  int i,j,k,l,f,ic,m,n,*indx,l0=0,t;
  double **a,*b,d,**poly,x,y,ix,iy,sum;
  for(i=order_poly+1;i>0;i--)l0+=i;
  m=order_PC*l0;
  a=dmatrix(0,m,0,m);
  b=dvector(0,m);
  indx=ivector(0,m);
  poly=dmatrix(0,l0,0,Nstar);
  FILE *fp;
  
  for(i=0;i<m;i++){
    b[i]=0;
    for(j=0;j<m;j++){
      a[i][j]=0;
    }
  }
  //creat polunomials
    for(ic=0;ic<Nstar;ic++){
      x=pos[ic][0];y=pos[ic][1];k=0;
      //printf("%e\t%e\t%d\n",x,y,ic );
      for(l=0;l<=order_poly;l++){
        for(f=0;f<=order_poly;f++){
          ix=l;iy=f;
          if(l+f<=order_poly){
            poly[k][ic]=pow(x,ix)*pow(y,iy);k++;
          }
        }
      }
    }
  printf("creat_matrix_coeff\n");
  creat_matrix_coeff(a,b,order_PC,l0,Nstar,poly,Coeff,weight,PCs,ngp);
  ludcmp(a, m, indx, &d);
  lubksb(a, m, indx, b);
    for(l=0;l<order_PC;l++){
      for(f=0;f<l0;f++){
        k=l*l0+f;
        dcoeff[l][f]=b[k];
      }
    }

  free_dmatrix(a,0,m,0,m);
  free_dvector(b,0,m);
  free_ivector(indx,0,m);
  free_dmatrix(poly,0,l0,0,Nstar);

}
void creat_matrix_coeff(double **a,double *b,int order_PC,int l0,int Nstar,
  double **poly,float **Coeff,float **weight,float ***PCs,int ngp){

    int get_fits_size(char *argv,int *dim);
  int write_fits_2D( char *argv,float **stars,int *dim);
  int write_fits_3D(char *argv,float ***stars,int *dim);
  int read_fits_2D(char *argv,float *galaxy,int imagesize);

    int i,j,k,l,f,m,n,ic,r,q,dim[3];
    m=order_PC*l0;
    for(n=0;n<order_PC;n++){
        for(q=0;q<l0;q++){
          i=n*l0+q;
          for(f=0;f<l0;f++){
            j=n*l0+f;
            for(ic=0;ic<Nstar;ic++){
              a[i][j]+=poly[q][ic]*poly[f][ic]*weight[ic][n];
            }
          }
        }
      }

    char *fname;
    FILE *fp;
    float **aa;
    r=order_PC*l0;
    aa=matrix(0,r,0,r);
    for(i=0;i<r;i++)for(j=0;j<r;j++)aa[i][j]=a[i][j];
    fname="fits/a.fits";
    dim[0]=dim[1]=order_PC*l0;
    //write_fits_2D(fname,aa,dim);
      for(n=0;n<order_PC;n++){
        for(q=0;q<l0;q++){
          i=n*l0+q;
          for(ic=0;ic<Nstar;ic++){
            b[i]+=poly[q][ic]*Coeff[ic][n]*weight[ic][n];
          }
        }
      }

}


void psf_interp(int Nobj,int order_poly,int order_PC,int ngp,float ***rPSFs,
  float **gpos,double **dcoeff,float ***PCs){
  char* itoa(int num,char *str,int radix);

  int i,j,t,ic,l,f,l0=0,k;
  for(i=order_poly+1;i>0;i--)l0+=i;printf("l=%d in psf_interp\n",l0 );
  double **xy,**Coeff,x,y,ix,iy,sum;
  FILE *fp;
  char fname[200],name[200]; 
  xy=dmatrix(0,Nobj,0,l0);
  Coeff=dmatrix(0,Nobj,0,order_PC);

      for(ic=0;ic<Nobj;ic++){
        x=gpos[ic][0];y=gpos[ic][1];
        k=0;
        for(l=0;l<=order_poly;l++){
          for(f=0;f<=order_poly;f++){
            ix=l;iy=f;
            if(l+f<=order_poly){
              xy[ic][k]=pow(x,ix)*pow(y,iy);
              k++;
            }
          }
        }
      }printf("k=%d\n",k );
    

      for(ic=0;ic<Nobj;ic++){
        for(l=0;l<order_PC;l++){
          Coeff[ic][l]=0;
          for(f=0;f<l0;f++)Coeff[ic][l]+=dcoeff[l][f]*xy[ic][f];
        }
      }printf("Coeff\n");
    
    
    for(ic=0;ic<Nobj;ic++){
      sum=0;
      for(i=0;i<ngp;i++)for(j=0;j<ngp;j++){
        rPSFs[ic][i][j]=0;
        for(l=0;l<order_PC;l++){
          rPSFs[ic][i][j]+=Coeff[ic][l]*PCs[l][i][j];//construct PSFs
        }sum+=rPSFs[ic][i][j];
      }//printf("sum=%e\n",sum );
      for(i=0;i<ngp;i++)for(j=0;j<ngp;j++){
        rPSFs[ic][i][j]/=sum;
      }
    }printf("rPSFs\n");

    free_dmatrix(xy,0,Nobj,0,l0);
    free_dmatrix(Coeff,0,Nobj,0,order_PC);
}
/*
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nrutil.h"
void npsf(int order_PC,int order_poly,int Nstar,float **pos,float **Coeff,
  float **weight,double **dcoeff,float ***PCs,int ngp){
  void creat_matrix_coeff(double **a,double *b,int order_PC,int l0,int Nstar,
    double **poly,float **Coeff,float **weight,float ***PCs,int ngp);
  void lubksb(double **a, int n, int *indx, double b[]);
  void ludcmp(double **a, int n, int *indx, double *d);

  int i,j,k,l,f,ic,m,n,*indx,l0=0,t;
  double **a,*b,d,**poly,x,y,ix,iy,sum,**aa,*c,*bb;
  float **wa,**wb;
  for(i=order_poly+1;i>0;i--)l0+=i;printf("l0=%d\n",l0 );
  m=order_PC*l0+2;
  a=dmatrix(0,m,0,m);wa=matrix(0,m,0,m);wb=matrix(0,1,0,m);
  b=dvector(0,m);
  indx=ivector(0,m);
  poly=dmatrix(0,l0,0,Nstar);
  FILE *fp;
  for(i=0;i<m;i++){
    b[i]=0;
    for(j=0;j<m;j++){
      a[i][j]=0;
    }
  }
  //creat polunomials
    for(ic=0;ic<Nstar;ic++){
      x=pos[ic][0];y=pos[ic][1];k=0;
      //printf("%e\t%e\t%d\n",x,y,ic );
      for(l=0;l<=order_poly;l++){
        for(f=0;f<=order_poly;f++){
          ix=l;iy=f;
          if(l+f<=order_poly){
            poly[k][ic]=pow(x,ix)*pow(y,iy);
            k++;
          }
        }
      }
    }printf("\n");
  printf("k=%d\n",k );
  creat_matrix_coeff(a,b,order_PC,l0,Nstar,poly,Coeff,weight,PCs,ngp);

  for(i=0;i<m;i++){
    for(j=0;j<m;j++){
      wa[i][j]=a[i][j];
    }wb[0][i]=b[i];
  }
  ludcmp(a, m, indx, &d);
  lubksb(a, m, indx, b);



    for(l=0;l<order_PC;l++){
      for(f=0;f<l0;f++){
        k=l*l0+f;
        dcoeff[l][f]=b[k];
      }
    }

    char *fname;int dim[5];
    int write_fits_2D( char *argv,float **stars,int *dim);
    fname="fits/a.fits";
    dim[0]=dim[1]=m;
    //write_fits_2D(fname,wa,dim);
    fname="fits/b.fits";
    dim[0]=1;dim[1]=m;
    //write_fits_2D(fname,wb,dim);


  free_dmatrix(a,0,m,0,m);
  free_dvector(b,0,m);
  free_ivector(indx,0,m);
  free_dmatrix(poly,0,l0,0,Nstar);
  free_matrix(0,m,0,m);
  free_matrix(0,1,0,m);

}
void creat_matrix_coeff(double **a,double *b,int order_PC,int l0,int Nstar,
  double **poly,float **Coeff,float **weight,float ***PCs,int ngp){

    int i,j,k,l,f,t,m,n,q,ic,r,s,dim[3];

      for(n=0;n<order_PC;n++){
        for(q=0;q<l0;q++){
          i=n*l0+q;
          for(f=0;f<l0;f++){
            j=n*l0+f;
            for(ic=0;ic<Nstar;ic++){
              a[i][j]+=2*poly[q][ic]*poly[f][ic]*weight[ic][n];
            }
          }
        }
      }//left uppper
      printf("left uppper\n");

      for(r=0;r<2;r++){
        j=order_PC*l0+r;
        if(r==0){
          for(n=0;n<order_PC;n++){
            for(q=0;q<1;q++){
              i=n*l0+q;
              for(l=0;l<ngp;l++)for(f=0;f<ngp;f++){
                a[i][j]+=PCs[n][l][f];
              }
            }
          }
        }
        else{
          for(n=0;n<order_PC;n++){
            for(q=1;q<l0;q++){
              i=n*l0+q;
              for(l=0;l<ngp;l++)for(f=0;f<ngp;f++){
                a[i][j]+=PCs[n][l][f];
              }
            }
          }
        }
      }//upper right
      printf("upper right\n");

      for(n=0;n<order_PC;n++){
        for(q=0;q<l0;q++){
          i=n*l0+q;
          for(ic=0;ic<Nstar;ic++){
            b[i]+=2*poly[q][ic]*Coeff[ic][n]*weight[ic][n];
          }
        }
      }printf("upper b\n");

      for(q=0;q<2;q++){
        i=order_PC*l0+q;
        if(q==0){
          for(n=0;n<order_PC;n++)for(f=0;f<1;f++){
            j=n*l0+f;
            for(r=0;r<ngp;r++)for(s=0;s<ngp;s++){
              a[i][j]+=PCs[n][r][s];
            }
          }
        }
        else{
          for(n=0;n<order_PC;n++)for(f=1;f<l0;f++){
            j=n*l0+f;
            for(r=0;r<ngp;r++)for(s=0;s<ngp;s++){
              a[i][j]+=PCs[n][r][s];
            }
          }
        }
      }//bottom left
      printf("bottom left\n");
      for(q=0;q<2;q++){
        k=order_PC*l0+q;
        if(q==0)b[k]=1.;
        else{
          b[k]=0;
        }
      }printf("down\n");

}


void psf_interp(int Nobj,int order_poly,int order_PC,int ngp,float ***rPSFs,
  float **gpos,double **dcoeff,float ***PCs){
  char* itoa(int num,char *str,int radix);

  int i,j,t,ic,l,f,l0=0,k;
  for(i=order_poly+1;i>0;i--)l0+=i;printf("l=%d in psf_interp\n",l0 );
  double **xy,**Coeff,x,y,ix,iy,sum;
  FILE *fp;
  char fname[200],name[200]; 
  xy=dmatrix(0,Nobj,0,l0);
  Coeff=dmatrix(0,Nobj,0,order_PC);

      for(ic=0;ic<Nobj;ic++){
        x=gpos[ic][0];y=gpos[ic][1];
        k=0;
        for(l=0;l<=order_poly;l++){
          for(f=0;f<=order_poly;f++){
            ix=l;iy=f;
            if(l+f<=order_poly){
              xy[ic][k]=pow(x,ix)*pow(y,iy);
              k++;
            }
          }
        }
      }printf("k=%d\n",k );
    

      for(ic=0;ic<Nobj;ic++){
        for(l=0;l<order_PC;l++){
          Coeff[ic][l]=0;
          for(f=0;f<l0;f++)Coeff[ic][l]+=dcoeff[l][f]*xy[ic][f];
        }
      }printf("Coeff\n");
    
    
    for(ic=0;ic<Nobj;ic++){
      sum=0;
      for(i=0;i<ngp;i++)for(j=0;j<ngp;j++){
        rPSFs[ic][i][j]=0;
        for(l=0;l<order_PC;l++){
          rPSFs[ic][i][j]+=Coeff[ic][l]*PCs[l][i][j];//construct PSFs
        }sum+=rPSFs[ic][i][j];
      }printf("sum=%e\n",sum );
      for(i=0;i<ngp;i++)for(j=0;j<ngp;j++){
        //rPSFs[ic][i][j]/=sum;
      }
    }printf("rPSFs\n");

    free_dmatrix(xy,0,Nobj,0,l0);
    free_dmatrix(Coeff,0,Nobj,0,order_PC);
}*/
/*
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nrutil.h"
#include "function.h"
void npsf(int order_PC,int order_poly,int Nstar,float **pos,float **Coeff,
  float **weight,double **dcoeff,float ***PCs,int ngp){
  void creat_matrix_coeff(double **a,double *b,int order_PC,int l0,int Nstar,
    double **poly,float **Coeff,float **weight,float ***PCs,int ngp);
  void lubksb(double **a, int n, int *indx, double b[]);
  void ludcmp(double **a, int n, int *indx, double *d);
  void matrix_vector(double **a,int m, int n,double *b,int k,double *c);

  int i,j,k,l,f,ic,m,n,*indx,l0=0,t,dim[5];
  double **a,*b,d,**poly,x,y,ix,iy,sum;
  float **wa,*c;
  for(i=order_poly+1;i>0;i--)l0+=i;printf("l0=%d\n",l0 );
  m=order_PC*l0+l0;
  a=dmatrix(0,m,0,m);wa=matrix(0,m,0,m);
  b=dvector(0,m);c=vector(0,m);
  indx=ivector(0,m);
  poly=dmatrix(0,l0,0,Nstar);
  FILE *fp;
  for(i=0;i<m;i++){
    b[i]=0;
    for(j=0;j<m;j++){
      a[i][j]=0;
    }
  }
  //creat polunomials
    for(ic=0;ic<Nstar;ic++){
      x=pos[ic][0];y=pos[ic][1];k=0;
      //printf("%e\t%e\t%d\n",x,y,ic );
      for(l=0;l<=order_poly;l++){
        for(f=0;f<=order_poly;f++){
          ix=l;iy=f;
          if(l+f<=order_poly){
            poly[k][ic]=pow(x,ix)*pow(y,iy);
            k++;
          }
        }
      }
    }printf("\n");
  printf("k=%d\n",k );
  creat_matrix_coeff(a,b,order_PC,l0,Nstar,poly,Coeff,weight,PCs,ngp);
  printf("creat_matrix_coeff\n");
  for(i=0;i<m;i++){
    for(j=0;j<m;j++){
      wa[i][j]=a[i][j];
    }
  }
  dim[0]=dim[1]=m;
  //write_fits_2D("fits/a.fits",wa,dim);
  for(i=0;i<1;i++){
    for(j=0;j<m;j++){
      wa[i][j]=b[j];
    }
  }
  dim[0]=1;dim[1]=m;
  //write_fits_2D("fits/b.fits",wa,dim);

  ludcmp(a, m, indx, &d);
  lubksb(a, m, indx, b);

  



    for(l=0;l<order_PC;l++){
      for(f=0;f<l0;f++){
        k=l*l0+f;
        dcoeff[l][f]=b[k];
      }
    }
    /*for(l=0;l<order_PC;l++){
      for(f=0;f<l0;f++){
        k=l*l0+f;
        printf("%f\t",dcoeff[l][f] );
      }printf("\n");
    }*/


  /*free_dmatrix(a,0,m,0,m);
  free_dvector(b,0,m);
  free_ivector(indx,0,m);
  free_dmatrix(poly,0,l0,0,Nstar);
  free_matrix(wa,0,m,0,m);
  free_vector(c,0,m);

}
void creat_matrix_coeff(double **a,double *b,int order_PC,int l0,int Nstar,
  double **poly,float **Coeff,float **weight,float ***PCs,int ngp){

    int i,j,k,l,f,t,m,n,q,ic,r,s;

      for(n=0;n<order_PC;n++){
        for(q=0;q<l0;q++){
          i=n*l0+q;
          for(f=0;f<l0;f++){
            j=n*l0+f;
            for(ic=0;ic<Nstar;ic++){
              a[i][j]+=2*poly[q][ic]*poly[f][ic]*weight[ic][n];
            }
          }
        }
      }
      for(n=0;n<order_PC;n++){
        for(q=0;q<l0;q++){
          i=n*l0+q;
          j=order_PC*l0+q;
          for(l=0;l<ngp;l++)for(f=0;f<ngp;f++){
            a[i][j]+=PCs[n][l][f];
          }
        }
      }

      for(n=0;n<order_PC;n++){
        for(q=0;q<l0;q++){
          i=n*l0+q;
          for(ic=0;ic<Nstar;ic++){
            b[i]+=2*poly[q][ic]*Coeff[ic][n]*weight[ic][n];
          }
        }
      }
      for(q=0;q<l0;q++){
        i=order_PC*l0+q;
        for(l=0;l<order_PC;l++){
          j=l*l0+q;
          for(r=0;r<ngp;r++)for(s=0;s<ngp;s++){
            a[i][j]+=PCs[l][r][s];
          }
        }
      }
      for(q=0;q<l0;q++){
        k=order_PC*l0+q;
        if(q==0)b[k]=1;
        else{
          b[k]=0;
        }
      }



}


void psf_interp(int Nobj,int order_poly,int order_PC,int ngp,float ***rPSFs,
  float **gpos,double **dcoeff,float ***PCs){
  char* itoa(int num,char *str,int radix);

  int i,j,t,ic,l,f,l0=0,k;
  for(i=order_poly+1;i>0;i--)l0+=i;printf("l=%d in psf_interp\n",l0 );
  double **xy,**Coeff,x,y,ix,iy,sum;
  FILE *fp;
  char fname[200],name[200]; 
  xy=dmatrix(0,Nobj,0,l0);
  Coeff=dmatrix(0,Nobj,0,order_PC);
  for(l=0;l<order_PC;l++){
      for(f=0;f<l0;f++){
        k=l*l0+f;
        //printf("%f\t",dcoeff[l][f] );
      }//printf("\n");
    }
    int dim[5];
    dim[0]=dim[1]=ngp;dim[2]=order_PC;
    //write_fits_3D("fits/get_PCs.fits",PCs,dim);

      for(ic=0;ic<Nobj;ic++){
        x=gpos[ic][0];y=gpos[ic][1];//printf("x=%f,y=%f\n",x,y );
        k=0;
        for(l=0;l<=order_poly;l++){
          for(f=0;f<=order_poly;f++){
            ix=l;iy=f;
            if(l+f<=order_poly){
              xy[ic][k]=pow(x,ix)*pow(y,iy);
              //printf("%f\t",xy[ic][k] );
              k++;
            }
          }
        }//printf("\n");
      }printf("k=%d\n",k );
    

      for(ic=0;ic<Nobj;ic++){
        for(l=0;l<order_PC;l++){
          Coeff[ic][l]=0;
          for(f=0;f<l0;f++)Coeff[ic][l]+=dcoeff[l][f]*xy[ic][f];
            //printf("%f\t",Coeff[ic][l] );
        }//printf("%d\n",ic);
      }printf("Coeff\n");
    
    
    for(ic=0;ic<Nobj;ic++){
      sum=0;
      for(i=0;i<ngp;i++)for(j=0;j<ngp;j++){
        rPSFs[ic][i][j]=0;
        for(l=0;l<order_PC;l++){
          rPSFs[ic][i][j]+=(float) (Coeff[ic][l]*PCs[l][i][j]);//construct PSFs
        }sum+=rPSFs[ic][i][j];
      }//printf("sum=%e\n",sum );
      for(i=0;i<ngp;i++)for(j=0;j<ngp;j++){
        rPSFs[ic][i][j]/=(float)sum;
      }
    }printf("rPSFs\n");

    free_dmatrix(xy,0,Nobj,0,l0);
    free_dmatrix(Coeff,0,Nobj,0,order_PC);
}*/
