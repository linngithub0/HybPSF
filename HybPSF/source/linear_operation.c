
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"
#include <string.h>
void matrix_matrix(double **a,int m, int n,double **b,int k, int l,double **c)
{
  double sum;
  int i,j,kt;
  if(n!=k){printf( "matrix production error n!=k stop");getchar();}
  for(i=0;i<m;i++){
    for(j=0;j<l;j++){
      sum=0;
      for(kt=0;kt<n;kt++)sum+=a[i][kt]*b[kt][j];
      c[i][j]=sum;
    }
  }
}
void matrix_vector(double **a,int m, int n,double *b,int k,double *c)
{
  double sum;
  int i,j,kt;
  if(n!=k){printf( "matrix_vector production error n!=k stop");getchar();}
  for(i=0;i<m;i++){
      sum=0;
      for(kt=0;kt<n;kt++)sum+=a[i][kt]*b[kt];
      c[i]=sum;
  }
}
void vector_matrix(double *a,int n,double **b,int k, int l,double *c)
{
  double sum;
  int i,j,kt;
  if(n!=k){printf( "vector_matrix production error n!=k stop");getchar();}
  for(j=0;j<l;j++){
    sum=0;
    for(kt=0;kt<n;kt++)sum+=a[kt]*b[kt][j];
    c[j]=sum;
  }
}
void vector_vector(double *a,int n,double *b,int k,double *c)
{
  double sum;
  int i,j,kt;
  if(n!=k){printf( "vector_vector production error n!=k stop");getchar();}
  sum=0;
  for(kt=0;kt<n;kt++)sum+=a[kt]*b[kt];
  *c=sum;
}
void matrix_transpose(double **a,int m, int n, double **b){
  int i,j;
  for(i=0;i<m;i++)for(j=0;j<n;j++)b[j][i]=a[i][j];
}


void matrix_invers(double **a,int N,double **b){
  void ludcmp(double **a, int n, int *indx, double *d);
  void lubksb(double **a, int n, int *indx, double b[]);
  double d,*col,**atmp;
  int i,j,*indx;
  atmp=dmatrix(0,N,0,N);
  col=dvector(0,N);
  indx=ivector(0,N);
  for(j=0;j<N;j++)for(i=0;i<N;i++){atmp[i][j]=a[i][j];}
  ludcmp(atmp,N,indx,&d);
  for(j=0;j<N;j++){
    for(i=0;i<N;i++)col[i]=0;
    col[j]=1.0;
    lubksb(atmp,N,indx,col);
    for(i=0;i<N;i++)b[i][j]=col[i];
  }
  free_dvector(col,0,N);
  free_ivector(indx,0,N);
  free_dmatrix(atmp,0,N,0,N);
}
double matrix_determinant(double **a, int N)
{
 void ludcmp(double **a, int n, int *indx, double *d);
 double d=1,**atmp;
 int i,j,*indx;
  atmp=dmatrix(0,N,0,N);
  indx=ivector(0,N);
  for(j=0;j<N;j++)for(i=0;i<N;i++){atmp[i][j]=a[i][j];}
  ludcmp(atmp,N,indx,&d);
  for(j=0;j<N;j++){
    d*=a[j][j];
  }
  free_ivector(indx,0,N);
  free_dmatrix(atmp,0,N,0,N);
  return d;
}
