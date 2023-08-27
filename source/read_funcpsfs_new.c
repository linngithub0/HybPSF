#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"
#include <string.h>
#include "dirent.h" 
void funcpsfs_size(char *fname,int *Nf,int *khead)
{
  char ftmp[500]={"#"};
  char comp[2]={"#"};
  FILE *fp;
  int k,px,py,i,j;
  float beta,fwhm,e1,e2;
  fp=fopen(fname,"r");
  k=0;
  do{
    if(fgets(ftmp,500,fp)!=NULL){ 
    k++;}
  }while(comp[0]==ftmp[0]);
  *khead=k-1;i=k-1;
  fclose(fp);
  fp=fopen(fname,"r");
  for(k=0;k<i;k++)(fgets(ftmp,500,fp));
  k=0;
  do{
    if(fgets(ftmp,500,fp)!=NULL){ 
    k++;}
  }while(!feof(fp));
  *Nf=k;
  //printf("Nf=%d\n",k );
}

void funcpsfs_read(char *fname,int Nf,int khead,float **para,int npara,char *delim){ 
  void exchange(float *a, float *b, int n);
  void read_character(FILE *fp,float *para,int colum);
  float ran1(long *idum);
  char ftmp[500]={"#"},*sget;
  char comp[2]={"#"};
  FILE *fp;
  int k,i,j;
  long iseed;
  //float beta,fwhm,e1,e2,**test;test=dmatrix(0,Nf,0,5);
  fp=fopen(fname,"r");
  for(k=0;k<khead;k++)(fgets(ftmp,500,fp));
  for(k=0;k<Nf;k++){//printf("k=%d\n",k );
    i=0;
    if(fgets(ftmp,500,fp)!=NULL){
      //printf("%s\n",ftmp );
      for(i=0;i<npara;i++){
        if(i==0){
          sget=strtok(ftmp,delim);
          para[k][i]=atof(sget);//printf("%f ", para[k][i]);
        }
        else{
          sget=strtok(NULL,delim);//printf("%s ", sget);
          para[k][i]=atof(sget);
        }
      }//printf("\n");
    }
    /*for(i=0;i<npara;i++){
      if(i==0){
        sget=strtok(ftmp," ");//printf("%s ", sget);
        para[k][i]=atof(sget);
      }
      else{
        sget=strtok(NULL," ");//printf("%s ", sget);
        para[k][i]=atof(sget);
      }
    }//printf("\n");
    */
  }
  for(k=0;k<Nf;k++){
    j=ran1(&iseed)*(Nf+1.);
    if(j>=Nf || j==k)    
      j=ran1(&iseed)*Nf;
    if(j>=Nf || j==k)
      j=ran1(&iseed)*Nf;
    exchange(para[k],para[j],npara);
  }

}


void exchange(float *a, float *b, int n)
{
  float tmp;
  int i;
  for(i=0;i<n;i++){tmp=a[i];a[i]=b[i];b[i]=tmp;}
}



/*****************************************************************************/
/*********************read file names for a given path in C*******************/
/*****************************************************************************/

void read_files(char *fpath,char *fname,int *f_num,char *appendix,int charnum){
  //  fpath      ;the given path
  //  fname      ;return the name in a string with f_num*charnum
  //  f_num      ;return the number of specified type files
  //  appendix   ;the specified file type
  //  charnum    ;the lenth for a single file name
  int i,j,k,l,m,n;
  char name[10000][charnum],a,tmp[charnum],appd[charnum];
  int fsize=0;
  DIR *dir=NULL;
  struct dirent *entry;//printf("%s\n",fpath );
  if((dir=opendir(fpath))==NULL){
    printf("open failed,path error %s\n",fpath);
    exit(EXIT_FAILURE);
  }
  else{
    while((entry=readdir(dir))!=NULL){
      //printf("%s\n",dir );
      //printf("filename = %s\t",entry->d_name);
      //printf("filetype = %d\n",entry->d_type);
      if(entry->d_type==DT_REG){
        strcpy(name[fsize],entry->d_name);
        fsize++;
      }
    }
  }
  closedir(dir);
  if(fsize>10000){
    printf("too many files\n");
    exit(EXIT_FAILURE);
  }//printf("fsize=%d\n",fsize );
  m=0;
  for(i=0;i<fsize;i++){
    k=strlen(name[i]);//printf("appendix=%d\n",strlen(appendix) );
    j=k;
    if(k>strlen(appendix)){
    for(n=0;n<strlen(appendix);n++){
      appd[n]=name[i][k-strlen(appendix)+n];
    }appd[n]='\0';//suffix the end char '\0'
    strcpy(tmp,appd);//printf("tmp=%s\n",tmp );
    if(strcmp(tmp,appendix)==0){
      //printf("appendix=%s\n",name[i] );
      for(l=0;l<charnum;l++){
        fname[m*charnum+l]=name[i][l];
      }
      m++;
    }
  }}
  *f_num=m;//printf("m=%d\n", m);
}

void read_list(char *fpath,char *fname,int *f_num,char *appendix,int charnum){
  //  fpath      ;the given path  for a text file
  //  fname      ;return the name in a string with f_num*charnum
  //  f_num      ;return the number of specified type files
  //  appendix   ;the specified file type
  //  charnum    ;the lenth for a single file name
  int i,j,k,l,m,n;
  char name[10000][charnum],a,tmp[charnum],appd[charnum];
  int fsize=0,Nstar;
  FILE *fp;
  if((fp=fopen(fpath,"r"))==NULL){
    printf("open failed,path error %s\n",fpath);
    exit(EXIT_FAILURE);
  }
  else{
    i=0;
    do{
      if(fgets(tmp,charnum,fp)!=NULL){
      strcpy(name[i],tmp);
      i++;}
      else{printf("get wrong in file %s",fpath);}
    }while(!feof(fp));
    i--;*f_num=i;
    Nstar=i;
    for(i=0;i<Nstar;i++){
      m=strlen(name[i]);
      j=m;
      do{a=name[i][j];
        j--;
      }while(a!='/');
      j+=2;
      k=67-1-j-strlen(appendix);//67 is the colum of strings attentions!!!
      for(n=0;n<k;n++){
        appd[n]=name[i][j+n];
      }
      for(n=0;n<charnum;n++){
        fname[i*charnum+n]=appd[n];
      }
      strcpy(tmp,appd);
    }
  }
}

/*
    DT_REG
        A regular file.

        常规文件


    DT_DIR
        A directory.

        目录


    DT_FIFO
        A named pipe, or FIFO. See FIFO Special Files.

        一个命名管道，或FIFO。


    DT_SOCK
        A local-domain socket.

        套接字


    DT_CHR
        A character device.

         字符设备


    DT_BLK
        A block device.

         块设备


    DT_LNK
        A symbolic link. 

         符号链接

具体数值 
enum
{ 
    DT_UNKNOWN = 0, 
 # define DT_UNKNOWN DT_UNKNOWN 
     DT_FIFO = 1, 
 # define DT_FIFO DT_FIFO 
     DT_CHR = 2, 
 # define DT_CHR DT_CHR 
     DT_DIR = 4, 
 # define DT_DIR DT_DIR 
     DT_BLK = 6, 
 # define DT_BLK DT_BLK 
     DT_REG = 8, 
 # define DT_REG DT_REG 
     DT_LNK = 10, 
 # define DT_LNK DT_LNK 
     DT_SOCK = 12, 
 # define DT_SOCK DT_SOCK 
     DT_WHT = 14 
 # define DT_WHT DT_WHT 
}; 
*/