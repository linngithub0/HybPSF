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
  for(k=0;k<i;k++)fgets(ftmp,500,fp);
  k=0;
  do{
    if(fgets(ftmp,500,fp)!=NULL){ 
    k++;}
  }while(!feof(fp));
  *Nf=k;
  //printf("Nf=%d\n",k );
}

void funcpsfs_read(char *fname,int Nf,int khead,float **para,int npara)
  {void exchange(float *a, float *b, int n);
    float ran1(long *idum);
  char ftmp[500]={"#"};
  char comp[2]={"#"};
  FILE *fp;
  int k,i,j;
  long iseed;
  //float beta,fwhm,e1,e2,**test;test=dmatrix(0,Nf,0,5);
  fp=fopen(fname,"r");
  for(k=0;k<khead;k++)fgets(ftmp,500,fp); 
  for(k=0;k<Nf;k++){
    switch(npara){
      case 1:  (fscanf(fp,"%e\n",&para[k][0]));break;
      case 2:  (fscanf(fp,"%e %e\n",&para[k][0],&para[k][1]));break;
      case 3:  (fscanf(fp,"%e %e %e\n",&para[k][0],&para[k][1],&para[k][2]));break;
      case 4:  (fscanf(fp,"%e %e %e %e\n",&para[k][0],&para[k][1],&para[k][2],&para[k][3]));break;
      case 5:  (fscanf(fp,"%e %e %e %e %e\n",&para[k][0],&para[k][1],&para[k][2],&para[k][3],&para[k][4]));break;
      case 6:  (fscanf(fp,"%e %e %e %e %e %e\n",&para[k][0],&para[k][1],&para[k][2],&para[k][3],&para[k][4],&para[k][5]));break;
      case 7:  (fscanf(fp,"%e %e %e %e %e %e %e\n",&para[k][0],&para[k][1],&para[k][2],&para[k][3],&para[k][4],&para[k][5],&para[k][6]));break;
      case 8:  (fscanf(fp,"%e %e %e %e %e %e %e %e\n",&para[k][0],&para[k][1],&para[k][2],&para[k][3],&para[k][4],&para[k][5],&para[k][6],&para[k][7]));break;
      case 9:  (fscanf(fp,"%e %e %e %e %e %e %e %e %e\n",&para[k][0],&para[k][1],&para[k][2],&para[k][3],&para[k][4],&para[k][5],&para[k][6],&para[k][7],&para[k][8]));break;
      case 10: (fscanf(fp,"%e %e %e %e %e %e %e %e %e %e\n",&para[k][0],&para[k][1],&para[k][2],&para[k][3],&para[k][4],&para[k][5],&para[k][6],&para[k][7],&para[k][8],&para[k][9]));break;
      case 11: (fscanf(fp,"%e %e %e %e %e %e %e %e %e %e %e\n",&para[k][0],&para[k][1],&para[k][2],&para[k][3],&para[k][4],&para[k][5],&para[k][6],&para[k][7],&para[k][8],&para[k][9],&para[k][10]));break;
      case 12: (fscanf(fp,"%e %e %e %e %e %e %e %e %e %e %e %e\n",&para[k][0],&para[k][1],&para[k][2],&para[k][3],&para[k][4],&para[k][5],&para[k][6],&para[k][7],&para[k][8],&para[k][9],&para[k][10],&para[k][11]));break;
      case 13: (fscanf(fp,"%e %e %e %e %e %e %e %e %e %e　%e %e %e\n",&para[k][0],&para[k][1],&para[k][2],&para[k][3],&para[k][4],&para[k][5],&para[k][6],&para[k][7],&para[k][8],&para[k][9],&para[k][10],&para[k][11],&para[k][12]));break;
      case 14: (fscanf(fp,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",&para[k][0],&para[k][1],&para[k][2],&para[k][3],&para[k][4],&para[k][5],&para[k][6],&para[k][7],&para[k][8],&para[k][9],&para[k][10],&para[k][11],&para[k][12],&para[k][13]));break;
      case 15: (fscanf(fp,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",&para[k][0],&para[k][1],&para[k][2],&para[k][3],&para[k][4],&para[k][5],&para[k][6],&para[k][7],&para[k][8],&para[k][9],&para[k][10],&para[k][11],&para[k][12],&para[k][13],&para[k][14]));break;
      case 16: (fscanf(fp,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",&para[k][0],&para[k][1],&para[k][2],&para[k][3],&para[k][4],&para[k][5],&para[k][6],&para[k][7],&para[k][8],&para[k][9],&para[k][10],&para[k][11],&para[k][12],&para[k][13],&para[k][14],&para[k][15]));break;
      case 17: (fscanf(fp,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",&para[k][0],&para[k][1],&para[k][2],&para[k][3],&para[k][4],&para[k][5],&para[k][6],&para[k][7],&para[k][8],&para[k][9],&para[k][10],&para[k][11],&para[k][12],&para[k][13],&para[k][14],&para[k][15],&para[k][16]));break;
      case 18: (fscanf(fp,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",&para[k][0],&para[k][1],&para[k][2],&para[k][3],&para[k][4],&para[k][5],&para[k][6],&para[k][7],&para[k][8],&para[k][9],&para[k][10],&para[k][11],&para[k][12],&para[k][13],&para[k][14],&para[k][15],&para[k][16],&para[k][17]));break;
      case 19: (fscanf(fp,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",&para[k][0],&para[k][1],&para[k][2],&para[k][3],&para[k][4],&para[k][5],&para[k][6],&para[k][7],&para[k][8],&para[k][9],&para[k][10],&para[k][11],&para[k][12],&para[k][13],&para[k][14],&para[k][15],&para[k][16],&para[k][17],&para[k][18]));break;
      case 20: (fscanf(fp,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",&para[k][0],&para[k][1],&para[k][2],&para[k][3],&para[k][4],&para[k][5],&para[k][6],&para[k][7],&para[k][8],&para[k][9],&para[k][10],&para[k][11],&para[k][12],&para[k][13],&para[k][14],&para[k][15],&para[k][16],&para[k][17],&para[k][18],&para[k][19]));break;
      default :{printf("npara=%d\tnpara cross the limits 20,change the input data strctures!\n",npara);exit(EXIT_FAILURE);}
    
    }
  }
  /*for(k=0;k<Nf;k++){
    //fgets(ftmp,500,fp);printf("%s\n",ftmp );
    fscanf(fp,"%e %e %e %e\n",&para[k][0],&para[k][1],&para[k][2],&para[k][3]);
  }*/

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
  struct dirent *entry;
  if((dir=opendir(fpath))==NULL){
    printf("open failed,path error %s\n",fpath);
    exit(EXIT_FAILURE);
  }
  else{
    while(entry=readdir(dir)){
      //printf("filename = %s",entry->d_name);printf("filetype = %d\n",entry->d_type);
      //if(entry->d_name[0]!='.' && entry->d_type==DT_REG){
        strcpy(name[fsize],entry->d_name);
        fsize++;
      //}
    }
  }
  closedir(dir);
  if(fsize>10000){
    printf("too many files\n");
    exit(EXIT_FAILURE);
  }
  m=0;
  for(i=0;i<fsize;i++){
    k=strlen(name[i]);
    j=k;
    if(k>strlen(appendix)){
    for(n=0;n<5;n++){
      appd[n]=name[i][k-5+n];
    }
    strcpy(tmp,appd);//printf("tmp=%s\n",tmp );
    if(strcmp(tmp,appendix)==0){//printf("appendix=%s\n",appd );
      for(l=0;l<charnum;l++){
        fname[m*charnum+l]=name[i][l];
      }
      m++;
    }
  }}
  *f_num=m;
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