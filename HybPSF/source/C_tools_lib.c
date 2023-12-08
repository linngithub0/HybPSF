#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fitsio.h"

void centriod_psflib(float *image0,int nx,int ny,double *cent){
    void star_gaus(float **y,int nx,int ny,float **yt,float **residu,double *para);
    int i,j,k;
    float **image,**starf,**residu;
    double *para;
    starf=(float **)calloc(nx,sizeof(float *));
    residu=(float **)calloc(nx,sizeof(float *));
    image=(float **)calloc(nx,sizeof(float *));
    para=(double *)calloc(10,sizeof(double));
    for(i=0;i<nx;i++){
        image[i]=(float *)calloc(ny,sizeof(float));
        starf[i]=(float *)calloc(ny,sizeof(float));
        residu[i]=(float *)calloc(ny,sizeof(float));
        for(j=0;j<ny;j++){
            k=i*ny+j;
            image[i][j]=image0[k];
        }
    }
    star_gaus(image,nx,ny,starf,residu,para);
    cent[0]=para[3];cent[1]=para[4];cent[2]=para[1];
    for(i=0;i<nx;i++){
        free(image[i]);
        free(starf[i]);
        free(residu[i]);
    }free(image);free(starf);free(residu);
    free(para);
}
void centriod_galaxylib(float *image0,int nx,int ny,double *cent){
    void galaxy_gaus(float **y,int nx,int ny,float **yt,float **residu,double *para);
    int i,j,k;
    float **image,**starf,**residu;
    double *para;
    starf=(float **)calloc(nx,sizeof(float *));
    residu=(float **)calloc(nx,sizeof(float *));
    image=(float **)calloc(nx,sizeof(float *));
    para=(double *)calloc(10,sizeof(double));
    for(i=0;i<nx;i++){
        image[i]=(float *)calloc(ny,sizeof(float));
        starf[i]=(float *)calloc(ny,sizeof(float));
        residu[i]=(float *)calloc(ny,sizeof(float));
        for(j=0;j<ny;j++){
            k=i*ny+j;
            image[i][j]=image0[k];
        }
    }
    galaxy_gaus(image,nx,ny,starf,residu,para); 
    cent[0]=para[3];cent[1]=para[4];
    for(i=0;i<nx;i++){
        free(image[i]);
        free(starf[i]);
        free(residu[i]);
    }free(image);free(starf);free(residu);
    free(para);
}


void centriod_burstlib(float *image0,int nx,int ny,double *cent){
    void burst_gaus(float **y,int nx,int ny,float **yt,float **residu,double *para);
    int i,j,k;
    float **image,**starf,**residu;
    double *para;
    starf=(float **)calloc(nx,sizeof(float *));
    residu=(float **)calloc(nx,sizeof(float *));
    image=(float **)calloc(nx,sizeof(float *));
    para=(double *)calloc(10,sizeof(double));
    for(i=0;i<nx;i++){
        image[i]=(float *)calloc(ny,sizeof(float));
        starf[i]=(float *)calloc(ny,sizeof(float));
        residu[i]=(float *)calloc(ny,sizeof(float));
        for(j=0;j<ny;j++){
            k=i*ny+j;
            image[i][j]=image0[k];
        }
    }
    burst_gaus(image,nx,ny,starf,residu,para); 
    cent[0]=para[3];cent[1]=para[4];
    for(i=0;i<nx;i++){
        free(image[i]);
        free(starf[i]);
        free(residu[i]);
    }free(image);free(starf);free(residu);
    free(para);
}


void Interp_bicubiclib(float *image0,int Ng1,int Ng2,
	int nx,int ny,double *cent,float *returns){
	void Interp_bicubic(int nx,int ny,float **a0,int nbound,
		    int nxt,int nyt, float **at,double *xcen);
	int i,j,k,nbond=(Ng1-nx)/2+1;
	float **image,**interpolated,**residu;
	image=(float **)calloc(Ng1,sizeof(float *));
	for(i=0;i<Ng1;i++){
		image[i]=(float *)calloc(Ng2,sizeof(float));
	}
	interpolated=(float **)calloc(nx,sizeof(float *));
	for(i=0;i<nx;i++){
		interpolated[i]=(float *)calloc(ny,sizeof(float));
	}
	for(i=0;i<Ng1;i++){
		for(j=0;j<Ng2;j++){
			k=i*Ng2+j;
			image[i][j]=image0[k];
		}
	}
  Interp_bicubic(Ng1,Ng2,image,nbond,nx,ny,interpolated,cent);
	for(i=0;i<nx;i++){
		for(j=0;j<ny;j++){
			k=i*ny+j;
			returns[k]=interpolated[i][j];
		}
	}
	for(i=0;i<Ng1;i++){
		free(image[i]);
	}free(image);
	for(i=0;i<nx;i++){
		free(interpolated[i]);
	}free(interpolated);
}

void BInterp_bicubiclib(float *image0,int Ng1,int Ng2,
	int nx,int ny,double *cent,float *returns){
	void BInterp_bicubic(int nx,int ny,float **a0,int nbound,
		    int nxt,int nyt, float **at,double *xcen);
	int i,j,k,nbond=(Ng1-nx)/2+1;
	float **image,**interpolated;
	image=(float **)calloc(Ng1,sizeof(float *));
	for(i=0;i<Ng1;i++){
		image[i]=(float *)calloc(Ng2,sizeof(float));
	}
	interpolated=(float **)calloc(nx,sizeof(float *));
	for(i=0;i<nx;i++){
		interpolated[i]=(float *)calloc(ny,sizeof(float));
	}
	for(i=0;i<Ng1;i++){
		for(j=0;j<Ng2;j++){
			k=i*Ng2+j;
			image[i][j]=image0[k];
		}
	}
  BInterp_bicubic(Ng1,Ng2,image,nbond,nx,ny,interpolated,cent);
	for(i=0;i<nx;i++){
		for(j=0;j<ny;j++){
			k=i*ny+j;
			returns[k]=interpolated[i][j];
		}
	}
	for(i=0;i<Ng1;i++){
		free(image[i]);
	}free(image);
	for(i=0;i<nx;i++){
		free(interpolated[i]);
	}free(interpolated);
}


void deblending(float *in_image,float *out_image,int Ngx,int Ngy,float sigma,int *blendings){
	void fof_image(float **image,double sigma,int nimgx,int nimgy,int *groupnum,
		int *tails,int **pgroup);
	void deblend(float **image,int nimgx,int nimgy,int ngroups,int *tails,
    	int **pgroup,int *nsub,float sigma);
	void deblend1(float **image,int nimgx,int nimgy,int ngroups,int *tails,
    	int **pgroup,int *nsub,float sigma);
	void mask_groups(int count,int ngroup,int *tails,int **pgroup,float **image,
    	double sigma,int nimgx,int nimgy,float **mask);
	int groupnum,*tails,**pgroups,i,j,k,ngroup,nsub1,nsub2;
	float **stamp,**mask;
	stamp=(float **)calloc(Ngx,sizeof(float *));
	mask=(float **)calloc(Ngx,sizeof(float *));
	tails=(int *)calloc(50,sizeof(int));
	pgroups=(int **)calloc(50,sizeof(float *));
	for(i=0;i<Ngx;i++){
		stamp[i]=(float *)calloc(Ngy,sizeof(float));
		mask[i]=(float *)calloc(Ngy,sizeof(float));
		for(j=0;j<Ngy;j++){
			mask[i][j]=1;
		}
	}
	for(i=0;i<50;i++){
		pgroups[i]=(int *)calloc(Ngx*Ngy,sizeof(int));
	}
	for(i=0;i<Ngx;i++){
		for(j=0;j<Ngy;j++){
			k=i*Ngy+j;
			stamp[i][j]=in_image[k];
		}
	}
	fof_image(stamp,sigma,Ngx,Ngy,&ngroup,tails,pgroups);//printf("ngroup=%d\n",ngroup );
	if(ngroup>=1){
		mask_groups(0,ngroup,tails,pgroups,stamp,sigma,Ngx,Ngy,mask);
		fof_image(stamp,sigma,Ngx,Ngy,&ngroup,tails,pgroups);
		//mask_groups(0,ngroup,tails,pgroups,stamp,sigma,Ngx,Ngy,mask);
		//fof_image(stamp,sigma*3,Ngx,Ngy,&ngroup,tails,pgroups);
		//mask_groups(0,ngroup,tails,pgroups,stamp,sigma,Ngx,Ngy,mask);
		//fof_image(stamp,sigma,Ngx,Ngy,&ngroup,tails,pgroups);
		//deblend(stamp,Ngx,Ngy,ngroup,tails,pgroups,&nsub1,sigma);
		//deblend1(stamp,Ngx,Ngy,ngroup,tails,pgroups,&nsub2,sigma);
	}
	//else{
	//	ngroup=1;
	//}
	
	blendings[0]=ngroup;blendings[1]=nsub1;blendings[2]=nsub2;
	for(i=0;i<Ngx;i++){
		for(j=0;j<Ngy;j++){
			k=i*Ngy+j;
			out_image[k]=stamp[i][j];
		}
	}

}



















