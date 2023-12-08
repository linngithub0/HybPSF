
void iSPCA_entr(float ***stars,float **spos,int npc,int Nstar0,int Ng0,int Ng,
  float ***PCs,float ***coeff,float gain,int *Nobj,float snrs,int osam,float ***wgt);
void EMPCA_entr(float ***stars,float **spos,int npc,int Nstar0,int Ng0,int Ng,
  float ***PCs,float ***coeff,float gain,int *Nobj,float snrs,float ***wgt);
void SPCA_entr(float ***stars,float **spos,int npc,int Nstar0,int Ng0,int Ng,
  float ***PCs,float ***coeff,float gain,int *Nobj,float snrs,float ***wgt);
int get_fits_size(char *argv,int *dim);
int write_fits_2D(char *argv,float **stars,int *dim);
int read_fits_2D(char *argv,float *galaxy,int imagesize);
int write_fits_3D(char *argv,float ***stars,int *dim);
void funcpsfs_size(char *fname,int *Nf,int *Nhead);
void funcpsfs_read(char *fname,int Nf,int khead,float **para,int npara,char *delim);
char* itoa(int num,char *str,int radix);
int get_fits_header(char *argv,double **keyva,char *time_obs);
void read_files(char *fpath,char *fname,int *f_num,char *appendix,int charnum);
void  gaus_estimate(float **star0,int Ng1,int Ng2,double *mean,double *sigma);
void star_gaus(float **y,int nx,int ny,float **yt,float **residu,double *pa);
void Interp_bicubic(int nx,int ny,float **a0,int nbound,
		    int nxt,int nyt, float **at,double *xcen);
void BInterp_bicubic(int nx,int ny,float **a0,int nbound,
        int nxt,int nyt, float **at,double *xcen);
void iSPCA(double **images,double **weight,double **basisf,int Mp,int Nstar,int nbf,double **cent,
      double ***compbf,double **coeff,int Nb,int Nk,int num,double **compbfr,int Ng);
void estimate_para(double **star,double **center,int Nstar,int Ng, int nbound,
       double *para);
void creat_moffatelets(double rd,double beta,int l0,int msign,int *mx,double **basisf0,int Mp,int ng,int *nbf);
void creat_moffatelet(double rd,double beta,int l0,int msign,int *mx,double **basisf0,int Mp,int ng,int *nbf);
void EMPCA(double **star,int Nstar,int Mp,double **basisf,double **coeff,double **coefferr,
     double **weight,int nbf);
void SPCA(double **images,double **weight,double **basisf,int Mp,int Nstar,int nbf,
		  double **compbf,double **coeff,int Nb);
void PCA(double **star,int Nstar,int Mp,double **basef,double **coff);
void creat_gaussianlets(double sigma,int l0,int msign,int *mx,double **basisf0,int Mp,int ng,int *nbf);
void creat_gaussianlet(double sigma,int l0,int msign,int *mx,double **basisf0,int Mp,int ng,int *nbf);
float ran1(long *idum);
void size(float **I,int Ng,double *cen,double rg,double *se);
int get_header_key(char *argv,double *keyva,char *keyna);
void web_psf_fit(float *get_stars,float *get_spos,int npc,int Nstar,int Ng0,int Ng,
  float *get_PCs,float *get_coeff,float gain,int method,int *get_Nobj,int osam,
  float snrs,float *imodel,float *wgt,float *get_wcoff);
void web_method_branches(float ***stars,float **spos,int npc,int Nstar,int Ng0,int Ng,
  float ***PCs,float ***coeff,float gain,int method,int *Nobj,int osam,
  float snrs,float ***imdoel,float ***wgt,float **wcoff);
void webiSPCA_entr(float ***stars,float **spos,int npc,int Nstar0,int Ng0,int Ng,
  float ***PCs,float ***coeff,float gain,int *Nobj,float snrs,int osam,
  float ***models,float ***wgt, float **wcoff);
void webSPCA_entr(float ***stars,float **spos,int npc,int Nstar0,int Ng0,int Ng,
  float ***PCs,float ***coeff,float gain,int *Nobj,float snrs,float ***models,
  float ***wgt);
void webEMPCA_entr(float ***stars,float **spos,int npc,int Nstar0,int Ng0,int Ng,
  float ***PCs,float ***coeff,float gain,int *Nobj,float snrs,float ***models,
  float ***wgt);








