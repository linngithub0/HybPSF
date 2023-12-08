#include <math.h>
#define NRANSI
#include "nrutil.h"
#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void fit(float x[], float y[], int ndata, float sig[], int mwt, float *a,
	float *b, float *siga, float *sigb, float *chi2, float *q)
{
	float gammq(float a, float x);
	int i;
	float wt,t,sxoss,sx=0.0,sy=0.0,st2=0.0,ss,sigdat;

	*b=0.0;
	if (mwt) {
		ss=0.0;
		for (i=0;i<ndata;i++) {
			wt=1.0/SQR(sig[i]);
			ss += wt;
			sx += x[i]*wt;
			sy += y[i]*wt;
		}
	} else {
		for (i=0;i<ndata;i++) {
			sx += x[i];
			sy += y[i];
		}
		ss=ndata;
	}
	sxoss=sx/ss;
	if (mwt) {
		for (i=0;i<ndata;i++) {
			t=(x[i]-sxoss)/sig[i];
			st2 += t*t;
			*b += t*y[i]/sig[i];
		}
	} else {
		for (i=0;i<ndata;i++) {
			t=x[i]-sxoss;
			st2 += t*t;
			*b += t*y[i];
		}
	}
	*b /= st2;
	*a=(sy-sx*(*b))/ss;
	*siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
	*sigb=sqrt(1.0/st2);
	*chi2=0.0;
	if (mwt == 0) {
		for (i=0;i<ndata;i++)
			*chi2 += SQR(y[i]-(*a)-(*b)*x[i]);
		*q=1.0;
		sigdat=sqrt((*chi2)/(ndata-2));
		*siga *= sigdat;
		*sigb *= sigdat;
	} else {
		for (i=0;i<ndata;i++)
			*chi2 += SQR((y[i]-(*a)-(*b)*x[i])/sig[i]);
		*q=gammq(0.5*(ndata-2),0.5*(*chi2));
	}
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software )1!. */

float gammq(float a, float x)
//Returns the incomplete gamma function Q(a, x) ≡ 1 − P (a, x). 
{
	void gcf(float *gammcf, float a, float x, float *gln);
	void gser(float *gamser, float a, float x, float *gln);
	void nrerror(char error_text[]);
	float gamser,gammcf,gln;
	if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine gammq");
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln); 
		return 1.0-gamser;
	} 
	else { 
		gcf(&gammcf,a,x,&gln); 
		return gammcf;
	} 
}

void gcf(float *gammcf, float a, float x, float *gln)
//Returns the incomplete gamma function Q(a, x) evaluated by its continued fraction represen- tation as gammcf. 
//Also returns lnΓ(a) as gln.
{
	float gammln(float xx);
	void nrerror(char error_text[]); 
	int i;
	float an,b,c,d,del,h;
	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d; 
	for(i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > ITMAX) nrerror("a too large, ITMAX too small in gcf");
	*gammcf=exp(-x+a*log(x)-(*gln))*h; //Put factors in front.
}

void gser(float *gamser, float a, float x, float *gln)
//Returns the incomplete gamma function P (a, x) evaluated by its series representation as gamser.
//Also returns lnΓ(a) as gln.
{
	float gammln(float xx);
	void nrerror(char error_text[]); 
	int n;
	float sum,del,ap;
	*gln=gammln(a); 
	if (x <= 0.0) {
		if (x < 0.0) nrerror("x less than 0 in routine gser"); 
		*gamser=0.0;
		return;
	} 
	else { 
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return; 
			}
		}
	nrerror("a too large, ITMAX too small in routine gser"); 
	return;
 	} 
}


float gammln(float xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5}; int j;
	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp); ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y; return -tmp+log(2.5066282746310005*ser/x);
}
