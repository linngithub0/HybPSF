import numpy as np
from astropy.io import fits
import os
import matplotlib.pyplot as plt
def write_fits(fitsname,data):
    try:
        wdata=fits.open(fitsname)
        os.system('rm '+fitsname)
        wdata=data#-image_data
        hdu = fits.PrimaryHDU(wdata)
        hdul = fits.HDUList([hdu])
        hdul.writeto(fitsname)
    except IOError:
        wdata=data#-image_data
        hdu = fits.PrimaryHDU(wdata)
        hdul = fits.HDUList([hdu])
        hdul.writeto(fitsname)

def write_mult_fits(fitsname,data1,data2=None,
    data3=None,data4=None,data5=None,data6=None):
    if (data2.any()==None and data3.any()==None and data4.any()==None 
        and data5.any()==None and data6.any()==None):
        hdu = fits.PrimaryHDU(data1)
        hdul = fits.HDUList([hdu])
        hdul.writeto(fitsname,overwrite=True) 
    if (data2.any()!=None and data3.any()==None and data4.any()==None
        and data5.any()==None and data6.any()==None):
        hdu = fits.PrimaryHDU(data1)
        hdu1=fits.ImageHDU(data2)
        hdul = fits.HDUList([hdu,hdu1])
        hdul.writeto(fitsname,overwrite=True) 
    if (data2.any()!=None and data3.any()!=None and data4.any()==None 
        and data5.any()==None and data6.any()==None):
        hdu = fits.PrimaryHDU(data1)
        hdu1=fits.ImageHDU(data2)
        hdu2=fits.ImageHDU(data3)
        hdul = fits.HDUList([hdu,hdu1,hdu2])
        hdul.writeto(fitsname,overwrite=True) 
    if (data2.any()!=None and data3.any()!=None and data4.any()!=None 
        and data5.any()==None and data6.any()==None):
        hdu = fits.PrimaryHDU(data1)
        hdu1=fits.ImageHDU(data2)
        hdu2=fits.ImageHDU(data3)
        hdu3=fits.ImageHDU(data4)
        hdul = fits.HDUList([hdu,hdu1,hdu2,hdu3])
        hdul.writeto(fitsname,overwrite=True)
    if (data2.any()!=None and data3.any()!=None and data4.any()!=None 
        and data5.any()!=None and data6.any()==None):
        hdu = fits.PrimaryHDU(data1)
        hdu1=fits.ImageHDU(data2)
        hdu2=fits.ImageHDU(data3)
        hdu3=fits.ImageHDU(data4)
        hdu4=fits.ImageHDU(data5)
        hdul = fits.HDUList([hdu,hdu1,hdu2,hdu3,hdu4])
        hdul.writeto(fitsname,overwrite=True)
    if (data2.any()!=None and data3.any()!=None and data4.any()!=None 
        and data5.any()!=None and data6.any()!=None):
        hdu = fits.PrimaryHDU(data1)
        hdu1=fits.ImageHDU(data2)
        hdu2=fits.ImageHDU(data3)
        hdu3=fits.ImageHDU(data4)
        hdu4=fits.ImageHDU(data5)
        hdu5=fits.ImageHDU(data6)
        hdul = fits.HDUList([hdu,hdu1,hdu2,hdu3,hdu4,hdu5])
        hdul.writeto(fitsname,overwrite=True)

hudl="spec.fits"
data=fits.open(hudl)
data=data[0].data
Mp=data.shape[0]
Nobj=data.shape[1];print("Nobj=",Nobj)
Nobj=1000
for i in range(Nobj):
    data[:,i]=data[:,i]/(np.sum(data[:,i]))
#data=data.T
tdata=data[:,0:Nobj];print(tdata.shape)
#Nobj=tdata.shape[1]
npc=20
#print(tdata.shape)
PC,sigma,v=np.linalg.svd(tdata)
#print(PC.shape,sigma.shape,v.shape)
coeff=np.zeros((npc,Nobj))
for i in range(npc):
    coeff[i,:]=sigma[i]*v.T[:,i]

print(PC.shape)
recon=np.dot(PC[:,0:npc],coeff)
error=tdata-recon
#print(error.shape)

x=range(tdata.shape[0])
plt.figure(figsize=(10,30))
for i in range(50):
    plt.subplot(50,1,i+1)
    plt.plot(x,PC[:,i],'-')
plt.savefig("PC1.pdf")
rnum=100
indx=np.random.randint(0,Nobj,size=rnum,dtype='int')
plt.figure(figsize=(10,80))
for i in range(rnum):
    plt.subplot(rnum,1,i+1)
    k=indx[i]
    plt.plot(x,tdata[:,k],'-',c='black',lw=5)
    plt.plot(x,recon[:,k],'-',c='blue',lw=1)
    plt.plot(x,error[:,k],'-',c='red')

plt.savefig("compare1.eps")

plt.figure(figsize=(10,2))
residual=np.zeros(Mp)
onesigma=np.zeros(Mp)
for i in range(Mp):
    residual[i]=np.mean(error[i,:])
    onesigma[i]=np.std(error[i,:])
plt.errorbar(x,residual,yerr=onesigma,marker='.')
plt.savefig("error1.pdf")








