from scipy import linalg
import numpy as np;
import h5py;
import matplotlib.pyplot as plt
fname = "z0.00-0-0-mags-sd-interp.h5";
npc=10
with h5py.File(fname, "r") as fp:
    spec=np.asarray(fp['spectra'])
spec=spec.T;print(spec.shape)
Mp=spec.shape[0];Ntotal=spec.shape[1]
Nobj=1000;
indx=np.random.randint(0,Ntotal,size=Nobj,dtype='int')
sample=np.zeros((Mp,Nobj),dtype='float')
for i in range(Nobj):
    sample[:,i]=spec[:,indx[i]] #normalize
    sample[:,i]/=np.sum(sample[:,i])

PC,sigma,v=np.linalg.svd(sample)
#reconstruction check
coeff=np.zeros((npc,Nobj))
for i in range(npc):
    coeff[i,:]=sigma[i]*v.T[:,i]
recon=np.dot(PC[:,0:npc],coeff)
error=sample-recon
x=range(sample.shape[0])
plt.figure(figsize=(10,30))
for i in range(npc):
    plt.subplot(npc,1,i+1)
    plt.plot(x,PC[:,i],'-')
plt.savefig("PC1.pdf")
rnum=100
indx=np.random.randint(0,Nobj,size=rnum,dtype='int')
plt.figure(figsize=(10,80))
for i in range(rnum):
    plt.subplot(rnum,1,i+1)
    k=indx[i]
    plt.plot(x,sample[:,k],'-',c='black',lw=5)
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


#sampling check
Nsample=10000;print(spec.shape)
indx=np.random.randint(0,Ntotal,size=Nsample,dtype='int')
samples=np.zeros((Mp,Nsample),dtype='float')
for i in range(Nsample):
    samples[:,i]=spec[:,indx[i]] #normalize
    samples[:,i]/=np.sum(samples[:,i])

coeffs=np.zeros((npc,Nsample))
for i in range(npc):
    for ic in range(Nsample):
        coeffs[i][ic]=np.dot(samples[:,ic],PC[:,i])

recons=np.dot(PC[:,0:npc],coeffs)
error=samples-recons
x=range(samples.shape[0])

rnum=100
indx=np.random.randint(0,Nsample,size=rnum,dtype='int')
plt.figure(figsize=(10,80))
for i in range(rnum):
    plt.subplot(rnum,1,i+1)
    k=indx[i]
    plt.plot(x,samples[:,k],'-',c='black',lw=5)
    plt.plot(x,recons[:,k],'-',c='blue',lw=1)
    plt.plot(x,error[:,k],'-',c='red')
plt.savefig("compare2.eps")
plt.figure(figsize=(10,2))
residual=np.zeros(Mp)
onesigma=np.zeros(Mp)
for i in range(Mp):
    residual[i]=np.mean(error[i,:])
    onesigma[i]=np.std(error[i,:])
plt.errorbar(x,residual,yerr=onesigma,marker='.')
plt.savefig("error2.pdf")


