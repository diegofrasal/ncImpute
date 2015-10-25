# Noisy matrix completion under the Coherent Model

```
library(ncImpute)

#### Coherent Model ####
set.seed(0)
n=800
p=400
J=10
SNR=10
missfrac=0.95
np=n*p
eps=0.001
n.lambda=100

#### Block Parameters for U ####
br=160
bc=2
#Forming orthogonal left-singular vectors
U=matrix(0, nrow=n, ncol=J)
U[1:br,1:bc]=matrix(rnorm(br*bc),br,bc)
U[(br+1):(2*br),(bc+1):(2*bc)]=matrix(rnorm(br*bc),br,bc)
U[(2*br+1):(3*br),(2*bc+1):(3*bc)]=matrix(rnorm(br*bc),br,bc)
U[(3*br+1):(4*br),(3*bc+1):(4*bc)]=matrix(rnorm(br*bc),br,bc)
U[(4*br+1):(5*br),(4*bc+1):(5*bc)]=matrix(rnorm(br*bc),br,bc)
decomp=qr(U)
U=qr.Q(decomp)

#### Block Parameters for V ####
br=80
bc=2
#Forming orthogonal right-singular vectors
V=matrix(0, nrow=p, ncol=J)
V[1:br,1:bc]=matrix(rnorm(br*bc),br,bc)
V[(br+1):(2*br),(bc+1):(2*bc)]=matrix(rnorm(br*bc),br,bc)
V[(2*br+1):(3*br),(2*bc+1):(3*bc)]=matrix(rnorm(br*bc),br,bc)
V[(3*br+1):(4*br),(3*bc+1):(4*bc)]=matrix(rnorm(br*bc),br,bc)
V[(4*br+1):(5*br),(4*bc+1):(5*bc)]=matrix(rnorm(br*bc),br,bc)
decomp=qr(V)
V=qr.Q(decomp)

set.seed(101)
#Forming singular values and noise with SNR
d=sort(runif(J, 0, 100), decreasing=TRUE)
d

#Forming signal matrix
signal=U%*%(d*t(V))
num=var(as.vector(signal))   
#Old noise matrix
E.old=matrix(rnorm(np),n,p)
den=SNR*var(as.vector(E.old))
tau=num/den
#New noise matrix
E.real=sqrt(tau)*E.old

#Forming matrix x
x=U%*%(d*t(V)) + E.real
#Creating Missing Values in a Full Matrix
ix=seq(np)
imiss=sample(ix,np*missfrac,replace=FALSE)
xna=x
xna[imiss]=NA
U
V

#Smallest value of lambda such that ncImpute returns the ZERO solution and lambda grid
max.lam=lambda0(xna)
max.lam=max.lam-0.1
min.lam=max.lam*eps
lam.list=seq(from=max.lam, to=min.lam, length.out=n.lambda)

#Gamma grid
gam.list=c(5000,1000,500,100,80,70,50,30,20,15,10,9,8,7,6,5.5,5,4.5,4,3.5,3,2.5,2,1.5,1.1)
Train=matrix(0, nrow=length(gam.list),ncol=length(lam.list))
Test=matrix(0, nrow=length(gam.list),ncol=length(lam.list))
ranks=matrix(0, nrow=length(gam.list),ncol=length(lam.list))

#Start Exploring Grid
fit.prev = vector("list",length(lam.list))
for(i in 1:length(gam.list)){
    for(j in 1:length(lam.list)){
        if(i==1){
           if(j==1){ 
              fit.prev[[j]]=ncImpute(xna,rank=50,lambda=lam.list[j],g=gam.list[i],type="svd",type.thresh="MC+",thresh=1e-05,maxit=70,trace.it=FALSE)
              ximp=complete(xna,fit.prev[[j]]$fit)
              Train[i,j]=Training.Error(U,d,V,E.real,fit.prev[[j]]$fit,imiss)
              Test[i,j]=Test.Error(U,d,V,ximp,imiss)
              ranks[i,j]=fit.prev[[j]]$est.rank
           }
           else{
              fit.prev[[j]]=ncImpute(xna,rank=50,lambda=lam.list[j],g=gam.list[i],type="svd",type.thresh="MC+",thresh=1e-05,maxit=70,trace.it=FALSE,warm.start=fit.prev[[j-1]]$fit)
              ximp=complete(xna,fit.prev[[j]]$fit)
              Train[i,j]=Training.Error(U,d,V,E.real,fit.prev[[j]]$fit,imiss)         
              Test[i,j]=Test.Error(U,d,V,ximp,imiss)
              ranks[i,j]=fit.prev[[j]]$est.rank            
           }
        }    
        else{
           if(j==1){ 
              fit.prev[[j]]=ncImpute(xna,rank=50,lambda=lam.list[j],g=gam.list[i],type="svd",type.thresh="MC+",thresh=1e-05,maxit=70,trace.it=FALSE,warm.start=fit.prev[[1]]$fit)
              ximp=complete(xna,fit.prev[[j]]$fit)
              Train[i,j]=Training.Error(U,d,V,E.real,fit.prev[[j]]$fit,imiss)
              Test[i,j]=Test.Error(U,d,V,ximp,imiss)     
              ranks[i,j]=fit.prev[[j]]$est.rank      
           }
           else{
              aux1=ncImpute(xna,rank=50,lambda=lam.list[j],g=gam.list[i],type="svd",type.thresh="MC+",thresh=1e-05,maxit=70,trace.it=FALSE,warm.start=fit.prev[[j]]$fit)
              aux2=ncImpute(xna,rank=50,lambda=lam.list[j],g=gam.list[i],type="svd",type.thresh="MC+",thresh=1e-05,maxit=70,trace.it=FALSE,warm.start=fit.prev[[j-1]]$fit)
              if(aux1$obj<aux2$obj) fit.prev[[j]]=aux1
              if(aux2$obj<aux1$obj) fit.prev[[j]]=aux2
              ximp=complete(xna,fit.prev[[j]]$fit)
              Train[i,j]=Training.Error(U,d,V,E.real,fit.prev[[j]]$fit,imiss)        
              Test[i,j]=Test.Error(U,d,V,ximp,imiss)       
              ranks[i,j]=fit.prev[[j]]$est.rank       
           }        
        }
    }
}
Train
Test
ranks
```
