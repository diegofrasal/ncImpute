# Matrix completion for the MovieLens100k data set

We provide within our ```ncImpute``` package the famous ```MovieLens 100k``` and ```1m``` data sets.

```
library(ncImpute)

### MovieLens 100k ###
set.seed(0)
J=200
missfrac=0.2
eps=0.01
n.lambda=100

data("ML100k", package = "ncImpute")
ML100k.m = Incomplete(i=ML100k[,"user_id"],j=ML100k[,"item_id"],x=ML100k[,"rating"])
Omega = which(ML100k.m!=0)
length(Omega)

# Verifying every row and column has at least one element
which(rowSums(ML100k.m)==0)
which(colSums(ML100k.m)==0)

#Bi-scaling of ML100k
ML100krc = biScale(ML100k.m, row.center=TRUE, row.scale=FALSE, col.center=TRUE, col.scale=FALSE, trace=TRUE)

itest=sample(seq(length(Omega)),20000,replace=FALSE)
print(length(itest))

testData=ML100k[itest,]
print(dim(testData)[1])

trainData=ML100k[-itest,]
print(dim(trainData)[1])

# Row/Column Centering & Assignment of Attributes
x=trainData[,"rating"]
x=x-(ML100krc@`biScale:row`$center)[trainData[,"user_id"]]
x=x-(ML100krc@`biScale:column`$center)[trainData[,"item_id"]]

ML100knew=Incomplete(i=trainData[,"user_id"],j=trainData[,"item_id"],x=x)
attributes(ML100knew)=c(attributes(ML100knew), list("biScale:row"=ML100krc@`biScale:row`,"biScale:column"=ML100krc@`biScale:column`,critmat=ML100krc@critmat))
names(attributes(ML100knew))
norm(ML100knew)
1 - length(which(ML100knew==0))/(dim(ML100knew)[1]*dim(ML100knew)[2])

#Lambda grid
max.lam = 35
min.lam = 2
n.lambda = 100
lam.list = exp(seq(from=log(max.lam), to=log(min.lam), length.out=n.lambda))

#Gamma grid
gam.list=c(5000,1000,500,100,80,70,50,30,20,15,10,9,8,7,6,5.5,5,4.5,4,3.5,3,2.5,2,1.5,1.1)
Train=matrix(0, nrow=length(gam.list),ncol=length(lam.list))
Test=matrix(0, nrow=length(gam.list),ncol=length(lam.list))
ranks=matrix(0, nrow=length(gam.list),ncol=length(lam.list))

#Start Exploring Grid
fit.prev = vector("list",length(lam.list))
for(i in 1:length(gam.list)){
    rank.max=3
    for(j in 1:length(lam.list)){
        if(rank.max>185) rank.max=J
        if(i==1){
           if(j==1){ 
              fit.prev[[j]]=ncImpute(ML100knew,rank=rank.max,lambda=lam.list[j],g=gam.list[i],type="als",type.thresh="MC+",thresh=1e-05,maxit=70,trace.it=TRUE)
              ximp=Round.Ratings(impute(fit.prev[[j]]$fit, trainData[,"user_id"], trainData[,"item_id"]))  
              Train[i,j]=RMSE(ximp,trainData[,"rating"])              
              ximp=Round.Ratings(impute(fit.prev[[j]]$fit, testData[,"user_id"], testData[,"item_id"]))
              Test[i,j]=RMSE(ximp,testData[,"rating"])
              ranks[i,j]=fit.prev[[j]]$est.rank             
           }
           else{
              fit.prev[[j]]=ncImpute(ML100knew,rank=rank.max,lambda=lam.list[j],g=gam.list[i],type="als",type.thresh="MC+",thresh=1e-05,maxit=70,trace.it=TRUE,warm.start=fit.prev[[j-1]]$fit)
              ximp=Round.Ratings(impute(fit.prev[[j]]$fit, trainData[,"user_id"], trainData[,"item_id"]))  
              Train[i,j]=RMSE(ximp,trainData[,"rating"])              
              ximp=Round.Ratings(impute(fit.prev[[j]]$fit, testData[,"user_id"], testData[,"item_id"]))
              Test[i,j]=RMSE(ximp,testData[,"rating"])
              ranks[i,j]=fit.prev[[j]]$est.rank     
           }
        }    
        else{
           fit.prev[[j]]=ncImpute(ML100knew,rank=rank.max,lambda=lam.list[j],g=gam.list[i],type="als",type.thresh="MC+",thresh=1e-05,maxit=70,trace.it=TRUE,warm.start=fit.prev[[j]]$fit)
              ximp=Round.Ratings(impute(fit.prev[[j]]$fit, trainData[,"user_id"], trainData[,"item_id"]))  
              Train[i,j]=RMSE(ximp,trainData[,"rating"])              
              ximp=Round.Ratings(impute(fit.prev[[j]]$fit, testData[,"user_id"], testData[,"item_id"]))
              Test[i,j]=RMSE(ximp,testData[,"rating"])
              ranks[i,j]=fit.prev[[j]]$est.rank      
        }
        cat(i,gam.list[i],j,lam.list[j],"train",Train[i,j],"test",Test[i,j],"rank",ranks[i,j],"rank.max",rank.max,"\n")        
        rank.max=ranks[i,j] + 15
    }
}
Train
Test
ranks
```
