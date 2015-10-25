# Matrix completion for the MovieLens1m data set

We provide within our ```ncImpute``` package the famous ```MovieLens 100k``` and ```1m``` data sets.

This example carries out the corresponding matrix completion for the particular value of gamma ```=20``` in the MC+ penalty.

```
library(ncImpute)

### MovieLens 1m ###
set.seed(0)
J=250
missfrac=0.2
eps=0.18
n.lambda=100

data("ML1m", package = "ncImpute")
ML1m.m = Incomplete(i=ML1m[,"user_id"],j=ML1m[,"item_id"],x=ML1m[,"rating"])
Omega = which(ML1m.m!=0)
length(Omega)

#Verifying every row and column has at least one element
which(rowSums(ML1m.m)==0)
which(colSums(ML1m.m)==0)

#Bi-scaling of ML1m
ML1mrc = biScale(ML1m.m, row.center=TRUE, row.scale=FALSE, col.center=TRUE, col.scale=FALSE, trace=TRUE)

itest=sample(seq(length(Omega)),200209,replace=FALSE)
print(length(itest))

testData=ML1m[itest,]
print(dim(testData)[1])

trainData=ML1m[-itest,]
print(dim(trainData)[1])

#Row/Column Centering & Assignment of Attributes
x=trainData[,"rating"]
x=x-(ML1mrc@`biScale:row`$center)[trainData[,"user_id"]]
x=x-(ML1mrc@`biScale:column`$center)[trainData[,"item_id"]]

ML1mnew=Incomplete(i=trainData[,"user_id"],j=trainData[,"item_id"],x=x)
attributes(ML1mnew)=c(attributes(ML1mnew), list("biScale:row"=ML1mrc@`biScale:row`,"biScale:column"=ML1mrc@`biScale:column`,critmat=ML1mrc@critmat))
names(attributes(ML1mnew))
norm(ML1mnew)
1 - length(which(ML1mnew==0))/(dim(ML1mnew)[1]*dim(ML1mnew)[2])

#Lambda grid
max.lam = 75
min.lam = 15
n.lambda = 100
lam.list = exp(seq(from=log(max.lam), to=log(min.lam), length.out=n.lambda))

#Gamma grid
gamma.index=9
gam.list=c(5000,1000,500,100,80,70,50,30,20,15,10,9,8,7,6,5.5,5,4.5,4,3.5,3,2.5,2,1.5,1.1)
gam.list=gam.list[gamma.index]
Train=matrix(0, nrow=length(gam.list),ncol=length(lam.list))
Test=matrix(0, nrow=length(gam.list),ncol=length(lam.list))
ranks=matrix(0, nrow=length(gam.list),ncol=length(lam.list))
inspect(gam.list)

#Start Exploring Grid
fit.prev = vector("list",length(lam.list))
for(i in 1:length(gam.list)){
	rank.max=3
    for(j in 1:length(lam.list)){
    	if(rank.max>235) rank.max=J
        if(i==1){
           if(j==1){ 
              fit.prev[[j]]=ncImpute(ML1mnew,rank=rank.max,lambda=lam.list[j],g=gam.list[i],type="als",type.thresh="MC+",thresh=1e-05,maxit=70,trace.it=TRUE)
              ximp=Round.Ratings(impute(fit.prev[[j]]$fit, trainData[,"user_id"], trainData[,"item_id"]))  
              Train[i,j]=RMSE(ximp,trainData[,"rating"])              
              ximp=Round.Ratings(impute(fit.prev[[j]]$fit, testData[,"user_id"], testData[,"item_id"]))
              Test[i,j]=RMSE(ximp,testData[,"rating"])
              ranks[i,j]=fit.prev[[j]]$est.rank  
           }
           else{
              fit.prev[[j]]=ncImpute(ML1mnew,rank=rank.max,lambda=lam.list[j],g=gam.list[i],type="als",type.thresh="MC+",thresh=1e-05,maxit=70,trace.it=TRUE,warm.start=fit.prev[[j-1]]$fit)
              ximp=Round.Ratings(impute(fit.prev[[j]]$fit, trainData[,"user_id"], trainData[,"item_id"]))  
              Train[i,j]=RMSE(ximp,trainData[,"rating"])              
              ximp=Round.Ratings(impute(fit.prev[[j]]$fit, testData[,"user_id"], testData[,"item_id"]))
              Test[i,j]=RMSE(ximp,testData[,"rating"])
              ranks[i,j]=fit.prev[[j]]$est.rank           
           }
        }    
        else{
              fit.prev[[j]]=ncImpute(ML1mnew,rank=rank.max,lambda=lam.list[j],g=gam.list[i],type="als",type.thresh="MC+",thresh=1e-05,maxit=70,trace.it=TRUE,warm.start=fit.prev[[j]]$fit)
              ximp=Round.Ratings(impute(fit.prev[[j]]$fit, trainData[,"user_id"], trainData[,"item_id"]))  
              Train[i,j]=RMSE(ximp,trainData[,"rating"])              
              ximp=Round.Ratings(impute(fit.prev[[j]]$fit, testData[,"user_id"], testData[,"item_id"]))
              Test[i,j]=RMSE(ximp,testData[,"rating"])
              ranks[i,j]=fit.prev[[j]]$est.rank          
        }
        cat(gamma.index,gam.list[i],j,lam.list[j],"train",Train[i,j],"test",Test[i,j],"rank",ranks[i,j],"rank.max",rank.max,"\n")        
        rank.max=ranks[i,j] + 15        
    }
}
Train
Test
ranks
```
