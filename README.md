# ncImpute
Matrix Completion via Non-Convex Regularization:

This ```R``` package provides iterative algorithms for matrix completion based on non-convex regularization of the singular values. The main approach uses iterative MC+ thresholded singular value decompositions to impute the missing values, and has an "EM" flavor, in that at each iteration the matrix is completed with the current estimate. For large matrices there is a special sparse-matrix class named "Incomplete" that efficiently handles all computations. The package includes procedures for centering and scaling rows, columns or both.

To install the current source version of our package in a ```Windows 10``` operating system, simply type the follwing instructions within your ```R``` session:

```
install.packages("/your_file_path/ncImpute_1.0.tar.gz", lib="/your_R_packages_library", repos=NULL, type="source")
```

For users with ```Mac OS X 10.4``` and higher, you must first make sure you have an updated version of [GNU Fortran](https://cran.r-project.org/bin/macosx/tools/){:target="_blank"} properly installed.

A series of worked examples follow below. More documentation to come. This is [a link](http://example.com){:target="_blank"} that opens in a new window. <a href="http://example.com" target="_blank"></a>

```
library(ncImpute)
set.seed(101)
n=200
p=100
J=50
np=n*p
missfrac=0.3
x=matrix(rnorm(n*J),n,J) %*% matrix(rnorm(J*p),J,p) + matrix(rnorm(np),n,p)/5
ix=seq(np)
imiss=sample(ix,np*missfrac,replace=FALSE)
xna=x
xna[imiss]=NA
###uses regular matrix method for matrices with NAs
fit1=ncImpute(xna,rank=50,lambda=30,g=20)
###uses sparse matrix method for matrices of class "Incomplete"
xnaC=as(xna,"Incomplete")
fit2=ncImpute(xnaC,rank=50,lambda=30,g=20)
###uses "als" algorithm
fit3=ncImpute(xnaC,rank=50,lambda=30,g=20,type="als")
fit4=ncImpute(xnaC,rank=50,lambda=30,type="als",type.thresh="SOFT")
ximp=complete(xna,fit1$fit)
### first scale xna
xnas=biScale(xna)
fit5=ncImpute(xnas,rank=50,lambda=10)
ximp=complete(xna,fit4$fit)
impute(fit5$fit,i=c(1,3,7),j=c(2,5,10))
impute(fit5$fit,i=c(1,3,7),j=c(2,5,10),unscale=FALSE)#ignore scaling and centering
```
