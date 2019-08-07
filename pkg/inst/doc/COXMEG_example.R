## ---- message=FALSE, warning=FALSE,echo=FALSE----------------------------
library(knitcitations)
cleanbib()
options("citation_format" = "pandoc")
r<-citep("10.3389/fpubh.2016.00003")
write.bibtex(file="references.bib")

## ----eval=FALSE----------------------------------------------------------
#  install.packages("coxmeg", repos="http://R-Forge.R-project.org")

## ----echo=TRUE-----------------------------------------------------------
library(coxmeg)
library(MASS)
library(Matrix)
n_f <- 200
mat_list <- list()
size <- rep(5,n_f)
offd <- 0.5
for(i in 1:n_f)
{
  mat_list[[i]] <- matrix(offd,size[i],size[i])
  diag(mat_list[[i]]) <- 1
}
sigma <- as.matrix(bdiag(mat_list))
sigma = as(sigma,'dgCMatrix')


## ----echo=TRUE-----------------------------------------------------------
n = nrow(sigma)
tau_var <- 0.2
x <- mvrnorm(1, rep(0,n), tau_var*sigma)
pred = rnorm(n,0,1)
myrates <- exp(x+0.1*pred-1)
y <- rexp(n, rate = myrates)
cen <- rexp(n, rate = 0.02 )
ycen <- pmin(y, cen)
outcome <- cbind(ycen,as.numeric(y <= cen))
head(outcome)
sigma[1:5,1:5]


## ----echo=TRUE-----------------------------------------------------------
re = coxmeg(outcome,sigma,pred,order=1,dense=FALSE)
re

## ----echo=TRUE-----------------------------------------------------------
library(coxme)
bls <- c(1)
for(i in (size[1]-1):1)
{bls <- c(bls, c(rep(offd,i),1))}
tmat <- bdsmatrix(blocksize=size, blocks=rep(bls,n_f),dimnames=list(as.character(1:n),as.character(1:n)))
re_coxme = coxme(Surv(outcome[,1],outcome[,2])~as.matrix(pred)+(1|as.character(1:n)), varlist=list(tmat),ties='breslow')
re_coxme

