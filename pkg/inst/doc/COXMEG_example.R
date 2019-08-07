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
n_f <- 600
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

## ----echo=TRUE-----------------------------------------------------------
re = coxmeg(outcome,sigma,pred,order=1,eigen=TRUE,dense=FALSE)
re

