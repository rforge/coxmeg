## ---- message=FALSE, echo=FALSE------------------------------------------
library(knitcitations)
cleanbib()
options("citation_format" = "pandoc")
r<-citep("10.1016/0040-5809(77)90005-3") 
r<-citep("10.1016/j.mbs.2006.11.006")
r<-citep("10.1080/08898480590932296")
r<-citep("10.1007/s10522-006-9073-3")
r<-citep("10.1016/j.jtbi.2009.01.023")
r<-citep("10.3389/fpubh.2014.00228")
r<-citep("10.1002/gepi.22058")
r<-citep("10.3389/fpubh.2016.00003")
write.bibtex(file="references.bib")

## ----eval=FALSE----------------------------------------------------------
#  install.packages("coxmeg", repos="http://R-Forge.R-project.org")

