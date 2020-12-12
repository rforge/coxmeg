setOldClass("SNPGDSFileClass")
setGeneric(".gdsSelectSNP", function(gdsobj, ...) standardGeneric(".gdsSelectSNP"))
setGeneric(".gdsGetGeno", function(gdsobj, ...) standardGeneric(".gdsGetGeno"))
setGeneric(".gdsSNPList", function(gdsobj, ...) standardGeneric(".gdsSNPList"))

setMethod(".gdsSelectSNP",
          "SNPGDSFileClass",
          function(gdsobj, sample.id=NULL, maf=NaN, missing.rate=NaN,
                   remove.monosnp=TRUE, verbose=TRUE){
              snpgdsSelectSNP(gdsobj, sample.id=sample.id,
                              maf=maf, missing.rate=missing.rate,
                              remove.monosnp=remove.monosnp,
                              verbose=verbose,
                              autosome.only=FALSE)
          })

setMethod(".gdsGetGeno",
          "SNPGDSFileClass",
          function(gdsobj, sample.id=NULL, snp.id=NULL, 
                   with.id=TRUE, verbose=TRUE){
              snpgdsGetGeno(gdsobj, sample.id=sample.id, snp.id=snp.id,
                            with.id=with.id, verbose=verbose)
          })


setMethod(".gdsSNPList",
          "SNPGDSFileClass",
          function(gdsobj, sample.id=NULL){
              snpgdsSNPList(gdsobj, sample.id=sample.id)
          })
