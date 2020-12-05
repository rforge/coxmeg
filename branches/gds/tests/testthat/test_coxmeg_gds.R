context("test gds function")

test_that("coxmeg_gds matches coxmeg_plink", {
    bed = system.file("extdata", "example_null.bed", package = "coxmeg")
    bed = substr(bed,1,nchar(bed)-4)
    pheno.file = system.file("extdata", "ex_pheno.txt", package = "coxmeg")
    cov.file = system.file("extdata", "ex_cov.txt", package = "coxmeg")
    
    ## building a relatedness matrix
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
    
    re.plink = coxmeg_plink(pheno.file,sigma,type='bd',bed=bed,tmp_dir=tempdir(),cov_file=cov.file)
    
    gdsfile <- tempfile()
    SNPRelate::snpgdsBED2GDS(bed.fn=paste0(bed,".bed"), fam.fn=paste0(bed,".fam"), bim.fn=paste0(bed,".bim"),
                             out.gdsfn=gdsfile, verbose=FALSE)
    gds <- SNPRelate::snpgdsOpen(gdsfile)
    pheno <- read.table(pheno.file, header=FALSE, as.is=TRUE, na.strings="-9")
    cov <- read.table(cov.file, header=FALSE, as.is=TRUE)
    
    re.gds <- coxmeg_gds(gds,pheno,sigma,type='bd',cov=cov)
    expect_equal(re.plink, re.gds, tolerance=1e-5)
    
    SNPRelate::snpgdsClose(gds)
    unlink(gdsfile)
})
