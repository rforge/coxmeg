context("test gds function")

.testSparseMatrix <- function(n, n_blocks) {
    n_f <- ceiling(n/n_blocks)
    mat_list <- list()
    size <- rep(n_blocks,n_f)
    offd <- 0.5
    for(i in 1:n_f){
        mat_list[[i]] <- matrix(offd,size[i],size[i])
        diag(mat_list[[i]]) <- 1
    }
    sigma <- as.matrix(bdiag(mat_list))
    sigma <- sigma[1:n, 1:n]
    sigma <- as(sigma,'dgCMatrix')
    return(sigma)
}

test_that("coxmeg_gds matches coxmeg_plink", {
    bed = system.file("extdata", "example_null.bed", package = "coxmeg")
    bed = substr(bed,1,nchar(bed)-4)
    pheno.file = system.file("extdata", "ex_pheno.txt", package = "coxmeg")
    cov.file = system.file("extdata", "ex_cov.txt", package = "coxmeg")
    
    ## building a relatedness matrix
    sigma <- .testSparseMatrix(n=3000, n_blocks=5)
    
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


test_that("SNPRelate and SeqArray methods match", {
    snpfile <- SNPRelate::snpgdsExampleFileName()
    seqfile <- tempfile()
    SeqArray::seqSNP2GDS(snpfile, seqfile, verbose=FALSE)
    
    snp <- SNPRelate::snpgdsOpen(snpfile)
    seq <- SeqArray::seqOpen(seqfile)
    
    snpsel <- .gdsSelectSNP(snp, maf=0.01, missing.rate=0.01, verbose=FALSE)
    seqsel <- .gdsSelectSNP(seq, maf=0.01, missing.rate=0.01, verbose=FALSE)
    expect_equal(snpsel, seqsel)
    
    seqResetFilter(seq, verbose=FALSE)
    sample.id <- seqGetData(seq, "sample.id")[1:50]
    snp.id <- seqGetData(seq, "variant.id")[1:100]
    snpref <- substr(gdsfmt::read.gdsn(gdsfmt::index.gdsn(snp, "snp.allele")), 1, 1)
    seqref <- seqGetData(seq, "$ref")
    allele.swap <- (snpref != seqref)
    swap.sel <- allele.swap[1:100]
    snpgeno <- .gdsGetGeno(snp, sample.id=sample.id, snp.id=snp.id, verbose=FALSE)
    seqgeno <- .gdsGetGeno(seq, sample.id=sample.id, snp.id=snp.id, verbose=FALSE)
    expect_equivalent(snpgeno[,!swap.sel], seqgeno[,!swap.sel])
    expect_equivalent(snpgeno[,swap.sel], 2-seqgeno[,swap.sel])
    
    snplist <- .gdsSNPList(snp)
    seqlist <- .gdsSNPList(seq)
    snplist$chromosome <- as.character(snplist$chromosome)
    expect_equivalent(snplist[!allele.swap,], seqlist[!allele.swap,])
    expect_equivalent(snplist[allele.swap,1:3], seqlist[allele.swap,1:3])
    expect_equal(snplist$afreq[allele.swap], 1-seqlist$afreq[allele.swap])
    
    SNPRelate::snpgdsClose(snp)
    SeqArray::seqClose(seq)
    unlink(seqfile)
})


test_that("SNPRelate and SeqArray coxmeg_gds match", {
    snpfile <- SNPRelate::snpgdsExampleFileName()
    seqfile <- tempfile()
    SeqArray::seqSNP2GDS(snpfile, seqfile, verbose=FALSE)
    
    snp <- SNPRelate::snpgdsOpen(snpfile)
    seq <- SeqArray::seqOpen(seqfile)
    
    sample.id <- seqGetData(seq, "sample.id")
    n <- length(sample.id)
    pheno <- data.frame(family=sample.id,
                        sample.id=sample.id,
                        time=rnorm(n, mean=100, sd=10),
                        status=rbinom(n, 1, 0.4),
                        stringsAsFactors=FALSE)
    
    # covariance matrix
    #covmat <- SNPRelate::snpgdsGRM(snp, verbose=FALSE)
    #sigma <- covmat$grm
    sigma <- .testSparseMatrix(n=n, n_blocks=5)
    
    # set high MAF so test runs faster
    re.snp <- coxmeg_gds(snp, pheno, sigma, type='bd', maf=0.47, verbose=FALSE)
    re.seq <- coxmeg_gds(seq, pheno, sigma, type='bd', maf=0.47, verbose=FALSE)
    allele.swap <- re.snp$summary$allele != re.seq$summary$allele
    expect_equal(re.snp$summary$beta[!allele.swap], re.seq$summary$beta[!allele.swap], tolerance=1e-4)
    expect_equal(re.snp$summary$beta[allele.swap], -re.seq$summary$beta[allele.swap], tolerance=1e-4)
    expect_equal(re.snp$summary$HR[!allele.swap], re.seq$summary$HR[!allele.swap], tolerance=1e-4)
    expect_equal(re.snp$summary$HR[allele.swap], 1/re.seq$summary$HR[allele.swap], tolerance=1e-4)
    expect_equal(re.snp$summary$sd_beta, re.seq$summary$sd_beta, tolerance=1e-4)
    expect_equal(re.snp$summary$p, re.seq$summary$p, tolerance=1e-4)
    
    SNPRelate::snpgdsClose(snp)
    SeqArray::seqClose(seq)
    unlink(seqfile)
})
