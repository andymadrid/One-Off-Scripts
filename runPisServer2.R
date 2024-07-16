runPis <- function(dss.formula, dmlFit, df, dmlTest, meth.mat) {
    X <- model.matrix(dss.formula, df)
    beta.mat <- data.frame(dmlFit$fit$beta)
    colnames(beta.mat) <- colnames(dmlFit$X)
    beta.mat <- as(beta.mat, "matrix")
    compute_pis <- function(X, beta.mat, test.result){ 
    N <- nrow(beta.mat)
 
    out <- mclapply(FUN = function(i) {
        compute_pi_for_cpg_i(i, X, beta.mat)
        }, 
        1:N,mc.cores=12
     )
 
    pi.mat <- do.call(rbind, out)
 
    # Assign column names based on sample id
    colnames(pi.mat) <- rownames(X)
 
    # Add chr, start, end as leading columns
    class(test.result) <- "data.frame" # remove DSS class
    bed.cols <- dplyr::transmute(test.result, chr, start = pos, end = pos + 1)
 
     pi.df <- data.table(cbind(bed.cols, pi.mat), rownames = NULL)
     pi.df
     }
    compute_pi_for_cpg_i <- function(i, X, beta.mat){
        # helper function that calls KMP's matrix multiply script
        #i is CpG index
        b <- beta.mat[i, ]

        # Inverse link function
        # y = arcsin(2*X\beta - 1)
        tmp <- (sin(multiply(as(X,'sparseMatrix'), as(b, 'sparseMatrix'))) + 1) / 2
        as.numeric(tmp)
    }

    Rcpp::sourceCpp("/media/data/WGBS/Datasets/LOAD_MCI_CONTROL/matrix_multiply.cpp")
    suppressPackageStartupMessages({
         library(data.table)
         library(magrittr)
         library(argparse)
         library(DSS)
         library(fdrtool)
         library(GenomicRanges)
         library(Rcpp)
     })
    pi.df <- compute_pis(X, beta.mat, dmlTest)
    pi.df <- as.data.frame(pi.df)
    rownames(pi.df) <- rownames(meth.mat)
    pi.df <- pi.df[,-c(1:3)]
    return(pi.df)
}
