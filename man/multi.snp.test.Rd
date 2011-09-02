%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% library(HapEstXXR)
% Created: November 30, 2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{multi.snp.test}
\alias{multi.snp.test}
\title{
Internal function used for multi-locus associations tests.
}
\description{
This function is used for internal computations. You should not used use it, but you can.
}
\usage{
multi.snp.test ( y , x , x.adj=NULL ,
        type = c("gaussian", "binomial") )
}
\arguments{
  \item{y}{
response
}
  \item{x}{
Matrix including SNPs or haplotypes
}
  \item{x.adj}{
Matrix of covariates
}
  \item{type}{
type of response
}
}
\author{
Sven Knueppel
}
