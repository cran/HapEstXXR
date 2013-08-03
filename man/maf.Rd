%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% library(HapEstXXR)
% Created: November 30, 2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{maf}
\alias{maf}
\title{
  Minor alle frequencies
}
\description{
  Calculation of minor allele frequencies (MAF),  call rate and asymptotic chisquare hardy-weinberg test
}
\usage{
maf ( geno , marker.label=NA )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geno}{
 (m,n)-genotype matrix

                m=number of individuals

                type R: n=number of snps
}
  \item{marker.label}{
Labels for the markers.
}
}
\details{
Call rate is defined by number of missing genotypes divided by sample size.

Testing deviation of the hardy-weinberg equilibirum is done by the usual
goodness-of-fit chisquare test: chi2 ~ sum( (observed-expected)^2/expected ).
}
\value{
This function returns an matrix with 8 coulmns.
}
\author{
Sven Knueppel
}

