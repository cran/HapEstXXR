%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% library(HapEstXXR)
% Created: November 30, 2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{itegeppXXR}
\alias{itegeppXXR}
\title{ Haplotype estimation routine for single individual data }
\description{
  itegeppXXR is haplotype estimation routine for samples of independent individual genotypes (EM-algorithm).
}
\usage{
itegeppXXR(geno, des = 0, lim = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geno}{ (n,m)-Matrix; n=No of Individuals, m=No of SNPs; R-Code: 1-column genotype matrix - code  1 = 1/1, 3 = 1/2, 2 = 2/2 }
  \item{des}{ des=1 haplotype pairs,   des=0 single haplotypes }
  \item{lim}{ Threshold for combining rare haplotypes }
}
\details{
Inferring haplotypes by EM-Algorithm
}
\value{

      \item{\code{hap.id} }{Id. of haplotypes}
      \item{\code{hap}    }{estimated haplotypes}
      \item{\code{freq}   }{haplotype frequencies}
      \item{\code{hapres} }{individual haplotypes}
      \item{\code{likres} }{Likelihood value}
      \item{\code{desres}}{Design matrix for the model (des=1 => Haplotype pairs, des=0 => single haplotypes)}

}
\references{ Excoffier L, Slatkin M (1995) Maximum-likelihood estimation of molecular haplotype
frequencies in a diploid population. Mol Biol Evol 12:921-927 }
\author{ Sven Knueppel and Klaus Rohde}
\note{
This function works only up to 15 SNP haplotypes
}
\examples{
  set.seed(123456)
  ns <- 4  # Number of SNPs
  N <- 2000   # Number of individuals
  patid <- N:1
  geno <- matrix(sample(c(1, 2, 3), ns * N, replace = TRUE), ncol = ns)
  iteHAP <- itegeppXXR(geno, des = 1, lim = 0.01)
}

