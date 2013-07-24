%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% library(HapEstXXR)
% Created: November 30, 2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{single.haplotype.test}
\alias{single.haplotype.test}
\alias{single.haplotype.test.binomial}
\alias{single.haplotype.test.gaussian}
\alias{single.haplotype.test.families}
\alias{single.haplotype.test.survival}
\alias{hapest.caco}
\title{
single haplotype test
}
\description{
  Association test based on haplotypes. Haplotypes are estimated by EM algorithm.
}
\usage{
 single.haplotype.test ( snps, trait, famid , patid , fid , mid ,
    adj.var=NULL , type = c("gaussian", "binomial", "families") ,
    prt=T , lim=0.05 , min.count=10 , alpha=0.05 ,
    sort=F )

}
\arguments{
  \item{snps}{(n.m)-Matrix; n=No. of individuals; m=no. of SNPs; Rohde-Code}
  \item{trait}{numeric; Outcome, phenotype }
  \item{famid}{vector; Identifier for every family; only need in case of type=families}
  \item{patid}{vector; Identifier for every individuals; only need in case of type=families}
  \item{fid}{vector; Identifier for father (0=unkown); only need in case of type=families}
  \item{mid}{vector; Identifier for mother (0=unkown); only need in case of type=families}
  \item{adj.var}{(n,m)-Matrix;   n=No. of individuals; m=no. of covariates; variables for adjustment;
     in case of type=families not available.}
  \item{type}{type of depending variable}
  \item{lim}{numeric; threshold for pooling of haplotypes and declare as rare.}
  \item{min.count}{Minimal count for using pooled rare haplotypes in the analysis.}
  \item{prt}{A logical value (TRUE or FALSE). If TRUE, an overview is printed.}
  \item{alpha}{In case of type=binomial the (1-alpha/2)-confidence intervals are computed. }
  \item{sort}{A logical value (TRUE or FALSE). Only usable with family data.
       If TRUE, families are sorted by famid and generation which is a condition of wTDT.}
}
\details{
  Haplotypes are infered by EM algorithm (Excoffier and Slatkin 1995).
  
  For normal distributed phenotypes from independent individuals we prefer an F
  test and for case control data we prefer the likelihood ratio test
  (logistic regression) in comparison of full model with genetic and non-genetic
  factors to a reduced model, which includes only non-genetic variables. In the
  case of no specified non-genetic variable only the intercept is used. If one of
  these tests are significance we assume a genetic effect.
  In case of family data the weigthed TDT statistic is used.

}
\value{
   \item{hap}{Haplotypes}
   \item{freq }{Estimated haplotype frequencies}
   \item{global.test}{Result of global test statistic.}
   \item{haplotype.i}{Result of haplotype specific tests}
}
\references{
   Excoffier L, Slatkin M. Mol Biol Evol. 1995 Sep;12(5):921-7.
   
   Rohde K, Fuerst R. Hum Hered. 2003;56(1-3):41-7.
   
   Rohde K, Fuerst R. Hum Mutat. 2001 Apr;17(4):289-95.

Knueppel S, Esparza-Gordillo J, Marenholz I, Holzhuetter HG,

Bauerfeind A, Ruether A, Weidinger S, Lee Y-A, Rohde K.

Multi-locus stepwise regression: a haplotype-based algorithm

for finding genetic associations applied to atopic dermatitis.

BMC Med Genet 2012;13(1):8.
}
\author{
  Sven Knueppel and Klaus Rohde
}
\seealso{
  \code{\link{single.snp.test}}
}
\keyword{ snp }
\keyword{ haplotypes }
\keyword{ logistic regression }
\keyword{ weigthed TDT }
\keyword{ survival analysis }

