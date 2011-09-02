%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% library(HapEstXXR)
% Created: November 30, 2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{msr}
\alias{msr}
\alias{msr.binomial.forward.adjusted}
\alias{msr.binomial.forward.unadjusted}
\alias{msr.binomial.haplotype.test.adjusted}
\alias{msr.binomial.haplotype.test.unadjusted}
\alias{msr.gaussian.haplotype.test.adjusted}
\alias{msr.gaussian.haplotype.test.unadjusted}
\alias{msr.families.unadjusted}
\alias{msr.gaussian.forward.adjusted}
\alias{msr.gaussian.forward.unadjusted} 
\title{ Multi-locus stepwise regression (MSR)}
\description{
  Stepwise regression for snp selection and haplotype testing
}
\usage{
msr ( famid, patid, fid, mid, snps, trait, adj.var = NA,
      lim = 0.05, maxSNP = 3, nt = 10, pair.begin = FALSE,
      pattern.begin.mat = NA, 
      type = "gaussian", baseline.hap = "max",
      min.count = 10, sort=F )
}
\arguments{
  \item{famid}{vector; Identifier for every family; only need in case of type=families}
  \item{patid}{vector; Identifier for every individuals; only need in case of type=families}
  \item{fid}{vector; Identifier for father (0=unkown); only need in case of type=families}
  \item{mid}{vector; Identifier for mother (0=unkown); only need in case of type=families}
  \item{snps}{(n.m)-Matrix; n=No. of individuals; m=no. of SNPs; Rohde-Code}
  \item{trait}{numeric; Outcome, phenotype }
  \item{adj.var}{(n,m)-Matrix;   n=No. of individuals; m=no. of covariates; variables for adjustment }
  \item{lim}{numeric; threshold for skipping haplotypes from analysis }
  \item{maxSNP}{integer; Number of SNPs maximal group to multilocus genotypes }
  \item{nt}{integer; Number of notice best hits (for every step) }
  \item{pair.begin}{If true then will be begin with first 2 SNP genotypes.
                                      Attention: k SNP lead to choose(k,2)=k*(k-1)/2 possible pairs }
  \item{pattern.begin.mat}{if begin.pattern.mat is not NA then is this starting point of \code{msr}
                                 n=No. of snp pattern, m=No. of SNPs}
  \item{type}{type of depending variable}
  \item{baseline.hap}{Choose baseline haplotype for statistical test to avoid singularity. "max" for most frequent haplotype and "min" for less frequent haplotype }
  \item{min.count}{minimal count of rare haplotypes}
  \item{sort}{A logical value (TRUE or FALSE). If TRUE, family data will be sorted.}
}
\details{

Haplotypes are infered by EM algorithm (Excoffier and Slatkin 1995).
Family haplotypes are inferred by modified EM algorithm proposed by Rohde (2001, 2003).

For normal distributed phenotypes from independent individuals we prefer an F test and for case
control data we prefer the likelihood ratio test (logistic regression) in comparison of full model with
genetic and non-genetic factors to a reduced model, which includes only non-genetic variables. In
the case of no specified non-genetic variable only the intercept is used. If one of these tests are
significance we assume a genetic effect. In case of family data the weigthed TDT statistic is used.

 The procedure of multi-locus stepwise regression could be time consuming.
}
\value{
  \code{msr} provides a list with maxSNP components.
  \item{list}{for every step one component:
              SNP numbers and test details like p value.}
}
\author{ Sven Knueppel and Klaus Rohde }
\references{
Excoffier L, Slatkin M. Mol Biol Evol. 1995 Sep;12(5):921-7.

Rohde K, Fuerst R. Hum Hered. 2003;56(1-3):41-7.

Rohde K, Fuerst R. Hum Mutat. 2001 Apr;17(4):289-95.

Knueppel S, Rohde K et al. Multi-locus stepwise regression: A haplotype-based
stepwise algorithm to find haplotypes associated with genetic phenotypes.
[in praparation]
                                          }
\keyword{ stepwise regression }
\keyword{ haplotypes }

