\name{msr.meta}
\alias{msr.meta}
\alias{msr.meta.binomial.unadjusted.families.unadjusted}
\alias{msr.binomial.haplotype.test.unadjusted msr.families.unadjusted}
\alias{msr.gaussian.forward.adjusted msr.gaussian.forward.unadjusted}
\alias{msr.meta.binomial.unadjusted.families.unadjusted msr.meta.error.check}
\alias{msr.meta.error.check}
\alias{meta.func}
\alias{meta.fisher}
\alias{meta.rank}
\title{
Multi-locus stepwise regression over two studies (as meta analysis)
}
\description{
\code{msr.meta} can be used as multi-locus stepwise regression over two studies.
}
\usage{
msr.meta ( data1 , control.options1 , data2 , control.options2 ,
    global.options = list( maxSNP=3 , lim=0.05 ,
    meta.method=c("fisher","rank") ) )
}
\arguments{
  \item{data1}{list; containing at least trait, genotypes, and covariates}
  \item{control.options1}{containing at least study type.}
  \item{data2}{list; containing at least trait, genotypes, and covariates}
  \item{control.options2}{containing at least study type.}
  \item{global.options}{list; Defining maxSNP,lim, nt, meta.method }
}
\details{
 The options you can define are listed in \code{\link{msr}}.
 
 
 Different methods to combine p values from independent samples
 
 
 (1) Fisher method: T = -2 * sum_i (  log(p_i)  ) ~ chi2 (2*k) ,
                    where k=sum of p values
                    
                    
 (2) Rank: sum of rank
 
 So far only combination of a case control study (binomial) and a family study
 (families) is supported.
}
\value{
  \code{msr} provides a list with maxSNP components.
  \item{list}{for every step one component:
              SNP numbers and test details like p value, rank of p values.}
}
\references{
Knueppel S, Rohde K et al. Multi-locus stepwise regression: A haplotype-based
stepwise algorithm to find haplotypes associated with genetic phenotypes.
[in praparation]

}
\author{
  Sven Knueppel and Klaus Rohde
}

\seealso{ \code{\link{msr}} }

\keyword{ meta analysis }
\keyword{ multi-locus stepwise regression }
\keyword{ msr }

