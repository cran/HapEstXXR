%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% library(HapEstXXR)
% Created: November 30, 2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{catt}
\alias{catt}
\title{Cochrane armitage trend test (CATT) for SNP genotypes}
\description{
Performs chi-squared test for SNP genotypes. By default, score is chosen as the number of alleles (0,1,2).
}
\usage{
 catt ( y , x , score=c(0,1,2) )
}
\arguments{
  \item{y}{Vector of trait values. \code{y} must have values of 1 for event, 0 for no event. }
  \item{x}{Vector of SNP genotypes, 1-column coding (SNP allele dosis: 0,1,2).}
  \item{score}{Group score.}
}
\details{
The Cochran-Armitage trend test is typically used in categorical data analysis when some categories are ordered. Here it is used as a genotype-based test for candidate gene association.
}
\value{
  \item{2x3-table}{Genotype distribution.}
  \item{chisq}{The value for the test statistic.}
  \item{df}{Degrees of freedom.}
  \item{p.value}{The p-value for the test.}
  \item{n.miss}{Number of individuals with missing values.}
}
\references{Sasieni PD. From genotypes to genes: doubling the sample size.Biometrics. 1997 Dec;53(4):1253-61.}
\author{Sven Knueppel}
\examples{
y <- sample ( c(0,1),100, replace=TRUE)
x <- sample ( c(0,1,2),100, replace=TRUE)
catt ( y , x )
}
\seealso{\code{\link{prop.trend.test}}}
\keyword{ catt }
\keyword{ cochrane }
\keyword{ armitage }
\keyword{ trend }
\keyword{ test }
\keyword{ SNP }
\keyword{ SNPs }

