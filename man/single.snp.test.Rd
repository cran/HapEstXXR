%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% library(HapEstXXR)
% Created: November 30, 2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{single.snp.test}
\alias{single.snp.test}
\alias{single.snp.test.binomial}
\alias{single.snp.test.gaussian}
\alias{single.snp.test.families}
\alias{single.snp.test.casecohort}
\alias{single.snp.test.casecohort.prentice}
\title{
Regression analysis with single SNP genotypes as independent variable
}
\description{
 This function fits a generalized linear model with quantative, dichotomous or survival
 trait as dependent variable and one or more potential covariates.
 In case of family data the weighted TDT statistic is used.
}
\usage{
single.snp.test(snps, trait, adj.var = NULL, 
           type = c("gaussian", "binomial", "families", "casecohort"), 
           famid, patid, fid, mid, 
           start.time, stop.time, subcohort, stratvar = NA, robust = FALSE, 
           marker.label = NA, 
           prt = TRUE, ties = "efron")
}
\arguments{
  \item{snps}{
  Matrix of alleles, such that each locus has one column of alleles
  (R code: 1 = 1/1, 3 = 1/2, 2 = 2/2, 0 = missing). Rows contains alleles for each subject.
}
  \item{trait}{
  Vector of trait values. For case control data use type= "binomial",
  trait must have values of 1 for event, 0 for no event.
}
  \item{adj.var}{
  Matrix of (non-genetic) covariates used to adjust the regression model.
}
  \item{type}{
Character string defining type of trait, with values of gaussian, binomial, families, survival, and casecohort.}
  \item{famid}{vector; Identifier for every family; needed by type="families".}
  \item{patid}{vector; Identifier for every individual; needed by type="families" and type="casecohort". }
  \item{fid}{vector; Identifier for father (0=unkown); needed by type="families".}
  \item{mid}{vector; Identifier for mother (0=unkown); needed by type="families".}
 \item{start.time}{vector; age at the start of the follow-up.}
 \item{stop.time}{vector; age at the end of the follow-up.}
 \item{subcohort}{A logical value (TRUE or FALSE). If TRUE, the individual is in the subcohort.}
 \item{stratvar}{vector; names the variables that determine the stratification.ss}
 \item{robust}{A logical value (TRUE or FALSE). If TRUE, request the robust sandwich estimate.}
\item{marker.label}{
  Vector of labels for marker.
}
  \item{prt}{
A logical value (TRUE or FALSE). If TRUE, an overview is printed.}
\item{ties}{defines the handling of ties in case-cohort design: "efron" (default),"breslow","exact". }
}
\details{
For normal distributed phenotypes from independent individuals we prefer an F
test and for case control data we prefer the likelihood ratio test
(logistic regression) in comparison of full model with genetic and non-genetic
factors to a reduced model, which includes only non-genetic variables. In the
case of no specified non-genetic variable only the intercept is used. If one of
these tests are significance we assume a genetic effect.
In case of family data the weigthed TDT statistic is used.

So far SURVIVAL data is not supported.

Cox proportional hazards regression modified for case cohort designs according
to the Prentice method will be used by type="casecohort".


As genetic effect the allele dosis (0, 1, 2) is modelled.
}
\value{
\code{single.snp.test} returns an object of class data.frame containing the following components:
 \item{snp}{snp number}
 \item{N}{number of individuals}
 \item{type}{type of depending variable}
 \item{beta}{estimation of beta coefficient out of full regression model}
 \item{se(beta)}{estimation of standard error of beta coefficient out of full regression model}
 \item{exp(beta)}{Odds ratio=exp(beta.estimate) are calculated, if type = "binomial". In case of type = "casecohort" hazard ratio is calculated.}
 \item{lower.95}{lower limit of 95 \% confidence intervall for exp(beta).}
 \item{upper.95}{upper limit of 95 \% confidence intervall for exp(beta).}
 \item{aic}{Akaike's An Information Criterion (AIC) of full model}
}
\references{
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
  \code{\link{single.haplotype.test}}
}
\examples{
N <- 2000
nloci <- 14
set.seed(1234)
y <- sample(c(0, 1), N, replace = TRUE)
snp <- matrix(sample(c(1, 2, 3), N * nloci, replace = TRUE), 
              ncol = nloci)
colnames(snp) <- paste("SNP", 1:nloci, sep = "")

adj.var <- matrix(rnorm(N * 3), ncol = 3)
colnames(adj.var) <- paste("A", 1:3, sep = "")

sst <- single.snp.test(snps = snp, trait = y, adj.var = adj.var,
    type = "binomial", prt = TRUE)
}
\keyword{ snp }
\keyword{ haplotypes }
\keyword{ logistic regression }
\keyword{ weigthed TDT }
\keyword{ survival analysis }




