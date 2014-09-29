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
\title{Multi-locus stepwise regression}
\description{
  Stepwise regression for snp selection and haplotype testing
}
\usage{
msr(snps, trait, famid, patid, fid, mid, 
    adj.var = NA, lim = 0.05, maxSNP = 3, 
    nt = 10, sort.by = "AICc", selection = 0,
    p.threshold = NA, 
    pair.begin = FALSE, pattern.begin.mat = NA, 
    type = "gaussian", 
    baseline.hap = "max", min.count = 10, sort = FALSE)
}
\arguments{
  \item{snps}{(n, m)-Matrix; n=No. of individuals; m=no. of SNPs; Rohde-Code}
  \item{trait}{numeric; Outcome, phenotype }
\item{famid}{vector; Identifier for every family; only need in case of type=families}
  \item{patid}{vector; Identifier for every individuals; only need in case of type=families}
  \item{fid}{vector; Identifier for father (0=unkown); only need in case of type=families}
  \item{mid}{vector; Identifier for mother (0=unkown); only need in case of type=families}
  \item{adj.var}{(n, m)-Matrix;   n = No. of individuals; m = no. of covariates; variables for adjustment }
  \item{lim}{numeric; threshold for skipping haplotypes from analysis }
  \item{maxSNP}{integer; Number of SNPs maximal group to multilocus genotypes }
  \item{nt}{integer; Number of notice best hits (for every step) }
  \item{sort.by}{the results in each step were sorted by "AIC", corrected ("AICc""), or p value ("p.value"). default = "AICc". }
  \item{selection}{0 = none, 1 = improve of the lowest corrected AIC (AICc) of the step before, 2 = improve of the lowest AIC of the step before, 3 = improve of p value, 4 = improve of best ten log10(p values), 5 = improve of the single AICc by adding one SNP to the noticed pattern}
  \item{p.threshold}{numeric vector; if global p value is lower than p.threshold[i], then the pattern will be stored for further processing. I indicates the number of SNPs. If your calculation should start with all pairwise SNPs, then p.threshold[1] will be not used but should be included.}
\item{pair.begin}{If true then will be begin with first 2 SNP genotypes.
                                      Attention: k SNP lead to choose(k, 2) = k * (k - 1) / 2 possible pairs }
  \item{pattern.begin.mat}{if begin.pattern.mat is not NA then is this starting point of \code{msr}
                                 n = No. of snp pattern, m = No. of SNPs}                                 
  \item{type}{type of depending variable}
  \item{baseline.hap}{Choose baseline haplotype for statistical test to avoid singularity.
                      "max" for most frequent haplotype and "min" for less frequent haplotype }
  \item{min.count}{minimal count of rare haplotypes. If the count of estimated haplotypes < min.count, then the combined rare haplotypes were excluded from the analysis of that specific pattern.}
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
              SNP numbers and test details. }
}
\author{Sven Knueppel and Klaus Rohde}
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
\keyword{stepwise regression}
