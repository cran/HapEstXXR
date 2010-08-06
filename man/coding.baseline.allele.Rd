\name{coding.baseline.allele}
\alias{coding.baseline.allele}
\title{
  Standardization of coding alleles
}
\description{
  Dependend on minor allele frequency the coding of the allele will be update.
}
\usage{
coding.baseline.allele(geno, coding = c("minor", "major"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geno}{
 (m,n)-genotype matrix

                m=number of individuals

                type R: n=number of snps
}
  \item{coding}{
which type of coding should be used.
}
}
\details{
Allele 1 is coded as the minor allele , if coding type "minor" is used.
Otherwise allele 1 is coded as major allele.
}
\value{
This function returns the updated genotype matrix.
}
\author{
Sven Knueppel
}
