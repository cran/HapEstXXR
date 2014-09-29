%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% library(HapEstXXR)
% Created: November 30, 2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{read.haploview}
\alias{read.haploview}

\title{
Read a haploview dataset
}
\description{
Data can be loaded in haploview format (linkage format) with columns of family,
individual, father, mother, gender  (1 = male, 2 = male), affected status
(0 = unkown, 1 = unaffected, 2 = affected), and genotypes(2 columns alleles).
}
\usage{
read.haploview(ped.file, map.file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ped.file}{Localizion of a pedigree file.}
  \item{map.file}{Localizion of a marker information file.}
}
\details{
 The marker information file should contain in the first column the marker name
 and the second column the physical position on the chromosome.
}
\value{
  \item{famid}{family identifier  }
  \item{patid }{individual identifier }
  \item{dad }{father identifier (0=unkown)}
  \item{mom }{mother identifier (0=unkown)}
  \item{sex}{sex (0=unkown,1=male,2=female) }
  \item{genotypes }{(n,m)-matrix; n=No. of individuals; m=No. of SNPs; 1-column allele dosis}
  \item{trait}{phenotype values}
  \item{marker.names}{marker.names}
  \item{marker.position}{Base-pair position on chromosome (base pair units)}
}
\references{
Barrett JC, Fry B, Maller J, Daly MJ. Haploview: analysis and visualization


of LD and haplotype maps. Bioinformatics. 2005; 21(2):263-265. [PubMed ID: 15297300]

Haploview/URL: http://www.broadinstitute.org/mpg/haploview}
\author{
Sven Knueppel
}
\seealso{
\code{\link{read.data}}
}
\keyword{ read.data }
\keyword{ haploview }

