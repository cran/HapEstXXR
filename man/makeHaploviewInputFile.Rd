%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% library(HapEstXXR)
% Created: November 30, 2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{makeHaploviewInputFile}
\alias{makeHaploviewInputFile}
\title{
Make Haploview input files
}
\description{
Create two data sets (*.ped and *.info) as input files for Haploview
}
\usage{
makeHaploviewInputFile(famid, patid, fid, mid, sex,
    aff, geno, marker.name,  marker.position,
    haploview.pedfile, haploview.infofile)
}
\arguments{
  \item{famid}{Family ID}
  \item{patid}{Individual ID}
  \item{fid}{Paternal ID}
  \item{mid}{Maternal ID}
  \item{sex}{1=male, 2=female, other=unknown}
  \item{aff}{disease phenotype (1=unaff, 2=aff, 0=missing/unkown)}
  \item{geno}{(n,m) genotype matrix (n=number of individuals,
  m=number of marker, 1-column for every marker, R-code:
  1 = 1/1, 3 = 1/2, 2 = 2/2); All markers should be biallelic.}
  \item{marker.name}{marker name}
  \item{marker.position}{marker position}
  \item{haploview.pedfile}{specify target of linkage file}
  \item{haploview.infofile}{specify target of marker Information file}
}
\details{
This function provides only limited options for creating Haploview input files.
For more details see Haploview/URL: http://www.broadinstitute.org/mpg/haploview.
}
\value{
no return values.
}
\references{
Barrett JC, Fry B, Maller J, Daly MJ. Haploview: analysis and visualization


of LD and haplotype maps. Bioinformatics. 2005; 21(2):263-265. [PubMed ID: 15297300]


Haploview/URL: http://www.broadinstitute.org/mpg/haploview
}
\author{
Sven Knueppel
}

\seealso{
\code{\link{makePlinkInputFile}}, \code{\link{allele1to2}}
}
