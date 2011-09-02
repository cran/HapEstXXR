%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% library(HapEstXXR)
% Created: November 30, 2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{makePlinkInputFile}
\alias{makePlinkInputFile}
\title{
Make PLINK input files
}
\description{
Create two data sets (*.ped and *.map) as input files for PLINK
}
\usage{
makePlinkInputFile ( famid , patid , fid , mid , sex , trait ,
    CHR , SNP , POS , geno.matrix , linkage.file , map.file ,
    cov.file )
}
\arguments{
  \item{famid}{Family ID}
  \item{patid}{Individual ID}
  \item{fid}{Paternal ID}
  \item{mid}{Maternal ID}
  \item{sex}{1=male, 2=female, other=unknown}
  \item{trait}{disease phenotype (1=unaff, 2=aff, -9 or 0=missing/unkown)}
  \item{CHR}{chromosome}
  \item{SNP}{marker name}
  \item{POS}{marker position}
  \item{geno.matrix}{(n,m) genotype matrix (n=number of individuals, m=number
  of marker, 1-column for every marker, R-code: 1 = 1/1, 3 = 1/2, 2 = 2/2);
  All markers should be biallelic.}
  \item{linkage.file}{specify target of linkage file}
  \item{map.file}{specify target of map file}
  \item{cov.file}{specify target of cov file}
}
\details{
This function provides only limited options for creating PLINK input files.
For more details see PLINK/URL: http://pngu.mgh.harvard.edu/~purcell/plink/.
}
\value{
no return values.
}
\references{
Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR,
Bender D, Maller J, Sklar P, de Bakker PIW, Daly MJ & Sham PC (2007)
PLINK: a toolset for whole-genome association and population-based
linkage analysis. American Journal of Human Genetics, 81.
(PLINK/URL: http://pngu.mgh.harvard.edu/~purcell/plink)
}
\author{
Sven Knueppel
}
\seealso{
\code{\link{makeHaploviewInputFile}}, \code{\link{allele1to2}}
}
