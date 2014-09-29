%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% library(HapEstXXR)
% Created: November 30, 2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{read.data}
\alias{read.data}
\title{ Read data from different input files }
\description{
   Data can be loaded in different formats.
}
\usage{
read.data(filename, linkage = TRUE, map = NA)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{ the name of the file which the data are to be read from. }
  \item{linkage}{a logical value indicating whether the file is in linkage format.}
  \item{map}{Localizion of a map file.}
}
\details{
  1) single individuals (3-columns)\cr
  expected columns\cr
  Individual identifier \cr
  genotype STRING (1=homozygot (wildtype) 2=homozygot (variant) 3=heterozygote
             0=missing value) >> Example: "1223" "3023" \cr
  phenotype  \cr
  
  2) family data (4-columns)\cr
  expected columns\cr
  Family identifier \cr
  Individual identifier\cr
  genotype (1=homozygot (wildtype) 2=homozygot (variant) 3=heterozygote
             0=missing value) >> Example: "1223" "3023" \cr
  phenotype \cr
  Remark 1: patid should not be 0 because 0 is unkown value for fid and mid.\cr
  Remark 2: Families are sorted. First two person in a family are adults
  (father and then mother) and after that all children. \cr
  
  3) Linkage format is expected, if \code{linkage=TRUE} :\cr
  Family identifier\cr
  Individual identifier  \cr
  Father identifier (0=unkown) \cr
  Mother identifier (0=unkown) \cr
  Sex (0=unkown,1=male,2=female)\cr
  Affectation_status (0=unkown,1=unaffected,2=affected) or trait_value\cr
  Marker_genotypes (M1_A1 M1_A2 M2_A1 ...)
  >> only 1, 2, or 0 for missing values \cr
  
  4) map file (4-columns), if specified:\cr
  chromosome (1-22, only autosomes)  \cr
  snp identifier  \cr
  Genetic distance (morgans)  \cr
  Base-pair position (bp units)
}
\value{
  \item{famid}{family identifier  }
  \item{patid }{individual identifier }
  \item{fid }{father identifier (0=unkown)}
  \item{mid }{mother identifier (0=unkown)}
  \item{sex}{sex (0=unkown,1=male,2=female) }
  \item{genotypes }{(n,m)-matrix; n=No. of individuals; m=No. of SNPs; Klaus format}
  \item{trait}{phenotype values}
  \item{chr}{chromosome}
  \item{snp}{snp identifier or rs id}
  \item{pos}{Base-pair position on chromosome (base pair units)}
}
%\references{ ~put references to the literature/web site here ~ }
\author{Sven Knueppel}
\keyword{ read data }
\keyword{ input files }% __ONLY ONE__ keyword per line
