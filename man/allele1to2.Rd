 \name{allele1to2}
\alias{allele1to2}
\alias{allele2to1}
\alias{alleleRto1}
\alias{alleleRto2}
\alias{allele1toR}
\alias{allele2toR}
\title{ Convert genotype matrix from two different types }
\description{
(not supported for x-linked markers)
}
\usage{
allele1to2( geno , marker.label=NULL , miss.val=NA )
allele2to1( geno , marker.label=NULL , miss.val=NA )
alleleRto1( geno , marker.label=NULL , miss.val=NA)
alleleRto2( geno , marker.label=NULL , miss.val=NA )
allele1toR( geno , marker.label=NULL , miss.val=c(-1,NA) )
allele2toR( geno , marker.label=NULL , miss.val=NA )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geno}{ (m,n)-genotype matrix
  
                m=number of individuals
                
                type 1 and R: n=number of snps
                
                type 2: n=2*number of snps }
  \item{marker.label}{ Vector of labels for marker, If a marker name is "SNP", its columns will be "SNP.1" and "SNP.2"  }
  \item{miss.val}{Vector of specified missing values.}
}
\details{
  3 different types of genotype matrices:
  
  Type 1 : 1-column genotype matrix : minor allele count (0,1,2)
  
  Type 2 : 2-column genotype matrix : each marker has a pair of two columns  (1/1, 1/2, 2/2)
  
  Type R : 1-column genotype matrix : code  (1 = 1/1, 3 = 1/2, 2 = 2/2)
}
\value{
converted genotype matrix
}
\author{ Sven Knueppel }

\examples{
 ## [A] allele1to2
 N <- 10
 ns <- 4
 (geno <- matrix(sample(c(NA,0:2), N*ns,replace=TRUE),nc=ns))
 allele1to2 (geno)

 ## [B] allele2to1
 (geno <- matrix ( c(0,0,1,1,2,1,1,2,
                     1,1,2,2,2,2,1,2,
                     0,0,1,1,2,1,0,0) ,nc=4,byrow=TRUE))
 allele2to1 (geno)

 ## [C] alleleRto1
 N <- 10
 ns <- 4
 (geno <- matrix(sample(c(NA,1:3), N*ns,replace=TRUE),nc=ns))
 alleleRto1 (geno)

 ## [D] alleleRto2
 N <- 10
 ns <- 4
 (geno <- matrix(sample(c(0,1:3), N*ns,replace=TRUE),nc=ns))
 alleleRto2 (geno)

 ## [E] allele1toR
 N <- 10
 ns <- 4
 (geno <- matrix(sample(c(NA,0:2), N*ns,replace=TRUE),nc=ns))
 allele1toR (geno)

 ## [F] allele2toR
 (geno <- matrix ( c(0,0,1,1,2,1,1,2,1,1,2,2,2,2,1,2,0,0,1,1,2,1,0,0) ,
     nc=4,byrow=TRUE))
 allele2toR (geno)
}
\keyword{ convert }
