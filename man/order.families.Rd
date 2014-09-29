%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% library(HapEstXXR)
% Created: November 30, 2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{order.families}
\alias{order.families}
\title{
Ordering of nuclear family data
}
\description{
\code{order.families} returns a permutation which rearranges the families into
ascending famid, generation, and sex, if given. }
\usage{
order.families(famid, patid, fid, mid, sex = NA)
}
\arguments{
  \item{famid}{vector; Identifier for every family}
  \item{patid}{vector; Identifier for every individual}
  \item{fid}{vector; Identifier for father (0 = unkown)}
  \item{mid}{vector; Identifier for mother (0 = unkown)}
  \item{sex}{vector; Individuals' gender (1 = male, 2 = female, 0 = unkown }
}
\author{
Sven Knueppel
}
\examples{

fam <- as.character(c(c(1, 1, 1, 1), c(0, 0, 0, 0, 0)))
pid <- as.character(c(c(1, 2, 3, 4), c(7, 8, 9, 10, 11)))
mid <- as.character(c(c(3, 3, 0, 0), c(10, 10, 10, 0, 0)))
fid <- as.character(c(c(4, 4, 0, 0), c(11, 11, 11, 0, 0)))
sex <- as.character(c(c(0, 2, 2, 1), c(1, 1, 2, 2, 1)))

ordfam <- order.families (fam, pid, fid, mid, sex)
print((cbind(fam, pid, fid, mid, sex))[ordfam, ])

    
}
\keyword{ families }
\keyword{ sort }

