%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% library(HapEstXXR)
% Created: November 30, 2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{powerset}
\alias{powerset}
\title{
Generating power set of a set
}
\description{
Generates the power set of a given set of values.
}
\usage{
 powerset(x, fileout = NA, only.file = FALSE)
}
\arguments{
  \item{x}{ a vector.}
  \item{fileout}{ a character string which contains name of the target file.}
  \item{only.file}{ a logical. Only a file is created, if true (default=FALSE).}
}
\details{
Suppose you have a set S. The power set is the set off all subsets of S, including empty set and S itself.
The number of elements of the power set is 2^(number of elements of S). You can save the powerset in a file,
if a filename \code{fileout} is specified.

Empty set will be excluded.
}
\value{
   \code{powerset} generates a list of all subsets of x, excluding empty set, if only.file=F.
}
\author{
  Sven Knueppel
}
\note{
   \code{powerset} is resctricted to vectors with maximum number of 15 elements.
   Using \code{only.file=T} you can create bigger powersets.
}
\examples{
ps <- powerset(1:10)
ps
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ powerset }
\keyword{ power set }
\keyword{ subset }% __ONLY ONE__ keyword per line

