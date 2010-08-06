\name{dec2bin}
\alias{dec2bin}
\title{
Decimal To Binary Conversion
}
\description{
\code{dec2bin} function converts a decimal number to a binary number.
}
\usage{
dec2bin ( vec, npos = NA )
}
\arguments{
  \item{vec}{
     	a numeric vector of positive values.
}
  \item{npos}{
      an optional number of length of the generating binary number.
}
}
\details{
This is a function to converting from decimal to binary number.
}
\value{
  \code{dec2bin} returns a matrix.
}
\examples{
binary <- dec2bin ( zz <- sample(0:100,10) )
print(zz)
print(binary)
}
\keyword{conversion}
\keyword{decimal}
\keyword{binary}
