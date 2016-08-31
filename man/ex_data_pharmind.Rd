\name{ex_data_pharmind}
\alias{ex_data_pharmind}
\docType{data}
\title{
Example dataset from Hoffelder (2016)
}
\description{
Multivariate example dataset of dissolution profiles. Dataset consists of two three-dimensional samples. The names of the three variables are "Diss_10_min","Diss_20_min" and "Diss_30_min". Variable "Group" discriminates between first sample (\code{Group == "REF"}) and second sample (\code{Group == "Test"}). Sample size is 12 per group.  
}
\usage{data("ex_data_pharmind")}
\format{
  A data frame with 24 observations on the following 4 variables.
  \describe{
    \item{\code{Diss_10_min}}{a numeric vector}
    \item{\code{Diss_20_min}}{a numeric vector}
    \item{\code{Diss_30_min}}{a numeric vector}
    \item{\code{Group}}{a character vector}
  }
}
\details{
Example dataset underlying Figure 1 in Hoffelder (2016). 
}
\source{
Hoffelder, T. (2016), "Highly Variable Dissolution Profiles: Comparison of \eqn{T^2}-Test for Equivalence and \eqn{f_2} Based Methods", \emph{pharmind}, 78:4, 587-592.
}
\references{
URL: \url{http://www.ecv.de/suse_item.php?suseId=Z|pi|8430}
}
\examples{
data(ex_data_pharmind)
}
\keyword{datasets}
