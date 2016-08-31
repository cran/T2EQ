\name{ex_data_JoBS}
\alias{ex_data_JoBS}
\docType{data}
\title{
Example dataset from Hoffelder et al. (2015)
}
\description{
Multivariate example dataset of dissolution profiles. Dataset consists of two three-dimensional samples. The names of the three variables are "Diss_15_min","Diss_20_min" and "Diss_25_min". Variable "Group" discriminates between first sample (\code{Group == "REF"}) and second sample (\code{Group == "Test"}). Sample size is 12 per group.   
}
\usage{data("ex_data_JoBS")}
\format{
  A data frame with 24 observations on the following 4 variables.
  \describe{
    \item{\code{Group}}{a factor with levels \code{REF} \code{TEST}}
    \item{\code{Diss_15_min}}{a numeric vector}
    \item{\code{Diss_20_min}}{a numeric vector}
    \item{\code{Diss_25_min}}{a numeric vector}
  }
}
\details{
Example dataset from Hoffelder et al. (2015).
}

\source{
Hoffelder, T., Goessl, R., Wellek, S. (2015), "Multivariate Equivalence Tests for Use in Pharmaceutical Development", \emph{Journal of Biopharmaceutical Statistics}, 25:3, 417-437.
}
\references{
URL: \url{http://dx.doi.org/10.1080/10543406.2014.920344}
}
\examples{
data(ex_data_JoBS)
}
\keyword{datasets}
