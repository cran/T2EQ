\name{T2EQ}
\alias{T2EQ}
\concept{multivariate statistics}
\concept{equivalence}
\concept{two-sample design}
\concept{Mahalanobis distance}
\concept{Hotelling's T^2}
\concept{multivariate normal distribution}

\title{
Function for applying the \eqn{T^2}-test for equivalence
}
\description{
The function \code{T2EQ()} implements the \eqn{T^2}-test for equivalence (see Wellek,2010 or Hoffelder et al., 2015). The \eqn{T^2}-test for equivalence is a multivariate two-sample equivalence test. 
Distance measure of the test is the Mahalanobis distance.
}
\usage{
T2EQ(X, Y, eq_margin, alpha = 0.05, print.results = TRUE)
}
\arguments{
  \item{X}{
numeric data matrix of the first sample. The rows of \code{X} contain the individual observations of the sample, the columns contain the variables/components of the multivariate sample. 
}
  \item{Y}{
numeric data matrix of the second sample. The rows of \code{X} contain the individual observations of the sample, the columns contain the variables/components of the multivariate sample. 
}
  \item{eq_margin}{
numeric (>0). The equivalence margin of the test. 
}
  \item{alpha}{
numeric (0<\code{alpha}<1). The significance level of the \eqn{T^2}-test for equivalence. Usually set to 0.05 which is the default.
}
  \item{print.results}{
logical; if TRUE (default) summary statistics and test results are printed in the output. If NO no output is created 
}
}
\details{
For multivariate normally distributed data the \eqn{T^2}-test for equivalence is exact and UMPI.
}
\value{
a data frame; three columns containing the results of the test
  \item{p.value }{numeric; the p-value of the \eqn{T^2}-test for equivalence}
  \item{testresult.num }{numeric; 0 (null hypothesis of nonequivalence not rejected) or 1 (null hypothesis of nonequivalence rejected, decision in favor of equivalence)}
  \item{testresult.text }{character; test result of the \eqn{T^2}-test for equivalence in text mode}
}
\references{
Wellek, S. (2010), \emph{Testing Statistical Hypotheses of Equivalence and Noninferiority. Second edition.} Boca Raton: Chapman & Hall/CRC.

Hoffelder, T., Goessl, R., Wellek, S. (2015). Multivariate Equivalence Tests for Use in Pharmaceutical Development. \emph{Journal of Biopharmaceutical Statistics}, 25:3, 417-437.
URL: \url{http://dx.doi.org/10.1080/10543406.2014.920344}
}
\author{
Thomas Hoffelder <thomas.hoffelder at boehringer-ingelheim.com>
}

\examples{
\dontrun{A recalculation of the example evaluation in Hoffelder et al. (2015) 
can be done with the following code:}

data(ex_data_JoBS)
REF_JoBS <- cbind(ex_data_JoBS[ which(ex_data_JoBS$Group=='REF'), ]
            [c("Diss_15_min","Diss_20_min","Diss_25_min")])
TEST_JoBS <- cbind(ex_data_JoBS[ which(ex_data_JoBS$Group=='TEST'), ]
            [c("Diss_15_min","Diss_20_min","Diss_25_min")])
equivalence_margin_JoBS <- 0.74^2
test_T2EQ_JoBS <- T2EQ(X=REF_JoBS,Y=TEST_JoBS,eq_margin = equivalence_margin_JoBS)
}

\keyword{multivariate}
\keyword{robust}

