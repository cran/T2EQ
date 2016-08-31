\name{T2EQ.dissolution.profiles.hoffelder}
\alias{T2EQ.dissolution.profiles.hoffelder}
\concept{multivariate statistics}
\concept{equivalence}
\concept{two-sample design}
\concept{Mahalanobis distance}
\concept{Hotelling's T^2}
\concept{dissolution profiles}

\title{
The \eqn{T^2}-test for equivalence for dissolution data
}
\description{
The function \code{T2EQ.dissolution.profiles.hoffelder()} implements a variant of the \eqn{T^2}-test for equivalence analyses of highly variable dissolution profiles (see Hoffelder,2016). It is a multivariate two-sample equivalence procedure. Distance measure of the test is the Mahalanobis distance.
}
\usage{
T2EQ.dissolution.profiles.hoffelder(X, Y, alpha = 0.05, print.results = TRUE)
}
\arguments{
  \item{X}{
numeric data matrix of the first sample (REF). The rows of \code{X} contain the individual observations of the REF sample, the columns contain the variables/components of the multivariate sample. More precisely, the variables are the measured dissolution time points and the rows contain the individual dissolution profiles.
}
  \item{Y}{
numeric data matrix of the second sample (TEST). The rows of \code{Y} contain the individual observations of the TEST sample, the columns contain the variables/components of the multivariate sample. More precisely, the variables are the measured dissolution time points and the rows contain the individual dissolution profiles.
}
  \item{alpha}{
numeric (0<\code{alpha}<1). The significance level of the test. Usually set to 0.05 which is the default.
}
  \item{print.results}{
logical; if TRUE (default) summary statistics and test results are printed in the output. If NO no output is created 
}
}
\details{
This function implements a variant of the \eqn{T^2}-test for equivalence suggested in Hoffelder (2016): The equivalence margin of the test is a compromise between the suggestions of Tsong et al. (1996) and EMA (2010) requirements. See Hoffelder (2016) for a discussion on that equivalence margin.   
}
\value{
a data frame; three columns containing the results of the test
  \item{p.value}{numeric; the p-value of the equivalence test according to Hoffelder (2016)}
  \item{testresult.num}{numeric; 0 (null hypothesis of nonequivalence not rejected) or 1 (null hypothesis of nonequivalence rejected, decision in favor of equivalence)}
  \item{testresult.text}{character; test result of the test in text mode}
}

\references{
Hoffelder, T. (2016). Highly Variable Dissolution Profiles: Comparison of \eqn{T^2}-Test for Equivalence and \eqn{f_2} Based Methods. \emph{pharmind}, 78:4, 587-592.
URL: \url{http://www.ecv.de/suse_item.php?suseId=Z|pi|8430}

Wellek, S. (2010), \emph{Testing Statistical Hypotheses of Equivalence and Noninferiority. Second edition.} Boca Raton: Chapman & Hall/CRC.

Tsong, Y., Hammerstrom, T., Sathe, P., Shah, V.P. (1996). Statistical Assessment of Mean Differences between two Dissolution Data Sets. \emph{Drug Information Journal}, 30:4, 1105-1112. 
URL: \url{http://dx.doi.org/10.1177/009286159603000427}

EMA (2010). Guidance on the Investigation of Bioequivalence.
URL: \url{http://www.ema.europa.eu/docs/en_GB/document_library/Scientific_guideline/2010/01/WC500070039.pdf}
}
\author{
Thomas Hoffelder <thomas.hoffelder at boehringer-ingelheim.com>
}
\examples{
\dontrun{A recalculation of the results underlying Figure 1 in Hoffelder (2016) 
can be done with the following code:}

data(ex_data_pharmind)
REF_pharmind <- cbind(ex_data_pharmind[ which(ex_data_pharmind$Group=='REF'), ]
                  [c("Diss_10_min","Diss_20_min","Diss_30_min")])
TEST_pharmind <- cbind(ex_data_pharmind[ which(ex_data_pharmind$Group=='TEST'), ]
                  [c("Diss_10_min","Diss_20_min","Diss_30_min")])
test_T2EQ.dissolution.profiles.hoffelder_pharmind <- 
      T2EQ.dissolution.profiles.hoffelder(X=REF_pharmind,Y=TEST_pharmind)
}

\keyword{multivariate}
\keyword{robust}