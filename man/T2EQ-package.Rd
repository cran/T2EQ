\name{T2EQ-package}
\alias{T2EQ-package}
\concept{multivariate statistics}
\concept{equivalence}
\concept{two-sample design}
\concept{Mahalanobis distance}
\concept{Hotelling's T^2}
\concept{dissolution profiles}
\concept{multivariate normal distribution}

\docType{package}
\title{
Functions for Applying the \eqn{T^2}-Test for Equivalence
}
\description{
Contains functions for applying the \eqn{T^2}-test for equivalence.
The \eqn{T^2}-test for equivalence is a multivariate two-sample equivalence test. 
Distance measure of the test is the Mahalanobis distance.
For multivariate normally distributed data the \eqn{T^2}-test for equivalence is exact and UMPI.
The function T2EQ() implements the \eqn{T^2}-test for equivalence according to Wellek (2010).
The function T2EQ.dissolution.profiles.hoffelder() implements a variant of the \eqn{T^2}-test for equivalence according to Hoffelder (2016) for the equivalence comparison of highly variable dissolution profiles.
}
\details{
\packageIndices{T2EQ}
}

\author{
Thomas Hoffelder

Maintainer: Thomas Hoffelder <thomas.hoffelder@boehringer-ingelheim.com>
}
\references{
Wellek, S. (2010), \emph{Testing Statistical Hypotheses of Equivalence and Noninferiority. Second edition.} Boca Raton: Chapman & Hall/CRC.

Hoffelder, T., Goessl, R., Wellek, S. (2015). Multivariate Equivalence Tests for Use in Pharmaceutical Development. \emph{Journal of Biopharmaceutical Statistics}, 25:3, 417-437.
URL: \url{http://dx.doi.org/10.1080/10543406.2014.920344}

Hoffelder, T. (2016). Highly Variable Dissolution Profiles: Comparison of \eqn{T^2}-Test for Equivalence and \eqn{f_2} Based Methods. \emph{pharmind}, 78:4, 587-592.
URL: \url{http://www.ecv.de/suse_item.php?suseId=Z|pi|8430}

Tsong, Y., Hammerstrom, T., Sathe, P., Shah, V.P. (1996). Statistical Assessment of Mean Differences between two Dissolution Data Sets. \emph{Drug Information Journal}, 30:4, 1105-1112. 
URL: \url{http://dx.doi.org/10.1177/009286159603000427}

EMA (2010). Guidance on the Investigation of Bioequivalence.
URL: \url{http://www.ema.europa.eu/docs/en_GB/document_library/Scientific_guideline/2010/01/WC500070039.pdf}
}

\keyword{ package }
\keyword{multivariate}
\keyword{robust}

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
