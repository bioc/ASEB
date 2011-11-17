\name{drawEScurve}
\alias{drawEScurve}
\alias{drawEScurve,data.frame-method}
\title{draw Enriment Score curves for specific sites}
\description{
   This function is used to draw Enriment Score curves for specific sites.
}
\usage{
drawEScurve(curveInfoDataFrame, sites=NULL, max_p_value=0.1, min_es=0.2, outputDir=NULL, figKind=c("pdf","jpeg"))
\S4method{drawEScurve}{data.frame}(curveInfoDataFrame, sites=NULL, max_p_value=0.1, min_es=0.2, outputDir=NULL, figKind=c("pdf","jpeg"))
}
\arguments{
\item{curveInfoDataFrame}{data.frame object: contains curve information for sites, see example.}
\item{sites}{character vector: only draw curves for sites with ids appear in this vector. \cr
             default: draw curves for all the sites appear in curveInfoDataFrame.}
\item{max_p_value}{numeric(1), only draw curves for sites with p-value less than this value.}
\item{min_es}{numeric(1), only draw curves for sites with Enriment Score more than this value.}
\item{outputDir}{character(1), output directory name for all the figures.}
\item{figKind}{character(1), fig format: \code{"pdf"} or \code{"jpeg"}.}
}
\details{
  This function is used to draw Enrichment Score curves for specific sites.
  These curves show running-sum process for calculating enrichment score.
  The data.frame object contains curve information is given by 
  \code{\link{asebSites}} or \code{\link{asebProteins}}. 
}
\seealso{
 \code{\link{SequenceInfo}},
 \code{\link{readSequence}},
 \code{\link{asebSites}},
 \code{\link{asebProteins}},
 \code{\link{drawStat}}.
}
\examples{
    backgroundSites <- readSequence(system.file("extdata", "background_sites.fa", package="ASEB")) 
    prodefinedSites <- readSequence(system.file("extdata", "predefined_sites.fa", package="ASEB"))
    testSites <- readSequence(system.file("extdata", "sites_to_test.fa", package="ASEB"))
    resultList <- asebSites(backgroundSites, prodefinedSites, testSites, permutationTimes=100)
    drawEScurve(resultList$curveInfo, max_p_value=0.1, min_es=0.1, outputDir=tempdir(), figKind="jpeg")
    cat("see figures in output dir:", tempdir(),"\n")
}
\keyword{methods}
