\name{drawStat}
\alias{drawStat}
\alias{drawStat,data.frame-method}
\title{draw P-values and enrichment scores for all lysine sites on a specific protein}
\description{
   This function is used to show P-values and enrichment scores for all lysine sites on a specific protein.
}
\usage{
drawStat(curveInfoDataFrame, proteinIds=NULL, outputDir=NULL, figKind=c("pdf","jpeg"))
\S4method{drawStat}{data.frame}(curveInfoDataFrame, proteinIds=NULL, outputDir=NULL, figKind=c("pdf","jpeg"))
}
\arguments{
\item{curveInfoDataFrame}{data.frame object: contains curve information for all proteins, see example.}
\item{proteinIds}{character vector: only draw curves for proteins with ids appear in this vector. \cr
                  default: draw curves for all the proteins appear in curveInfoDataFrame.}
\item{outputDir}{character(1), output directory name for all the figures.}
\item{figKind}{character(1), fig format: \code{"pdf"} or \code{"jpeg"}.}
}
\details{
  This function is used to draw P-values and enrichment scores for all lysine sites on a specific protein.
  The X-axis shows positions of all lysine sites on a specific protein, and Y-axis shows the enrichment scores (0~1) and P-values (0~1)
  for each lysine site. The data.frame object contains curve information is given by \code{\link{asebProteins}}. 
}
\seealso{
 \code{\link{SequenceInfo}},
 \code{\link{readSequence}},
 \code{\link{asebSites}},
 \code{\link{asebProteins}},
 \code{\link{drawEScurve}}.
}
\examples{
    backgroundSites <- readSequence(system.file("extdata", "background_sites.fa", package="ASEB")) 
    prodefinedSites <- readSequence(system.file("extdata", "predefined_sites.fa", package="ASEB"))
    testProteins <- readSequence(system.file("extdata", "proteins_to_test.fa", package="ASEB"))
    resultList <- asebProteins(backgroundSites, prodefinedSites, testProteins, permutationTimes=100)
    #drawEScurve(resultList$curveInfo, max_p_value=0.5, min_es=0, outputDir=tempdir(), figKind="jpeg")
    drawStat(resultList$curveInfo, outputDir=tempdir(), figKind="jpeg");
    cat("see figures in output dir:", tempdir(),"\n")
}
\keyword{methods}
