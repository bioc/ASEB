\name{asebSites}
\alias{asebSites}
\alias{asebSites,SequenceInfo,SequenceInfo,SequenceInfo-method}
\alias{asebSites,character,character,character-method}
\title{prediction of lysine sites that can be acetylated}
\description{
   This function is used to predict lysine sites that can be acetylated by a specific KAT-family.  
}
\usage{
asebSites(backgroundSites, prodefinedSites, testSites, outputFile=NULL, permutationTimes=10000)
\S4method{asebSites}{character,character,character}(backgroundSites, prodefinedSites, testSites, outputFile=NULL, permutationTimes=10000)
\S4method{asebSites}{SequenceInfo,SequenceInfo,SequenceInfo}(backgroundSites, prodefinedSites, testSites, outputFile=NULL, permutationTimes=10000)
}
\arguments{
 \item{backgroundSites}{\link{SequenceInfo} object or file name (character(1)) for background peptides set.}
 \item{prodefinedSites}{\link{SequenceInfo} object or file name (character(1)) for KAT special peptides set.}
 \item{testSites}{\link{SequenceInfo} object or file name (character(1)) for query peptides set.}
 \item{outputFile}{file name for output (character(1)).}
 \item{permutationTimes}{permutation times (integer(1)), default and recommended: 10000.}
}
\details{
  This function is used to predict lysine sites that can be acetylated by a specific KAT-family.
  It will give an enrichment score and a P-value for each query lysine site.
  The whole process is similar with the GSEA method (permuting gene sets). Please see the references for details. \cr
  
  The first three arguments of method asebSites can be \link{SequenceInfo} objects or file names.
  If these arguments are \link{SequenceInfo} objects, this method returns a list to the users besides an output file.
  Otherwise, this method processes the FASTA format files directly and outputs all results to a file. 
  In this case, this method can process huge number of sites each time without loading any sequences to R workspace.
}
\value{
  The output file contains enrichment scores and P-values for each query site.
  The \code{asebSites,SequenceInfo,SequenceInfo,SequenceInfo-method} also returns a list contains two \code{data.frame} objects: \code{results} and \code{curveInfo}.
  \item{results}{contains enrichment scores and P-values for each query site.}
  \item{curveInfo}{contains information for enrichment score curves.}
}
\note{
  The acetylated lysine sites and their surrounding amino acids (8 on each side) are treated as acetylated peptides. \cr
  Example for peptides sequence :  "KEHDDIFDKLKEAVKEE". \cr
  All input file should follow FASTA format. 
}
\references{
 Subramanian, A. et al. (2005) Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. 
 \emph{Proc Natl Acad Sci U S A}, \bold{102}, 15545-15550.
 
 Mootha, V.K. et al. (2003) PGC-1alpha-responsive genes involved in oxidative phosphorylation are coordinately downregulated in human diabetes.
 \emph{Nat Genet}, \bold{34}, 267-273.
 
Guttman, M. et al. (2009) Chromatin signature reveals over a thousand highly conserved large non-coding RNAs in mammals. 
\emph{Nature},  \bold{458}, 223-227.

Li, T.T. et al. Characterization and prediction of lysine (K)-acetyl-transferase (KAT) specific acetylation sites.
\emph{Mol Cell Proteomics}, \bold{in press}.
}
\seealso{
 \code{\link{SequenceInfo}},
 \code{\link{readSequence}},
 \code{\link{asebProteins}},
 \code{\link{drawStat}},
 \code{\link{drawEScurve}}.
}
\examples{
    backgroundSites <- readSequence(system.file("extdata", "background_sites.fa", package="ASEB")) 
    prodefinedSites <- readSequence(system.file("extdata", "predefined_sites.fa", package="ASEB"))
    testSites <- readSequence(system.file("extdata", "sites_to_test.fa", package="ASEB"))
    resultList <- asebSites(backgroundSites, prodefinedSites, testSites, permutationTimes=100)
    resultList$results[1:2,]
}
\keyword{methods}
