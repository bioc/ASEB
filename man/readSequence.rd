\name{readSequence}
\alias{readSequence}
\title{read sequences from file}
\description{
   This function is used to read sequences from FASTA format file.  
}
\usage{
readSequence(file)
}
\arguments{
 \item{file}{character(1), file name for input. \cr The input file should follow FASTA format.}
}
\details{
  This function return an object of \code{\link{SequenceInfo}} that 
  contains sequences and identifiers from FASTA format input file.
}
\value{
  A \code{\link{SequenceInfo}} object containing sequences and identifiers from input file. 
}
\seealso{
 \code{\link{SequenceInfo}},
 \code{\link{asebSites}},
 \code{\link{asebProteins}},
 \code{\link{drawStat}},
 \code{\link{drawEScurve}}.
}
\examples{
   ff <- system.file("extdata", "background_sites.fa", package="ASEB")
   readSequence(ff)
}
\keyword{methods}
