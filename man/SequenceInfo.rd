\name{SequenceInfo}
\docType{class}
\alias{SequenceInfo}
\title{"SequenceInfo" class}

\description{
  This class is used to store sequences and identifiers for lysine sites or proteins.
}

\section{Objects from the Class}{
  Objects from this class are created by constructor \code{SequenceInfo}, as outlined below.
}

\section{Slots}{
  \describe{
    \item{\code{sequences}:}{\code{"character"}
      containing sequences.}
    \item{\code{id}:}{\code{"character"} containing identifiers.}
  }
}

\section{Methods}{
  Constructor:
  \describe{
    \item{SequenceInfo}{\code{signature(sequences = "character", id = "character")}:
      Create a \code{SequenceInfo} object from sequences and their
      identifiers. The length of \code{id} must match that of \code{sequences}.}
  }
}
\seealso{
 \code{\link{readSequence}},
 \code{\link{asebSites}},
 \code{\link{asebProteins}},
 \code{\link{drawStat}},
 \code{\link{drawEScurve}}.
}
\examples{
showClass("SequenceInfo")
}
\keyword{classes}
