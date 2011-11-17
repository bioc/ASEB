readSequence <- function(file){
   ff <- Biostrings::readFASTA(file,strip.descs=TRUE, checkComments=TRUE)
   sequences <- sapply(ff, function(x) x$seq)
   ids <- sapply(ff, function(x) x$desc)
   SequenceInfo(sequences, ids)
}
