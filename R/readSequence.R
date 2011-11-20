readSequence <- function(file){
   rd <- read.table(file, fill=TRUE,sep="")
   rd <- as.matrix(rd)
   sequences <- ""
   ids <- ""
   id <- ""
   sequence <- ""
   count <- 0
   for(i in (1:nrow(rd))){
       line <- as.vector(rd[i,])
       if(length(grep(">", line[1])) > 0){
          if(id != ""){
             count <- count+1
             ids[count] <- id
             sequences[count] <- sequence
          }
          if(nchar(line[1]) > 1){
             id <- substring(line[1],2) 
          }else{
             id <- line[2]
          }
          sequence <- ""
       }else{
          for(j in (1:length(line))){
              if(is.na(line[j])){
                 break
              }
              sequence <- paste(sequence, line[j], sep="")
              sequence <- gsub("\\s+","",sequence)
          }
       }
   }
   if(id != ""){
      count <- count+1
      ids[count] <- id
      sequences[count] <- sequence
   }
   SequenceInfo(sequences, ids)
}