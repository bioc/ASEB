SequenceInfo <- function(sequences, ids = character(0)){
   if(length(ids) == 0){
      ids <- sequences
   }
   if(length(ids) != length(sequences)){
      ids <- sequences
      cat("Waring: length of ids is not equal to sequences!\n")
   }
   object <- new("SequenceInfo")
   sequences(object) <- sequences
   ids(object) <- ids
   object
}

setMethod(show, "SequenceInfo", function(object) {
   cat("object of SequenceInfo\n")
   cat("total", length(sequences(object)), "sequences\n")
   showRows <- length(sequences(object))
   if(showRows > 10){
      showRows = 10
   }
   cat("first", showRows, "sequences\n")
   cat("Slot \"ids\":\n")
   show(ids(object)[1:showRows])
   cat("Slot \"sequences\":\n")
   show(sequences(object)[1:showRows])
})


setMethod("sequences", "SequenceInfo", function(object) slot(object, "sequences"))
setMethod("ids", "SequenceInfo", function(object) slot(object, "ids"))

setReplaceMethod("sequences", signature("SequenceInfo", "character"),
  function(object, value) {
     slot(object,"sequences")<- value
     object
  })
  
setReplaceMethod("ids", signature("SequenceInfo", "character"),
  function(object, value) {
     slot(object,"ids")<- value
     object
  })
