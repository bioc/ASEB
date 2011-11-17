setMethod(asebSites, c("character","character","character"),
function(backgroundSites, prodefinedSites, testSites, outputFile, permutationTimes=10000)
{   
    if(is.null(outputFile[1])){
       outputFile <- paste(tempdir(),"/","outputFile",sep="")
    }
    cat("Please wait patiently!\n")
    kk <- .C(".asebC",as.character(backgroundSites),as.character(prodefinedSites), as.character(testSites),as.character(outputFile),as.integer(permutationTimes),as.integer(0))
})

split2ShortMers <- function(longString, merLen=50, sep="\n"){
    i <- 1
    len <- nchar(longString)
    mers <- "";
    while(i <= len){
        mer50 <- substr(longString, i, i+merLen-1)
        if(i==1){
           mers <- mer50
        }else{
           mers <- paste(mers, mer50, sep=sep)
        }
        i <- i+merLen
    }
    mers
}
split2ShortMers2 <- function(longStrings, merLen=50, sep="\n"){
    sapply(longStrings, function(x) split2ShortMers(x,merLen=merLen, sep=sep))
}

setMethod(asebSites, c("SequenceInfo", "SequenceInfo", "SequenceInfo"),
function(backgroundSites, prodefinedSites, testSites, outputFile=NULL, permutationTimes=10000)
{
    backgroundSitesSeqs <- sequences(backgroundSites)
    backgroundSitesIds <- ids(backgroundSites)
    prodefinedSitesSeqs <- sequences(prodefinedSites)
    prodefinedSitesIds <- ids(prodefinedSites)
    testSitesSeqs <- sequences(testSites)
    testSitesIds <- ids(testSites)
    backgroundSitesFile <- paste(tempdir(),"/","backgroundSitesFile",sep="")
    prodefinedSitesFile<- paste(tempdir(),"/","prodefinedSitesFile",sep="")
    testSitesFile<- paste(tempdir(),"/","testSitesFile",sep="")
    if(is.null(outputFile[1])){
       outputFile <- paste(tempdir(),"/","outputFile",sep="")
    }
    
    backgroundSitesIds <- paste(">", backgroundSitesIds, sep = "")
    prodefinedSitesIds <- paste(">", prodefinedSitesIds, sep = "")
    testSitesIds <- paste(">", testSitesIds, sep = "")
        
    write.table(file=backgroundSitesFile, cbind(backgroundSitesIds,backgroundSitesSeqs),row.names = FALSE,col.names = FALSE,quote = FALSE,sep="\n")
    write.table(file=prodefinedSitesFile, cbind(prodefinedSitesIds,prodefinedSitesSeqs),row.names = FALSE,col.names = FALSE,quote = FALSE,sep="\n")
    write.table(file=testSitesFile, cbind(testSitesIds,testSitesSeqs),row.names = FALSE,col.names = FALSE,quote = FALSE,sep="\n")
    asebSites(backgroundSitesFile, prodefinedSitesFile, testSitesFile, outputFile, permutationTimes)
    curveInfoFile <- paste(outputFile, ".curves",sep="")
    results <- read.table(outputFile,header=FALSE)
    colnames(results)[1:3] <- c("site", "ES", "p-value")
    curveInfo <- read.table(curveInfoFile,header=FALSE,sep="\t")
    colnames(curveInfo)[1] <- c("site")
    return(list(results=results, curveInfo=curveInfo))
})
