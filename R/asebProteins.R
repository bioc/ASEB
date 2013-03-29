setMethod(asebProteins, c("character","character","character"),
function(backgroundSites, prodefinedSites, testProteins, outputFile, permutationTimes=10000)
{   
    if(is.null(outputFile[1])){
       outputFile <- paste(tempdir(),"/","outputFile",sep="")
    }
    cat("Please wait patiently!\n")
    kk <- .C(".asebC",as.character(backgroundSites),as.character(prodefinedSites), as.character(testProteins), as.character(outputFile),as.integer(permutationTimes),as.integer(1))
})

setMethod(asebProteins, c("SequenceInfo", "SequenceInfo", "SequenceInfo"),
function(backgroundSites, prodefinedSites, testProteins, outputFile=NULL, permutationTimes=10000)
{
    backgroundSitesSeqs <- sequences(backgroundSites)
    backgroundSitesIds <- ids(backgroundSites)
    prodefinedSitesSeqs <- sequences(prodefinedSites)
    prodefinedSitesIds <- ids(prodefinedSites)
    testProteinsSeqs <- sequences(testProteins)
    testProteinsIds <- ids(testProteins)
    backgroundSitesFile <- paste(tempdir(),"/","backgroundSitesFile",sep="")
    prodefinedSitesFile<- paste(tempdir(),"/","prodefinedSitesFile",sep="")
    testProteinsFile<- paste(tempdir(),"/","testProteinsFile",sep="")
    if(is.null(outputFile[1])){
       outputFile <- paste(tempdir(),"/","outputFile",sep="")
    }
    
    backgroundSitesIds <- paste(">", backgroundSitesIds, sep = "")
    prodefinedSitesIds <- paste(">", prodefinedSitesIds, sep = "")
    testProteinsIds <- paste(">", testProteinsIds, sep = "")
    testProteinsSeqs <- split2ShortMers2(testProteinsSeqs)
    
    write.table(file=backgroundSitesFile, cbind(backgroundSitesIds,backgroundSitesSeqs),row.names = FALSE,col.names = FALSE,quote = FALSE,sep="\n")
    write.table(file=prodefinedSitesFile, cbind(prodefinedSitesIds,prodefinedSitesSeqs),row.names = FALSE,col.names = FALSE,quote = FALSE,sep="\n")
    write.table(file=testProteinsFile, cbind(testProteinsIds,testProteinsSeqs),row.names = FALSE,col.names = FALSE,quote = FALSE,sep="\n")
    asebProteins(backgroundSitesFile, prodefinedSitesFile, testProteinsFile, outputFile, permutationTimes)
    curveInfoFile <- paste(outputFile, ".curves",sep="")
    results <- read.table(outputFile,header=FALSE,sep="\t")
    colnames(results)[1:3] <- c("site", "ES", "p-value")
    curveInfo <- read.table(curveInfoFile,header=FALSE)
    colnames(curveInfo)[1] <- c("site")
    return(list(results=results, curveInfo=curveInfo))
})
