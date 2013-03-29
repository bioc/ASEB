setGeneric("asebSites",
           function(backgroundSites, prodefinedSites, testSites, outputFile=NULL, permutationTimes=10000) standardGeneric("asebSites"), signature=c("backgroundSites", "prodefinedSites", "testSites"))
setGeneric("asebProteins",
           function(backgroundSites, prodefinedSites, testProteins, outputFile=NULL, permutationTimes=10000) standardGeneric("asebProteins"), signature=c("backgroundSites", "prodefinedSites", "testProteins"))
setGeneric("drawStat",
           function(curveInfoDataFrame, proteinIds=NULL, outputDir=NULL, figKind=c("pdf","jpeg")) standardGeneric("drawStat"), signature="curveInfoDataFrame")
setGeneric("drawEScurve",
           function(curveInfoDataFrame, sites=NULL, max_p_value=0.1, min_es=0.2, outputDir=NULL, figKind=c("pdf","jpeg")) standardGeneric("drawEScurve"), signature="curveInfoDataFrame")

setGeneric("sequences", function(object) standardGeneric("sequences"))

setGeneric("ids", function(object) standardGeneric("ids"))

setGeneric("sequences<-", function(object, value) standardGeneric("sequences<-"))

setGeneric("ids<-", function(object, value) standardGeneric("ids<-"))

