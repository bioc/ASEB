setMethod(drawStat, "data.frame",
function(curveInfoDataFrame, proteinIds=NULL, outputDir=NULL, figKind=c("pdf","jpeg"))
{
  valid_rows <- (1:(nrow(curveInfoDataFrame)/2))*2-1
  rt <- as.matrix(curveInfoDataFrame[valid_rows,c(1,ncol(curveInfoDataFrame)-1,ncol(curveInfoDataFrame))])
  positions <- as.vector(rt[,2])
  if((!is.null(outputDir[1]))||(length(outputDir) > 1)){
      figKind <- match.arg(figKind)     
      dir.create(outputDir[1], showWarnings = FALSE, recursive = TRUE)
      #cat("create ", outputDir[1], "\n")
  }else{
      figKind <- "none"
  }
  for(i in (1:(nrow(rt)))){
      proteinName <- as.vector(rt[i,1])
      stokens <- unlist(strsplit(proteinName, "_"))
      site <- stokens[length(stokens)]
      stop <- nchar(proteinName)-nchar(site)-1
      rt[i,1]<-substr(proteinName, 1, stop)
      positions[i]<- as.numeric(site)
  }
  if((is.null(proteinIds[1]))&&(length(proteinIds) == 0)){
      proteinIds <- as.vector(rt[,1])
  }
  name2name <- c()
  name2name[proteinIds] <- proteinIds
  for(i in (1:length(name2name))){
      valid_rows <- rt[,1]==name2name[i]
      es_scores <- rt[valid_rows,2]
      p_values <- rt[valid_rows,3]
      x_values <- as.numeric(positions[valid_rows])
      if(figKind == "none"){
         #dev.new() ;
      }
      if(figKind == "pdf"){
         pdf(paste(outputDir,"/",name2name[i],".pdf",sep=""))
      }
      if(figKind == "jpeg"){
         jpeg(paste(outputDir,"/",name2name[i],".jpeg",sep=""))
      }
      plot(x_values, es_scores,pch=24,col=2,ylim=c(0,1),lwd=3,xlab=paste("sites of protein ", name2name[i],sep=""),ylab="ES & p-value")
      lines(x_values, es_scores,pch=24,col=2, lwd=2)
      points(x_values,p_values,pch=23,col=3, ylim=c(0,1),lwd=3)
      lines(x_values,p_values,pch=23,col=3, lwd=2)
      points((max(x_values)-min(x_values))/10+min(x_values), 1.0, pch=23, lwd=3, col=3)
      text((max(x_values)-min(x_values))/10*1.3+min(x_values), 1.0,label="nominal p-value",adj=0)
      points((max(x_values)-min(x_values))/10+min(x_values), 0.95, pch=24, lwd=3, col=2)
      text((max(x_values)-min(x_values))/10*1.3+min(x_values), 0.95,label="Enrichment Score",adj=0)
      if(figKind == "pdf"){
         dev.off()
      }
      if(figKind == "jpeg"){
         dev.off()
      }
  }
})
