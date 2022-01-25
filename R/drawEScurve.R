setMethod(drawEScurve, "data.frame",
function(curveInfoDataFrame, sites=NULL, max_p_value=0.1, min_es=0.2, outputDir=NULL, figKind=c("pdf","jpeg"))
{ 
  rt <- curveInfoDataFrame
  nrow(rt)/2
  name2name <- c()
  name2name[sites]<- sites
  if((!is.null(outputDir[1]))||(length(outputDir) > 1)){
      figKind <- match.arg(figKind)
      dir.create(outputDir[1], showWarnings = FALSE, recursive = TRUE)
      #cat("create ", outputDir[1], "\n")
      if(file.access(outputDir, mode = 0) != 0){
         cat("warning: can not creat ", outputDir)
      }else{
         #cat("create down", outputDir[1], "\n")
      }
  }else{
      figKind <- "none"
  }
  #cat(figKind,"\t")
  for(i in (1:(nrow(rt)/2))){
      ES <- rt[i*2-1,ncol(rt)-1]
      p_value <- rt[i*2-1,ncol(rt)]
      if((!is.null(sites[1]))||(length(sites) > 1)){
         proteinName <- as.vector(rt[i*2-1,1])
         if(is.null(name2name[proteinName])){
            next
         }
      }
      #cat(as.vector(rt[i*2-1,1]), " ES:",ES," p_value:",p_value,"\n")
      if(ES < min_es){
         next
      }
      if(p_value > max_p_value){
         next
      }
      indexes <- as.vector(rt[i*2-1,])
      values <- as.vector(rt[i*2,])
      potins_num <- length(indexes)-2
      max_x <- as.numeric(indexes[potins_num])
      max_y_index <- which.max(as.numeric(values[2:potins_num]))[1]
      max_y <- max(as.numeric(values[2:potins_num]))
      if(figKind == "none"){
         #dev.new() ;
      }
      if(figKind == "pdf"){
         proteinName <- as.vector(rt[i*2-1,1])
         pdf(paste(outputDir,"/",proteinName,".pdf",sep=""))
      }
      if(figKind == "jpeg"){
         proteinName <- as.vector(rt[i*2-1,1])
         jpeg(paste(outputDir,"/",proteinName,".jpeg",sep=""))
      }
      plot(as.numeric(indexes[2:potins_num]),y=as.numeric(values[2:potins_num]),type="l",ylim=c(-1,1),xlab=paste("site: ", indexes[1],sep=""),ylab="ES Score",col=2,lwd=3)
      #lines(c(0,max_x),c(0.9,0.9),lwd=3)
      lines(c(0,max_x),c(0,0),lty=2,col=8)
      lines(c(-1000,as.numeric(indexes[max_y_index+1])),c(max_y,max_y),lwd=1,lty=2)
      lines(c(as.numeric(indexes[max_y_index+1]),as.numeric(indexes[max_y_index+1])),c(max(-0.7,max_y-0.7),min(max_y+0.7,0.7)),lwd=1,lty=2)
      text(as.numeric(indexes[max_y_index+1])+0.05, max_y, max_y,pos=4)
      
      for(j in (2:potins_num)){
          if(j%%2==1){
             next
          }
          if(j <= max_y_index+1){
             lines(c(indexes[j],indexes[j]),c(0.8,0.95),col=rgb(0,0,200,maxColorValue = 255))
          }else{
             lines(c(indexes[j],indexes[j]),c(0.8,0.95),col=rgb(160,160,255,maxColorValue = 255))
          }
      }
      if(figKind == "pdf"){
         dev.off()
      }
      if(figKind == "jpeg"){
         dev.off()
      }
  }
})
