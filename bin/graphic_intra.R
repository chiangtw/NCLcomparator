#! /usr/bin/env Rscript

# Usage:
# ./graphic.R Merged.result OUTPUTname

args <- commandArgs(TRUE)
in_file <- args[1]
out_file <- args[2]

input <- read.table(in_file, sep="\t", header=TRUE)
outputName <- paste(out_file, "pdf", sep=".")
###########################################################################

if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
  library(RColorBrewer)
}


inputEventNum <- dim(input)[1]
inputLastcol <- dim(input)[2]
input2 <- data.matrix(input[,c(8:inputLastcol)])
inputToolsNum <- inputLastcol-7
inputToolsNames <- colnames(input2)

if(inputToolsNum >1)
{  
  input3 <- 1*(input2 > 0)
  input3.RS <- rowSums(input3)
  inputTools.p <- c()
  inputTools.len <- c()
  supportNum <- matrix(rep(0, inputToolsNum), nrow=1)
  for (i in 1:inputToolsNum)
  {
    toolone <- input3[input3[,i]==1,]
    toolone.len <- dim(toolone)[1]
    inputTools.len <- c(inputTools.len, toolone.len)
    
    p.one <- c()
    for (j in 1:inputToolsNum)
    { 
      p.1 <- sum(rowSums(toolone) ==j)/toolone.len
      p.one <- c(p.one, p.1)
    }
    inputTools.p <- rbind(inputTools.p, p.one)
    supportNum[,i]<- sum(1*(input3.RS==i))
   }
  
  rownames(inputTools.p) <- inputToolsNames
  input3.RS <- rowSums(input3)
  
  
  colorIn <- brewer.pal(inputToolsNum, "Paired")
  maxY <- 2000+ max(inputTools.len)
  maxY2 <- 1000+max(c(supportNum))
  
  pdf(outputName, width=12, height=5) 
  
  par(mfrow=c(1,2), mar=c(5,5,2,5))
  barplot(t(inputTools.p), col=colorIn , las=2, beside = FALSE, cex.lab=0.8, cex.axis=0.8, cex.names=0.7, yaxs="i", ylab="Proportion")
  par(new=T)
  plot(c(1:inputToolsNum), inputTools.len, type="o" ,axes=F, xlab=NA, ylab=NA, ylim=c(0,maxY),yaxs="i" )
  axis(side=4, cex.axis=0.8, las=2)
  mtext(side=4, line=3, "Number of identified intragenic NCL events (circle point)", cex=0.8)
  barplot(c(supportNum), names.arg=c(1:inputToolsNum) ,col=colorIn, cex.lab=0.8,cex.axis=0.8, cex.names=1, cex.lab=0.8, yaxs="i", ylim=c(0, maxY2), ylab="Number of identified intragenic NCL events")
  legend(inputToolsNum/2,(maxY2-1000), c(1:inputToolsNum), col=colorIn, pch=15, bty = "n", cex=0.8, pt.cex = 1.2)
  text(3*(inputToolsNum/4), (maxY2-900), "Number of supporting tools:", cex=0.8)
  dev.off()
}  









