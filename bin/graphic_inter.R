#! /usr/bin/env Rscript

# Usage:
# ./graphic.R interMerged.result OUTPUTname
#
# Output files: interMerged_junction.result and OUTPUTname.pdf

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
input2 <- data.matrix(input[,c(9:inputLastcol)])
inputToolsNum <- inputLastcol-8
inputToolsNames <- colnames(input2)

## generate interMerged_junction.result
input.out <- c()
input.RS <- rowSums(1*(input[,-c(1:8)] >0))
for (i in 1:inputEventNum)
{
  E.one <- as.numeric(input[i,-c(1:8)])
  Jmed <- median(E.one)
  Jtau <- sum(1-log(E.one+1)/log(max(E.one)+1))/(length(E.one)-1)
  NCLscore <- -1*log10((Jtau+0.01)/(((Jmed)^2)+0.01))
  NumSupportedTools <- input.RS[i]
  one.out <- cbind(input[i,], NumSupportedTools, Jmed, Jtau, NCLscore)
  input.out <- rbind(input.out, one.out) 
}
write.table(input.out, file="interMerged_junction.result", append=FALSE, 
            quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)


## plot OUTPUTname.pdf 
if(inputToolsNum >1)
{  
  input3 <- 1*(input2 > 0)
  input3.RS <- rowSums(input3)
  supportNum <- matrix(rep(0, inputToolsNum), nrow=1)
  inputTools.p <- c()
  inputTools.len <- c()
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
    supportNum[1,i]<- sum(1*(input3.RS==i))
    
  }
  
  rownames(inputTools.p) <- inputToolsNames

  
  colorIn <- brewer.pal(inputToolsNum, "Paired")
  maxY <- 100+max(inputTools.len)
  maxY2 <- 100+max(c(supportNum))
  
  pdf(outputName, width=12, height=5) 
  par(mfrow=c(1,2), mar=c(5,5,2,5))
  barplot(t(inputTools.p), col=colorIn , las=2, beside = FALSE, cex.lab=0.8, cex.axis=0.8, cex.names=0.7, yaxs="i", ylab="Proportion")
  par(new=T)
  plot(c(1:inputToolsNum), inputTools.len, type="o" ,axes=F, xlab=NA, ylab=NA, ylim=c(0,maxY),yaxs="i" )
  axis(side=4, cex.axis=0.8, las=2)
  mtext(side=4, line=3, "Number of identified intergeinc NCL events (circle point)", cex=0.8)
  barplot(c(supportNum), names.arg=c(1:inputToolsNum), col=colorIn ,cex.lab=0.8,cex.axis=0.8, cex.names=1, cex.lab=0.8, yaxs="i", ylim=c(0, maxY2), ylab="Number of identified intergenic NCL events")
  legend(2*(inputToolsNum/3),(maxY2-50), c(1:inputToolsNum), col=colorIn, pch=15, bty = "n", cex=0.8, pt.cex = 1.2)
  text(2*(inputToolsNum/3), (maxY2-30), "Number of supporting tools:", cex=0.8)
  dev.off()
}  






