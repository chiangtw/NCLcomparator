#! /usr/bin/env Rscript

# Usage:
# ./graphic.R Merged.result OUTPUTname

args <- commandArgs(TRUE)
in_file <- args[1]
out_file <- args[2]

intra <- read.table(in_file, sep="\t", header=TRUE)
outputName <- paste(out_file, "pdf", sep=".")

########################################################################

if (!require("ggplot2")) {
  install.packages("ggplot2")
  library(ggplot2)
}

if (!require("gridExtra")) {
  install.packages("gridExtra")
  library(gridExtra)
}

intraEventNum <- dim(intra)[1]
intraLastcol <- dim(intra)[2]
intra2 <- data.matrix(intra[,c(8:intraLastcol)])
intraToolsNum <- intraLastcol-7
intraToolsNames <- colnames(intra2)

intraTools.p <- c()
intraTools.len <- c()
for (i in 1:intraToolsNum)
{
  toolone <- intra2[intra2[,i]==1,]
  
  toolone.len <- dim(toolone)[1]
  intraTools.len <- c(intraTools.len, toolone.len)
  
  p.one <- c()
  for (j in 1:intraToolsNum)
  { 
    p.1 <- sum(rowSums(toolone) ==j)/toolone.len
    p.one <- c(p.one, p.1)
  }
  intraTools.p <- rbind(intraTools.p, p.one)
}


intraToolNameV <- as.factor(rep(intraToolsNames,each=intraToolsNum))
intraCountTool <- as.numeric(rep(c(1:intraToolsNum), intraToolsNum))
intraProportion <- as.numeric((c(t(intraTools.p))))
intra3 <- data.frame(intraToolNameV, intraCountTool, intraProportion)
intra4 <- data.frame(intraToolsNames, intraTools.len)

intra4$intraToolsNames <- factor((intra4$intraToolsNames), levels=as.character(intra4$intraToolsNames))
intra3$intraToolNameV <- factor(intra3$intraToolNameV, levels=as.character(intra4$intraToolsNames))

pdf(outputName, width=12, height=5) 
 
intraPlot1 <- ggplot(intra4, aes(x=intraToolsNames, y=intraTools.len, fill=intraToolsNames)) + geom_bar(stat="identity") + theme(legend.position="none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + labs(y ="Number of identified events", title="Reported number by each tools  ")
intraPlot2 <- ggplot(intra3, aes(x=intraCountTool, y=intraProportion, colour=intraToolNameV)) + geom_line() + scale_x_discrete(limit=c(1:intraToolsNum)) + labs(x="Number of supported tools", y="Proportion", colour="Tools", title="Supported by tools")
grid.arrange(intraPlot1, intraPlot2, nrow=1, ncol=2)

dev.off()









