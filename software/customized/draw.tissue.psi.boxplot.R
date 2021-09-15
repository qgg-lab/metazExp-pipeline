# _*_ coding: utf-8 _*_
Sys.setenv("DISPLAY"=":0")
library(getopt)
library(ggplot2)
library(Cairo)
library(ggsignif)
library(grid)

command=matrix(c(
        'psiDataFile', 'p', 1, 'character', 'e.g Tissue\tPsi',
        'tissueNum','n', 1, 'character', 'tissue num',
	'experimentNum','e', 1, 'character', 'at least experiment num',
	'imageFile','i', 1, 'character', 'Image file name',
        'format','f', 1, 'character', 'Image format, png or pdf',
        'help',     'H', 0, 'logical', 'help document'
), byrow=T, ncol=5)

args <- getopt(command)

if(is.null(args$psiDataFile) || is.null(args$imageFile) || is.null(args$format) || is.null(args$tissueNum)){
        cat(paste(getopt(command, usage = T), "\n"))
        q(status=1)
}

if(args$format=="png"){
	CairoPNG(file=paste(args$imageFile, ".", args$format, sep=""), width=800,height=600)
}else{
	CairoPDF(file=paste(args$imageFile, ".", args$format, sep=""))
}

psi <- read.table(args$psiDataFile, sep='\t', header=TRUE)


line <- "#1F3552"

p <- ggplot(psi, aes(Tissue, Psi)) + geom_boxplot(colour = line, fill="red", alpha = 0.8)

p + theme(panel.background = element_blank(), axis.line.x = element_line(size = 0.5, colour = "black"), axis.line.y = element_line(size = 0.5, colour = "black"), axis.title = element_text(size=18), axis.text.y = element_text(size=14), axis.text.x = element_text(size=14, angle=65, hjust=1)) + scale_y_continuous(name="Percentage Spliced-In(PSI)", limits=c(0, 1)) + scale_x_discrete(name = "") + ggtitle(paste("Boxplot of PSI in ", args$tissueNum , " tissues (>=", args$experimentNum , "exp.)", sep="")) + theme(plot.title = element_text(size=24, hjust = 0.5), legend.position = "none")


dev.off()

