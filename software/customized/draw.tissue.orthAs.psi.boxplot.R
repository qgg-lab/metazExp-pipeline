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
        'speciesNum','e', 1, 'character', 'species num',
	'imageFile','i', 1, 'character', 'Image file name',
        'format','f', 1, 'character', 'Image format, png or pdf',
	'nncol','c', 1, 'integer', 'ncol for facet ncol',
	'fontSize', 's', 1, 'integer', 'size for facet ncol',
        'help',     'H', 0, 'logical', 'help document'
), byrow=T, ncol=5)

args <- getopt(command)

if(is.null(args$psiDataFile) || is.null(args$imageFile) || is.null(args$format) || is.null(args$tissueNum) || is.null(args$nncol) || is.null(args$fontSize)|| is.null(args$speciesNum)){
        cat(paste(getopt(command, usage = T), "\n"))
        q(status=1)
}

if(args$format=="png"){
	CairoPNG(file=paste(args$imageFile, ".", args$format, sep=""), width=1000,height=1000)
}else{
	CairoPDF(file=paste(args$imageFile, ".", args$format, sep=""))
}

psi <- read.table(args$psiDataFile, sep='\t', header=TRUE)


line <- "#1F3552"

p <- ggplot(psi, aes(y=Psi, x=Species, fill=Species)) + facet_wrap(~Tissue, ncol = args$nncol) +  geom_boxplot(colour = line, alpha = 0.8)

p + theme(panel.background = element_blank(), axis.line.x = element_line(size = 0.5, colour = "black"), axis.line.y = element_line(size = 0.5, colour = "black"), axis.title = element_text(size=args$fontSize + 6), axis.text.y = element_text(size= args$fontSize + 4 ), axis.text.x = element_text(size= args$fontSize + 6, angle=75, hjust=1)) + scale_y_continuous(name="Percentage Spliced-In(PSI)", limits=c(0, 1)) + scale_x_discrete(name = "") + ggtitle(paste("Boxplot of PSI in ", args$tissueNum , " tissues in ", args$speciesNum, " species", sep="")) + theme(plot.title = element_text(size= args$fontSize + 12, hjust = 0.5), legend.position = "none", strip.text.x = element_text(size = args$fontSize + 6, colour = "black"))


dev.off()
