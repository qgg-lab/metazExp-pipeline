# _*_ coding: utf-8 _*_
library(ggplot2)
library(getopt)
command=matrix(c(
	'dataFile', 'd', 1, 'character', 'Data file of 4 classes of numbers(number and class label)',
	'xLabel', 'X', 1, 'character', 'The description of number(e.g Exon_length(no space))',
	'xlimMin', 'm', 1, 'integer', 'Min value of x-axis',
	'xlimMax', 'x', 1, 'integer', 'Max value of x-axis',
	'imageFile', 'i', 1, 'character', 'Image file name',
	'mageType', 't', 1, 'character', 'Image type include png',
	'help', 'h', 0, 'logical', 'help document'
), byrow=T, ncol=5)

args <- getopt(command)
if(is.null(args$dataFile)||is.null(args$xLabel) || is.null(args$xlimMin) || is.null(args$xlimMax)){
	cat(paste(getopt(command, usage = T), "\n"))
	q(status=1)
}

dataSrc<-read.table(args$dataFile, sep="\t", header=FALSE)

#names(dataSrc)<-c(args$numberDesp, args$classDesp)

png(file=args$imageFile)

ggplot(dataSrc, aes(V1, colour = V2)) + geom_density(alpha = 0.1) + xlim(args$xlimMin, args$xlimMax) + xlab(args$xLabel) + theme(legend.position="right") +  theme(legend.title=element_blank())

dev.off()
