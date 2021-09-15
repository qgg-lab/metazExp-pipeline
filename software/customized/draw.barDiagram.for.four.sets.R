# _*_ coding: utf-8 _*_
library(ggplot2)
library(getopt)
command=matrix(c(
	'dataFile', 'd', 1, 'character', 'Data file of 4 classes of numbers(number and class label)',
	'xAxisLabel','X', 1, 'character', 'The label of X-axis, e.g Annotation',
	'yAxisLabel','Y', 1, 'character', 'The label of Y-axis, e.g Permill',
	'titleLegend','L', 1, 'character', 'The title of legend, e.g exon num',
	'imageFile', 'i', 1, 'character', 'Image file name',
	'mageType', 't', 1, 'character', 'Image type include png',
	'help', 'h', 0, 'logical', 'help document'
), byrow=T, ncol=5)

args <- getopt(command)
if(is.null(args$dataFile) || is.null(args$xAxisLabel) || is.null(args$xAxisLabel) || is.null(args$titleLegend)){
	cat(paste(getopt(command, usage = T), "\n"))
	q(status=1)
}

dataSrc<-read.table(args$dataFile, sep="\t", header=FALSE)

png(file=args$imageFile)

g <- ggplot(dataSrc, aes(V2))

g + geom_bar(aes(fill = V1)) + xlab(args$xAxisLabel) + ylab(args$yAxisLabel) + theme(legend.position="bottom") + guides(fill=guide_legend(title=args$titleLegend))

dev.off()
