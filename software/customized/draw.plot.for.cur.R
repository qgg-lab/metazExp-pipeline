# _*_ coding: utf-8 _*_
library(getopt)
command=matrix(c(
	'dataFile', 'd', 1, 'character', 'e.g ExperimentNum, Gene, transcript, exon num',
	'title',    't', 1, 'character', 'Title of diagram',
	'xLabel',   'x', 1, 'character', 'Label for x-axis',
	'yLabel',   'y', 1, 'character', 'Label for y-axis',
	'width',    'w', 1, 'integer', 'width of image(px)',
	'height',   'h', 1, 'integer', 'height of image(px)',
	'imageFile','i', 1, 'character', 'Image file name',
	'help',     'H', 0, 'logical', 'help document'
), byrow=T, ncol=5)

args <- getopt(command)

colList <- c("red", "blue", "purple", "yellow", "orange", "black", "cyan", "green")

if(is.null(args$dataFile)||is.null(args$imageFile)||is.null(args$xLabel)||is.null(args$yLabel)||is.null(args$title)||is.null(args$width)||is.null(args$height)){
	cat(paste(getopt(command, usage = T), "\n"))
	q(status=1)
}

#create png file
png(args$imageFile, width=args$width, height=args$height, units="px")

#
# par(oma=c(1,1,1,1), mar=c(2,2,2,2))
#
# figure include plot. xLabel and yLabel of axis out of plot and in figure.
# oma is the region out of figure
# mar is the region in figure and out of plot
# box is used to point which is plot, figure and outter
# box(which = "plot",  col = "red",    lwd = 2)
# box(which = "figure",col = "blue",   lwd = 4)
# box(which = "outer", col = "black",  lty = 8)
#
#

#load data into data frame
dataSrc <- read.table(args$dataFile, sep=",", header=TRUE)

#obtain colummn name from dataSrc
cName <- names(dataSrc)

#obtain column name
cNum <- length(cName)

#calculate max Y value
maxY <- max(dataSrc)

#pch
pchValue <- 15

#plot first

x <- data.frame(dataSrc[1], dataSrc[2])

plot(x, type = "o", xlab = args$xLabel, ylab = args$yLabel, ylim=c(0,maxY), col=colList[1], main = args$title, pch=c(15))

colNum <- 1
#loop to create curve for(i in c(3:5))
for(i in c(3: cNum)){

	pchValue <- pchValue + 1

	colNum <- colNum + 1

	x <- data.frame(dataSrc[1], dataSrc[i])	

	lines(x, type="o", col = colList[colNum], pch=c(pchValue))

}


legend("topleft", cName[2:cNum], col=colList[1:cNum-1], pch=15:pchValue)

dev.off()
