# _*_ coding: utf-8 _*_
library(grid)
library(futile.logger)
library(VennDiagram)
library(getopt)
command=matrix(c(
	'aSetFile', 'A', 1, 'character', 'Data file of 1st set',
	'aSetLabel', 'a', 1, 'character', 'Label of 1st set',
	'bSetFile', 'B', 1, 'character', 'Data file of 2nd set',
	'bSetLabel', 'b', 1, 'character', 'Label of 2nd set',
	'cSetFile', 'C', 1, 'character', 'Data file of 3rd set',
	'cSetLabel', 'c', 1, 'character', 'Label of 3rd set',
	'dSetFile', 'D', 1, 'character', 'Data file of 4th set',
	'dSetLabel', 'd', 1, 'character', 'Label of 4th set',
	'eSetFile', 'E', 1, 'character', 'Data file of 5th set',
	'eSetLabel', 'e', 1, 'character', 'Label of 5th set',
	'imageFile', 'i', 1, 'character', 'Image file name',
	'mageType', 't', 1, 'character', 'Image type include png, pdf',
	'help', 'h', 0, 'logical', 'help document'
), byrow=T, ncol=5)

args <- getopt(command)
if(is.null(args$aSetFile) || is.null(args$bSetFile) || is.null(args$cSetFile) || is.null(args$dSetFile) || is.null(args$imageFile) || is.null(args$eSetFile) || is.null(args$mageType)){
	cat(paste(getopt(command, usage = T), "\n"))
	q(status=1)
}

A<-read.table(args$aSetFile, sep="\t", header=FALSE)
B<-read.table(args$bSetFile, sep="\t", header=FALSE)
C<-read.table(args$cSetFile, sep="\t", header=FALSE)
D<-read.table(args$dSetFile, sep="\t", header=FALSE)
E<-read.table(args$eSetFile, sep="\t", header=FALSE)

venn.diagram(x=list(first=A[,1], second=B[,1], third=C[,1], fourth=D[,1], fifth=E[,1]), category.names = c(args$aSetLabel, args$bSetLabel, args$cSetLabel, args$dSetLabel, args$eSetLabel), col="transparent",fill=c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),alpha = 0.50,cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,1, 0.8, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),cat.col=c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),cat.cex=1,cat.fontface="bold",margin = 0.08, filename = args$imageFile, imagetype = args$mageType)

