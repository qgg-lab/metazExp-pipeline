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
	'imageFile', 'i', 1, 'character', 'Image file name',
	'mageType', 't', 1, 'character', 'Image type include png, pdf',
	'help', 'h', 0, 'logical', 'help document'
), byrow=T, ncol=5)

args <- getopt(command)
if(is.null(args$aSetFile) || is.null(args$bSetFile) || is.null(args$cSetFile)){
	cat(paste(getopt(command, usage = T), "\n"))
	q(status=1)
}

A<-read.table(args$aSetFile, sep="\t", header=FALSE)
B<-read.table(args$bSetFile, sep="\t", header=FALSE)
C<-read.table(args$cSetFile, sep="\t", header=FALSE)

venn.diagram(x=list(first1=A[,1], second2=B[,1], third3=C[,1]), alpha=0.5, cex=1.2, category.names = c(args$aSetLabel, args$bSetLabel, args$cSetLabel), fill = c("tomato2", "skyblue3","yellow2"), euler.d = TRUE, cat.pos = c(-20,20,180), cat.cex = 1.5, reverse = TRUE, col = "transparent", cat.col = c("tomato2", "skyblue3", "yellow2"), fontfamily = "serif", margin = 0.05, filename = args$imageFile, imagetype = args$mageType)

