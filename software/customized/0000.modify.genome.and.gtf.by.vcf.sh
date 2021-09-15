#!/bin/bash
#SBATCH --job-name=replace
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=110G
#SBATCH --time=24:00:00

export currDir=`pwd`
export scriptDir=$currDir/script
export genomefa=/mnt/research/qgg/share/jianping/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa
export annotationGtf=/mnt/research/qgg/share/jianping/Sus_scrofa.Sscrofa11.1.98.gtf
export replaceBed=/mnt/research/qgg/share/jianping/test.bed




grep -vP "\tstart_codon|stop_codon\t" $annotationGtf |grep -vP "\tSelenocysteine\t" |grep -vP "#" > $currDir/orig.gtf

#:<<blockGenome
perl $scriptDir/003.singleLine.genomeSeq.with.onlyId.pl \
	$genomefa \
	$currDir/ori.genome.fa
#blockGenome

#:<<blockGene
grep -P "\tgene\t" $currDir/orig.gtf > $currDir/gene.gtf
perl $scriptDir/001.replace.pl \
        $currDir/ori.genome.fa \
        $currDir/replaced.genome.fa \
        $currDir/gene.gtf \
        $currDir/modified.gene.gtf \
        $replaceBed

#blockGene

#:<<block5utr
grep -P "\tfive_prime_utr\t" $currDir/orig.gtf > $currDir/5utr.gtf 
perl $scriptDir/001.replace.pl \
        $currDir/ori.genome.fa \
        $currDir/replaced.genome.fa \
        $currDir/5utr.gtf \
        $currDir/modified.5utr.gtf \
	$replaceBed
#block5utr

#:<<blockTrspt
grep -P "\ttranscript\t" $currDir/orig.gtf > $currDir/transcript.gtf
perl $scriptDir/001.replace.pl \
        $currDir/ori.genome.fa \
        $currDir/replaced.genome.fa \
        $currDir/transcript.gtf \
        $currDir/modified.transcript.gtf \
	$replaceBed
#blockTrspt

#:<<blockExon
grep -P "\texon\t" $currDir/orig.gtf > $currDir/exon.gtf
perl $scriptDir/001.replace.pl \
        $currDir/test.ori.genome.fa \
        $currDir/replaced.genome.fa \
        $currDir/test.exon.gtf \
        $currDir/modified.test.exon.gtf \
	$replaceBed
#blockExon

#:<<blockStartCodon
grep -P "\tstart_codon\t" $currDir/orig.gtf > $currDir/start_codon.gtf 
perl $scriptDir/001.replace.pl \
        $currDir/ori.genome.fa \
        $currDir/replaced.genome.fa \
        $currDir/start_codon.gtf \
        $currDir/modified.start_codon.gtf \
	$replaceBed
#blockBlockCDS


#:<<blockCDS
grep -P "\tCDS\t" $currDir/orig.gtf > $currDir/CDS.gtf 
perl $scriptDir/001.replace.pl \
        $currDir/ori.genome.fa \
        $currDir/replaced.genome.fa \
        $currDir/CDS.gtf \
        $currDir/modified.CDS.gtf \
	$replaceBed
#blockBlockCDS

#:<<blockStopCodon
grep -P "\tstop_codon\t" $currDir/orig.gtf > $currDir/stop_codon.gtf 
perl $scriptDir/001.replace.pl \
        $currDir/ori.genome.fa \
        $currDir/replaced.genome.fa \
        $currDir/stop_codon.gtf \
        $currDir/modified.stop_codon.gtf \
	$replaceBed
#block3utr

#:<<block3utr
grep -P "\tthree_prime_utr\t" $currDir/orig.gtf > $currDir/3utr.gtf 
perl $scriptDir/001.replace.pl \
        $currDir/ori.genome.fa \
        $currDir/replaced.genome.fa \
        $currDir/3utr.gtf \
        $currDir/modified.3utr.gtf \
	$replaceBed
#block3utr

# 将modified.gene.gtf, modified.transcript.gtf, modified.exon.gtf, modified.CDS.gtf重新组装到一起
#:<<blockCombine
perl $scriptDir/002.combine.gene.transcript.exon.CDS.together.pl \
	$currDir/modified.gene.gtf,$currDir/modified.transcript.gtf,$currDir/modified.exon.gtf,$currDir/modified.5utr.gtf,$currDir/modified.start_codon.gtf,$currDir/modified.CDS.gtf,$currDir/modified.stop_codon.gtf,$currDir/modified.3utr.gtf \
	gene,transcript,exon,5utr,start_codon,CDS,stop_codon,3utr \
	$currDir/modified.annotation.gtf
#blockCombine
