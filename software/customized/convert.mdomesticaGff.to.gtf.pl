#!/usr/bin/perl
my $gff=$ARGV[0];
my $newGff = $ARGV[1];
# 将苹果的gff3注释转换为gtf注释
# GFF3:
# Chr00   phytozomev12    gene    25312   25566   .       -       .       ID=MD00G1000100.v1.1;Name=MD00G1000100
# Chr00   phytozomev12    mRNA    25312   25566   .       -       .       ID=MD00G1000100.v1.1.491;Name=MD00G1000100;pacid=40123409;longest=1;Parent=MD00G1000100.v1.1
# Chr00   phytozomev12    exon    25312   25566   .       -       .       ID=MD00G1000100.v1.1.491.exon.1;Parent=MD00G1000100.v1.1.491;pacid=40123409
# Chr00   phytozomev12    CDS     25312   25566   .       -       0       ID=MD00G1000100.v1.1.491.CDS.1;Parent=MD00G1000100.v1.1.491;pacid=40123409
# Chr00   phytozomev12    five_prime_UTR  38005   38090   .       +       .       ID=MD00G1000200.v1.1.491.five_prime_UTR.1;Parent=MD00G1000200.v1.1.491;pacid=40123412
 
# GTF:
# （1）mRNA改成transcript
# （2）Id=* => gene_id "xxxx"; transcript_id "yyyy"
#  (3) geneId 从MD00G1000100.v1.1 => MD00G1000100
#  (4) transcriptId 从MD00G1000100.v1.1.491 => MD00G1000100.491
my ($line, @field, $geneId, $trsptId);
open FF, "<$gff";
open WW, ">$newGff";
while($line=<FF>){
	chomp($line);
	@field = split(/\t/, $line);
	if($field[2] eq "gene"){
		if($field[8]=~/ID=(.*?)\.v1\.1;Name=(.*)/){
			$geneId = $1 . "g";
			pop(@field);
			print WW join("\t", @field, "gene_id \"" . $geneId . "\";") . "\n";
		}
	}elsif($field[2] eq "mRNA"){
		if($field[8]=~/ID=(.*?)\.v1\.1\.491;Name.*/){
			pop(@field);
			$trsptId = $1;
			$field[2] = "transcript";
			print WW join("\t", @field, "gene_id \"" . $geneId . "\"; transcript_id \"" . $trsptId . "\";") . "\n";
		}
	}elsif($field[2] eq "exon"){
		if($field[8]=~/ID=(.*?)\.v1\.1\.491\.exon\.(\d+);Parent=.*/){
			pop(@field);
			$exonNum = $2;
			print WW join("\t", @field, "gene_id \"" . $geneId . "\"; transcript_id \"" . $trsptId . "\"; exon_number \"" . $exonNum . "\";") . "\n";
		}
	}elsif($field[2] eq "CDS"){
		if($field[8]=~/ID=(.*?)\.v1\.1\.491\.CDS\.(\d+);Parent=.*/){
			pop(@field);
			print WW join("\t", @field, "gene_id \"" . $geneId . "\"; transcript_id \"" . $trsptId . "\"; exon_number \"" . $exonNum . "\";") . "\n";
		}
	}elsif($field[2] eq "five_prime_UTR"){
		# ID=MD00G1000200.v1.1.491.five_prime_UTR.1;Parent=MD00G1000200.v1.1.491;pacid=40123412
		if($field[8]=~/ID=(.*?)\.v1\.1\.491\.five_prime_UTR\.(\d+);Parent=.*/){
			pop(@field);
			$exonNum = $2;
			print WW join("\t", @field, "gene_id \"" . $geneId . "\"; transcript_id \"" . $trsptId . "\";") . "\n";
		}
	}elsif($field[2] eq "three_prime_UTR"){
		# ID=MD00G1000200.v1.1.491.three_prime_UTR.1;Parent=MD00G1000200.v1.1.491;pacid=40123412
		if($field[8]=~/ID=(.*?)\.v1\.1\.491\.three_prime_UTR\.(\d+);Parent=.*/){
			pop(@field);
			$exonNum = $2;
			print WW join("\t", @field, "gene_id \"" . $geneId . "\"; transcript_id \"" . $trsptId . "\";") . "\n";
		}
	}

}
close FF;
close WW;
