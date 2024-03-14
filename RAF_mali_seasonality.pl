#!/usr/bin/perl -w                                                                                                                                             \
                                                                                                                                                                            
use strict;
use warnings;

open (SAMPLES, "/local/projects-t3/SerreDLab-3/kieran.tebben/mali/sample_lists/redosample.txt");
my @sample_names;

while(<SAMPLES>) {
	chomp;
	push @sample_names, $_;
}

my $mpileup;
my @read;
my $position;
my $allele;
my $reference;
my $RAF; 

for my $sample (@sample_names) {
	print "$sample\n";
	
	my $out = $sample."_RAF.txt";
	
	open(OUT, ">", "/local/projects-t3/SerreDLab-3/kieran.tebben/mali/hisat2_outputs/pf_rmdup_bams/mpileup/subset_mpileup/$out");
	print OUT "chromosome\tposition\tcoverage\treads_with_ref_allele\tRAF\n";
	
	chomp $sample;
	my $mpileup_path = "/local/projects-t3/SerreDLab-3/kieran.tebben/mali/hisat2_outputs/pf_rmdup_bams/mpileup/subset_mpileup/${sample}.pf.subset.rmdup_subset.mpileup";
	open(IN, $mpileup_path);
	while(my $line = <IN>){
		chomp $line;
#		print $line;
		$reference = 0;
		@read = split(/\t+/, $line);  ## splitting line into array
#		print "$read[3]\n";
#		sleep(2)
		if($read[3] >= 50){
#			print "$read[0]\n";
			for($position = 0; $position < length $read[4]; $position++){
				$allele = substr($read[4], $position, 1);
#				print "$allele\n";
				if($allele eq "\." or $allele eq "\," or $allele eq "\>" or $allele eq "\<"){
#					print "$allele\n";
#					print "reference\n";
					$reference++;
#					print "$reference\n";
				}
			}	
#			print "$reference\n";
#			print "$read[3]\n";		
			$RAF = $reference/$read[3];
#			print "$RAF\n";
			print OUT "$read[0]\t";
			print OUT "$read[1]\t";
			print OUT "$read[3]\t";
			print OUT "$reference\t";
			print OUT "$RAF\n";
		}
	}
	close(OUT);
}

exit;