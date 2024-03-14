#!/usr/bin/perl -w                                                                                                                                             \
                                                                                                                                                                            
use strict;
use warnings;

my %keep;

open(GENES, "/local/projects-t3/SerreDLab-3/kieran.tebben/mali/genomes/pf_genes.txt");
#print "opened\n";

print "generating hash\n";
while(my $line = <GENES>){
#	print $line;
	my @info = split(/\t+/, $line);
#	print "$info[1]\n";
	my $chr = $info[0];
#	print "$chr\n";
	my $start = $info[1];
	my $end = $info[2];
#	print "$chr\t$start\t$end\n";
	for(my $i = $start; $i <= $end; $i++){
		my $position = $chr."_$i";
		$keep{"$position"} = 1;
#		print "$position\n";
	}
}

close(GENES);

my @mpileup_files = </local/projects-t3/SerreDLab-3/kieran.tebben/mali/sample_lists/*.mpileup>;
my $mpileup;

for $mpileup (@mpileup_files){
#	print "$mpileup\n";
	
	my @files = split /\//, $mpileup;
	my $sample_name = $files[7];	
	(my $sample = $sample_name) =~ s/\.[^.]+$//;
	print "$sample\n";
	
	my $out = $sample."_subset.mpileup";
#	print "$out\n";
	open(OUT, ">", "/local/projects-t3/SerreDLab-3/kieran.tebben/mali/sample_lists/$out");
	
	open(IN, $mpileup);
		
	while(my $pileup_line = <IN>){
#		print "$pileup_line\n";
		my @info = split(/\t+/, $pileup_line);
		my $chr = $info[0];
		my $loc = $info[1];
		my $position = $chr."_$loc";
#		print "$position\n";
		if (defined $keep{$position}){
#			print "$position\n";
#			print "$pileup_line\n";
			print OUT "$pileup_line";
		}
	}
	close(OUT);
}


exit;