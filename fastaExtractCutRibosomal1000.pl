#!/usr/bin/perl


# script to cut ribosomal 16S RNA sequences as annotated by barrnap out of contigs (as long as the contigs are still at least 1000nt)
# takes 6 inputs: - the fasta file with the contigs;
#		  - the output from barrnap for eukaryotes, archaea, bacteria, mitochondria (in this order)
#		  - the names of the output files (first the fasta, then the table)
# writes 2 outputs: - a fasta file with the all contigs (whole contigs < 1000nt, whole contigs without rRNA, cut contigs with rRNA removed)
#		    - a table with the length of the removed sequences

# written by Anna Heintz-Buschart (February 2015)

use strict;
use warnings;
use Bio::DB::Fasta;

my $fastaFile = shift;
my $eukFile = shift;
my $arcFile = shift;
my $bacFile = shift;
my $mitoFile = shift;
my $outFile = shift;
my $cutFile = shift;

my %rRNAs = ();

open (IN, $eukFile);
while (my $line = <IN>){
    chomp $line;
    unless ($line =~ /Name=5S/) {
        my($contig, undef, undef, $start, $end, undef, undef, undef, undef) = split "\t", $line;
        push @{$rRNAs{$contig}[0]}, $start;
        push @{$rRNAs{$contig}[1]}, $end;
    }
}
close(IN);

open (IN, $arcFile);
while (my $line = <IN>){
    chomp $line;
    unless ($line =~ /Name=5S/) {
        my($contig, undef, undef, $start, $end, undef, undef, undef, undef) = split "\t", $line;
        push @{$rRNAs{$contig}[0]}, $start;
        push @{$rRNAs{$contig}[1]}, $end;
    }
}
close(IN);

open (IN, $bacFile);
while (my $line = <IN>){
    chomp $line;
    unless ($line =~ /Name=5S/) {
        my($contig, undef, undef, $start, $end, undef, undef, undef, undef) = split "\t", $line;
        push @{$rRNAs{$contig}[0]}, $start;
        push @{$rRNAs{$contig}[1]}, $end;
    }
}
close(IN);

open (IN, $mitoFile);
while (my $line = <IN>){
    chomp $line;
    unless ($line =~ /Name=5S/) {
        my($contig, undef, undef, $start, $end, undef, undef, undef, undef) = split "\t", $line;
        push @{$rRNAs{$contig}[0]}, $start;
        push @{$rRNAs{$contig}[1]}, $end;
    }
}
close(IN);


open(OUT, '>', $outFile) or die "could not open file '$outFile' \n";
open(CUT, '>', $cutFile) or die "could not open file '$cutFile' \n";
my $db = Bio::DB::Fasta->new( $fastaFile );
my @ids = $db->get_all_ids;
foreach my $contig (@ids){
    my $sequence = $db->seq($contig);
    if  (!defined( $sequence )) {
            print STDERR "Sequence $contig not found. \n";
            next;
    }
    my $lenSeq = length $sequence;
    unless ($lenSeq < 1000) {
	if (exists $rRNAs{$contig}) {
	    my @useStarts = ();
	    my @useEnds = ();
	    if ($#{$rRNAs{$contig}[0]} > 0) {
		my @starts = @{$rRNAs{$contig}[0]};
		my @ends = @{$rRNAs{$contig}[1]};
		my @idx = sort { $starts[$a]+0 <=> 0+$starts[$b] } 0 .. $#starts;
		@starts = @starts[@idx];
		@ends = @ends[@idx];
		my $currentStart = $starts[0];
		my $currentEnd = $ends[0];
		for (my $i = 0; $i <= $#starts; $i++) { 
		    if ($starts[$i] <= $currentEnd){
			$currentEnd = $ends[$i] if ($ends[$i] >= $currentEnd)
		    } else {
			unshift @useStarts, $currentStart;
			unshift @useEnds, $currentEnd;
			$currentStart = $starts[$i];
			$currentEnd = $ends[$i];
		    }
		}
		unshift @useStarts, $currentStart;
		unshift @useEnds, $currentEnd;
	    } else {
		push @useStarts, $rRNAs{$contig}[0][0];
		push @useEnds, $rRNAs{$contig}[1][0];
	    }
	    my $lenCutSeq = length $sequence;
	    for (my $j = 0 ; $j <= $#useStarts ; $j++ ) {
		my $startcut = $useStarts[$j]; 
		my $length = $useEnds[$j] +1 - $startcut;
		$lenCutSeq = (length $sequence) - $length;
		if ($lenCutSeq >= 1000) {
		    substr($sequence,$startcut,$length) = "";
		} else {
		    $length = (length $sequence) - 1000;
		    substr($sequence,$startcut,$length) = "";
		    $lenCutSeq = length $sequence;
		    last;
		}
	    }
	    print CUT "$contig\t", "$lenSeq\t", "$lenCutSeq\n";
	}
    print OUT ">$contig\n", "$sequence\n";
    }
}
close(OUT);
close(CUT);


sub get_all_ids {

 grep {!/^__/} keys %{shift->{offsets}}

}

