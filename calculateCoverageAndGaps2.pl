#!/usr/bin/perl


# script to calculate average coverage of genes or contigs from file with depth of coverage at every position (except at the positions with 0 depth of coverage; output of samtools depth);
# it replaces Mads Albertsen's calculateCoverageAndGaps.pl, which miscalculated the covered length
# takes 2 inputs: - the fasta file with genes or contigs
#		  - the depth file
# 1 output is written to standard out: a tab-separated table with the gene or contig IDs, lengths, average coverage and covered length

# written by Anna Heintz-Buschart (April 2015)


use strict;

my $contigs=$ARGV[0];
my $indepth=$ARGV[1];

my %contigs=();

my ($id);

#
my %coverage;
#my @order;
my %length;

open(CON,$contigs) or die $!;
while(my $str=<CON>){
    chomp($str);
    next if length($str)==0;
    #print STDERR "str: $str\n";
    if($str=~/>(\S+)/){
        $id=$1;
        #print STDERR "$id\n";
        die "$id already exists\n" if exists($contigs{$id});
    }else{
        #print STDERR "$id\n";
        $contigs{$id}{'seq'}.=$str;
    }
}
close(CON);

open(INdepth, "$indepth") or die("Cannot read file: $indepth\n");

while ( my $line = <INdepth> ) {
        chomp $line;
        my @splitline = split(/\t/,$line);
        if (exists($coverage{$splitline[0]})){
                $coverage{$splitline[0]} = $coverage{$splitline[0]} + $splitline[2];
                $length{$splitline[0]}++;
        }
        else{
                $coverage{$splitline[0]} = $splitline[2];
                $length{$splitline[0]} = 1;
 #               push (@order , $splitline[0]);
        }
}

close INdepth;

#open(OUT, ">$outputfile") or die("Cannot create file: $outputfile\n");

print "SequenceID\tReference.length\tAverage.coverage\tCovered.length\n";
#

#print "sequenceID\tlength\tGCperc\n"; 


foreach $id (keys(%contigs)){
    next if length($id)==0;
    my $seq=$contigs{$id}{'seq'};
    my $reflength=length($seq); 
    if ($reflength eq 0){
	print STDERR "ERROR: $id $reflength $seq\n";
	next;
    }
    my $cov = 0;
    my $covlength = 0;
    if (exists $coverage{$id}){
	$cov = sprintf("%.3f", $coverage{$id} / $reflength);
	$covlength = $length{$id}
    } 
    $id=~s/\s+/_/g;
    print "$id\t$reflength\t$cov\t$covlength\n"; 
}


