#!/usr/bin/perl

# this script takes a fasta file and a space-delimited file with IDs present in the fasta headers and start- and end- coordinates (with the first letter counted as position 1)
# 2 inputs: - the fasta file (eg. mOTU.v1.padded)
#           - a space-separated table with IDs and coordinates (eg. mOTU.v1.padded.coord)
# 1 output to stdout: a fasta file with sequences indicated by the table

# Anna Heintz-Buschart, February 2015



use strict;
use warnings;
use Bio::DB::Fasta;


my $fastaFile = shift;
my $queryFile = shift;

my $db = Bio::DB::Fasta->new( $fastaFile );

open (IN, $queryFile);
while (<IN>){
    chomp;
    my @fields = split /\t/;
    my $seq = $fields[0];
    my $sequence = $db->seq($seq);
    if  (!defined( $sequence )) {
            print STDERR "Sequence $seq not found. \n";
            next;
    }
    my $start=0;
    my $stop=0;
    if ($fields[1] < $fields[2]) {
        $sequence = substr($sequence,$fields[1]-1,($fields[2]-$fields[1]+1));
        $start = $fields[1];
        $stop = $fields[2];
    } else {
        $sequence = substr($sequence,$fields[2]-1,($fields[1]-$fields[2]+1));
        $start = $fields[2];
        $stop = $fields[1];
        $sequence = join('', reverse(split(//, $sequence)));
	$sequence =~ tr/ACGT/TGCA/;
    }
    print ">$seq\n", "$sequence\n";
}
close(IN);
