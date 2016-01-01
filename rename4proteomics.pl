#!/usr/bin/perl

# script to rename all genes in a fasta file to >sp|xxxxxxx| genename
# takes 1 input - the fasta file with the genes;
# writes new fasta file to standard out

# written by Anna Heintz-Buschart (February 2015)


use strict;
use warnings;
use Bio::DB::Fasta;

my $fastaFile = shift;
my $batch = shift;

my $cnt = 1;
my $db = Bio::DB::Fasta->new( $fastaFile );
my @ids = $db->get_all_ids;
foreach my $gene (@ids){
    my $sequence = $db->seq($gene);
    if  (!defined( $sequence )) {
            print STDERR "Sequence $gene not found. \n";
            next;
    }
    my $geneNo = sprintf("%07d",$cnt);
    $cnt++;
    $gene = "sp|$batch$geneNo| $gene";
    print ">$gene\n", "$sequence\n";
}


sub get_all_ids {

 grep {!/^__/} keys %{shift->{offsets}}

}

