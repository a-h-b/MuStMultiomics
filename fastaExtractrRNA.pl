#!/usr/bin/perl

# this script takes a fasta file with contigs and the output of one Barrnap run to produce a fasta file with rRNA sequences and a .tab file in the style of prodigal
# 4 inputs: - the fasta file of the contigs
#           - the .gff file from Barrnap
#           - the name of the .tab output
#           - the name of the fasta output
# 2 outputs: - a fasta file with the rRNA gene sequences, naming is the contig name appended with _r and a continuous number per contig
#	     - a table with contig name rRNA gene name (see above), the sense, length, start position, end position, completeness and kind (16S, 5S etc)

# Anna Heintz-Buschart, November 2014


use strict;
use Bio::DB::Fasta;

my $fastaFile = shift;
my $bacFile = shift;
my $listFile = shift;
my $geneFile = shift;

my %allContigs = ();

my $db = Bio::DB::Fasta->new( $fastaFile );
open (IN, $bacFile);
open(GEN, ">", $geneFile) or die "cannot open $geneFile \n";
open(TAB, ">", $listFile) or die "cannot open $listFile \n";
print TAB "contig\tgene\tsense\tlength\tstart\tend\tcompleteness\tkind\n";
while (my $line = <IN>){
    unless ($line =~ /^#/) { 
        chomp $line;
        my($contig, undef, undef, $start, $end, undef, $sense, undef, $attribute) = split "\t", $line;
        my $completeness = "";
	if ($attribute =~ /partial/){
            $completeness = "incomplete"
        } else {
            $completeness = "complete"
            }
        if (exists $allContigs{$contig}) {
            $allContigs{$contig}++
        } else {
            $allContigs{$contig} = 1
        }
        my $gene = join("", $contig, "_r", $allContigs{$contig});
        my $length = 1 + $end - $start;
        my @attributes = split ";", $attribute,2;
        my $kind = $attributes[0];
        substr($kind, 0 , 5) = "";
        print TAB join("\t",$contig,$gene,$sense,$length,$start,$end,$completeness,$kind), "\n";
        my $sequence = "";
        if ($sense == "+") {
            $sequence = $db->seq($contig, $start, $end );
        } else {
            $sequence = $db->seq($contig, $end, $start);
        }
        if  (!defined( $sequence )) {
            print STDERR "Sequence $contig not found. \n";
            next;
        } else {
            print GEN ">$gene\n", "$sequence\n";
        }
    }
}
close(IN);
close(GEN);
close(TAB);
