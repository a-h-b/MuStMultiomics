#!/usr/bin/perl

# script to remove non-tryptic peptides from incomplete gene predictions
# takes 2 inputs: - the fasta file with the genes;
#                 - a tab-separated table listing the presence of the start and stop codon in the prediction in the 6th and 7th column (like in the prodigal or "variant_annotateRepairedTab.pl" .tab outputs)
# writes new fasta file to standard out
# rule for tryptic peptides: "cut behind K or R unless followed by P"

# written by Anna Heintz-Buschart (January 2015)


use Bio::DB::Fasta;

my $fastaFile = shift;
my $tabFile = shift;

my $db = Bio::DB::Fasta->new( $fastaFile );


open (IN, $tabFile);
while (<IN>) {
    next if (/^\#/);
    chomp;
    my @f = split(/\t/);
    my $seq = $f[0];
    my $sequence = $db->seq($seq);
        if  (!defined( $sequence )) {
            print STDERR "Sequence $seq not found. \n";
            next;
        }
    if ($f[5] ne "yes" && $f[6] eq "yes") {
        my $ncut = 0;
        while ($sequence ne "" && $ncut == 0) {
            my $k = index($sequence, 'K') >= 0? index($sequence, 'K') : length($sequence);
            my $r = index($sequence, 'R') >= 0? index($sequence, 'R') : length($sequence);
            my $cut = $k < $r ? $k : $r;
            $sequence = substr($sequence, $cut+1);
            if ($sequence eq "" || substr($sequence,0,1) ne "P") {
                $ncut = 1;
            } else {
                $sequence = substr($sequence,1);              
            }
        }
    } elsif ($f[6] ne "yes" && $f[5] eq "yes") {
        my $sub = '';
        my $ccut = 0; 
        while ($sequence ne "" && $ccut == 0) {
            my $k = rindex($sequence, 'K');
            my $r = rindex($sequence, 'R');
            my $cut = $k > $r ? $k : $r;
            $cut = $cut > 0? $cut : -1; 
            my $sub = substr($sequence, $cut+1);
            $sequence = substr($sequence,0,$cut+1);
            if ($sequence eq "" || substr($sub,0,1) ne "P") {
                $ccut = 1;
            } else {
                $sequence = substr($sequence,0,$cut);
            }
        }
    } elsif ($f[5] ne "yes" && $f[6] ne "yes") {
        my $sub = '';
        my $ccut = 0; 
        while ($sequence ne "" && $ccut == 0) {
            my $k = rindex($sequence, 'K');
            my $r = rindex($sequence, 'R');
            my $cut = $k > $r ? $k : $r;
            $cut = $cut > 0? $cut : -1; 
            my $sub = substr($sequence, $cut+1);
            $sequence = substr($sequence,0,$cut+1);
            if ($sequence eq "" || substr($sub,0,1) ne "P") {
                $ccut = 1;
            } else {
                $sequence = substr($sequence,0,$cut);
            }
        }
        my $ncut = 0;
        while ($sequence ne "" && $ncut == 0) {
            my $k = index($sequence, 'K') >= 0? index($sequence, 'K') : length($sequence);
            my $r = index($sequence, 'R') >= 0? index($sequence, 'R') : length($sequence);
            my $cut = $k < $r ? $k : $r;
            $sequence = substr($sequence, $cut+1);
            if ($sequence eq "" || substr($sequence,0,1) ne "P") {
                $ncut = 1;
            } else {
                $sequence = substr($sequence,1);              
            }
        }
    }
    print ">$seq\n", "$sequence\n" unless $sequence eq  "";
}
close(IN)
