#!/usr/bin/perl -w

# script to make a tabular overview of the variant positions and coverages ofthese positions by DNA and RNA reads
# takes 3 inputs: -a, a table with the gene positions, starts and stops (like the .tab output of prodigal);
#                 -v, a .vcf file of the same contigs
#		  -o, a name for the output file (optional, is the same as the vcf followed by .type if not given)
# 1 output - a table with the positions of the variants
#
#written by Nicolás Pinel (October 12, 2013); changed by Anna Heintz-Buschart (February 6, 2015)


use strict;
use Getopt::Long;

my ($annotation,$variants,$output);
GetOptions('a=s' => \$annotation,
	   'v=s' => \$variants,
	   'o=s' => \$output);
$output = $output ? $output : "$variants.type";


my $ann = &loadAnnotation($annotation);

print "annotation loaded\n";

my %features = ();
open my $out, '>', $output;
print $out qq(@{[join("\t", "contig","pos","type","diff","gene","frameshift","totalCov","totalVarCov","DNACov","DNAvarCov","RNACov","RNAvarCov")]}\n);
open my $vcf, '<', $variants;
while (<$vcf>) {
    next if (/^\#/);
    my @f = split(/\t/);
    my $ctg = $f[0];
    my $pos = $f[1];
    my $ref = $f[3];
    my $reflength = length($ref);
    my @alts = split ",",$f[4];
    my @DNAdetail = split(/:/,$f[9]);
    my @RNAdetail = split(/:/,$f[10]);
    my @DNAcov = split ",",$DNAdetail[4];
    my @RNAcov = split ",",$RNAdetail[4];
    my @DNAvarCov = split ",",$DNAdetail[5];
    chomp $RNAdetail[5];
    my @RNAvarCov = split ",",$RNAdetail[5];
    my $i = 0;
    foreach my $alt (@alts){
	my $type = "SNP";
	my $diff = 0;
	my $altlength = length($alt);
	unless ($reflength == 1 && $altlength == 1) {
	    if ($reflength == $altlength) {
		$type = "MNP";
	    } elsif ($reflength > $altlength) {
		$type = "del";
		$diff = $altlength - $reflength;
	    } else {
		$type = "insert";
		$diff = $altlength - $reflength;
	    }
	}
	my $DC = $DNAcov[0] == 0 ? 0 : $DNAcov[$i];
	my $RC = $RNAcov[0] == 0 ? 0 : $RNAcov[$i];
	my $DVC = $DNAvarCov[0] == 0 ? 0 : $DNAvarCov[$i];
	my $RVC = $RNAvarCov[0] == 0 ? 0 : $RNAvarCov[$i];
	my $cov = $DC + $RC;
	my $varcov = $DVC + $RVC;
	my $gene = "intergenic";
	my $frameshift = 0;
	my $cdiff = $diff;
	my $firstgene = 0;
	foreach my $nuc ($pos..($pos+$reflength-1)){
	    foreach my $feature (@{$$ann{$f[0]}}){
		if ($nuc >= $$feature{'start'} && $nuc <= $$feature{'stop'}) {
		    if ($gene eq "intergenic") {
			$gene = $$feature{'gene'};
			if ($cdiff > 0 && $firstgene != 0) {
			    $cdiff--;
			} elsif ($cdiff < 0&& $firstgene != 0) {
			    $cdiff++;
			}
			$firstgene = 1;
		    } else {
			$gene = "$gene\;$$feature{'gene'}";
		    }
		}
	    }
	}
	my @allgenes = split("\;",$gene);
	@allgenes = &uniq(@allgenes); 
	$gene = join("\;",@allgenes);
	if ($gene ne "intergenic" && $cdiff % 3 != 0) {
	    $frameshift = $cdiff % 3;
	}
	print $out qq(@{[join("\t", $ctg,$pos,$type,$diff,$gene,$frameshift,$cov,$varcov,$DC,$DVC,$RC,$RVC,"\n")]});
	$i++;
    }
}


exit;

#### subroutines ###

sub usage {

exit;
}

sub uniq {
    my %seen;
    return grep { !$seen{$_}++ } @_;
}

sub loadAnnotation {
    my $f = shift;
    open my $in, '<', $f;
    my @headers = ('gene','start','stop');
    my %feats = ();
    while (<$in>) {
	next if (/^\#/);
	chomp;
	my %hash = ();
	my @f = split(/\t/);
	my $ctg = $f[0];
	$ctg =~ s/_\d+$//g;
	$ctg =~ s/gene//;
	@hash{@headers} = ($f[0],$f[2],$f[3]);
	push(@{$feats{$ctg}}, \%hash);
    }
    return \%feats;
}

