#!/usr/bin/perl -w


# written by Venkata Satagopam, with a few corrections and limitation to KEGG by Anna Heintz-Buschart

# this is a version of "consolidata_hmmscan_results.pl" which only uses the KEGG results.
#takes 2 inputs: the name of the output, will be concatenated with the database name and .tsv; and the result of hmmscan for kegg


use strict;

my $result_file=$ARGV[0];
my $kegg_hmmscan_result_file=$ARGV[1];

my %hit_hash=();
my @all = ();
my $line = ();

### KEGG

print STDERR "# KEGG: $kegg_hmmscan_result_file\n";

my %kegg_hmmscan_result_hash = ();
close(FH);
open(FH, "$kegg_hmmscan_result_file") or die "Failed to open $kegg_hmmscan_result_file\n";
@all = <FH>;	
shift @all; # to get rid off header
while (@all) {
	$line = shift @all;
	chomp($line);
	next if $line=~/^#/;
	my ($hit_name, $acc1, $query_name, $acc2,$seqEvalue,$seqscore,$seqbias,$dom1Evalue,$dom1score,$dom1bias,$exp,$reg,$clu,$ov,$env,$dom,$rep,$inc,@desc_of_target) = split(/ +/, $line);
	push @{$kegg_hmmscan_result_hash{$query_name}}, $line;
	$hit_hash{$query_name} = 1;
} 
close(FH);

my $outfile=$result_file."kegg.tsv";
open RESULT,">$outfile";
print RESULT "QUERY_NAME\tHIT_NAME\tSCORE\tSIGNIFICANCE\n";

foreach my $query_name (keys %hit_hash) {
	# kegg
	foreach $line (@{$kegg_hmmscan_result_hash{$query_name}}){
		my ($hit_name, $acc1, $query_name, $acc2,$seqEvalue,$seqscore,$seqbias,$dom1Evalue,$dom1score,$dom1bias,$exp,$reg,$clu,$ov,$env,$dom,$rep,$inc,$desc_of_target) = split(/ +/, $line);
		print RESULT "$query_name\t$hit_name\t$seqscore\t$seqEvalue\n";
	} 
}
close(RESULT);
