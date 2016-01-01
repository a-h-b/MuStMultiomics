#!/usr/bin/perl -w


# written by Venkata Satagopam, with a few corrections by Anna Heintz-Buschart

#takes 6 inputs: the name of the output, will be concatenated with the database name and .tsv; the results of hmmscan for kegg, metacyc, pfam, swissprot and tigrpfam

use strict;

my $result_file=$ARGV[0];
my $kegg_hmmscan_result_file=$ARGV[1];
my $metacyc_hmmscan_result_file=$ARGV[2];
my $pfam_hmmscan_result_file=$ARGV[3];
my $swissprot_hmmscan_result_file=$ARGV[4];
my $tigrpfam_hmmscan_result_file=$ARGV[5];

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
	print STDERR "query_name : $query_name line: $line\n";		
	push @{$kegg_hmmscan_result_hash{$query_name}}, $line;
	$hit_hash{$query_name} = 1;
} 
close(FH);

my $outfile=$result_file."kegg.tsv";
open RESULT,">$outfile";
print RESULT "QUERY_NAME\tHIT_NAME\tSCORE\tSIGNIFICANCE\n";

#kegg
foreach my $query_name (keys %hit_hash) { 
	foreach $line (@{$kegg_hmmscan_result_hash{$query_name}}){
		my ($hit_name, $acc1, $query_name, $acc2,$seqEvalue,$seqscore,$seqbias,$dom1Evalue,$dom1score,$dom1bias,$exp,$reg,$clu,$ov,$env,$dom,$rep,$inc,$desc_of_target) = split(/ +/, $line);
		print RESULT "$query_name\t$hit_name\t$seqscore\t$seqEvalue\n";
	} 
}
close(RESULT);
%kegg_hmmscan_result_hash = ();
%hit_hash = ();

### Pfam
my %pfam_hmmscan_result_hash = ();
close(FH);
open(FH, "$pfam_hmmscan_result_file") or die "Failed to open $pfam_hmmscan_result_file\n";
@all = <FH>;	
shift @all; # to get rid off header
while (@all) {
	$line = shift @all;
	next if $line=~/^#/;
	chomp($line);
	my ($hit_name, $acc1, $query_name, $acc2,$seqEvalue,$seqscore,$seqbias,$dom1Evalue,$dom1score,$dom1bias,$exp,$reg,$clu,$ov,$env,$dom,$rep,$inc,$desc_of_target) = split(/ +/, $line);
	push @{$pfam_hmmscan_result_hash{$query_name}}, $line;
	$hit_hash{$query_name} = 1;
} 
close(FH);

my $outfile=$result_file."pfam.tsv";
open RESULT,">$outfile";
print RESULT "QUERY_NAME\tHIT_NAME\tSCORE\tSIGNIFICANCE\n";

foreach my $query_name (keys %hit_hash) {

    # pfam
    foreach $line (@{$pfam_hmmscan_result_hash{$query_name}}){
	 my ($hit_name, $acc1, $query_name, $acc2,$seqEvalue,$seqscore,$seqbias,$dom1Evalue,$dom1score,$dom1bias,$exp,$reg,$clu,$ov,$env,$dom,$rep,$inc,$desc_of_target) = split(/ +/, $line);
	 print RESULT "$query_name\t$hit_name\t$seqscore\t$seqEvalue\n";
    }				          
}
close(RESULT);

%pfam_hmmscan_result_hash = ();
%hit_hash = ();

### metacyc
my %metacyc_hmmscan_result_hash = ();
close(FH);
open(FH, "$metacyc_hmmscan_result_file") or die "Failed to open $metacyc_hmmscan_result_file\n";
@all = <FH>;	
shift @all; # to get rid off header
while (@all) {
	$line = shift @all;
	chomp($line);
	next if $line=~/^#/;
	my ($hit_name, $acc1, $query_name, $acc2,$seqEvalue,$seqscore,$seqbias,$dom1Evalue,$dom1score,$dom1bias,$exp,$reg,$clu,$ov,$env,$dom,$rep,$inc,$desc_of_target) = split(/ +/, $line);
	push @{$metacyc_hmmscan_result_hash{$query_name}}, $line;
	$hit_hash{$query_name} = 1;
} 
close(FH);

$outfile=$result_file."metacyc.tsv";
open RESULT,">$outfile";
print RESULT "QUERY_NAME\tHIT_NAME\tSCORE\tSIGNIFICANCE\n";

foreach my $query_name (sort keys %hit_hash) {
    # metacyc
    foreach $line (@{$metacyc_hmmscan_result_hash{$query_name}}){
	my ($hit_name, $acc1, $query_name, $acc2,$seqEvalue,$seqscore,$seqbias,$dom1Evalue,$dom1score,$dom1bias,$exp,$reg,$clu,$ov,$env,$dom,$rep,$inc,$desc_of_target) = split(/ +/, $line);
	print RESULT "$query_name\t$hit_name\t$seqscore\t$seqEvalue\n";
    }	  
}
close(RESULT);

%metacyc_hmmscan_result_hash = ();
%hit_hash = ();

### swissprot
my %swissprot_hmmscan_result_hash = ();
close(FH);
open(FH, "$swissprot_hmmscan_result_file") or die "Failed to open $swissprot_hmmscan_result_file\n";
@all = <FH>;
shift @all; # to get rid off header
while (@all) {
        $line = shift @all;
        chomp($line);
	next if $line=~/^#/;
	my ($hit_name, $acc1, $query_name, $acc2,$seqEvalue,$seqscore,$seqbias,$dom1Evalue,$dom1score,$dom1bias,$exp,$reg,$clu,$ov,$env,$dom,$rep,$inc,$desc_of_target) = split(/ +/, $line);
        push @{$swissprot_hmmscan_result_hash{$query_name}},$line;
        $hit_hash{$query_name} = 1;
} 
close(FH);

$outfile=$result_file."swissprot.tsv";
open RESULT,">$outfile";
print RESULT "QUERY_NAME\tHIT_NAME\tSCORE\tSIGNIFICANCE\n";

foreach my $query_name (sort keys %hit_hash) {

    # swissprot
    foreach $line (@{$swissprot_hmmscan_result_hash{$query_name}}){
	my ($hit_name, $acc1, $query_name, $acc2,$seqEvalue,$seqscore,$seqbias,$dom1Evalue,$dom1score,$dom1bias,$exp,$reg,$clu,$ov,$env,$dom,$rep,$inc,$desc_of_target) = split(/ +/, $line);
	print RESULT "$query_name\t$hit_name\t$seqscore\t$seqEvalue\n";
    }
}
close(RESULT);

%swissprot_hmmscan_result_hash = ();
%hit_hash = ();

### tigrpfam
my %tigrpfam_hmmscan_result_hash = ();
close(FH);
open(FH, "$tigrpfam_hmmscan_result_file") or die "Failed to open $tigrpfam_hmmscan_result_file\n";
@all = <FH>;
shift @all; # to get rid off header
while (@all) {
        $line = shift @all;
        chomp($line);
	next if $line=~/^#/;
	my ($hit_name, $acc1, $query_name, $acc2,$seqEvalue,$seqscore,$seqbias,$dom1Evalue,$dom1score,$dom1bias,$exp,$reg,$clu,$ov,$env,$dom,$rep,$inc,$desc_of_target) = split(/ +/, $line);
	push @{$tigrpfam_hmmscan_result_hash{$query_name}}, $line;
        $hit_hash{$query_name} = 1;
}  
close(FH);      

$outfile=$result_file."tigrpfam.tsv";
open RESULT,">$outfile";
print RESULT "QUERY_NAME\tHIT_NAME\tSCORE\tSIGNIFICANCE\n";

foreach my $query_name (sort keys %hit_hash) {

    # tigrpfam  
    foreach $line (@{$tigrpfam_hmmscan_result_hash{$query_name}}){
        my ($hit_name, $acc1, $query_name, $acc2,$seqEvalue,$seqscore,$seqbias,$dom1Evalue,$dom1score,$dom1bias,$exp,$reg,$clu,$ov,$env,$dom,$rep,$inc,$desc_of_target) = split(/ +/, $line);
	print RESULT "$query_name\t$hit_name\t$seqscore\t$seqEvalue\n";
    }	      
}
close RESULT;
%tigrpfam_hmmscan_result_hash = ();
%hit_hash = ();

