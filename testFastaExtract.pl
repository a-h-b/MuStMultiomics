
#!/usr/bin/perl

use Bio::DB::Fasta;

my $fastaFile = shift;
my $queryFile = shift;

my $db = Bio::DB::Fasta->new( $fastaFile );
open (IN, $queryFile);
while (<IN>){
    chomp; 
    $seq = $_;
    my $sequence = $db->seq($seq);
    if  (!defined( $sequence )) {
            print STDERR "Sequence $seq not found. \n";
            next;
    }   
    print ">$seq\n", "$sequence\n";
}
