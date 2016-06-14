#!/usr/bin/perl -w

# script to include variant proteins in protein database
# takes 4 inputs: -a, a table with the gene positions, starts and stops (like the .tab output of prodigal);
#		  -s, a fasta file with the contigs
#                 -v, a .vcf file of the same contigs
#		  -p, if not given, no fasta files will be generated
# 4 outputs: - a fasta file with the amino acid sequences of all proteins, including variants; <name>.variants.faa
#	     - a fasta file with the nucleotide sequences of all genes, including variants; <name>.variants.fna
# 	     - a table with the positions of the predicted variant genes, their start and stop codons <name>.variants.tab
#	     - a log file <name>.variants.log
#
#written by Nicolas Pinel (July 04, 2013), some errors repaired by Anna Heintz Buschart (Jan 29, 2015 - Jun 14, 2016)
# this version handles the mix up of start- and stop-codons in the prodigal .tab output and returns a .tab with correct start- and stop-codons



use strict;
use Getopt::Long;
use Bio::SeqIO;
#use Array::Utils qw(:all);

# Obtain inputs
my ($ann,$seq,$vcf,$print);
GetOptions('a|ann=s' => \$ann,
	   's|seq=s' => \$seq,
	   'v|vcf=s' => \$vcf,
	   'p|print' => \$print);

# Load sequences
my $seqbank = &loadSequences($seq);

# Load annotation (tab file)
my $annot = &loadAnnotation($ann);




# Load variant file (subroutine written)
&loadVariants($vcf);

# Load amino acid table (manually written subroutine)
my $aa_table = &loadAA();
my $syn = 0;
my $non = 0;
my $fix = 0;
(my $prov_out = $seq) =~ s/\.\w+$/\.variants.faa/;
(my $dnav_out = $seq) =~ s/\.\w+$/\.variants.fna/;
(my $tabv_out = $seq) =~ s/\.\w+$/\.variants.tab/;
(my $log_out = $seq) =~ s/\.\w+$/\.variants.log/;

# Print output files
my ($o_dnav,$o_prov,$o_tabv, $o_log);
if ($print) {
    open $o_dnav, '>', $dnav_out;
    open $o_prov, '>', $prov_out;
    open $o_tabv, '>', $tabv_out;
    open $o_log, '>', $log_out;
}
foreach my $k (sort keys %{$annot}) {
    foreach my $feat (@{$$annot{$k}}) {
	my $id = $$feat{'gene'};
	print qq($k\t$id\toriginal\n);
	if ($$feat{'variants'}) {
	    my ($prov,$idv);
	    my @nuc = &getSeq($k, $feat);
	    my $pro = &translate($nuc[0]);
	    print $o_dnav qq(\>$id\n$nuc[0]\n) if ($print);
	    print $o_prov qq(\>$id\n$pro\n) if ($print);
	    my @varcoll = &getAllVars($feat,$$feat{'variants'});
	    my $vcnt = 0;
	    my $vind = 0;
	    ## and now print variants by possibilites per feature
	    while ($varcoll[$vcnt]){
		my @nucvs = &getSeq($k, $feat, $varcoll[$vcnt]);
		foreach my $nucv (@nucvs) {
		    $idv = $id.'.variant'.$vind;
		    $prov = &translate($nucv);
		    print $o_dnav qq(\>$idv\n$nucv\n) if ($print);
		    my $stopc = $$feat{'stopcodon'};
		    my $compl = $$feat{'completeness'};
		    if ($pro eq $prov) {
			++$syn;
		    } else {
			++$non;
			if ($$feat{'startcodon'} eq "yes" && substr($prov,0,3) ne substr($pro,0,3)){
			    print $o_log qq($idv has variant on startcodon \n) if ($print);
			}
			if (substr($pro,-1) eq "*" && substr($prov,-1) ne "*") {
			    $stopc = "no";
			    $compl = "incomplete";
			    print $o_log qq($idv has lost stopcodon \n) if ($print);
			} elsif ($$feat{'stopcodon'} eq "no" && substr($prov,-1) eq "*"){
			    $stopc = "yes";
			    if ($$feat{'startcodon'} eq "yes") {
				$compl = "complete";
			    }
			    print $o_log qq($idv has gained stopcodon \n) if ($print);
			}
			print $o_prov qq(\>$idv\n$prov\n) if ($print);
			print $o_tabv qq(@{[join("\t", $idv,$$feat{'strand'},$$feat{'start'},$$feat{'stop'},"na",$$feat{'startcodon'},$stopc,$compl)]}\n) if ($print);
		    }
		    $vind++;
		}
		$vcnt++;
		$vind++;
	    }
	} else {
	    ++$fix;
	    my @nuc = &getSeq($k, $feat);
	    my $pro = &translate($nuc[0]);
	    print $o_dnav qq(\>$id\n$nuc[0]\n) if ($print);
	    print $o_prov qq(\>$id\n$pro\n) if ($print);
	}
    }
}
print qq(invariant\: $fix\nsynonymous\: $syn\nnonsynon\: $non\n);

exit;

#### subroutines ###

sub usage {
    
    exit;
}

sub array_minus(\@\@) {
          my %e = map{ $_ => undef } @{$_[1]};
	          return grep( ! exists( $e{$_} ), @{$_[0]} );
		}

sub loadSequences {
    my $f = shift;
    my $inseq = Bio::SeqIO->new(-file => $f,
				-format => 'fasta');
    my %seqs = ();
    while (my $seqobj = $inseq->next_seq()) {
	$seqs{$seqobj->display_id()} = $seqobj->seq();
    }
    return \%seqs;
}

sub loadAnnotation {
    my $f = shift;
    open my $in, '<', $f;
    my @headers = ('gene','strand','start','stop','length','startcodon','stopcodon','completeness','startcodonvar','stopcodonvar','completenessvar');
    my %feats = ();
    while (<$in>) {
	next if (/^\#/);
	chomp;
	my %hash = ();
	my @f = split(/\t/);
	next if ($f[1] =~ /CUFF/);
	my $ctg = $f[0];
	$ctg =~ s/_\d+$//g;
	$ctg =~ s/gene//;
	if ($f[1] eq "+") {
	    @hash{@headers} = (@f,@f[5 .. $#f]);
	} else{
	    @hash{@headers} = (@f[0 .. 4],$f[6],$f[5],$f[7],$f[6],$f[5],$f[7]);
	}
	push(@{$feats{$ctg}}, \%hash);
    }
    return \%feats;
}

sub loadVariants {
    my $f = shift;
    if ($f =~ /gz$/) {
	(my $nf = $f) =~ s/\.gz$//;
	system("bgzip -c -d $f > $nf");
	$f = $nf;
    }
    open my $in, '<', $f;
    while (<$in>) {
	next if (/^\#/);
	chomp;
	my @f = split(/\t/);
	foreach my $feat (@{$$annot{$f[0]}}) {
	    if ($f[1] >= $$feat{'start'} && $f[1] <= $$feat{'stop'}) {
		push(@{$$feat{'variants'}}, {'pos' => $f[1], 'ref' => $f[3], 'new' => $f[4]});
	    }
	}
    }
}

sub loadAA {
    my %aa = ('ATT' => 'I', 'ATC' => 'I', 'ATA' => 'I', 'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG' => 'L', 'TTA' => 'L',
	      'TTG' => 'L', 'GTT' => 'V', 'GTC' => 'V', 'GTA' => 'V', 'GTG' => 'V', 'TTT' => 'F', 'TTC' => 'F', 'ATG' => 'M',
	      'TGT' => 'C', 'TGC' => 'C', 'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A', 'GGT' => 'G', 'GGC' => 'G',
	      'GGA' => 'G', 'GGG' => 'G', 'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P', 'CCG' => 'P', 'ACT' => 'T', 'ACC' => 'T',
	      'ACA' => 'T', 'ACG' => 'T', 'TCT' => 'S', 'TCC' => 'S', 'TCA' => 'S', 'TCG' => 'S', 'AGT' => 'S', 'AGC' => 'S',
	      'TAT' => 'Y', 'TAC' => 'Y', 'TGG' => 'W', 'CAA' => 'Q', 'CAG' => 'Q', 'AAT' => 'N', 'AAC' => 'N', 'CAT' => 'H',
	      'CAC' => 'H', 'GAA' => 'E', 'GAG' => 'E', 'GAT' => 'D', 'GAC' => 'D', 'AAA' => 'K', 'AAG' => 'K', 'CGT' => 'R',
	      'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R', 'AGA' => 'R', 'AGG' => 'R', 'TAA' => '*', 'TAG' => '*', 'TGA' => '*');
    return \%aa;
}

sub translate {
    my $s = shift;
    my @codons = $s =~ /(.{3})/g;
    my $aa = '';
    foreach my $cod (@codons) {
	my $add = $$aa_table{$cod} ? $$aa_table{$cod} : 'X';
	$aa .= $add;
	if ($add eq "*") {
	    last
	}
    }
    return $aa;
}


sub getAllVars {
    my ($ref, $variants) = @_;
    my @varcoll = ();
    my $vnum = 0;
    foreach my $v (sort {$$a{'pos'} <=> $$b{'pos'}} @{$variants}){
	push @varcoll, $v;
	$vnum++;
    }
    my @newcoll = ();
    if ($vnum > 1) {
	my @testcol = [@varcoll];
	@newcoll = &getAdjVars(@testcol);
	my $i = 0;
	while ($newcoll[$i]){
	    my $ii = 0;
	    while ($newcoll[$ii]){
		my @var1 = ();
		my $j = 0;
		while ($newcoll[$i][$j]) {
		    push @var1, $newcoll[$i][$j]{'pos'};
		    $j++;
		}
		my $jj = 0;
		my @var2 = ();
		while ($newcoll[$ii][$jj]) {
			push @var2, $newcoll[$ii][$jj]{'pos'};
		    $jj++;
		}
		my @minus = array_minus(@var1,@var2);
		if ($i != $ii && $#minus == -1) {
		    #print "Yay";
		    splice @newcoll, $i, 1;
		    $ii = 0;
		} else {
		    $ii++;
		}
	    }
	    $i++;
	}
    } else {
	@newcoll = [@varcoll];
    }
    return @newcoll;
}


sub getAdjVars{
    my (@variants) = @_;
    my @varcoll = ();
    foreach my $v (@variants){
	push @varcoll,$v;
    }
    my @newcoll= ();
    my $i=0;
    while ($varcoll[$i]) {
	my $j = 0;
	while ($varcoll[$i][$j]) {
	    push @newcoll,$varcoll[$i][$j];
	    $j++;
	}
	$i++;
    }
    my @newcollA = [@newcoll];
    my $vid =0;
    while ($newcoll[$vid+1]){
	if ($newcoll[$vid]{'pos'}+length($newcoll[$vid]{'ref'})>$newcoll[$vid+1]{'pos'} ) {
	    my @idx1 = (0..$#newcoll);
	    my @idx2 = (0..$#newcoll);
	    splice @idx1,$vid,1;
	    splice @idx2,$vid+1,1;
	    my @test1 = [@newcoll[@idx1]];
	    my @newcoll1 = &getAdjVars(@test1);
	    my @test2 = [@newcoll[@idx2]];
	    my @newcoll2 = &getAdjVars(@test2);
	    shift @newcollA;
	    push @newcollA, @newcoll1;
	    push @newcollA, @newcoll2;
	    last;
	} else {
	    $vid++;
	}
    }
    return @newcollA;
}
    

sub getSeq {
    my ($ctg,$ref, $variants) = @_;
    my $length = $$ref{'stop'} - $$ref{'start'} + 1;
    my $seq = substr($$seqbank{$ctg}, $$ref{'start'} - 1, $length);
    my @seqs = ();
    if ($variants && $$ref{'strand'} eq "+") {
	push @seqs,$seq;
	foreach my $v (sort {$$b{'pos'} <=> $$a{'pos'}} @{$variants}) {
	    my $seqnum = 0;
	    while ($seqnum <= $#seqs){
		my $ss1 = substr($seqs[$seqnum], 0, ($$v{'pos'} - $$ref{'start'}));
		my $secstart = ($$v{'pos'} - $$ref{'start'} + length($$v{'ref'}));
		my $ss2 ="";
		unless ($secstart > length($seqs[$seqnum])) {
		    $ss2 = substr($seqs[$seqnum], ($$v{'pos'} - $$ref{'start'} + length($$v{'ref'})));
		}
		my $c = $$v{'new'} =~ tr/,//;
		if ($c >0) {
		    splice @seqs,$seqnum,1;
		    my @allAlts = split(",",$$v{'new'});
		    foreach my $alt (@allAlts){
			unshift @seqs, $ss1.$alt.$ss2;
			if ((length($alt)-length($$v{'ref'}))%3 !=0 ){
			    print qq(frameshift in $$ref{'gene'}\n);
			}
			$seqnum++;
		    }
		} else {
		    $seqs[$seqnum] = $ss1.$$v{'new'}.$ss2;
		    if ((length($$v{'new'})-length($$v{'ref'}))%3 !=0 ){
			print qq(frameshift in $$ref{'gene'}\n);
		    }
		    $seqnum++;
		}
	    }
	}
    } elsif ($variants && $$ref{'strand'} eq "-") {
	push @seqs,$seq;
	foreach my $v (sort {$$b{'pos'} <=> $$a{'pos'}} @{$variants}) {
	    my $seqnum = 0;
	    while ($seqnum <= $#seqs){
		my $ss1 = substr($seqs[$seqnum], 0, ($$v{'pos'} - $$ref{'start'}));
		my $secstart = ($$v{'pos'} - $$ref{'start'} + length($$v{'ref'}));
		my $ss2 ="";
		unless ($secstart > length($seqs[$seqnum])) {
		    $ss2 = substr($seqs[$seqnum], ($$v{'pos'} - $$ref{'start'} + length($$v{'ref'})));
		}
		my $c = $$v{'new'} =~ tr/,//;
		if ($c >0) {
		    splice @seqs,$seqnum,1;
		    my @allAlts = split(",",$$v{'new'});
		    foreach my $alt (@allAlts){
			unshift @seqs, $ss1.$alt.$ss2;
			if ((length($alt)-length($$v{'ref'}))%3 !=0 ){
			    print qq(frameshift in $$ref{'gene'}\n);
			}
			$seqnum++;
		    }
		} else {
		    $seqs[$seqnum] = $ss1.$$v{'new'}.$ss2;
		    if ((length($$v{'new'})-length($$v{'ref'}))%3 !=0 ){
			print qq(frameshift in $$ref{'gene'}\n);
		    }
		    $seqnum++;
		}
	    }
	}
	my $i = 0;
	while ($seqs[$i]) {
	    $seqs[$i] = join('', reverse(split(//, $seqs[$i])));
	    $seqs[$i] =~ tr/ACGT/TGCA/;
	    $i++;
	}
    } elsif ($$ref{'strand'} eq '-') {
	$seq = join('', reverse(split(//, $seq)));
	$seq =~ tr/ACGT/TGCA/;
	push @seqs, $seq;
    } else {
	push @seqs, $seq;
    }
    return @seqs;
}
