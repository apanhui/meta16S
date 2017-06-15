#!/Bio/bin/perl
#-----------------------------------------------+
#    [APM] This script was created by amp.pl    |
#    [APM] Created time: 2017-02-15 18:32:20    |
#-----------------------------------------------+
# name: fetch_non_chimera.pl
# func: 
# version:: 0.1

use strict;
use warnings;
use Bio::SeqIO;

my ($fasta,$chimera) = @ARGV;

# hash to save chimera ids
my %hash;
open IN,$chimera or die $!;
while(<IN>)
{
	chomp;
	my ($id,$flag) = (split /\t/)[1,-1];
	$hash{$id} = $flag if ($flag eq "Y");
}
close IN;

# read the total sequence and fetch the chimera
my $seqio = Bio::SeqIO->new(-file=>$fasta,-format=>"fasta");

while(my$seq = $seqio->next_seq)
{
	my $id = $seq->id;
	next if $hash{$id};
	
	my $seqstr = $seq->seq;
	print ">$id\n$seqstr\n";
}
