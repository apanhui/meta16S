#!/Bio/bin/perl
#-----------------------------------------------+
#    [APM] This script was created by amp.pl    |
#    [APM] Created time: 2017-02-22 15:45:42    |
#-----------------------------------------------+
# name: uparse2mothur.pl
# func: turn uparse result to Muthor format 
# version: 2.0
# update: change shared input file to table file 

use strict;
use warnings;

use List::Util qw/sum/;
use Getopt::Std;
use lib "/home/aipeng/work/develepment/lib";
use SEQUENCE qw/read_fasta/;

my %opts = (n=>2);
getopts('n:s:',\%opts);

die qq(
Usage: perl $0 [-n 2] [-s samples] <raw otu sequence file> <raw uparse otu table> <new otu sequence file> <new otu shared file>

Note:  -n is the minmum tags numner allowed as a true OTU, default is 2.
\n) unless @ARGV == 4;

my $otus = shift @ARGV;
my $tab = shift @ARGV;
my $newotus = shift @ARGV;
my $newtab = shift @ARGV;

my @samples;
my %otus = read_out_table($tab);
my %seqs = read_fasta($otus);

# reset the samples order 
@samples = split /,/ , $opts{s} if ($opts{s});

open my $ofh_tab , ">$newtab" or die "can't open file, $newtab $!";
open my $ofh_seq , ">$newotus" or die "can't open file, $newotus $!";

my $count = "000000";
my %shared;
my @new_otuids;
my $otu_num = scalar keys %otus;
foreach my $otuid (sort {$otus{$b}{total} <=> $otus{$a}{total}} keys %otus)
{
	$count ++;
	my $newid = "Otu$count";
	push @new_otuids , $newid;
	
	for my $i (0 .. $#samples)
	{
		push @{$shared{$samples[$i]}} , $otus{$otuid}{tags}->[$i];
	}

	# print new otu represstive sequence
	my $seq = $seqs{$otuid};
	print $ofh_seq ">$newid\n$seq\n";
}

# print new shared file 
my $header = join "\t" , ("label","Group","numOtus",@new_otuids);
print $ofh_tab $header , "\n";
foreach my $sample (@samples)
{
	my $line = join "\t" , ("0.03",$sample,$otu_num,@{$shared{$sample}});
	print $ofh_tab $line , "\n";
}

close $ofh_seq;
close $ofh_tab;

#-------------------------------------------------------------------------------
#  sub function
#-------------------------------------------------------------------------------
# read the uparse otu tab file 
sub read_out_table
{
	my $tab = shift;

	open TAB,$tab or die "can't open file, $tab $!";
	my $header = <TAB>;
	chomp $header;
	@samples = split /\t/ , $header;
	shift @samples;

	my %otus;
	while(<TAB>)
	{
		chomp;
		my ($otuid,@tags) = split /\t/;
		my $sum = sum(@tags);

		next if ($sum < $opts{n});
		$otus{$otuid}{total} = $sum;
		$otus{$otuid}{tags} = \@tags;
	}
	close TAB;
	
	return %otus;
}
