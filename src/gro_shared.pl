#!/Bio/bin/perl
#-----------------------------------------------+
#    [APM] This script was created by amp.pl    |
#    [APM] Created time: 2017-05-25 21:59:00    |
#-----------------------------------------------+
# name: gro_shared.pl
# func: 
# version: 0.1

use strict;
use warnings;
use List::Util qw/sum/;

die qq(
Usage: perl $0 <samples shared file> <group file>
\n) unless @ARGV == 2;

my $shared_file = shift @ARGV;
my $group_file  = shift @ARGV;

# read shared file 
open my $fh_shared , $shared_file or die $!;

my $header = <$fh_shared>;
my %tags;
my ($label,$sample,$num);

while(<$fh_shared>)
{
	chomp;
	my @tags;
	($label,$sample,$num,@tags) = split /\t/,$_;
	$tags{$sample} = \@tags;
}

close $fh_shared;

# create group shared file 
open my $fh_group , $group_file or die $!;
open my $ofh , ">all.otus.group.shared" or die $!;
print $ofh $header;

while(<$fh_group>)
{
	chomp;
	my ($group,$samples) = split / = / , $_;
	my @samples = split /\,/,$samples;
	
	my @tags2;
	for my $i (0 .. $num-1)
	{
		my @vals = map { $tags{$_}->[$i] } @samples;
		my $mean = int ( sum(@vals)/($#samples+1) );
		push @tags2 , $mean;
	}

	my $line = join "\t" , ($label,$group,$num,@tags2);
	print $ofh $line , "\n";
}

close $ofh;
close $fh_group;
