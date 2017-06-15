#!/Bio/bin/perl
#-----------------------------------------------+
#    [APM] This script was created by amp.pl    |
#    [APM] Created time: 2017-06-04 18:21:46    |
#-----------------------------------------------+
# name: qiime2ggplot2.pl
# func: turn qiime result file of alpha diversity to ggplot2 readable format
# version: 1.0

use strict;
use warnings;
use utf8;
use File::Basename qw/basename/;

die qq(
Usage: perl $0 <qiime alpha diversity average table>
\n) unless @ARGV == 1;

my $file = shift @ARGV;
my $file_name = basename($file);
$file_name =~ s/\.txt$//;

open my $fh , $file or die "can't open the input file, $file, $!";
open my $ofh , ">" , "$file_name.forGGPLOT2.txt" or die $!; 

$/ = ">> ";

# fetch the name of values ($name) and is Group or Sample info ($flag)
# fetch the xaxis vals and x axis max val 
my $head = <$fh>;
chomp $head;

my ($name,$flag,$xaxis,$xmax) = split /\n/,$head;

$name  = $name =~ /^\# (\w+)\.txt/ ? $1 : die "input file's format is error!";
$flag  = $flag =~ /^\# (\w+)/ ? $1 : die "input file's format is error!";
$xaxis = $xaxis =~ /^xaxis: (.+)\t/ ? $1 : die "input file's format is error!";
$xmax  = $xmax =~ /xmax: ([\w\.]+)/ ? $1 : die "input file's format is error!";

my @xvals = split /\t/ , $xaxis;

print $ofh "$flag\tcolor\tx\tval\tse\n";
while(<$fh>)
{
	chomp;
	my ($sample,$color,$vals,$errors) = split /\n/;
	
	$color  = (split /\s+/,$color)[1];
	$vals   = (split /\s+/,$vals,2)[1];
	$errors = (split /\s+/,$errors,2)[1];
	$errors =~ s/nan/0/g;

	my @vals   = split /\t/,$vals;
	my @errors = split /\t/,$errors;

	# remove the last <TAB>
	pop @vals;
	pop @errors;

	for my$i (0 .. $#vals)
	{
		print $ofh "$sample\t$color\t$xvals[$i]\t$vals[$i]\t$errors[$i]\n";
	}
}

close $fh;
close $ofh;

$/ = "\n";
