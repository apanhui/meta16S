#!/Bio/bin/perl
#-----------------------------------------------+
#    [APM] This script was created by amp.pl    |
#    [APM] Created time: 2017-06-14 16:06:58    |
#-----------------------------------------------+
# name: diversity_diff_summary.pl
# func: 
# version: 1.0

use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use File::Basename qw/basename/;

my %OPTS = ();
GetOptions(\%OPTS,'alpha','beta');

&usage unless @ARGV == 1;

my $dir = abs_path(shift @ARGV);

my @alpha = qw/chao1 ace goods_coverage observed_species shannon simpson/;
my @beta  = qw/weighted_unifrac unweighted_unifrac/;

my (@subdirs,$index,$outname);

if ($OPTS{alpha})
{
	@subdirs = @alpha;
	$index   = "Alpha diversity";
	$outname = "all.alpha_diversity.groups_diff.summary.xls";
}
elsif ($OPTS{beta})
{
	@subdirs = @beta;
	$index   = "Beta diversity";
	$outname = "all.beta_diversity.groups_diff.summary.xls";
}
else 
{
	die "either --alpha or --beta must be defined";
}

my %stat;
my @methods = qw/ttest Wilcox Kruskal Tukey/;
my %method_fullnames = (ttest=>"T test",Wilcox=>"Wilcoxon",Kruskal=>"Kruskal-Wallis",Tukey=>"Tukey HSD");

open my $ofh_summary , ">" , "$outname" or die $!;

my $head = "$index\tdiff_groups\ttest method\tpvalue\tsignificant\n";
print $ofh_summary $head;

foreach my $subdir (@subdirs)
{
	die "[$dir/$subdir] is not exists, please check!" unless -d "$dir/$subdir";
	
	my @two_groups = read_normal_test  ("$dir/$subdir","ttest");
	read_normal_test ("$dir/$subdir","Wilcox");

	my @methods = qw/ttest Wilcox/;
	# output the summary result 
	foreach my $group (@two_groups)
	{
		foreach my $method (@methods)
		{
			my $pvalue = $stat{$group}{$method};
			my $sig = fetch_significant($pvalue);
			print $ofh_summary "$subdir\t$group\t$method_fullnames{$method}\t$pvalue\t$sig\n";
		}
	}

	my @multi_groups = read_normal_test("$dir/$subdir","Kruskal");
	read_Tukey  ("$dir/$subdir") if (@multi_groups);

	@methods = qw/Kruskal Tukey/;
	# output the summary result 
	foreach my $group (@multi_groups)
	{
		foreach my $method (@methods)
		{
			my $pvalue = $stat{$group}{$method};
			my $sig = fetch_significant($pvalue);
			print $ofh_summary "$subdir\t$group\t$method_fullnames{$method}\t$pvalue\t$sig\n";
		}
	}
}

close $ofh_summary;

#-------------------------------------------------------------------------------
#  sub function
#-------------------------------------------------------------------------------
sub read_normal_test
{
	my $odir   = shift;
	my $method = shift;
	
	my %files  = ("ttest"=>"Ttest.txt",Wilcox=>"Wilcox.txt",Kruskal=>"Kruskal.txt");
	
	my $res = "$odir/$files{$method}";

	return () unless -e $res;

	my @groups;
	open my $fh , $res or die $!;
	<$fh>;
	while(<$fh>)
	{
		chomp;
		my ($group,$pvalue) = split /\t/;
		$stat{$group}{$method} = $pvalue;
		push @groups , $group;
	}
	close $fh;

	return @groups;
}

sub read_Tukey
{
	my $odir = shift;

	my @files = glob("$odir/*.TukeyHSD.txt");
	
	foreach my $file (@files)
	{
		my $group = basename($file,".TukeyHSD.txt");
		open my $fh , $file or die $!;
		<$fh>;
		my $min = 1;
		while(<$fh>)
		{
			chomp;
			my ($group,$pvalue) = (split /\t/)[0,-1];
			$min = $pvalue if ($min > $pvalue);
		}
		close $fh;

		$stat{$group}{Tukey} = $min;
	}
}

sub fetch_significant
{
	my $pvalue = shift;

	if ($pvalue <= 0.01)
	{
		return "**";
	}
	elsif ($pvalue <= 0.05)
	{
		return "*";
	}
	else 
	{
		return "";
	}
}

sub usage
{
	my $help = <<HELP;

Usage: perl $0 <--alpha|--beta> <dir>

HELP
	print $help;
	exit;
}
