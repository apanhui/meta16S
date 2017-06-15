#!/Bio/bin/perl
#-----------------------------------------------+
#    [APM] This script was created by amp.pl    |
#    [APM] Created time: 2017-06-09 19:44:24    |
#-----------------------------------------------+
# name: tags_otus_count.pl
# func: 
# version: 1.0

use strict;
use warnings;
use List::Util qw/sum/;

my $dir = shift @ARGV;

my %stat;

# stat the total tags 
stat_tags("$dir/01.data_process/5.tags_stat/stat_tag.xls","Total_Tags");


# stat the unique tags 
stat_tags("$dir/01.data_process/5.tags_stat/stat_uniqTag.xls","Unique_Tags");

# stat the otu number 
my @samples = stat_otus("$dir/02.otus/all.otus.shared");

# stat the taxon and uniclassfied tags 
stat_taxon("$dir/03.taxonomy/2.taxa_profiling/all.otus.profiling.xls");

my @attrs = ("Total_Tags","Unique_Tags","Taxon_Tags","Unclassified_Tags","Singleton_Tags","OTUs");

open my $ofh , ">" , "all.tags_otus.stat.xls" or die $!;
open my $ofh_fig , ">" , "all.tags_otus.stat.forFig.xls" or die $!;

my $header = join "\t" , ("SampleID",@attrs);
$header =~ s/_/ /g;
print $ofh $header , "\n";

my $fig_head = "SampleID\tattr\tValue\n";
print $ofh_fig $fig_head;

foreach my $sample (@samples)
{
	my @vals = map { $stat{$sample}{$_} || 0} @attrs;
	my $line = join "\t" , ($sample,@vals);
	print $ofh $line , "\n";
	
	foreach (0 .. $#attrs)
	{
		print $ofh_fig "$sample\t$attrs[$_]\t$vals[$_]\n";
	}

}

# print avg vals 
my @avg_vals = map {
	my $attr = $_;
	my @vals = map { $stat{$_}{$attr} } @samples;
	my $avg  = int sum(@vals)/($#vals+1);
	$avg;
} @attrs;

my $line = join "\t" , ("Avg",@avg_vals);
print $ofh $line , "\n";

for my$i (0 .. $#attrs)
{
	print $ofh_fig "Avg\t$attrs[$i]\t$avg_vals[$i]\n";
}

close $ofh;
close $ofh_fig;

#-------------------------------------------------------------------------------
#  sub function
#-------------------------------------------------------------------------------
sub stat_tags
{
	my $file = shift;
	my $attr = shift;

	open my$fh , $file or die $!;
	<$fh>;
	while(<$fh>)
	{
		chomp;
		my ($sample,$num) = (split /\t/)[0,1];
		$stat{$sample}{$attr} = $num;
	}
	close $fh;
}

sub stat_otus
{
	my $file = shift;
	my @samples ;
	open my$fh , $file or die $!;
	<$fh>;
	while(<$fh>)
	{
		chomp;
		my ($label,$sample,$total,@numbers) = split /\t/;
		@numbers = grep { $_ != 0 } @numbers;
		$stat{$sample}{OTUs} = scalar @numbers;
		push @samples , $sample;
	}
	close $fh;
	return @samples;
}

sub stat_taxon
{
	my $file = shift ;
	open my$fh , $file or die $!;
	my $header = <$fh>;
	
	my $sample_num = scalar @samples;
	
	while(<$fh>)
	{
		chomp;
		my @arr = split /\t/;
		my @taxons = @arr[$sample_num*2+2 .. $sample_num*2+8];
		my @tags = @arr[2 .. $sample_num+1];
		
		for my $i (0 .. $#samples)
		{
			if ($taxons[0] ne "")
			{
				$stat{$samples[$i]}{Taxon_Tags} += $tags[$i]
			}
			else 
			{
				$stat{$samples[$i]}{Unclassified_Tags} += $tags[$i]
			}
		}
	}
	
	map { $stat{$_}{Unclassified_Tags} ||= 0 } @samples;
	map { $stat{$_}{Singleton_Tags} = $stat{$_}{Total_Tags} - $stat{$_}{Taxon_Tags} - $stat{$_}{Unclassified_Tags} } @samples;
}
