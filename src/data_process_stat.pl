#!/Bio/bin/perl
#-----------------------------------------------+
#    [APM] This script was created by amp.pl    |
#    [APM] Created time: 2017-06-02 12:17:51    |
#-----------------------------------------------+
# name: data_process_stat.pl
# func: stat the data process of meta16/ITS
# version: 1.0

use strict;
use warnings;

use File::Basename qw/basename/;

die qq(
Usage: perl $0 <data process dir> <outdir>
\n) unless @ARGV == 2;

my $dir = shift @ARGV;
my $outdir = shift @ARGV;

my %stat;
my $index = 1;

my ($readsQC,$tagsQC);
# read reads QC result 
if (-d "$dir/$index.reads_filter")
{
	read_reads_filter("$dir/$index.reads_filter");
	$readsQC = 1;
	$index ++;
}

read_overlap_flash("$dir/$index.overlap");
$index ++;

if (-d "$dir/$index.tags_filter")
{
	read_tags_filter("$dir/$index.tags_filter");
	$tagsQC = 1;
	$index ++;
}

# the remove Chimera dir 
$index ++;

my @samples = read_effective_tags("$dir/$index.tags_stat");

# print the stat result 
open my $ofh_stat , ">" , "$outdir/all.data_process_stat.xls" or die $!;
open my $ofh_fig ,  ">" , "$outdir/all.data_process_stat.forFig.xls" or die $!;

my @attrs = grep { $stat{$samples[0]}{$_} } qw/Raw_PE Clean_PE Raw_Tags Clean_Tags Effective_Tags/;
my @attrs_head = @attrs;
@attrs_head = map { s/_/ /; $_; } @attrs_head;
my $header = join "\t" , ("Sample Name",@attrs_head,"Effective Ratio (%)");
print $ofh_stat $header , "\n";
print $ofh_fig ${ [join "\t" , qw/Sample Reads Number Ratio/ ] }[0] , "\n";
foreach my $sample (@samples)
{
	my @vals = map { $stat{$sample}{$_}; } @attrs;
	my $effective_ratio = sprintf "%.2f" , $vals[-1]*100/$vals[0];
	my $line = join "\t" , ($sample,@vals,$effective_ratio);
	print $ofh_stat "$line\n";
	
	my (@elements,@evals);
	
	my $index = 1;
	if ($readsQC)
	{
		push @elements , "Reads QC filter";
		push @evals , $vals[$index-1] - $vals[$index];
		$index ++;
	}

	push @elements , "Non-overlap";
	push @evals , $vals[$index-1] - $vals[$index];
	$index ++;

	if ($tagsQC)
	{
		push @elements , "Tags QC filter";
		push @evals , $vals[$index-1] - $vals[$index];
		$index ++;
	}

	push @elements , "Chimera";
	my $chimera = $vals[$index-1] - $vals[$index];
	push @evals , $vals[$index-1] - $vals[$index];

	push @elements , "Effective Tags";
	push @evals , $vals[-1];

	for my$i (0 .. $#evals)
	{
		my $ratio = sprintf "%.2f" , $evals[$i]*100/$vals[0];
		print $ofh_fig "$sample\t$elements[$i]\t$evals[$i]\t$ratio\n";
	}
}

close $ofh_stat;
close $ofh_fig;

#-------------------------------------------------------------------------------
#  sub function
#-------------------------------------------------------------------------------
sub read_reads_filter
{
	my $dir = shift;
	
	my @files = glob("$dir/*.stat.xls");
	die "There is no reads filter result files" unless @files > 0;

	foreach my $file (@files)
	{
		open my $fh , $file or die "can't open file, $file $!";
		<$fh>;
		my $vals = <$fh>;
		my ($sample,$rawPE,$cleanPE) = (split /\t/,$vals)[0,5,11];
		close $fh;

		$stat{$sample}{Raw_PE}   = $rawPE/2;
		$stat{$sample}{Clean_PE} = $cleanPE/2;
	}

	return 1;
}

sub read_overlap
{
	my $dir = shift;

	my @files = glob("$dir/*.log.o");
	die "There is no reads overlap result files" unless @files > 0;

	foreach my $file (@files)
	{
		my $sample = basename($file,".log.o");
		open my $fh , $file or die $!;
		my @array = <$fh>;
		close $fh;
		my ($raw_tags)  = $array[0] =~ /overlap pairs:\s*(\d+)/;
		my ($nooverlap) = $array[1] =~ /left pairs:\s*(\d+)/;
		my $rawPE = $raw_tags + $nooverlap;
		$stat{$sample}{Raw_Tags} = $raw_tags;
		$stat{$sample}{Raw_PE} = $rawPE unless $stat{$sample}{Raw_PE};
	}

	return 1;
}

sub read_overlap_flash
{
	my $dir = shift;

	my @files = glob("$dir/*.log.o");
	die "There is no reads overlap result files" unless @files > 0;

	foreach my $file (@files)
	{
		my $sample = basename($file,".log.o");
		$/ = undef;
		open my $fh , $file or die $!;
		my $str = <$fh>;
		close $fh;
		$/ = "\n";
		my ($raw_tags)  = $str =~ /Combined pairs:\s*(\d+)/;
		my ($nooverlap) = $str =~ /Uncombined pairs:\s*(\d+)/;
		my $rawPE = $raw_tags + $nooverlap;
		$stat{$sample}{Raw_Tags} = $raw_tags;
		$stat{$sample}{Raw_PE} = $rawPE unless $stat{$sample}{Raw_PE};
	}

	return 1;
}
sub read_tags_filter
{
	my $dir = shift;
	
	my @files = glob("$dir/*.QC.tags.fasta");
	die "There is no tags QC result files" unless @files > 0;

	foreach my $file (@files)
	{
		my $sample = basename($file,".QC.tags.fasta");
		my $clean_tags = `grep ">" -c $file`;
		chomp $clean_tags;
		$stat{$sample}{Clean_Tags} = $clean_tags;
	}
	
	return 1;
}

sub read_effective_tags
{
	my $dir = shift;
	
	die "the stat_tag.xls is not exists" unless -e "$dir/stat_tag.xls";
	
	my @samples;
	
	open my $fh , "$dir/stat_tag.xls";
	<$fh>;
	while(<$fh>)
	{
		my ($sample,$number) = (split /\t/)[0,1];
		$stat{$sample}{Effective_Tags} = $number;
		push @samples , $sample;
	}
	close $fh;

	return @samples;
}
