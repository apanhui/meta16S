#!/Bio/bin/perl
#-----------------------------------------------+
#    [APM] This script was created by amp.pl    |
#    [APM] Created time: 2017-04-01 15:37:04    |
#-----------------------------------------------+
# name: rdp2profiling.pl
# func: turn rdp taxon assign result to OTU profiling file and taxonomy profiling files 
# version: 0.1

use strict;
use warnings;
use List::Util qw(sum);

use Getopt::Long;
my %opts = ('mtfp'=>2,'mohp'=>0.1);

GetOptions(\%opts,
'mtfp:f',
'mohp:f'
);

die &usage unless @ARGV == 3;

my $otu_table  = shift @ARGV;
my $rdp_result = shift @ARGV;
my $outdir     = shift @ARGV;

mkdir $outdir unless -d $outdir;
my @taxon_levels  = qw/Domain Phylum Class Order Family Genus Species/;

#-------------------------------------------------------------------------------
#  prepare: read input data 
#-------------------------------------------------------------------------------
my @samples;
my %sums;
my %abundance = read_otu_shared($otu_table,\@samples,\%sums);
my %taxonomy  = read_rdp_res($rdp_result);

# used to count the annot otu number of each taxonomy level
my %taxon_count;

# used to save the taxonomy detail info
my %taxon_detail;

#-------------------------------------------------------------------------------
# create OTU profiling file 
#-------------------------------------------------------------------------------
open my$ofh_otu_pro     , ">$outdir/all.otus.profiling.xls" or die $!;
open my$ofh_otu_tab     , ">$outdir/all.otus.tab" or die $!;
open my$ofh_otu_heatmap , ">$outdir/all.otus.exp_for_heatmap" or die $!;
open my$ofh_otu_unifrac , ">$outdir/all.otus.exp_for_unifrac" or die $!;
open my$ofh_otu_FUNGuild , ">$outdir/all.otus.exp_for_FUNGuild" or die $!;

# print the header
my @tags_head     = map { "${_}_tags" } @samples;
my @relative_head = map { "${_}_relative_abundance" } @samples;
my $otu_pro_head  = join "\t" , ("Otu_id","Total_tags",@tags_head,@relative_head,@taxon_levels);
print $ofh_otu_pro "$otu_pro_head\n";

my $otu_tab_head  = join "\t" , ("Otu_id",@samples,"Taxonomy");
print $ofh_otu_tab "$otu_tab_head\n";

my $otu_heatmap_head = join "\t" , ("Otu_id",@samples);
print $ofh_otu_heatmap $otu_heatmap_head , "\n";

my $otu_unifrac_head = join "\t" , ("Otu",@samples);
print $ofh_otu_unifrac $otu_unifrac_head , "\n";

my $otu_total_num = scalar keys %abundance;
foreach my $otuid (sort keys %abundance)
{
	my @taxons     = parse_taxon($taxonomy{$otuid});
	my @tags       = map { $abundance{$otuid}{$_}{tags} } @samples;
	my @relative   = map { $abundance{$otuid}{$_}{relative} } @samples;
	my $total_tags = sum(@tags);
	my $line       = join "\t" , ($otuid , $total_tags , @tags , @relative , @taxons);
	print $ofh_otu_pro $line , "\n";
	
	my $unifrac_line = join "\t" , ($otuid, @tags);
	print $ofh_otu_unifrac $unifrac_line , "\n";

	my $tab_line = join "\t" , ($otuid , @tags , $taxonomy{$otuid});
	print $ofh_otu_tab "$tab_line\n";
	
	my $new_taxon = join ";" , @taxons;
	my $FunGuild_line = join "\t" , ($otuid , @tags, $new_taxon);
	print $ofh_otu_FUNGuild $FunGuild_line , "\n";

	for (@relative)
	{
		if ($_ >= $opts{mohp})
		{
			my $heatmap_line = join "\t" , ($otuid , @relative);
			print $ofh_otu_heatmap "$heatmap_line\n";
			last;
		}
	}

	foreach my $i (0 .. $#taxons)
	{
		my $level = $i+1;
		next unless $taxons[$i];
		
		my $total = 0;
		foreach my $sample (@samples)
		{
			$taxon_count{$taxon_levels[$i]}{$sample}{accum} += $abundance{$otuid}{$sample}{tags};
			$taxon_count{$taxon_levels[$i]}{$sample}{leaf}  += $abundance{$otuid}{$sample}{tags} unless ($taxons[$i+1]);
			
			$taxon_detail{$level}{$taxons[$i]}{$sample}{tags} += $abundance{$otuid}{$sample}{tags};
			$taxon_detail{$level}{$taxons[$i]}{$sample}{relative} += $abundance{$otuid}{$sample}{relative};
			
			$total += $abundance{$otuid}{$sample}{tags};
		}
		$taxon_detail{$level}{$taxons[$i]}{total}{tags} += $total;
	}
}
close $ofh_otu_pro;
close $ofh_otu_tab;
close $ofh_otu_heatmap;
close $ofh_otu_unifrac;
close $ofh_otu_FUNGuild;

#-------------------------------------------------------------------------------
#  create taxonomy annot OTU number stat file
#-------------------------------------------------------------------------------
open my $ofh_taxon_ratio , ">$outdir/all.taxonomy.stat.forFig.xls" or die $!;
print $ofh_taxon_ratio "Level\tSample\tSum\tNumber\n";
foreach my $level (@taxon_levels)
{
	foreach my $sample (@samples)
	{
		my $num   = $taxon_count{$level} && $taxon_count{$level}{$sample} ? $taxon_count{$level}{$sample}{accum} : 0;
		my $leaf  = $taxon_count{$level} && $taxon_count{$level}{$sample} && $taxon_count{$level}{$sample}{leaf} ? $taxon_count{$level}{$sample}{leaf} : 0;
		print $ofh_taxon_ratio "$level\t$sample\t$num\t$leaf\n";
	}
}
close $ofh_taxon_ratio;

open my $ofh_taxon_stat , ">" , "$outdir/all.taxonomy.stat.xls" or die $!;
print $ofh_taxon_stat "\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n";
foreach my $sample (@samples)
{
	my @tags = map { $taxon_count{$_}{$sample}{accum} } @taxon_levels;
	my $line = join "\t" , ($sample,@tags);
	print $ofh_taxon_stat $line , "\n";
}
close $ofh_taxon_stat;

#-------------------------------------------------------------------------------
#  create taxonomy profiling info
#-------------------------------------------------------------------------------
my @taxon_init = map { "" } 0 .. $#taxon_levels;
foreach my $i (1 .. @taxon_levels)
{
	# data file which contain the profiling info of all taxonomy info with each taxonomy level
	my $fname1 = "profiling.L$i.$taxon_levels[$i-1].xls";
	open my $ofh_taxon_pro , ">" , "$outdir/$fname1" or die "can't open $fname1, $!";
	my $header = join "\t" , ("#Level",@taxon_levels,@tags_head,@relative_head);
	print $ofh_taxon_pro $header , "\n";
	
	# data file for draw fill bar fig (whcih contain others and unclassified value)
	my $fname2 = "profiling.L$i.$taxon_levels[$i-1].min$opts{mtfp}.xls";
	my @other   = ();
	my @accumu  = ();
	my @unclass = ();
	open my $ofh_taxon_pro_fig , ">" , "$outdir/$fname2" or die "can't open $fname2, $!";
	$header = join "\t" , ("Taxonomy",@samples);
	print $ofh_taxon_pro_fig $header , "\n";
	
	# data file for draw heatmap, do PCA and cluster
	my $fname3 = "profiling.L$i.$taxon_levels[$i-1].min$opts{mohp}.xls";
	open my $ofh_taxon_pro_pca , ">" , "$outdir/$fname3" or die "can't open $fname3, $!";
	$header = join "\t" , ("Taxonomy",@samples);
	print $ofh_taxon_pro_pca $header , "\n";

	foreach my $name (sort { $taxon_detail{$i}{$b}{total}{tags} <=> $taxon_detail{$i}{$a}{total}{tags} } keys %{$taxon_detail{$i}} )
	{
		$taxon_init[$i-1] = $name;
		my @taxon_tags     = map { $taxon_detail{$i}{$name}{$_} ? $taxon_detail{$i}{$name}{$_}{tags} : 0 } @samples;
		my @taxon_relative = map { $taxon_detail{$i}{$name}{$_} ? sprintf "%.4f" , $taxon_detail{$i}{$name}{$_}{relative} : 0 } @samples;
		@accumu            = map { $accumu[$_] += $taxon_relative[$_] } 0 .. $#samples;

		my $line = join "\t" , ($i,@taxon_init,@taxon_tags,@taxon_relative);
		print $ofh_taxon_pro $line , "\n";
		
		# check the taxon is others or not (at least one sample relative > $opts{mtfp})
		my $flag = check_ratio($opts{mtfp},@taxon_relative);
		if ($flag)
		{
			my $fig_line = join "\t" , ($name,@taxon_relative);
			print $ofh_taxon_pro_fig $fig_line , "\n";
		}
		else 
		{
			@other = map { $other[$_] += $taxon_relative[$_] } 0 .. $#samples;
		}

		# check the taxon is need to be displayed or analysis or not 
		$flag = check_ratio($opts{mohp},@taxon_relative);
		if ($flag)
		{
			my $fig_line = join "\t" , ($name,@taxon_relative);
			print $ofh_taxon_pro_pca $fig_line , "\n";
		}
	}
	
	# print others and unclassified
	@unclass = map { sprintf "%.3f" , abs(100 - $accumu[$_]) } 0 .. $#samples;
	my $other_line   = join "\t" , ("Other",@other);
	my $unclass_line = join "\t" , ("Unclassified",@unclass);
	print $ofh_taxon_pro_fig "$other_line\n$unclass_line\n";

	close $ofh_taxon_pro;
	close $ofh_taxon_pro_fig;
	close $ofh_taxon_pro_pca;
	
	$taxon_init[$i-1] = "";
}

#-------------------------------------------------------------------------------
#  sub function
#-------------------------------------------------------------------------------
sub read_otu_table
{
	my $file = shift ;
	my $samples = shift;
	my %otu_abundance;

	my @sums;
	open IN,$file or die "can't open OTU table file, $file $!";

	my $header = <IN>;
	chomp $header;
	my @temp = split /\t/,$header;
	shift @temp;
	@$samples = @temp;

	while(<IN>)
	{
		chomp;
		my ($otuid,@vals) = split /\t/;

		for my$i ( 0 .. $#vals)
		{
			$sums[$i] += $vals[$i];
		}

		my ($otu_order) = $otuid =~ /Otu(\d+)/;
		$otu_abundance{$otuid}{order} = $otu_order;
		$otu_abundance{$otuid}{abundance} = \@vals;
	}
	close IN;
	
	# calc the relative abundance 
	foreach my $otuid (keys %otu_abundance) 
	{
		my @abundance = @{$otu_abundance{$otuid}{abundance}};
		my @relative  = map { sprintf ( "%.6f" , $abundance[$_]/$sums[$_]*100 ) } 0 .. $#abundance;
		$otu_abundance{$otuid}{relative} = \@relative;
	}

	return %otu_abundance;
}

sub read_otu_shared
{
	my $file    = shift;
	my $samples = shift;
	my $sums    = shift;

	my %otu_abundance;
	my @otus;
	
	my ($label,$flag,$num);
	
	open IN,$file or die "can't open all shared file, $file, $!";
	while(<IN>)
	{
		chomp;

		if (1 == $.)
		{
			($label,$flag,$num,@otus) = split /\t/;
		}
		else 
		{
			my ($label,$sample,$num,@vals) = split /\t/;
			push @$samples , $sample;
			
			my $sum = sum(@vals);
			$sums->{$sample} = $sum;

			for my $i (0 .. $#otus)
			{
				my $otuid = $otus[$i];
				my $relative = sprintf( "%.6f" , $vals[$i]/$sum*100 );
				$otu_abundance{$otuid}{$sample}{tags} = $vals[$i];
				$otu_abundance{$otuid}{$sample}{relative} = $relative;
			}
		}
	}
	close IN;
	
	return %otu_abundance;
}

sub read_rdp_res
{
	my $file = shift;
	my %taxon;
	
	open IN,$file or die "can't open RDP Classifier result file, $file";
	while(<IN>)
	{
		my ($otuid,$taxon) = (split /\t/,$_)[0,1];
		$taxon{$otuid} = $taxon;
	}

	return %taxon;
	close IN;
}

sub parse_taxon
{
	my $taxon = shift;
	
	# init the taxons array
	my @last_taxons = map { "" } 0 .. 6;
	
	return @last_taxons unless $taxon;

	my @taxons = map { $_ = substr($_,3); $_} split /;/,$taxon;
	shift @taxons; # skip the Root
	for my $i (0 .. 6)
	{
		$last_taxons[$i] = $taxons[$i] if ($taxons[$i]);
	}

	return @last_taxons;
}

sub check_ratio
{
	my $cutoff = shift;
	for (@_)
	{
		return 1 if ($_ >= $cutoff);
	}
	return 0;
}

sub usage 
{
	print <<HELP;
Usage:   perl $0 <all.otus.shared> <rdp assign result> <outdir>

Options: --mtfp <FLOAT>   the minmum value of taxonomy relative abundance allowed to 
                          be displayed in the taxonomy fill bar fig, default is 2, means 2%
         --mohp <FLOAT>   the minmum value of OTU relative abundance allowed to be 
                          displayed in the otu heatmap fig, default is 0.1, means 0.1%
HELP
	exit;
}
