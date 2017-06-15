#!usr/bin/perl -w
use strict;
use Cwd 'abs_path';

#-------------------------------------------------------------------------------
#  update log 
#-------------------------------------------------------------------------------
#  version: v6
#-------------------------------------------------------------------------------
#  date:    2016-05-26
#  Author:  aipeng
#  update:  replace the report create script and fix some bugs, details is in the UPDATE.md
#-------------------------------------------------------------------------------
#  version: v6.1
#  date:    2016-09-07
#  Author:  aipeng
#  update:  1. replace the Miseq 2500 with Hiseq 2500
#           2. modify the rawdata filter options in the report 
#           3. add the description of chao1, ACE and coverage in the alpha diversity 
#           4. add the log2FC value as the last field in the metastat result and add its deacription
#           5. change the GeneNumber bug in the figure of Pathway enrichment result display (to KoNumber)
#-------------------------------------------------------------------------------
#  version: v6.1.1
#  date:    2017-01-19
#  Author:  aipeng
#  update:  1. add NSTI score to PICRUSt result
#-------------------------------------------------------------------------------

our $version = 6.1.1;

&usage unless @ARGV == 2;

# set the important vars
my $bin_dir="/Bio/Bin/pipe/meta/16s";
my $komap="/Bio/Database/Database/kegg/latest_kegg/map_class/meta_old.tab";

# the input 
my $config_file=shift @ARGV;
my $output_dir=shift @ARGV;

$output_dir=abs_path($output_dir);
my $group_config=abs_path($config_file);
`mkdir -p $output_dir` unless -d $output_dir;

#-------------------------------------------------------------------------------
# read the config file by function 
my %conf = read_conf($config_file);

#-------------------------------------------------------------------------------
# read the config file 
my $SequenceTitle;
my $SequenceRegion;
my $FixTagNum=-1;
my $CutTagNum=-1;
my $Phread;
my $FilterOptions;
my $OverlapOptions;
my $TaxonamyFasta;
my $TaxonamyTxt;
my $AlignFasta;
my $TaxonamyOptions;
my $OtuLabel;
my @SamplesOrder;
my $SamplesOrderLst;
my $len_min=0;
my $len_max=0;
my $project;
my %samples;
my $Group="no";
my %fq1;
my %fq2;
my $type;
my $ExpPercentCutOff=0.01;  #mean  0.01%
my $StackCutOff=0.01;      #mean  1%
my $lefse="F";
my $metastat="F";
my $cpu = 4;

open(IN,$config_file)||die"cannot open:$!";
while(<IN>)
{
	chomp;
	next if(/^#/);
	if (/^Group_diff_Lefse\s*:/)
	{
		$lefse="T";
	}
	if (/^Group_diff_metastat\s*:/)
	{
		$metastat="T";
	}
	if (/^komap\s*=\s(.*)/)
	{
		$komap=$1;
	}
	if (/^Group_diff_.*?:.*/)
	{
		$Group="yes";
	}
	if (/^SequenceTitle\s*=\s*(.*)/)
	{
		$type=$1;
	}
	if(/^SequenceTitle\s*=\s*(.*)/)
	{
		$SequenceTitle=$1;
		print "$SequenceTitle\n";
	}
	if(/^Project\s*=\s*(.*)/)
	{
		$project=$1;
		print "$project\n";
	}
	if(/^SequenceRegion\s*=\s*(.*)/)
	{
		$SequenceRegion=$1;
		print "$SequenceRegion\n";
	}
	if(/^FixTagNum\s*=\s*(.*)/)
	{
		$FixTagNum=$1;
		print "$FixTagNum\n";
	}
	if(/^CutTagNum\s*=\s*(.*)/)
	{
		$CutTagNum=$1;
		print "$CutTagNum\n";
	}
	if(/^LengthRange\s*=\s*(.*)/)
	{
		my @a=split(/-/,$1);
		$len_min=$a[0];
		$len_max=$a[1];
	}
	if(/^Phread\s*=\s*(.*)/)
	{
		$Phread=$1;
		print "$Phread\n";
	}
	if(/^FilterOptions\s*=\s*(.*)/)
	{
		$FilterOptions=$1;
		print "$FilterOptions\n";
	}
	if(/^OverlapOptions\s*=\s*(.*)/)
	{
		$OverlapOptions=$1;
		print "$OverlapOptions\n";
	}
	if(/^TaxonamyFasta\s*=\s*(.*)/)
	{
		$TaxonamyFasta=$1;
		print "$TaxonamyFasta\n";
	}
	if(/^TaxonamyTxt\s*=\s*(.*)/)
	{
		$TaxonamyTxt=$1;
		print "$TaxonamyTxt\n";
	}
	if(/^AlignFasta\s*=\s*(.*)/)
	{
		$AlignFasta=$1;
		print "$AlignFasta\n";
	}
	if(/^TaxonamyOptions\s*=\s*(.*)/)
	{
		$TaxonamyOptions=$1;
		print "$TaxonamyOptions\n";
	}
	if(/^OtuLabel\s*=\s*(.*)/)
	{
		$OtuLabel=$1;
		print "$OtuLabel\n";
	}
	if(/^ExpPercentCutOff\s*=\s*(.*)/)
	{
		$ExpPercentCutOff=$1;
		print "$ExpPercentCutOff\n";
	}
	if(/^StackCutOff\s*=\s*(.*)/)
	{
		$StackCutOff=$1;
		print "$StackCutOff\n";
	}
	if(/^SamplesOrder\s*=\s*(.*)/)
	{
		my $tmp=$1;
		$SamplesOrderLst=$tmp;
		@SamplesOrder=split(/,/,$tmp);
	}
	if(/^Rawdata\s*=\s*(.*)/)
	{
		my @rawdata=`ls $1/*1.fq`;
		my $dir=$1;
		for my $nname(@rawdata)
		{
			$nname=`basename $nname`;
			chomp $nname;
			$nname=~s/_1\.fq//;
			$samples{$nname}++;
			$fq1{$nname}="$dir/$nname\_1.fq";
			$fq2{$nname}="$dir/$nname\_2.fq";
			print "$nname\n$fq1{$nname}\n$fq2{$nname}\n";
		}
	}
}

if(! defined $SequenceTitle || ! defined $SequenceRegion || ! defined $OtuLabel || ! defined $SamplesOrderLst || ! defined $TaxonamyFasta || ! defined $TaxonamyTxt || ! defined $AlignFasta || ! defined $Phread ||! defined $SequenceRegion)
{
	print STDERR "some important parameter not found\nplease check again\n";
	exit;
}

if($FixTagNum>0 && $CutTagNum>0)
{
	print STDERR "FixTagNum and CutTagNum cannot existance in same time\n";
	exit
}

if(scalar @SamplesOrder != scalar keys %samples)
{
	print STDERR "Fatal Error: your config file have different sample count between SamplesOrder and samples , please check again\n";
	for (@SamplesOrder)
	{
		print "$_\n" unless defined$samples{$_};
	}
	exit;
}

foreach my $out(sort @SamplesOrder)
{
	if(!exists $samples{$out})
	{
		print STDERR "Fatal Error: your config file have wrongs samples names , please check again\n";
		exit;
	}
	
	if($samples{$out}>1)
	{
		print STDERR "Fatal Error: your config file has same name samples, sampls:$out show $samples{$out} times ,please check again\n";
		exit;
	}
	
	if($fq1{$out} eq $fq2{$out})
	{
		print STDERR "Fatal Error: your config file sample:$out  fq1 and fq2 were same, please check again\n";
		exit;
	}
}
close IN;
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# main 

`mkdir -p $output_dir/shell` unless -d "$output_dir/shell";

#### step1  filter_fq
`mkdir -p $output_dir/1.filter` unless -d "$output_dir/1.filter";

#my $filter_fq_software;
open(OUT,">$output_dir/shell/1.filter_fq.sh")||die"cannot open:$!";
foreach my $out(@SamplesOrder)
{
	print OUT "perl $bin_dir/software/filter_fq/filter_fq_quality-1.06.pl $fq1{$out} $fq2{$out} $output_dir/1.filter/ $out $Phread $FilterOptions\n";
}
close OUT;

#### step2  reads_overlap  and   filter length
`mkdir -p $output_dir/2.overlap` unless -d "$output_dir/2.overlap";
open(OUT,">$output_dir/shell/2.overlap.sh")||die"cannot open:$!";
my @tag_stat_lst;
my @uniq_stat_lst;
foreach my $out(@SamplesOrder)
{
	print OUT "$bin_dir/software/reads_overlap/find_reads_overlap_modify_by_chgao -a $output_dir/1.filter/$out\_1.fq -b $output_dir/1.filter/$out\_2.fq -c $output_dir/2.overlap/$out.left -d $output_dir/2.overlap/$out.right -o $output_dir/2.overlap/$out.overlap -q $output_dir/2.overlap/$out.qual $OverlapOptions > $output_dir/2.overlap/$out.log.o 2> $output_dir/2.overlap/$out.log.e\n";
	print OUT "perl $bin_dir/bin/rename_and_filter.pl $output_dir/2.overlap/$out.overlap $output_dir/2.overlap/ $out $len_min $len_max $FixTagNum $CutTagNum\n";
	print OUT "$bin_dir/software/mothur/mothur/mothur \"#unique.seqs(fasta=$output_dir/2.overlap/$out.fasta);\"\n";
	
	push @tag_stat_lst,"$output_dir/2.overlap/$out.fasta";
	push @uniq_stat_lst,"$output_dir/2.overlap/$out.unique.fasta";
}

print OUT "perl $bin_dir/bin/tags_stat.pl Tags ",(join",",@tag_stat_lst)," > $output_dir/2.overlap/stat_tag.xls\n";
print OUT "perl $bin_dir/bin/tags_stat.pl \"Unique Tags\" ",(join",",@uniq_stat_lst)," > $output_dir/2.overlap/stat_uniqTag.xls\n";
close OUT;

#### step3 cluster analysis
`mkdir -p $output_dir/3.cluster` unless -d "$output_dir/3.cluster";
open(OUT,">$output_dir/shell/3.cluster.sh")||die"cannot open:$!";
print OUT "PYTHONPATH_old=\$PYTHONPATH_old\n";
print OUT "export PYTHONPATH=$bin_dir/software/python/lib/python2.7/site-packages/\n";
print OUT "export RDP_JAR_PATH=$bin_dir/software/rdp_classifier/rdp_classifier_2.2/rdp_classifier-2.2.jar\n";
print OUT "cat";
foreach my $out(@SamplesOrder)
{
	print OUT " $output_dir/2.overlap/$out.fasta";
}
print OUT " > $output_dir/3.cluster/all.fasta\n";
print OUT "cat";
foreach my $out(@SamplesOrder)
{
	print OUT " $output_dir/2.overlap/$out.groups";
}
print OUT " > $output_dir/3.cluster/all.groups\n";
print OUT "$bin_dir/software/mothur/mothur/mothur \"#unique.seqs(fasta=$output_dir/3.cluster/all.fasta);count.seqs(name=$output_dir/3.cluster/all.names,group=$output_dir/3.cluster/all.groups);align.seqs(fasta=$output_dir/3.cluster/all.unique.fasta,reference=$AlignFasta,flip=T,processors=8);summary.seqs(fasta=$output_dir/3.cluster/all.unique.align,count=$output_dir/3.cluster/all.count_table)\"\n";

print OUT "python /Bio/bin/assign_taxonomy.py -i $output_dir/3.cluster/all.unique.fasta -t $TaxonamyTxt -r $TaxonamyFasta -o $output_dir/3.cluster/ $TaxonamyOptions\n";
print OUT "perl $bin_dir/bin/result_change.pl $output_dir/3.cluster/all.unique_tax_assignments.txt > $output_dir/3.cluster/all.unique.tax\n";
print OUT "export PYTHONPATH=\$PYTHONPATH_old\n";
my $diffs=int((($len_min+$len_max)/2)*($OtuLabel-0.01));
#my $diffs=0;
#print OUT "$bin_dir/software/mothur/mothur/mothur \"#filter.seqs(fasta=$output_dir/3.cluster/all.unique.align);pre.cluster(fasta=$output_dir/3.cluster/all.unique.filter.fasta, count=$output_dir/3.cluster/all.count_table, diffs=$diffs)\"\n";
print OUT "$bin_dir/software/mothur/mothur/mothur \"#filter.seqs(fasta=$output_dir/3.cluster/all.unique.align);pre.cluster(fasta=$output_dir/3.cluster/all.unique.filter.fasta, name=$output_dir/3.cluster/all.names, diffs=$diffs, processors=8);count.seqs(name=$output_dir/3.cluster/all.unique.filter.precluster.names,group=$output_dir/3.cluster/all.groups)\"\n";
my $dist_cutoff=$OtuLabel+0.03;
print OUT "$bin_dir/software/mothur/mothur/mothur \"#dist.seqs(fasta=$output_dir/3.cluster/all.unique.filter.precluster.fasta,cutoff=$dist_cutoff,processors=8);cluster(column=$output_dir/3.cluster/all.unique.filter.precluster.dist,count=$output_dir/3.cluster/all.unique.filter.precluster.count_table,method=furthest,cutoff=$OtuLabel)\"\n";
print OUT "$bin_dir/software/mothur/mothur/mothur \"#make.shared(list=$output_dir/3.cluster/all.unique.filter.precluster.fn.unique_list.list, count=$output_dir/3.cluster/all.unique.filter.precluster.count_table,label=$OtuLabel);rarefaction.single(shared=$output_dir/3.cluster/all.unique.filter.precluster.fn.unique_list.shared,processors=8);rarefaction.single(shared=$output_dir/3.cluster/all.unique.filter.precluster.fn.unique_list.shared,calc=shannon,processors=8);summary.single(shared=$output_dir/3.cluster/all.unique.filter.precluster.fn.unique_list.shared,calc=sobs-chao-ace-jack-shannon-npshannon-simpson-coverage)\"\n";
if(@SamplesOrder==1)
{
	print OUT "sed -i 's/\\t\\t/\\t$SamplesOrder[0]\\t/' $output_dir/3.cluster/all.unique.filter.precluster.fn.unique_list.groups.summary\n";
	print OUT "sed -i 's/-\\t/-$SamplesOrder[0]\\t/g' $output_dir/3.cluster/all.unique.filter.precluster.fn.unique_list.groups.r_shannon\n";
	print OUT "sed -i 's/-\\t/-$SamplesOrder[0]\\t/g' $output_dir/3.cluster/all.unique.filter.precluster.fn.unique_list.groups.rarefaction\n";
}
close OUT;

#### step4 stat
`mkdir -p $output_dir/4.stat` unless -d "$output_dir/4.stat";
open(OUT,">$output_dir/shell/4.stat.sh")||die"cannot open:$!";
`mkdir -p $output_dir/4.stat/tag_len_distribution` unless -d "$output_dir/4.stat/tag_len_distribution";
`mkdir -p $output_dir/4.stat/tag_abd_distribution` unless -d "$output_dir/4.stat/tag_abd_distribution";
`mkdir -p $output_dir/4.stat/otu_abd_distribution` unless -d "$output_dir/4.stat/otu_abd_distribution";
`mkdir -p $output_dir/4.stat/otu_rarefaction` unless -d "$output_dir/4.stat/otu_rarefaction";
`mkdir -p $output_dir/4.stat/shannon_rarefaction` unless -d "$output_dir/4.stat/shannon_rarefaction";
`mkdir -p $output_dir/4.stat/otu_alpha_dirvisty` unless -d "$output_dir/4.stat/otu_alpha_dirvisty";
`mkdir -p $output_dir/4.stat/taxa_stat` unless -d "$output_dir/4.stat/taxa_stat";
`mkdir -p $output_dir/4.stat/otu_pca` unless -d "$output_dir/4.stat/otu_pca";
`mkdir -p $output_dir/4.stat/taxa_pca` unless -d "$output_dir/4.stat/taxa_pca";
`mkdir -p $output_dir/4.stat/otu_beta_diversity` unless -d "$output_dir/4.stat/otu_beta_diversity";
`mkdir -p $output_dir/4.stat/taxa_beta_diversity` unless -d "$output_dir/4.stat/taxa_beta_diversity";
`mkdir -p $output_dir/4.stat/otu_bray` unless -d "$output_dir/4.stat/otu_bray";
`mkdir -p $output_dir/4.stat/taxa_bray` unless -d "$output_dir/4.stat/taxa_bray";
`mkdir -p $output_dir/4.stat/otu_heatmap` unless -d "$output_dir/4.stat/otu_heatmap";
`mkdir -p $output_dir/4.stat/taxa_heatmap` unless -d "$output_dir/4.stat/taxa_heatmap";
`mkdir -p $output_dir/4.stat/otu_sample_cluster` unless -d "$output_dir/4.stat/otu_sample_cluster";
`mkdir -p $output_dir/4.stat/taxa_sample_cluster` unless -d "$output_dir/4.stat/taxa_sample_cluster";
`mkdir -p $output_dir/4.stat/otu_rank_abd_distribution` unless -d "$output_dir/4.stat/otu_rank_abd_distribution";
`mkdir -p $output_dir/4.stat/otu_tree_stack` unless -d "$output_dir/4.stat/otu_tree_stack";

my $sample_num=scalar @SamplesOrder;

open(OUT1,">$output_dir/4.stat/tag_len_distribution/tag.list")||die"cannot open:$!";
foreach my $out(@SamplesOrder)
{
	print OUT "ln -s -f -t $output_dir/4.stat/tag_len_distribution $output_dir/2.overlap/$out.fasta \n";
	print OUT1 "$output_dir/4.stat/tag_len_distribution/$out.fasta\n";
	print OUT "ln -s -f -t $output_dir/4.stat/tag_abd_distribution $output_dir/2.overlap/$out.names \n";
}
close OUT1;
#print OUT "perl $bin_dir/bin/muti_sample_len_dist.pl $output_dir/4.stat/tag_len_distribution/*\n";
print OUT "\n### tag length distribution\n";
print OUT "perl $bin_dir/bin/customizeLenToPlot.pl $output_dir/4.stat/tag_len_distribution/tag.list $sample_num $output_dir/4.stat/tag_len_distribution/tagLenDistr\n";
print OUT "perl $bin_dir/bin/taglen_v1.pl $output_dir/4.stat/tag_len_distribution/tagLenDistr $output_dir/4.stat/tag_len_distribution/tagLenDistr.forPlot 0.1\n";
print OUT "perl $bin_dir/bin/lenBarPlot.pl $output_dir/4.stat/tag_len_distribution/tagLenDistr.forPlot $sample_num $SamplesOrderLst $output_dir/4.stat/tag_len_distribution/tagLenDistr\n\n";

#print OUT "perl $bin_dir/bin/";

print OUT "\n### tag abundance distribution\n";
print OUT "perl $bin_dir/bin/hebing_abd_distri.pl $output_dir/4.stat/tag_abd_distribution/*\n";
print OUT "perl $bin_dir/bin/customizeAbdToPlot.pl $output_dir/4.stat/tag_abd_distribution/abd_distri_hunhe $sample_num $output_dir/4.stat/tag_abd_distribution/tagAbdDistr\n";
print OUT "perl $bin_dir/bin/tagAbdBarPlot.pl $output_dir/4.stat/tag_abd_distribution/tagAbdDistr $sample_num $SamplesOrderLst $output_dir/4.stat/tag_abd_distribution/tagAbdDistr\n";

print OUT "\n### otu abundance distribution\n";
print OUT "perl $bin_dir/bin/otu_stat.pl $output_dir/3.cluster/all.unique.filter.precluster.fn.unique_list.shared $SamplesOrderLst $output_dir/4.stat/otu_abd_distribution/\n";
print OUT "perl $bin_dir/bin/customizeAbdToPlot.pl $output_dir/4.stat/otu_abd_distribution/stat_otu_distribution.xls $sample_num $output_dir/4.stat/otu_abd_distribution/otuAbdDistr\n";
print OUT "perl $bin_dir/bin/otuAbdBarPlot.pl $output_dir/4.stat/otu_abd_distribution/otuAbdDistr $sample_num $SamplesOrderLst $output_dir/4.stat/otu_abd_distribution/otuAbdDistr\n";

print OUT "\n### rarefaction\n";
print OUT "perl $bin_dir/bin/rarefactionPrepare.pl $output_dir/3.cluster/all.unique.filter.precluster.fn.unique_list.groups.rarefaction $OtuLabel > $output_dir/4.stat/otu_rarefaction/otu.rarefaction\n";
print OUT "perl $bin_dir/bin/otu_rarefaction_plot.pl $output_dir/4.stat/otu_rarefaction/otu.rarefaction $sample_num $output_dir/4.stat/otu_rarefaction/otu.rarefaction\n";
print OUT "perl $bin_dir/bin/rarefactionPrepare.pl $output_dir/3.cluster/all.unique.filter.precluster.fn.unique_list.groups.r_shannon $OtuLabel > $output_dir/4.stat/shannon_rarefaction/shannon.rarefaction\n";
print OUT "perl $bin_dir/bin/shannon_rarefaction_plot.pl $output_dir/4.stat/shannon_rarefaction/shannon.rarefaction $sample_num $output_dir/4.stat/shannon_rarefaction/shannon.rarefaction\n";

print OUT "\n### OTU_Rank_abd_distribution\n";
print OUT "perl $bin_dir/bin/draw_rank_abund.pl $output_dir/3.cluster/all.unique.filter.precluster.fn.unique_list.shared $output_dir/4.stat/otu_rank_abd_distribution\n";

print OUT "\n### alpha diversity\n";
print OUT "perl $bin_dir/bin/get_alpha_diversity_file.pl $output_dir/3.cluster/all.unique.filter.precluster.fn.unique_list.groups.summary $output_dir/4.stat/otu_alpha_dirvisty\n";


print OUT "\n### taxa stat\n";
print OUT "perl $bin_dir/bin/all_otu_reformat.pl $output_dir/3.cluster/all.unique.filter.precluster.fn.unique_list.list $OtuLabel $output_dir/3.cluster/all.unique.filter.precluster.map $output_dir/3.cluster/all.names $output_dir/3.cluster/all.count_table $output_dir/3.cluster/all.fasta $output_dir/4.stat/taxa_stat/all_sample.reformat $output_dir/4.stat/taxa_stat/otu.representative.fasta\n";
print OUT "perl $bin_dir/bin/Mothur_reformated_find_tax_and_abundance_for_Mothur_Names_file_2.pl $output_dir/3.cluster/all.unique.tax $output_dir/3.cluster/all.unique.filter.precluster.names $output_dir/4.stat/taxa_stat/all_sample.reformat $output_dir/4.stat/taxa_stat/all_sample.reformat.tax\n";
print OUT "perl $bin_dir/bin/Mothur_form_OTU_tax_table_2.11.pl $output_dir/4.stat/taxa_stat/all_sample.reformat.tax $output_dir/4.stat/taxa_stat/all_sample.reformat.tax_out\n";
print OUT "perl $bin_dir/bin/Mothur_form_tax_abund_table.pl $output_dir/4.stat/taxa_stat/all_sample.reformat.tax_out $output_dir/4.stat/taxa_stat/all_sample.reformat.tax_out_table\n";
print OUT "perl $bin_dir/bin/get_otu_profile.pl $output_dir/4.stat/taxa_stat/all_sample.reformat.tax_out $output_dir/3.cluster/all.unique.filter.precluster.fn.unique_list.shared $SamplesOrderLst > $output_dir/4.stat/taxa_stat/otu.profile.xls\n";
print OUT "perl $bin_dir/bin/otu_level_stat.pl $output_dir/4.stat/taxa_stat/otu.profile.xls $SamplesOrderLst $output_dir/4.stat/taxa_stat/\n";
print OUT "perl $bin_dir/bin/wholeTaxTrans.pl $output_dir/4.stat/taxa_stat/level_stat.xls $output_dir/4.stat/taxa_stat/wholeTaxStat\n";
print OUT "perl $bin_dir/bin/wholeTaxStat.pl $output_dir/4.stat/taxa_stat/wholeTaxStat $output_dir/2.overlap/stat_tag.xls\n";
print OUT "perl $bin_dir/bin/get_samples_tax_out.pl $output_dir/4.stat/taxa_stat/otu.profile.xls $output_dir/4.stat/taxa_stat/all_sample.reformat.tax_out $sample_num $output_dir/4.stat/taxa_stat/samples\n";
foreach my $out(@SamplesOrder)
{
	print OUT "perl $bin_dir/bin/Mothur_form_tax_abund_table.pl $output_dir/4.stat/taxa_stat/samples/$out.reformat.tax_out $output_dir/4.stat/taxa_stat/samples/$out.reformat.tax_out_table\n";
#	print OUT "awk 'NR>=4' $output_dir/4.stat/taxa_stat/samples/$out.reformat.tax_out_table|cut -f 1,3|awk -F \";\" '{print NF-1\"\\t\"$_}' > $output_dir/4.stat/taxa_stat/samples/$out.reformat.tax_out_table.out\n";
#	print OUT "perl $bin_dir/bin/Standard_Format.pl $output_dir/4.stat/taxa_stat/samples/$out.reformat.tax_out_table.out $output_dir/4.stat/taxa_stat/samples/$out.temporarymap\n";
#	print OUT "awk '{pritn \"$out\\t\"\$1\"\\t\"\$2}' $output_dir/4.stat/taxa_stat/samples/$out.temporarymap > $output_dir/4.stat/taxa_stat/samples/$out.temporary\n";
}
print OUT "perl $bin_dir/bin/profiling_select.pl $output_dir/4.stat/taxa_stat/samples/ 0 $SamplesOrderLst $StackCutOff $output_dir/4.stat/taxa_stat/samples/*.reformat.tax_out_table\n";
print OUT "ls  $output_dir/4.stat/taxa_stat/samples/*_StackMap|awk '{print \"perl $bin_dir/bin/taxStackMapPlot.pl \"\$1\" \"\$1}'|bash\n";
print OUT "ls $output_dir/4.stat/taxa_stat/samples/*.reformat.tax_out_table > $output_dir/4.stat/taxa_stat/samples/all.list\n";
print OUT "perl $bin_dir/bin/level_tongji.pl $output_dir/4.stat/taxa_stat/samples/all.list > $output_dir/4.stat/taxa_stat/samples/level.tongji\n";

print OUT "perl $bin_dir/bin/taxa_tree_plot.pl 2,",$sample_num+1," $output_dir/4.stat/taxa_stat/best_level $output_dir/2.overlap/stat_tag.xls $output_dir/4.stat/taxa_stat/samples/ $output_dir/4.stat/taxa_stat/samples/for_tree\n";
#print OUT "cat $output_dir/4.stat/taxa_stat/samples/*.temporary > $output_dir/4.stat/taxa_stat/samples/whole.temporaryout\n";
#print OUT "perl $bin_dir/bin/Preproceed.pl $output_dir/4.stat/taxa_stat/samples/whole.temporaryout $output_dir/4.stat/taxa_stat/samples/whole\n";
#print OUT "perl $bin_dir/bin/Phylomap.pl $output_dir/4.stat/taxa_stat/samples/whole $output_dir/4.stat/taxa_stat/samples/whole.phylomap $output_dir/4.stat/taxa_stat/samples/whole.map\n";

my $heatmap_cut_begin=$sample_num+2;
my $heatmap_cut_end=$heatmap_cut_begin+$sample_num-1;

print OUT "##########################tree_stack\n";
print OUT "perl $bin_dir/software/SBV/SBV_circos.pl $output_dir/4.stat/taxa_stat/otu.profile.xls $heatmap_cut_begin,$heatmap_cut_end $output_dir/4.stat/taxa_stat/otu.representative.fasta $output_dir/4.stat/otu_tree_stack/ 0.001\n";

if ($type=~/16s/i)
{
	`mkdir -p $output_dir/4.stat/otu_annot` unless -d "$output_dir/4.stat/otu_annot";
	print OUT "###################################otu annotation\n";
#	print OUT "export PYTHONPATH=/Bio/Bin/pipe/meta/16s/software/python/lib/python2.7/site-packages/\n";
	print OUT "perl $bin_dir/software/python/picrust/picrust-1.0.0/meta_annot.pl $output_dir/3.cluster/all.unique.align.report $output_dir/4.stat/taxa_stat/all_sample.reformat.tax $output_dir/4.stat/taxa_stat/otu.profile.xls $heatmap_cut_begin,$heatmap_cut_end $output_dir/4.stat/otu_annot\n";
	print OUT "mkdir -p $output_dir/4.stat/otu_annot/Sample/\n";
	print OUT "rm $output_dir/4.stat/otu_annot/Sample/*\n";
	print OUT "perl $bin_dir/software/koenrich/meta_annot_ko.pl $output_dir/4.stat/otu_annot/meta_ko_annot.xls $komap $output_dir/4.stat/otu_annot/Sample/\n";
}

if($sample_num>=2)
{
	print OUT "\n### OTU heatmap\n";
	print OUT "perl $bin_dir/bin/get_percent_cutoff_for_heatmap.pl $output_dir/4.stat/taxa_stat/otu.profile.xls $heatmap_cut_begin,$heatmap_cut_end $ExpPercentCutOff > $output_dir/4.stat/otu_heatmap/otu_exp.filter\n";
	print OUT "Rscript $bin_dir/bin/heatmap.R $output_dir/4.stat/otu_heatmap/otu_exp.filter $heatmap_cut_begin,$heatmap_cut_end \n";
	print OUT "\n### taxa heatmap\n";
	print OUT "perl $bin_dir/bin/level_exp_to_dir.pl $output_dir/4.stat/taxa_stat/samples $output_dir/4.stat/taxa_heatmap\n";
	print OUT "ls $output_dir/4.stat/taxa_heatmap/ | awk '{print \"perl $bin_dir/bin/get_percent_cutoff_for_heatmap.pl $output_dir/4.stat/taxa_heatmap/\"\$1\"/\"\$1\"_exp $heatmap_cut_begin,$heatmap_cut_end $ExpPercentCutOff > $output_dir/4.stat/taxa_heatmap/\"\$1\"/\"\$1\"_exp.filter\"}'|bash\n";
	print OUT "ls $output_dir/4.stat/taxa_heatmap/ | awk '{print \"Rscript $bin_dir/bin/heatmap.R $output_dir/4.stat/taxa_heatmap/\"\$1\"/\"\$1\"_exp.filter $heatmap_cut_begin,$heatmap_cut_end\"}'|bash\n";
}

if($sample_num>=3)
{
	print OUT "\n### OTU PCA\n";
	print OUT "ln -s -f -t $output_dir/4.stat/otu_pca $output_dir/4.stat/otu_heatmap/otu_exp.filter\n";
	print OUT "Rscript $bin_dir/bin/pca.R $output_dir/4.stat/otu_pca/otu_exp.filter $heatmap_cut_begin,$heatmap_cut_end $conf{group_list}\n";
	
	print OUT "\n### taxa PCA\n";
	print OUT "perl $bin_dir/bin/level_to_dir.pl $output_dir/4.stat/taxa_stat/samples $output_dir/4.stat/taxa_pca\n";
	print OUT "ls $output_dir/4.stat/taxa_pca/| awk '{print \"ln -s -f -t $output_dir/4.stat/taxa_pca/\"\$1\" $output_dir/4.stat/taxa_heatmap/\"\$1\"/\"\$1\"_exp.filter\"}'|bash\n";
	print OUT "ls $output_dir/4.stat/taxa_pca/ | awk '{print \"Rscript $bin_dir/bin/pca.R $output_dir/4.stat/taxa_pca/\"\$1\"/\"\$1\"_exp.filter $heatmap_cut_begin,$heatmap_cut_end $conf{group_list}\"}'|bash\n";
}

if($sample_num>=2)
{
	print OUT "\n### OTU beta diversity\n";
#	print OUT "cut -f 1,$otu_cut_begin-$otu_cut_end $output_dir/4.stat/taxa_stat/otu.profile.xls > $output_dir/4.stat/otu_beta_diversity/otu_exp\n";
	print OUT "ln -s -f -t $output_dir/4.stat/otu_beta_diversity $output_dir/4.stat/otu_heatmap/otu_exp.filter\n";
	print OUT "perl $bin_dir/bin/beta_flow.pl $output_dir/4.stat/otu_beta_diversity/otu_exp.filter $heatmap_cut_begin,$heatmap_cut_end $output_dir/4.stat/otu_beta_diversity/ 0 $sample_num\n";
	print OUT "\n### taxa beta diversity\n";
	print OUT "perl $bin_dir/bin/level_to_dir.pl $output_dir/4.stat/taxa_stat/samples $output_dir/4.stat/taxa_beta_diversity\n";
#	print OUT "ls $output_dir/4.stat/taxa_beta_diversity/ | awk '{print \"perl $bin_dir/bin/beta_flow.pl $output_dir/4.stat/taxa_beta_diversity/\"\$1\"/\"\$1\"_exp $output_dir/4.stat/taxa_beta_diversity/\"\$1\"/ 0 $sample_num\"}'|bash\n";
	print OUT "ls $output_dir/4.stat/taxa_beta_diversity/| awk '{print \"ln -s -f -t $output_dir/4.stat/taxa_beta_diversity/\"\$1\" $output_dir/4.stat/taxa_heatmap/\"\$1\"/\"\$1\"_exp.filter\"}'|bash\n";
	print OUT "ls $output_dir/4.stat/taxa_beta_diversity/ | awk '{print \"perl $bin_dir/bin/beta_flow.pl $output_dir/4.stat/taxa_beta_diversity/\"\$1\"/\"\$1\"_exp.filter $heatmap_cut_begin,$heatmap_cut_end $output_dir/4.stat/taxa_beta_diversity/\"\$1\"/ 0 $sample_num\"}'|bash\n";
}

if($sample_num>=2)
{
	print OUT "\n### OTU bray\n";
	print OUT "ln -s -f -t $output_dir/4.stat/otu_bray $output_dir/4.stat/otu_beta_diversity/otu_exp.filter.colm\n";
	print OUT "Rscript $bin_dir/bin/bray.R $output_dir/4.stat/otu_bray/otu_exp.filter.colm\n";
	print OUT "\n### taxa bray stack\n";
	print OUT "perl $bin_dir/bin/level_to_dir.pl $output_dir/4.stat/taxa_stat/samples $output_dir/4.stat/taxa_bray\n";
	print OUT "ls $output_dir/4.stat/taxa_bray/| awk '{print \"Rscript $bin_dir/bin/bray_stackmap.R $output_dir/4.stat/taxa_beta_diversity/\"\$1\"/\"\$1\"_exp.filter.colm $output_dir/4.stat/taxa_stat/samples/\"\$1\"_StackMap $output_dir/4.stat/taxa_bray/\"\$1\"/\"\$1 }'|bash\n";
}

if($sample_num>=3)
{
	print OUT "\n### OTU sample cluster\n";
#	print OUT "cut -f 1,$otu_cut_begin-$otu_cut_end $output_dir/4.stat/taxa_stat/otu.profile.xls > $output_dir/4.stat/otu_sample_cluster/otu_exp\n";
	print OUT "ln -s -f -t $output_dir/4.stat/otu_sample_cluster $output_dir/4.stat/otu_heatmap/otu_exp.filter\n";
	print OUT "perl $bin_dir/bin/bootstrap.pl $output_dir/4.stat/otu_sample_cluster/otu_exp.filter $heatmap_cut_begin,$heatmap_cut_end 100 kld_jsd $output_dir/4.stat/otu_sample_cluster/\n";
#	print OUT "R CMD BATCH $output_dir/4.stat/otu_sample_cluster/otu_exp.filter.kld_jsd.R $output_dir/4.stat/otu_sample_cluster/otu_exp.filter.kld_jsd.Rout\n";
	print OUT "\n### taxa sample cluster\n";
	print OUT "perl $bin_dir/bin/level_to_dir.pl $output_dir/4.stat/taxa_stat/samples $output_dir/4.stat/taxa_sample_cluster\n";
	print OUT "ls $output_dir/4.stat/taxa_sample_cluster/| awk '{print \"ln -s -f -t $output_dir/4.stat/taxa_sample_cluster/\"\$1\" $output_dir/4.stat/taxa_heatmap/\"\$1\"/\"\$1\"_exp.filter\"}'|bash\n";
	print OUT "ls $output_dir/4.stat/taxa_sample_cluster/| awk '{print \"perl $bin_dir/bin/bootstrap.pl $output_dir/4.stat/taxa_sample_cluster/\"\$1\"/\"\$1\"_exp.filter $heatmap_cut_begin,$heatmap_cut_end 100 kld_jsd $output_dir/4.stat/taxa_sample_cluster/\"\$1\"/\"}'|bash \n";
#	print OUT "ls $output_dir/4.stat/taxa_sample_cluster/| awk '{print \"perl $bin_dir/bin/bootstrap.pl $output_dir/4.stat/taxa_sample_cluster/\"\$1\"/\"\$1\"_exp.filter $heatmap_cut_begin,$heatmap_cut_end 100 kld_jsd $output_dir/4.stat/taxa_sample_cluster/\"\$1\"/\\nR CMD BATCH $output_dir/4.stat/taxa_sample_cluster/\"\$1\"/\"\$1\"_exp.filter.kld_jsd.R $output_dir/4.stat/taxa_sample_cluster/\"\$1\"/\"\$1\"_exp.filter.kld_jsd.Rout\"}'|bash\n";
}
if ($Group eq "yes")
{
	`mkdir -p $output_dir/4.stat/group_diff` unless -d "$output_dir/4.stat/group_diff";
	print OUT "######################group_diff\n";
	print OUT "perl $bin_dir/software/group_diff/meta_diff.pl -ml -r multi -c 4 -o $output_dir/4.stat/group_diff/ $group_config $output_dir/4.stat/taxa_stat/otu.profile.xls\n";
	#print OUT "perl $bin_dir/software/group_diff/meta_group_diff.pl $group_config $output_dir/4.stat/taxa_stat/otu.profile.xls $heatmap_cut_begin,$heatmap_cut_end $output_dir/4.stat/group_diff/\n";
	if ($type=~/16s/i && $metastat eq "T")
	{
		print OUT "for i in $output_dir/4.stat/group_diff/metastat/*/*.result.xls;do dirname=`dirname \$i`;name=`basename \$i|sed 's/.metastat.result.xls//'`;perl $bin_dir/software/koenrich/metastat_ko.pl \$dirname/\$name.metastat.result.xls \$dirname/\$name.metastat.xls $output_dir/ \$dirname \$name;done\n";
	}
}
close OUT;

#### step5 unload
`mkdir -p $output_dir/upload` unless -d "$output_dir/upload";
`mkdir -p $output_dir/upload/1.tag` unless -d "$output_dir/upload/1.tag";
`mkdir -p $output_dir/upload/1.tag/fasta` unless -d "$output_dir/upload/1.tag/fasta";
`mkdir -p $output_dir/upload/1.tag/stat` unless -d "$output_dir/upload/1.tag/stat";
`mkdir -p $output_dir/upload/2.otu` unless -d "$output_dir/upload/2.otu";
`mkdir -p $output_dir/upload/2.otu/stat` unless -d "$output_dir/upload/2.otu/stat";
`mkdir -p $output_dir/upload/2.otu/expression_profile` unless -d "$output_dir/upload/2.otu/expression_profile";
`mkdir -p $output_dir/upload/2.otu/analysis` unless -d "$output_dir/upload/2.otu/analysis";
`mkdir -p $output_dir/upload/3.taxa` unless -d "$output_dir/upload/3.taxa";
`mkdir -p $output_dir/upload/3.taxa/stat` unless -d "$output_dir/upload/3.taxa/stat";
`mkdir -p $output_dir/upload/3.taxa/expression_profile` unless -d "$output_dir/upload/3.taxa/expression_profile";
`mkdir -p $output_dir/upload/3.taxa/expression_profile/tree_stack` unless -d "$output_dir/upload/3.taxa/expression_profile/tree_stack";
if($sample_num>=2)
{
	`mkdir -p $output_dir/upload/3.taxa/analysis` unless -d "$output_dir/upload/3.taxa/analysis";
}
open(OUT,">$output_dir/shell/5.upload.sh")||die"cannot open:$!";
print OUT "\n#### cp tag\n";
foreach my $out(@SamplesOrder)
{
	print OUT "cp $output_dir/2.overlap/$out.fasta $output_dir/upload/1.tag/fasta/\n";
}
print OUT "cp $output_dir/2.overlap/stat_tag.xls $output_dir/upload/1.tag/stat/\n";
print OUT "cp $output_dir/2.overlap/stat_uniqTag.xls $output_dir/upload/1.tag/stat/\n";
print OUT "cp $output_dir/4.stat/tag_len_distribution/tagLenDistr $output_dir/upload/1.tag/stat/\n";
print OUT "cp $output_dir/4.stat/tag_len_distribution/tagLenDistr.forPlot $output_dir/upload/1.tag/stat/\n";
print OUT "cp $output_dir/4.stat/tag_len_distribution/tagLenDistr.png $output_dir/upload/1.tag/stat/\n";
print OUT "cp $output_dir/4.stat/tag_len_distribution/tagLenDistr.pdf $output_dir/upload/1.tag/stat/\n";
print OUT "cp $output_dir/4.stat/tag_abd_distribution/tagAbdDistr $output_dir/upload/1.tag/stat/\n";
print OUT "cp $output_dir/4.stat/tag_abd_distribution/tagAbdDistr.png $output_dir/upload/1.tag/stat/\n";
print OUT "cp $output_dir/4.stat/tag_abd_distribution/tagAbdDistr.pdf $output_dir/upload/1.tag/stat/\n";

print OUT "\n#### cp otu\n";
print OUT "cp $output_dir/4.stat/otu_abd_distribution/stat_otu.xls $output_dir/upload/2.otu/stat/\n";
print OUT "cp $output_dir/4.stat/otu_abd_distribution/otuAbdDistr $output_dir/upload/2.otu/stat/\n";
print OUT "cp $output_dir/4.stat/otu_abd_distribution/otuAbdDistr.png $output_dir/upload/2.otu/stat/\n";
print OUT "cp $output_dir/4.stat/otu_abd_distribution/otuAbdDistr.pdf $output_dir/upload/2.otu/stat/\n";
print OUT "cp $output_dir/4.stat/taxa_stat/otu.profile.xls $output_dir/upload/2.otu/expression_profile/\n";
print OUT "cp $output_dir/4.stat/taxa_stat/otu.representative.fasta $output_dir/upload/2.otu/expression_profile/\n";
print OUT "cp $output_dir/4.stat/otu_rarefaction/otu.rarefaction $output_dir/upload/2.otu/analysis/\n";
print OUT "cp $output_dir/4.stat/otu_rarefaction/otu.rarefaction.png $output_dir/upload/2.otu/analysis/\n";
print OUT "cp $output_dir/4.stat/otu_rarefaction/otu.rarefaction.pdf $output_dir/upload/2.otu/analysis/\n";
print OUT "cp $output_dir/4.stat/shannon_rarefaction/shannon.rarefaction $output_dir/upload/2.otu/analysis/\n";
print OUT "cp $output_dir/4.stat/shannon_rarefaction/shannon.rarefaction.png $output_dir/upload/2.otu/analysis/\n";
print OUT "cp $output_dir/4.stat/shannon_rarefaction/shannon.rarefaction.pdf $output_dir/upload/2.otu/analysis/\n";
print OUT "cp $output_dir/4.stat/otu_alpha_dirvisty/*.xls $output_dir/upload/2.otu/analysis/\n";
print OUT "cp $output_dir/4.stat/otu_rank_abd_distribution/Rank_abd_Distr.png $output_dir/upload/2.otu/analysis/\n";
print OUT "cp $output_dir/4.stat/otu_rank_abd_distribution/Rank_abd_Distr.pdf $output_dir/upload/2.otu/analysis/\n";


if ($Group eq "yes")
{
	if ($metastat eq "T")
	{
		`mkdir -p $output_dir/upload/4.Group_diff/metastat` unless -d "$output_dir/upload/4.Group_diff/metastat";
		print OUT "rm -rf $output_dir/upload/4.Group_diff/metastat/*\n";
		print OUT "cp -r $output_dir/4.stat/group_diff/metastat/* $output_dir/upload/4.Group_diff/metastat/\n";
		print OUT "rm $output_dir/upload/4.Group_diff/metastat/*/*metastat.xls\n";
		print OUT "mv $output_dir/upload/4.Group_diff/metastat/Domain $output_dir/upload/4.Group_diff/metastat/1.Domain;mv $output_dir/upload/4.Group_diff/metastat/Phylum $output_dir/upload/4.Group_diff/metastat/2.Phylum;mv $output_dir/upload/4.Group_diff/metastat/Class $output_dir/upload/4.Group_diff/metastat/3.Class;mv $output_dir/upload/4.Group_diff/metastat/Order $output_dir/upload/4.Group_diff/metastat/4.Order;mv $output_dir/upload/4.Group_diff/metastat/Family $output_dir/upload/4.Group_diff/metastat/5.Family;mv $output_dir/upload/4.Group_diff/metastat/Genus $output_dir/upload/4.Group_diff/metastat/6.Genus;mv $output_dir/upload/4.Group_diff/metastat/Species $output_dir/upload/4.Group_diff/metastat/7.Species\n";
		print OUT "rm $output_dir/upload/4.Group_diff/metastat/*/*/*annot*\n";
	}
	
	if ($lefse eq "T")
	{
		`mkdir -p $output_dir/upload/4.Group_diff/LefSe` unless -d "$output_dir/upload/4.Group_diff/LefSe/";
		print OUT "cp -r $output_dir/4.stat/group_diff/LefSe/* $output_dir/upload/4.Group_diff/LefSe/\n";
		print OUT <<SH;
for i in $output_dir/upload/4.Group_diff/LefSe/*Lefse.svg
do
	name=`echo \$i|sed 's/.Lefse.svg//'`
	if [ -s \$i ];then
		echo \$i
	else
		rm \$i \$name.Lefse.png \$name.Lefse.cladogram.svg \$name.Lefse.cladogram.png \$name -rf
	fi
done
SH
		print OUT "rm $output_dir/upload/4.Group_diff/LefSe/*.xls $output_dir/upload/4.Group_diff/LefSe/*.in $output_dir/upload/4.Group_diff/LefSe/*.res\n";
	}
}

if ($type=~/16s/i)
{
	`mkdir -p $output_dir/upload/2.otu/annot` unless -d "$output_dir/upload/2.otu/annot";
	`mkdir -p $output_dir/upload/2.otu/annot/Sample/` unless -d "$output_dir/upload/2.otu/annot/Sample/";
	print OUT "cp $output_dir/4.stat/otu_annot/meta_ko_annot.xls $output_dir/upload/2.otu/annot/\n";
	print OUT "cp $output_dir/4.stat/otu_annot/Sample/* $output_dir/upload/2.otu/annot/Sample/\n";
}
if($sample_num>=3)
{
	print OUT "cp $output_dir/4.stat/otu_pca/otu_exp.filter.* $output_dir/upload/2.otu/analysis/\n";
	print OUT "cp $output_dir/4.stat/otu_sample_cluster/otu_exp.filter.kld_jsd.pdf $output_dir/upload/2.otu/analysis/\n";
	print OUT "cp $output_dir/4.stat/otu_sample_cluster/otu_exp.filter.kld_jsd.png $output_dir/upload/2.otu/analysis/\n";
	print OUT "cp $output_dir/4.stat/otu_sample_cluster/otu_exp.filter.kld_jsd.newick $output_dir/upload/2.otu/analysis/\n";
}
if($sample_num>=2)
{
	print OUT "cp $output_dir/4.stat/otu_heatmap/otu_exp.filter $output_dir/upload/2.otu/analysis/\n";
	print OUT "cp $output_dir/4.stat/otu_heatmap/otu_exp.filter.heatmap.* $output_dir/upload/2.otu/analysis/\n";
	print OUT "cp $output_dir/4.stat/otu_beta_diversity/otu_exp.filter.matrix $output_dir/upload/2.otu/analysis/\n";
	print OUT "cp $output_dir/4.stat/otu_beta_diversity/otu_exp.filter.betadiversity.pdf $output_dir/upload/2.otu/analysis/\n";
	print OUT "cp $output_dir/4.stat/otu_beta_diversity/otu_exp.filter.betadiversity.png $output_dir/upload/2.otu/analysis/\n";
	print OUT "cp $output_dir/4.stat/otu_bray/otu_exp.filter.colm.bray.pdf $output_dir/upload/2.otu/analysis/\n";
	print OUT "cp $output_dir/4.stat/otu_bray/otu_exp.filter.colm.bray.png $output_dir/upload/2.otu/analysis/\n";
}

print OUT "\n#### cp taxa\n";
print OUT "cp $output_dir/4.stat/taxa_stat/wholeTaxStat $output_dir/upload/3.taxa/stat/\n";
print OUT "cp $output_dir/4.stat/taxa_stat/wholeTaxStat.*.png $output_dir/upload/3.taxa/stat/\n";
print OUT "cp $output_dir/4.stat/taxa_stat/wholeTaxStat.*.pdf $output_dir/upload/3.taxa/stat/\n";
print OUT "cp $output_dir/4.stat/taxa_stat/best_level $output_dir/upload/3.taxa/stat/\n";
print OUT "cp $output_dir/4.stat/taxa_stat/best_level_stat.xls $output_dir/upload/3.taxa/stat/\n";
print OUT "cp $output_dir/4.stat/taxa_stat/samples/*.xls $output_dir/upload/3.taxa/expression_profile/\n";
print OUT "cp $output_dir/4.stat/taxa_stat/samples/*_StackMap $output_dir/upload/3.taxa/expression_profile/\n";
print OUT "cp $output_dir/4.stat/taxa_stat/samples/*_StackMap.pdf $output_dir/upload/3.taxa/expression_profile/\n";
print OUT "cp $output_dir/4.stat/taxa_stat/samples/*_StackMap.png $output_dir/upload/3.taxa/expression_profile/\n";
print OUT "cp $output_dir/4.stat/taxa_stat/samples/for_tree/taxtree.svg $output_dir/upload/3.taxa/expression_profile/\n";
print OUT "cp $output_dir/4.stat/taxa_stat/samples/for_tree/taxtree.png $output_dir/upload/3.taxa/expression_profile/\n";
print OUT "cp $output_dir/4.stat/otu_tree_stack/*/*legend* $output_dir/upload/3.taxa/expression_profile/tree_stack/\n";
print OUT "cp $output_dir/4.stat/otu_tree_stack/*/*range $output_dir/upload/3.taxa/expression_profile/tree_stack/\n";
print OUT "cp $output_dir/4.stat/otu_tree_stack/*/*bar_color $output_dir/upload/3.taxa/expression_profile/tree_stack/\n";
print OUT "cp $output_dir/4.stat/otu_tree_stack/*/*.fa $output_dir/upload/3.taxa/expression_profile/tree_stack/\n";
print OUT "rm -rf $output_dir/upload/3.taxa/analysis/*\n" if -d "$output_dir/upload/3.taxa/analysis";
if($sample_num>=3)
{
#	print OUT "rm $output_dir/upload/3.taxa/analysis/*/*.filter\n";
	print OUT "cp -a $output_dir/4.stat/taxa_pca/* $output_dir/upload/3.taxa/analysis/\n";
	print OUT "cp -a $output_dir/4.stat/taxa_pca/* $output_dir/upload/3.taxa/analysis/\n";
	print OUT "rm $output_dir/upload/3.taxa/analysis/*/*.filter\n";
	print OUT "cp -a $output_dir/4.stat/taxa_sample_cluster/* $output_dir/upload/3.taxa/analysis/\n";
	print OUT "rm $output_dir/upload/3.taxa/analysis/*/*.filter\n";
	print OUT "rm $output_dir/upload/3.taxa/analysis/*/*.colm\n";
}
if($sample_num>=2)
{
#	print OUT "rm $output_dir/upload/3.taxa/analysis/*/*.filter\n";
	print OUT "cp -a $output_dir/4.stat/taxa_beta_diversity/* $output_dir/upload/3.taxa/analysis/\n";
	print OUT "rm $output_dir/upload/3.taxa/analysis/*/*.filter\n";
	print OUT "cp -a $output_dir/4.stat/taxa_bray/* $output_dir/upload/3.taxa/analysis/\n";
	print OUT "cp -a $output_dir/4.stat/taxa_heatmap/* $output_dir/upload/3.taxa/analysis/\n";
	print OUT "rm $output_dir/upload/3.taxa/analysis/*/*.R $output_dir/upload/3.taxa/analysis/*/*.Rout $output_dir/upload/3.taxa/analysis/*/*.colm\n";
	print OUT "cp -a $output_dir/upload/3.taxa/analysis/Domain $output_dir/upload/3.taxa/analysis/1.Domain;cp -a $output_dir/upload/3.taxa/analysis/Phylum $output_dir/upload/3.taxa/analysis/2.phylum;cp -a $output_dir/upload/3.taxa/analysis/Class $output_dir/upload/3.taxa/analysis/3.Class;cp -a $output_dir/upload/3.taxa/analysis/Order $output_dir/upload/3.taxa/analysis/4.Order;cp -a $output_dir/upload/3.taxa/analysis/Family $output_dir/upload/3.taxa/analysis/5.Family;cp -a $output_dir/upload/3.taxa/analysis/Genus $output_dir/upload/3.taxa/analysis/6.Genus;cp -a $output_dir/upload/3.taxa/analysis/Species $output_dir/upload/3.taxa/analysis/7.Species\n";
	print OUT "rm -rf $output_dir/upload/3.taxa/analysis/Domain $output_dir/upload/3.taxa/analysis/Phylum $output_dir/upload/3.taxa/analysis/Class $output_dir/upload/3.taxa/analysis/Order $output_dir/upload/3.taxa/analysis/Species $output_dir/upload/3.taxa/analysis/Family  $output_dir/upload/3.taxa/analysis/Genus\n";
}
close OUT;

# create the report 
$SequenceTitle =~ s/\s+/_/;
my $out_pack_name = "${project}_${SequenceTitle}_result";
open(OUT,">$output_dir/shell/6.report.sh")||die"cannot open:$!";
print OUT "perl /Bio/User/aipeng/bin/meta16S/meta16S_report.pl -n $out_pack_name $group_config $output_dir/upload\n";
close OUT;

# crate the total shell 
open(OUT,">$output_dir/shell/total.sh")||die"cannot open:$!";
print OUT "#! /bin/bash\n#\$ -S /bin/bash\n";
print OUT "perl $bin_dir/bin/qsub-sge.pl --maxjob $cpu --convert no $output_dir/shell/1.filter_fq.sh\n";
print OUT "sh $output_dir/shell/2.overlap.sh\n";
print OUT "sh $output_dir/shell/3.cluster.sh\n";
print OUT "sh $output_dir/shell/4.stat.sh\n";
print OUT "sh $output_dir/shell/5.upload.sh\n";
print OUT "sh $output_dir/shell/6.report.sh\n";
close OUT;

#===============================================================================
#  sub function 
#-------------------------------------------------------------------------------
# print the help info.
#-------------------------------------------------------------------------------
sub usage 
{
	print <<HELP;

Intro: A pipeline for create shell for Meta 16S / ITS 

Usage: perl $0 <config file> <outdir>

Version: $version

HELP
	exit;
}

#-------------------------------------------------------------------------------
#  A function to read config file 
#-------------------------------------------------------------------------------
sub read_conf
{
	my $conf = shift;
	
	my %groups;
	my %conf;

	open IN,$conf or die $!;
	
	while(<IN>)
	{
		chomp;
		next if /^#/;
		next if ($_ eq "");
		
		# group define 
		if (/Group:\s*(\S+)\s*=\s*(\S+)/)
		{
			my @tmp = split /,/,$2;
			map { $groups{$_} = $1 } @tmp;
		}
		# normal pattern
		elsif (/(\S+)\s*[:=]\s*(\S+)/)
		{
			$conf{$1} = $2;
		}
		else 
		{
			die "Your config can't be recognized, [$_]";
		}
	}
	
	close IN;
	
	my @samples = split /,/ , $conf{SamplesOrder};
	my @groups = map { $groups{$_} } @samples;
	my $group_list = join "," , @groups;
	
	$conf{group_list} = $group_list;
	$conf{groups} = \@groups;
	$conf{samples} = \@samples;
	
	return %conf;
}
