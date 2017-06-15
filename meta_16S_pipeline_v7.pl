#!/Bio/bin/perl
#-----------------------------------------------+
#    [APM] This script was created by amp.pl    |
#    [APM] Created time: 2017-02-10 14:50:54    |
#-----------------------------------------------+
=pod

=head2 v7.0

Date: 2017-02-10 14:50:54

=head1 Name

meta_16S_pipeline_v7.pl

=head1 Synopsis



=head1 Feedback

Author: Peng Ai
Email:  pai@genedenovo.com

=head1 Version

Version history

=cut


use strict;
use warnings;

use Getopt::Long;
use File::Basename qw(basename dirname);
use FindBin qw($Bin);
use Cwd;
use List::Util qw/min sum/;

use lib "$Bin/lib";
use DEBUG;
use CONF;

our %OPTS = (run=>"no");
GetOptions(\%OPTS,
	'conf:s',
	'run:s',
	'steps:s',
	'help'
);

&usage if ($OPTS{help});

if (! $OPTS{conf})
{
	WARN("the conf file is not defined!");
	&usage;
}

my $conf_file = check_path($OPTS{conf});
my $conf = load_conf($OPTS{conf},-check=>1,-init=>1);
#-------------------------------------------------------------------------------
#  main 
#-------------------------------------------------------------------------------
# init some global vars
our @samples = split /[;,\s\t]/ , $conf->{samples};
our $sample_num = scalar @samples;
our @group_names = fetch_group_names($conf);
our %groups = fetch_groups($conf);
our ($min_len,$max_len) = split /-/ , $conf->{length_range};
our $split_line = '-' x 60;
our $src = "$Bin/src";
our $target = $conf->{sequence_target};

# load some path 
my $python     = check_path($conf->{soft}->{python});
my $python_lib = check_path($conf->{soft}->{python_lib});
my $Rscript    = check_path($conf->{soft}->{rscript});
my $sbv        = check_path($conf->{soft}->{sbv});
my $mothur     = check_path($conf->{soft}->{mothur});
my $qiime      = check_path($conf->{soft}->{qiime});

# create directory 
&create_dir($conf);
our $dirs = $conf->{dirs};

# run 
&run($conf);

#-------------------------------------------------------------------------------
#  sub functions
#-------------------------------------------------------------------------------
# create the directory
sub create_dir
{
	my $conf = shift;
	my $odir = $conf->{out_dir};
	
	my $upload_name = "$conf->{project_id}_$conf->{sequence_target}_rDNA_result";
	my $brief_name  = "$conf->{project_id}_$conf->{sequence_target}_rDNA_result_brief";

	# init dirs 
	$conf->{dirs} = {
		odir     => $odir,
		shell    => "$odir/shell",
		upload   => "$odir/$upload_name",
		brief    => "$odir/$brief_name",
		data     => "$odir/01.data_process",
		otus     => "$odir/02.otus",
		taxonomy => "$odir/03.taxonomy",
		alpha    => "$odir/04.alpha_diversity",
		beta     => "$odir/05.beta_diversity",
		function => "$odir/06.function_annot",
		compare  => "$odir/07.different_analysis"
	};

	my @order = qw /odir shell upload data otus taxonomy alpha beta function compare/;
	map {
		mkdir $conf->{dirs}->{$_} unless -d $conf->{dirs}->{$_};
		$conf->{dirs}->{$_} = check_path($conf->{dirs}->{$_})
	} @order;
	
	timeLOG($split_line);
	timeLOG("start the Meta $conf->{sequence_target} rDNA pipeline ...");
	timeLOG("create the directory done, main out dir is [$conf->{dirs}->{odir}] ... ");
}

sub run 
{
	my $conf = shift;
	
	&data_process($conf);
	&otus($conf);
	&taxonomy($conf);
	&alpha_diversity($conf);
	&beta_diversity($conf);
	&function_annot($conf);
	&different_analysis($conf);
	&pack_result($conf);
}

#-------------------------------------------------------------------------------
#  step1: data process 
#-------------------------------------------------------------------------------
#  1. raw data Quality Control, with Filter_fq
#  2. Overlap the pair reads, with FLASH
#  3. Tags Quality Control
#  4. remove the Chimera Sequence, with Usearch and GOLD databsse
#-------------------------------------------------------------------------------
sub data_process
{
	timeLOG($split_line);
	timeLOG("**** 01: data process ****");

	my $conf = shift;
	my $index = 0;
	
	$index = &rawdata_QC($conf,$index);
	$index = &overlap($conf,$index);
	$index = &tags_QC($conf,$index);
	$index = &remove_chimera($conf,$index);
	$index = &tags_stat($conf,$index);
}

sub rawdata_QC
{
	my ($conf,$index) = @_;
	return $index unless $conf->{reads_filter};
	
	my $filter_fq = check_path($conf->{soft}->{filter_fq});

	# create dir 
	$index ++;
	my $dir = "$dirs->{data}/${index}.reads_filter";
	mkdir $dir unless -d $dir;
	timeLOG("**** 01_${index}: rawdata Quality Control ****");
	
	# create shell 
	open SH,">$dirs->{shell}/SH01_${index}.reads_filter.sh" or die $!;
	foreach my $sample (@samples)
	{
		my($read1,$read2) = fetch_rawdata($conf,$sample);
		print SH "perl $filter_fq $read1 $read2 $dir $sample $conf->{phred} $conf->{filter_options}\n";
		
		$conf->{cleandata}->{$sample} = [ "$dir/${sample}_1.fq" , "$dir/${sample}_2.fq" ];
	}
	close SH;
	
	run_shell("$dirs->{shell}/SH01_${index}.reads_filter.sh",$conf->{cpus}->{reads_filter}) if ($OPTS{steps} =~ /1/);

	return $index;
}

# assemble one pair read to one tag if overlap 
sub overlap
{
	my ($conf,$index) = @_;
	my $flash = check_path($conf->{soft}->{flash});
	my $overlap = check_path($conf->{soft}->{overlap});
	
	# create dir 
	$index ++;
	my $dir = "$dirs->{data}/${index}.overlap";
	mkdir $dir unless -d $dir;
	$dirs->{data_overlap} = $dir;
	timeLOG("**** 01_${index}: reads overlap assemble ****");

	# create shell 
	open SH,">$dirs->{shell}/SH01_${index}.overlap.sh" or die $!;
	foreach my $sample (@samples)
	{
		my($read1,$read2) = $conf->{reads_filter} ? @{$conf->{cleandata}->{$sample}} : fetch_rawdata($conf,$sample);
		print SH "$flash $conf->{overlap_options} -z -d $dir -p $conf->{phred} -o $sample $read1 $read2 1>$dir/$sample.log.o \n";
		$conf->{overlap_reads}->{$sample} = "$dir/$sample.extendedFrags.fastq.gz";

		#print SH "$overlap $conf->{overlap_options} -a $read1 -b $read2 -c $dir/$sample.left -d $dir/$sample.right -o $dir/$sample.overlap -q $dir/$sample.qual 1>$dir/$sample.log.o 2>$dir/$sample.log.e\n";
		#$conf->{overlap_reads}->{$sample} = "$dir/$sample.overlap";
		#$conf->{overlap_qual}->{$sample} = "$dir/$sample.qual";
	}
	close SH;
	
	run_shell("$dirs->{shell}/SH01_${index}.overlap.sh",$conf->{cpus}->{overlap}) if ($OPTS{steps} =~ /1/);
	return $index;
}

# tags 
sub tags_QC
{
	my ($conf,$index) = @_;
	return $index unless $conf->{tags_filter};
	
	# create dir 
	$index ++;
	my $dir = "$dirs->{data}/${index}.tags_filter";
	mkdir $dir unless -d $dir;
	timeLOG("**** 01_${index}: Tags Quality Control ****");
	
	# create shell 
	open SH,">$dirs->{shell}/SH01_${index}.tags_filter.sh" or die $!;
	foreach my $sample (@samples)
	{
		#print SH "perl $src/Trim_seq.pl -minl $min_len -maxl $max_len -phread $conf->{phred} -name $sample -outfile $dir/$sample.QC.tags.fasta $conf->{tags_filter_options} $conf->{overlap_reads}->{$sample} $conf->{overlap_qual}->{$sample}\n";
		print SH "perl $src/Trim_seq.pl -minl $min_len -maxl $max_len -phread $conf->{phred} -name $sample -outfile $dir/$sample.QC.tags.fasta $conf->{tags_filter_options} $conf->{overlap_reads}->{$sample}\n";
		
		$conf->{overlap_reads}->{$sample} = "$dir/$sample.QC.tags.fasta";
	}
	close SH;

	run_shell("$dirs->{shell}/SH01_${index}.tags_filter.sh",$conf->{cpus}->{tags_filter}) if ($OPTS{steps} =~ /1/);
	return $index;
}

sub remove_chimera
{
	my ($conf,$index) = @_;
	return $index unless $conf->{remove_chimera};
	
	my $chimera_db = check_path($conf->{db}->{chimera_db});
	my $usearch8 = check_path($conf->{soft}->{usearch8});

	# create dir 
	$index ++;
	my $dir = "$dirs->{data}/${index}.remove_chimera";
	mkdir $dir unless -d $dir;
	timeLOG("**** 01_${index}: remove the chimera from the tags ****");
	
	# create shell 
	open SH,">$dirs->{shell}/SH01_${index}.remove_chimera.sh" or die $!;
	foreach my $sample (@samples)
	{
		print SH "$usearch8 -uchime_ref $conf->{overlap_reads}->{$sample} -db $chimera_db -uchimeout $dir/$sample.chimera.txt -strand plus; ";
		print SH "perl $src/fetch_non_chimera.pl $conf->{overlap_reads}->{$sample} $dir/$sample.chimera.txt > $dir/$sample.non_chimera.fasta\n";
		$conf->{overlap_reads}->{$sample} = "$dir/$sample.non_chimera.fasta";
	}
	close SH;
	
	run_shell("$dirs->{shell}/SH01_${index}.remove_chimera.sh",$conf->{cpus}->{remove_chimera}) if ($OPTS{steps} =~ /1/);
	return $index;
}

sub tags_stat
{
	my ($conf,$index) = @_;
	
	my $usearch9 = check_path($conf->{soft}->{usearch9});

	# create dir 
	$index ++;
	my $dir = "$dirs->{data}/${index}.tags_stat";
	mkdir $dir unless -d $dir;
	timeLOG("**** 01_${index}: stat the tags ****");
	
	my $fix_tag_num = $conf->{fix_tag_num};
	my $cut_tag_num = $conf->{cut_tag_num};
	
	# create shell 
	my @tag_list;
	my @unique_tag_list;
	open SH,">$dirs->{shell}/SH01_${index}.tags_stat.sh" or die $!;
	foreach my $sample (@samples)
	{
		# rename tags and filter length 
		print SH "perl $src/rename_and_filter.pl $conf->{overlap_reads}->{$sample} $dir $sample $min_len $max_len $fix_tag_num $cut_tag_num; ";

		# fetch the unique tags and stat the duplication size 
		print SH "$usearch9 -fastx_uniques $dir/$sample.fasta -fastaout $dir/$sample.unique.fasta -sizeout\n";
		$conf->{overlap_reads}->{$sample} = "$dir/$sample.unique.fasta";

		push @tag_list , "$dir/$sample.fasta";
		push @unique_tag_list , "$dir/$sample.unique.fasta";
		
		$conf->{total_tag}->{$sample} = "$dir/$sample.fasta";
		$conf->{unique_tag}->{$sample} = "$dir/$sample.unique.fasta";
		$conf->{tag_group}->{$sample} = "$dir/$sample.groups";
	}
	
	# stat tags and unique tags 
	print SH qq(perl $src/tags_stat.pl Tags ${[ join ",",@tag_list ]}[0] > $dir/stat_tag.xls\n);
	print SH qq(perl $src/tags_stat.pl "Unique Tags" ${[ join ",",@unique_tag_list ]}[0] > $dir/stat_uniqTag.xls\n);
	print SH qq(perl $src/data_process_stat.pl $dirs->{data} $dir\n);
	print SH qq($Rscript $src/reads_stack.R $dir/all.data_process_stat.forFig.xls $dir\n);
	close SH;

	run_shell("$dirs->{shell}/SH01_${index}.tags_stat.sh",1) if ($OPTS{steps} =~ /1/);
	
	# create group files for some analysis
	groups2file($conf,"$dir/all.groups");
	groups2file_qiime($conf,"$dir/all.qiime.groups");
	$conf->{groups_file}  = "$dir/all.groups";
	$conf->{groups_qiime} = "$dir/all.qiime.groups";
	$conf->{stat_tag} = "$dir/stat_tag.xls";
	$conf->{tags_num} = read_tags_num("$dir/stat_tag.xls") if ( -e "$dir/stat_tag.xls" );
	
	$dirs->{tags_stat} = $dir;
	return $index;
}

#-------------------------------------------------------------------------------
#  step2: OTUs analysis 
#-------------------------------------------------------------------------------
#  1. merge all unique tags with Uparse (97% identity)
#  2. compare the OTUs of different samples/groups and draw venn diagram
#-------------------------------------------------------------------------------
sub otus
{
	timeLOG($split_line);
	timeLOG("**** 02: OTUs analysis ****");
	
	my $conf = shift;
	
	my $usearch9 = check_path($conf->{soft}->{usearch9});
	my $dir = "$dirs->{otus}";
	
	## merge Tags to OTUs 
	my @tag_list = map { $conf->{total_tag}->{$_} } @samples;
	my @groups_list = map { $conf->{tag_group}->{$_} } @samples;
	my $cmd1;

	if ($conf->{cluster_soft} eq "uparse")
	{
		$cmd1 = <<SH;
#!/bin/sh
# create all tags sequence file
cat ${[ join " ",@tag_list ]}[0] | perl -ne 'if(\$.%2==1){chomp;\$sample=(split /\_/)[1]; \$_ .= ";sample=\$sample;\n"; } print \$_' > $dir/all.tags.fasta

# fetch unique tags with usearch9
$usearch9 -fastx_uniques $dir/all.tags.fasta -fastaout $dir/all.unique_tag.fasta -sizeout

# cluster unique tags to OTUs with usearch9 (uparse)
$usearch9 -cluster_otus $dir/all.unique_tag.fasta -otus $dir/all.otus.fasta -uparseout $dir/all.up -relabel Otu -minsize 1

# create OTU abundance table with usearch9 and get the representative sequence
$usearch9 -usearch_global $dir/all.tags.fasta -db $dir/all.otus.fasta -strand plus -id 0.97 -otutabout $dir/all.otus.table

# 1. turn uparse result file (table format) to shared format
# 2. rename the OTU representative sequence
perl $src/uparse2mothur.pl $dir/all.otus.fasta $dir/all.otus.table $dir/all.otus.representative.fasta $dir/all.otus.shared

# create group shared file 
cd $dir
perl $src/gro_shared.pl $dir/all.otus.shared $conf->{groups_file}
SH
	}
	elsif ($conf->{cluster_soft} eq "mothur")
	{
		my $sequence_target = lc $conf->{sequence_target};
		my $align_seq = check_path($conf->{db}->{mothur_align_db}->{$sequence_target});
		
		my $cluster_cutoff = $conf->{otu_label};
		my $dist_cutoff = $cluster_cutoff + 0.03;
		my $diffs = int (($max_len+$min_len)/2 * ($cluster_cutoff-0.01));
		my $mothur_method = $conf->{mothur_cluster_method};

		$cmd1 = <<SH;
#!/bin/sh
# set the Python lib path for Mothur
PYTHONPATH_old=\$PYTHONPATH
export PYTHONPATH=$python_lib

cd $dir

# prepare the input files for Mothur
cat ${[ join " ",@tag_list ]}[0] > $dir/all.fasta
cat ${[ join " ",@groups_list ]}[0] > $dir/all.groups

# fetch the unique tag and align to reference
$mothur "#unique.seqs(fasta=$dir/all.fasta);count.seqs(name=$dir/all.names,group=$dir/all.groups);align.seqs(fasta=$dir/all.unique.fasta,reference=$align_seq,flip=T,processors=8);summary.seqs(fasta=$dir/all.unique.align,count=$dir/all.count_table)"

# pre-cluster
$mothur "#filter.seqs(fasta=$dir/all.unique.align);pre.cluster(fasta=$dir/all.unique.filter.fasta,name=$dir/all.names, diffs=$diffs, processors=8);count.seqs(name=$dir/all.unique.filter.precluster.names,group=$dir/all.groups)"

# calc the distance
$mothur "#dist.seqs(fasta=$dir/all.unique.filter.precluster.fasta,cutoff=$dist_cutoff,processors=8)"

# cluster
$mothur "#cluster(column=$dir/all.unique.filter.precluster.dist,count=$dir/all.unique.filter.precluster.count_table,method=$mothur_method,cutoff=$cluster_cutoff)"

# create the shared format file (contain abundance info)
$mothur "#make.shared(list=$dir/all.unique.filter.precluster.fn.unique_list.list,count=$dir/all.unique.filter.precluster.count_table,label=$cluster_cutoff);"

# rename the shared and fetch OTU representative seq file 
ln -s $dir/all.unique.filter.precluster.fn.unique_list.shared $dir/all.otus.shared
perl $src/all_otu_reformat.pl $dir/all.unique.filter.precluster.fn.unique_list.list $cluster_cutoff $dir/all.unique.filter.precluster.map $dir/all.names $dir/all.count_table $dir/all.fasta $dir/all_sample.reformat $dir/all.otus.representative.fasta

# create group shared file 
cd $dir
perl $src/gro_shared.pl $dir/all.otus.shared $conf->{groups_file}

# recovery the raw Python lib path
export PYTHONPATH=\$PYTHONPATH_old
SH
	}
	else 
	{
		ERROR("'$conf->{cluster_method}' was wrong, cluster_method must be <uparse> or <mothur>");
	}
	
	$conf->{otus_seq} = "$dir/all.otus.representative.fasta";
	$conf->{otus_shared} = "$dir/all.otus.shared";
	$conf->{otus_group_shared} = "$dir/all.otus.group.shared";

	timeLOG("**** 02_1: cluster Tags into OTUs ****");
	create_shell("SH02_1.otus.sh",$cmd1);
	run_shell("$dirs->{shell}/SH02_1.otus.sh",1) if ($OPTS{steps} =~ /2/);

	## Compare OTUs (venn diagram)
	my $cmd2 = <<CMD;
# compare OTUs of different groups or samples with venn diagram
cd $dir; mkdir -p venn; cd venn

# draw venn diagram
CMD
	
	$cmd2 .= qq($Rscript $src/otus_venn.R $conf->{otus_group_shared} "$conf->{otus_samples_venn}"\n) if ($conf->{otus_samples_venn});
	$cmd2 .= qq($Rscript $src/otus_venn.R $conf->{otus_group_shared} "$conf->{otus_groups_venn}"\n) if ($conf->{otus_groups_venn});
	$cmd2 .= "rm -f *.log\n";

	timeLOG("**** 02_2: OTUs compare analysis (venn diagram) ****");
	create_shell("SH02_2.otus.sh",$cmd2);
	run_shell("$dirs->{shell}/SH02_2.otus.sh",1) if ($OPTS{steps} =~ /2/);
}

sub taxonomy
{
	timeLOG($split_line);
	timeLOG("**** 03: OTUs taxonomy analysis ****");
	
	my $conf = shift;
	my $dir = $dirs->{taxonomy};
	
	my $rdp          = check_path($conf->{soft}->{rdp});
	my $rdp_py_lib   = check_path($conf->{soft}->{rdp_python_lib});
	my $assign_taxon = check_path($conf->{soft}->{assign_taxon});
	my $biom         = check_path($conf->{soft}->{biom});
	my $taxon_db     = $conf->{taxon_db};
	my $taxon_seq    = check_path($conf->{db}->{$taxon_db}->{seq});
	my $taxon_tax    = check_path($conf->{db}->{$taxon_db}->{tax});
	
	# create shell 
	my $cmd1 = <<CMD;
#!/bin/sh
#----------------------------
# 3.1 do taxonomy annot 
#----------------------------
cd $dir; mkdir -p 1.taxa_annot; cd 1.taxa_annot

#----------------------------
# set the python path 
PYTHONPATH_old=\$PYTHONPATH
export PYTHONPATH=$rdp_py_lib

# set the rdp PATH
export RDP_JAR_PATH=$rdp

#----------------------------
# do taxonomy annot with RDP classifier
$assign_taxon -i $conf->{otus_seq} -t $taxon_tax -r $taxon_seq -o $dir/1.taxa_annot $conf->{taxon_options}

# reback to raw python path
PYTHONPATH=\$PYTHONPATH_old
CMD
	
	$conf->{otus_rdp_res} = "$dir/1.taxa_annot/all.otus.representative_tax_assignments.txt";
	my $mttp = $conf->{min_taxtree_percent};
	my $mtfp = $conf->{min_taxon_fill_percent};
	my $mtfn = $conf->{max_taxon_fill_number};
	my $mthp = $conf->{min_taxon_heatmap_percent};
	my $mthn = $conf->{max_taxon_heatmap_number};
	my $mohp = $conf->{min_otu_heatmap_percent};

	my $cmd2 = <<CMD;
#!/bin/sh
#----------------------------
# 3.2 fetch the taxonomy profiling info
#----------------------------
PYTHONPATH_old=\$PYTHONPATH
export PYTHONPATH=$python_lib

cd $dir; mkdir -p 2.taxa_profiling; cd 2.taxa_profiling
#----------------------------

# turn rdp annot result to taxonomy profiling file 
perl $src/rdp2profiling.pl --mohp $mohp --mtfp $mtfp $conf->{otus_shared} $conf->{otus_rdp_res} $dir/2.taxa_profiling
$biom convert -i $dir/2.taxa_profiling/all.otus.tab -o $dir/2.taxa_profiling/all.otus_tab.biom --table-type "OTU table" --to-hdf5

# create the total taxonomy profiling file 
perl $src/profiling2tree.pl $dir/2.taxa_profiling/all.otus.profiling.xls

#----------------------------
# 3.3 do some stat for taxonomy info 
#----------------------------
cd $dir; mkdir -p 3.taxa_stat; cd 3.taxa_stat
#----------------------------

# draw a figure to stat the taxonomy annot percentage of each level
$Rscript $src/taxon_stack.R $dir/2.taxa_profiling/all.taxonomy.stat.forFig.xls

# draw taxonomy fill bar figure
# draw OTU phylogenetic tree and taxonomy fill figure 
mkdir -p fill_bar; cd fill_bar
perl $src/profiling_draw.pl -n $mtfn -p $mtfp -b $dir/2.taxa_profiling/profiling.L2.Phylum.min$mtfp.xls Phylum
perl $src/profiling_draw.pl -n $mtfn -p $mtfp -b $dir/2.taxa_profiling/profiling.L3.Class.min$mtfp.xls Class
perl $src/profiling_draw.pl -n $mtfn -p $mtfp -b $dir/2.taxa_profiling/profiling.L4.Order.min$mtfp.xls Order
perl $src/profiling_draw.pl -n $mtfn -p $mtfp -b $dir/2.taxa_profiling/profiling.L5.Family.min$mtfp.xls Family
perl $src/profiling_draw.pl -n $mtfn -p $mtfp -b $dir/2.taxa_profiling/profiling.L6.Genus.min$mtfp.xls Genus
perl $src/profiling_draw.pl -n $mtfn -p $mtfp -b $dir/2.taxa_profiling/profiling.L7.Species.min$mtfp.xls Species
cd ..

# KRONA
mkdir -p krona; cd krona;
perl $src/run_krona.pl $dir/2.taxa_profiling/all.otus.profiling.xls 
cd ..

# draw taxonomy tree with abundance info 
mkdir -p taxtree; cd taxtree;
perl $src/taxtree.pl -m $mttp -k all -g $conf->{groups_file} $dir/2.taxa_profiling/profiling.all.taxonomy.v2.xls
cd ..

rm -f Rplots.pdf

# stat the number of tags (toatl, unique, singleton, taxon, unclassified) and OTUs 
perl $src/tags_otus_count.pl $conf->{dirs}->{odir}
$Rscript $src/tags_otus.count.R $dir/3.taxa_stat/all.tags_otus.stat.forFig.xls all.tags_otus_count

#----------------------------
# 3.4 best taxonomy level 
#----------------------------
cd $dir ; mkdir -p 4.best_level; cd 4.best_level
perl $src/wholeTaxStat.pl $dir/2.taxa_profiling/all.taxonomy.stat.xls $conf->{stat_tag}
ln -sf $dir/2.taxa_profiling/all.taxonomy.stat.xls

#----------------------------
# 3.5 taxonomy abundance clust 
#----------------------------
cd $dir; mkdir -p 5.taxa_heatmap; cd 5.taxa_heatmap
#----------------------------
# taxonomy abundance heatmap 
perl $src/profiling_draw.pl -n $mthn -p $mthp -h $dir/2.taxa_profiling/profiling.L2.Phylum.min$mohp.xls Phylum
perl $src/profiling_draw.pl -n $mthn -p $mthp -h $dir/2.taxa_profiling/profiling.L3.Class.min$mohp.xls Class
perl $src/profiling_draw.pl -n $mthn -p $mthp -h $dir/2.taxa_profiling/profiling.L4.Order.min$mohp.xls Order
perl $src/profiling_draw.pl -n $mthn -p $mthp -h $dir/2.taxa_profiling/profiling.L5.Family.min$mohp.xls Family
perl $src/profiling_draw.pl -n $mthn -p $mthp -h $dir/2.taxa_profiling/profiling.L6.Genus.min$mohp.xls Genus
perl $src/profiling_draw.pl -n $mthn -p $mthp -h $dir/2.taxa_profiling/profiling.L7.Species.min$mohp.xls Species

# OTU abundance heatmap (need revise)
$Rscript $src/heatmap.R $dir/2.taxa_profiling/all.otus.exp_for_heatmap $dir/5.taxa_heatmap/all.otus

# interactive OTU and taxonomy abundance heatmap (QIIME)
# underdevelopment

# reback to raw python path
PYTHONPATH=\$PYTHONPATH_old
CMD
	
	timeLOG("**** 03_1: OTUs taxonomy annot ****");
	create_shell("SH03_1.taxonomy_annot.sh",$cmd1);
	run_shell("$dirs->{shell}/SH03_1.taxonomy_annot.sh",1) if ($OPTS{steps} =~ /3/);

	timeLOG("**** 03_2: OTUs taxonomy stat ****");
	create_shell("SH03_2.taxonomy_stat.sh",$cmd2);
	run_shell("$dirs->{shell}/SH03_2.taxonomy_stat.sh",1) if ($OPTS{steps} =~ /3/);
	
	$conf->{otu_abu_file} = "$dir/2.taxa_profiling/all.otus.tab";
	$conf->{otu_pro_file} = "$dir/2.taxa_profiling/all.otus.profiling.xls";
	$conf->{tax_pro_file} = "$dir/2.taxa_profiling/profiling.all.taxonomy.v2.xls";
}

sub alpha_diversity
{
	timeLOG($split_line);
	timeLOG("**** 04: Alpha diversity analysis ****");
	
	my $conf = shift;
	my $dir = $dirs->{alpha};
	my $tax_dir = $dirs->{taxonomy};

	my $step = $conf->{alpha_diversity_step} || 1000;
	my $min_tags = $conf->{tags_num} ? min(values %{$conf->{tags_num}}) : 50010;

	my $cmd1 = <<SH;
#!/bin/sh
# set the Python lib path for QIIME
PYTHONPATH_old=\$PYTHONPATH
export PYTHONPATH=$python_lib

#-------------------------------------------------------------------------------
# 4.1 calc the alpha diversity with QIIME
#-------------------------------------------------------------------------------
cd $dir; mkdir -p 1.prepare; cd 1.prepare

# calc the multiple rarefaction alpha diversity
$python $qiime/multiple_rarefactions.py -i $tax_dir/2.taxa_profiling/all.otus_tab.biom -o $dir/1.prepare/multiple_rarefactions -m 10 -x $min_tags -s $step
$python $qiime/alpha_diversity.py -i $dir/1.prepare/multiple_rarefactions -m chao1,goods_coverage,observed_species,shannon,simpson -o $dir/1.prepare/div_alpha
$python $qiime/collate_alpha.py -i $dir/1.prepare/div_alpha -o $dir/1.prepare/collate_alpha
$python $qiime/make_rarefaction_plots.py -i $dir/1.prepare/collate_alpha -m $conf->{groups_qiime} -o $dir/1.prepare/alpha_plot --generate_average_tables --generate_per_sample_plots -d 300

# calc the summary alpha diversity
$python $qiime/alpha_diversity.py -i $tax_dir/2.taxa_profiling/all.otus_tab.biom -m chao1,ace,goods_coverage,observed_species,shannon,simpson -o $dir/1.prepare/summary.alpha_diversity

#-------------------------------------------------------------------------------
# 4.2 draw rarefaction curve
#-------------------------------------------------------------------------------
cd $dir; mkdir -p 2.rarefaction_curve; cd 2.rarefaction_curve
ln -sf $dir/1.prepare/alpha_plot/average_tables/observed_speciesGroup.txt
ln -sf $dir/1.prepare/alpha_plot/average_tables/observed_speciesSampleID.txt
perl $src/qiime2ggplot2.pl observed_speciesGroup.txt 
perl $src/qiime2ggplot2.pl observed_speciesSampleID.txt
$Rscript $src/fold_line.R observed_speciesGroup.forGGPLOT2.txt 
$Rscript $src/fold_line.R observed_speciesSampleID.forGGPLOT2.txt
rm Rplots.pdf

#-------------------------------------------------------------------------------
# 4.3 draw shannon rarefaction curve
#-------------------------------------------------------------------------------
cd $dir; mkdir -p 3.shannon_rarefaction_curve; cd 3.shannon_rarefaction_curve
ln -sf $dir/1.prepare/alpha_plot/average_tables/shannonSampleID.txt
ln -sf $dir/1.prepare/alpha_plot/average_tables/shannonGroup.txt
ln -sf $dir/1.prepare/alpha_plot/average_tables/chao1Group.txt
ln -sf $dir/1.prepare/alpha_plot/average_tables/chao1SampleID.txt

perl $src/qiime2ggplot2.pl shannonSampleID.txt
perl $src/qiime2ggplot2.pl shannonGroup.txt
perl $src/qiime2ggplot2.pl chao1SampleID.txt
perl $src/qiime2ggplot2.pl chao1Group.txt

$Rscript $src/fold_line.R shannonGroup.forGGPLOT2.txt
$Rscript $src/fold_line.R shannonSampleID.forGGPLOT2.txt
$Rscript $src/fold_line.R chao1Group.forGGPLOT2.txt
$Rscript $src/fold_line.R chao1SampleID.forGGPLOT2.txt
rm Rplots.pdf

#-------------------------------------------------------------------------------
# 4.4 draw OTU Rank Abundance
#-------------------------------------------------------------------------------
cd $dir ; mkdir -p 4.otu_rank_abundance; cd 4.otu_rank_abundance
perl $src/draw_rank_abund.pl $conf->{otus_shared} $dir/4.otu_rank_abundance

#-------------------------------------------------------------------------------
# 4.5 calc the summary alpha diversity with Mothur
#-------------------------------------------------------------------------------
# calc the alpha diversity
# $mothur "#rarefaction.single(shared=$conf->{otus_shared},processors=8);rarefaction.single(shared=$conf->{otus_shared},calc=shannon,processors=8);summary.single(shared=$conf->{otus_shared},calc=sobs-chao-ace-jack-shannon-npshannon-simpson-coverage)"

export PYTHONPATH=\$PYTHONPATH_old
SH
	
	my $cmd2 = <<CMD;
#-------------------------------------------------------------------------------
# diff alpha diversity analysis
#-------------------------------------------------------------------------------
cd $dir; mkdir -p 5.diff_alpha_diversity
perl $src/diffs/meta16s_diversity_diff.pl -conf $conf_file -alpha $dir/1.prepare/summary.alpha_diversity -outdir $dir/5.diff_alpha_diversity
perl $src/diversity_diff_summary.pl --alpha $dir/5.diff_alpha_diversity
CMD

	timeLOG("**** 04_1: Alpha diversity calc and visulization ****");
	create_shell("SH04_1.alpha_diversity.sh",$cmd1);
	run_shell("$dirs->{shell}/SH04_1.alpha_diversity.sh",1) if ($OPTS{steps} =~ /4/);
	
	timeLOG("**** 04_2: Alpha diversity diff analysis ****");
	create_shell("SH04_2.diff_alpha_diversity.sh",$cmd2);
	run_shell("$dirs->{shell}/SH04_2.diff_alpha_diversity.sh",1) if ($OPTS{steps} =~ /4/);
}

sub beta_diversity
{
	timeLOG($split_line);
	timeLOG("**** 05: Beta diversity analysis ****");
	
	my $conf = shift;
	my $dir = $dirs->{beta};
	my $tax_dir = $dirs->{taxonomy};
	
	if ($sample_num < 3)
	{
		WARN("The samples number is less than 3, skip the beta diversity analysis");
		return 1;
	}

	my $treebest   = check_path($conf->{soft}->{treebest});
	my $muscle     = check_path($conf->{soft}->{muscle});
	my $mohp       = $conf->{min_otu_heatmap_percent};

	my $cmd1 = <<CMD;
#!/bin/sh
PYTHONPATH_old=\$PYTHONPATH
export PYTHONPATH=$python_lib

#-------------------------------------------------------------------------------
#  5.1 
#  a. calc the (un)weighted unifrac distance
#  b. do PCoA and draw scatter plot 
#  c. do NMDS and draw scatter plot 
#  d. do UPGMA and draw tree plot 
#-------------------------------------------------------------------------------
cd $dir; mkdir -p 1.unifrac; cd 1.unifrac

# build the OTU phylogenetic tree for calc the unifrac distance
$muscle -in $conf->{otus_seq} -out $dir/1.unifrac/all.otus.aln.fasta
$treebest nj -b 10 $dir/1.unifrac/all.otus.aln.fasta > $dir/1.unifrac/all.otus.nhx
sed 's/\\[\\&\\&NH.*\\]//' $dir/1.unifrac/all.otus.nhx > $dir/1.unifrac/all.otus.nwk

# calc unifrac and do PCoA, NMDS
$Rscript $src/beta_diversity.R  $tax_dir/2.taxa_profiling/all.otus.exp_for_unifrac  $dir/1.unifrac/all.otus.nwk $conf->{groups_qiime}

# do UPGMA analysis and draw phylo tree and taxonomy fill bar
$python $qiime/upgma_cluster.py -i $dir/1.unifrac/weighted_unifrac.tsv -o $dir/1.unifrac/weighted_unifrac.UPGMA.tre
$python $qiime/upgma_cluster.py -i $dir/1.unifrac/unweighted_unifrac.tsv -o $dir/1.unifrac/unweighted_unifrac.UPGMA.tre

for i in $tax_dir/2.taxa_profiling/profiling.L*.*.min2.xls
do
	fname=`basename \$i`
	level=`echo \$fname | awk -F '.' '{print \$2}'`
	taxon=`echo \$fname | awk -F '.' '{print \$3}'`
	if [ \$level == "L1" ]
	then
		continue
	fi
	$Rscript $src/upgma_plot.R $dir/1.unifrac/weighted_unifrac.UPGMA.tre \$i \$taxon weighted 20
	$Rscript $src/upgma_plot.R $dir/1.unifrac/unweighted_unifrac.UPGMA.tre \$i \$taxon unweighted 20
done

#-------------------------------------------------------------------------------
#  5.2 do PCA analysis with OTU abundance table
#-------------------------------------------------------------------------------
cd $dir; mkdir -p 2.PCA; cd 2.PCA

# do PCA with OTU abundance table 
$Rscript $src/pca.v2.R $tax_dir/2.taxa_profiling/all.otus.exp_for_unifrac 2,${[ $sample_num+1 ]}[0] $conf->{groups_qiime}

# do PCA with Taxonomy abundance table 
$Rscript $src/pca.v2.R $tax_dir/2.taxa_profiling/profiling.L2.Phylum.min$mohp.xls 2,${[ $sample_num+1 ]}[0] $conf->{groups_qiime}
$Rscript $src/pca.v2.R $tax_dir/2.taxa_profiling/profiling.L3.Class.min$mohp.xls 2,${[ $sample_num+1 ]}[0] $conf->{groups_qiime}
$Rscript $src/pca.v2.R $tax_dir/2.taxa_profiling/profiling.L4.Order.min$mohp.xls 2,${[ $sample_num+1 ]}[0] $conf->{groups_qiime}
$Rscript $src/pca.v2.R $tax_dir/2.taxa_profiling/profiling.L5.Family.min$mohp.xls 2,${[ $sample_num+1 ]}[0] $conf->{groups_qiime}
$Rscript $src/pca.v2.R $tax_dir/2.taxa_profiling/profiling.L6.Genus.min$mohp.xls 2,${[ $sample_num+1 ]}[0] $conf->{groups_qiime}
$Rscript $src/pca.v2.R $tax_dir/2.taxa_profiling/profiling.L7.Species.min$mohp.xls 2,${[ $sample_num+1 ]}[0] $conf->{groups_qiime}

export PYTHONPATH=\$PYTHONPATH_old
CMD

	my $cmd2 = <<CMD;
#!/bin/sh
# for diff beta diversity
cd $dir; mkdir -p 3.diff_beta_diversity; cd 3.diff_beta_diversity
perl $src/diffs/meta16s_diversity_diff.pl -conf $conf_file -dist $dir/1.unifrac/unweighted_unifrac.tsv -outdir $dir/3.diff_beta_diversity/unweighted_unifrac
perl $src/diffs/meta16s_diversity_diff.pl -conf $conf_file -dist $dir/1.unifrac/weighted_unifrac.tsv -outdir $dir/3.diff_beta_diversity/weighted_unifrac
perl $src/diversity_diff_summary.pl --beta $dir/3.diff_beta_diversity

$Rscript $src/diffs/Dissimilarity_diff.R $dir/1.unifrac/weighted_unifrac.tsv $conf->{groups_file} "$conf->{adonis_groups_diff}"
$Rscript $src/diffs/Dissimilarity_diff.R $dir/1.unifrac/unweighted_unifrac.tsv $conf->{groups_file} "$conf->{adonis_groups_diff}"
CMD
	
	timeLOG("**** 05_1: Beta diversity calc and stat ****");
	create_shell("SH05_1.beta_diversity.sh",$cmd1);
	run_shell("$dirs->{shell}/SH05_1.beta_diversity.sh",1) if ($OPTS{steps} =~ /5/);
	
	timeLOG("**** 05_2: Beta diversity diff analysis ****");
	create_shell("SH05_2.diff_beta_diversity.sh",$cmd2);
	run_shell("$dirs->{shell}/SH05_2.diff_beta_diversity.sh",1) if ($OPTS{steps} =~ /5/);
}

sub function_annot
{
	timeLOG($split_line);
	timeLOG("**** 06: Function annotation ****");
	
	my $conf = shift;
	my $dir = $dirs->{function};
	my $otu_dir = $dirs->{otus};
	
	my $cmd1 = <<CMD;
#!/bin/sh
# do function annotation 
CMD

	# set the default annot software
	$conf->{anno_soft} = $target =~ /16s/i ? "PICRUSt"  :
						 $target =~ /its/i ? "FUNGuild" : "none" unless $conf->{anno_soft};

	if ($target =~ /16s/i && $conf->{anno_soft} eq "PICRUSt")
	{
		timeLOG("do pathway annotation with PICRUSt");
		
		my $greengene_align_seq = check_path($conf->{db}->{greengene_align_seq});
		$cmd1 .= <<CMD;
# set the Python lib path 
PYTHONPATH_old=\$PYTHONPATH
export PYTHONPATH=$python_lib

# align otu representative seq to <greengene aligned seq> with Mothur
$mothur "#align.seqs(fasta=$conf->{otus_seq},reference=$greengene_align_seq,flip=T,processors=8);"

# do pathway annotation with PICRUSt
# 1. create KEGG L4 abundance table 
# 2. create KEGG L3 (Pathway) abundance table
# 3. create spf file for STAMP which can do Function different analysis
perl $src/meta16S_pathway_analysis.v1.3.pl -k all -a $otu_dir/all.otus.representative.align.report -o $dir -s 12345 $conf->{otu_abu_file}

export PYTHONPATH=\$PYTHONPATH_old
CMD
	}
	elsif ($target =~ /16s/i && $conf->{anno_soft} eq "tax4fun")
	{
		$cmd1 .= <<CMD;

CMD
	}
	elsif ($target =~ /its/i && $conf->{anno_soft} eq "FUNGuild")
	{
		my $FUNGuild = check_path($conf->{soft}->{funguild});
		$cmd1 .= <<CMD;
$python $FUNGuild -otu $dirs->{taxonomy}/2.taxa_profiling/all.otus.expr_for_FUNGuild -db fungi
CMD
	}
	elsif ($conf->{anno_soft} eq "none")
	{
		WARN("skip function annotation");
		return 1;
	}
	else 
	{
		WARN("skip function annotation, $conf->{anno_soft} is not for $target");
		$conf->{anno_soft} = "none";
		return 1;
	}
	
	create_shell("SH06_1.function.sh",$cmd1);
	run_shell("$dirs->{shell}/SH06_1.function.sh",1) if ($OPTS{steps} =~ /6/);
}

sub different_analysis
{
	timeLOG($split_line);
	timeLOG("**** 07: different analysis ****");
	
	my $conf = shift;
	my $dir = $dirs->{compare};
	
	my $cmd1 = <<CMD;
#!/bin/sh
# do taxonomy diff analysis with Metastats and Lefse
perl $src/meta_diff.v2.pl -ml -r multi -c $conf->{cpus}->{diff} -o $dir $conf_file $dirs->{taxonomy}/2.taxa_profiling/all.otus.profiling.xls
CMD

	create_shell("SH07_1.diff.sh",$cmd1);
	run_shell("$dirs->{shell}/SH07_1.diff.sh",1) if ($OPTS{steps} =~ /7/);
}

sub pack_result
{
	timeLOG($split_line);
	timeLOG("**** 08: pack the result ****");

	my $conf = shift;
	my $pack_dir = $dirs->{upload};
	
	my $cmd = <<CMD;
#!/bin/sh 
# pack the the result 

#------------------------------------------
# prepare
# get into upload dir 
cd $pack_dir

# empty the upload dir
rm -rf $pack_dir/*

# define the push files function
# usage: push_files "\${files[*]} current_dir source_dir"
push_files (){
	files=\$1
	for i in \${files[*]}
	do
		ln -sf \$2/\$i \$3
	done
}

#-----------------------------------------
### create directory tree
## Tags
mkdir 01.tags
mkdir 01.tags/01.sequence
mkdir 01.tags/02.stat

# push unique tags and all tags sequence
for i in $dirs->{tags_stat}/*.fasta
do
	ln -sf \$i 01.tags/01.sequence
done

# push stat files 
files=(stat_tag.xls stat_uniqTag.xls all.data_process_stat.xls reads_fill.pdf reads_fill.png reads_stack.pdf reads_stack.png)
push_files "\${files[*]}" $dirs->{tags_stat} 01.tags/02.stat

## OTUs
mkdir 02.otus
mkdir 02.otus/01.sequence
#mkdir 02.otus/02.compare
ln -sf $dirs->{otus}/all.otus.representative.fasta 02.otus/01.sequence
ln -sf $dirs->{taxonomy}/2.taxa_profiling/all.otus.profiling.xls 02.otus/01.sequence
ln -sf $dirs->{taxonomy}/3.taxa_stat/all.tags_otus.stat.xls 02.otus/01.sequence
ln -sf $dirs->{taxonomy}/3.taxa_stat/all.tags_otus_count.png 02.otus/01.sequence
ln -sf $dirs->{taxonomy}/3.taxa_stat/all.tags_otus_count.pdf 02.otus/01.sequence

ln -sf $dirs->{otus}/venn 02.otus/02.compare

## Taxonomy
mkdir 03.taxonomy
mkdir 03.taxonomy/01.taxa_annot
mkdir 03.taxonomy/02.taxa_profiling
mkdir 03.taxonomy/03.taxa_tree
mkdir 03.taxonomy/04.fill_bar
mkdir 03.taxonomy/05.heatmap

# push rdp annot result and summary stat result
ln -sf $dirs->{taxonomy}/1.taxa_annot/all.otus.representative_tax_assignments.txt 03.taxonomy/01.taxa_annot
ln -sf $dirs->{taxonomy}/2.taxa_profiling/all.taxonomy.stat.xls 03.taxonomy/01.taxa_annot
ln -sf $dirs->{taxonomy}/3.taxa_stat/krona/all.taxonomy.krona.html 03.taxonomy/01.taxa_annot
files=(taxonomy_fill.png taxonomy_fill.pdf taxonomy_stack.png taxonomy_stack.pdf)
push_files "\${files[*]}" $dirs->{taxonomy}/3.taxa_stat 03.taxonomy/01.taxa_annot

# push taxonomy profiling files
files=(profiling.all.taxonomy.xls profiling.all.taxonomy.v2.xls)
push_files "\${files[*]}" $dirs->{taxonomy}/2.taxa_profiling 03.taxonomy/02.taxa_profiling
for i in $dirs->{taxonomy}/2.taxa_profiling/profiling.L?.*.xls
do
	ln -sf \$i 03.taxonomy/02.taxa_profiling
done

# taxonomy tree
files=(all.groups_taxtree.png all.groups_taxtree.svg all.groups.percent all.samples_taxtree.png all.samples_taxtree.svg all.samples.percent)
push_files "\${files[*]}" $dirs->{taxonomy}/3.taxa_stat/taxtree 03.taxonomy/03.taxa_tree

# taxonomy fill bar 
for i in $dirs->{taxonomy}/3.taxa_stat/fill_bar/*.bar.p??
do
	ln -sf \$i 03.taxonomy/04.fill_bar
done

# best level result
#files=(best_level_stat.xls all.TagRatio.pdf all.TagRatio.png all.TagNumber.pdf all.TagNumber.png)
#push_files "\${files[*]}" $dirs->{taxonomy}/4.best_level 03.taxonomy/05.best_level

# heatmap 
for i in $dirs->{taxonomy}/5.taxa_heatmap/*.heatmap.p?? 
do
	ln -sf \$i 03.taxonomy/05.heatmap
done

## alpha diversity
mkdir 04.alpha_diversity
#mkdir 04.alpha_diversity/01.alpha_plot
mkdir 04.alpha_diversity/02.rarefaction_curve
mkdir 04.alpha_diversity/03.shannon_rarefaction_curve
mkdir 04.alpha_diversity/04.otu_rank_abundance

ln -sf $dirs->{alpha}/1.prepare/summary.alpha_diversity 04.alpha_diversity/summary.alpha_diversity.xls
ln -sf $dirs->{alpha}/1.prepare/alpha_plot 04.alpha_diversity/01.alpha_plot

files=(observed_speciesGroup.png observed_speciesGroup.pdf observed_speciesSampleID.png observed_speciesSampleID.pdf)
push_files "\${files[*]}" $dirs->{alpha}/2.rarefaction_curve 04.alpha_diversity/02.rarefaction_curve

files=(chao1Group.pdf chao1Group.png chao1SampleID.pdf chao1SampleID.png shannonGroup.pdf shannonGroup.png shannonSampleID.pdf shannonSampleID.png)
push_files "\${files[*]}" $dirs->{alpha}/3.shannon_rarefaction_curve 04.alpha_diversity/03.shannon_rarefaction_curve

ln -sf $dirs->{alpha}/4.otu_rank_abundance/Rank_abd_Distr.png 04.alpha_diversity/04.otu_rank_abundance
ln -sf $dirs->{alpha}/4.otu_rank_abundance/Rank_abd_Distr.pdf 04.alpha_diversity/04.otu_rank_abundance
CMD
	
	if ($conf->{two_pass_diffs} || $conf->{multi_pass_diffs})
	{
		$cmd .= <<CMD;
# for diff alpha diversity
mkdir 04.alpha_diversity/05.diff_alpha_diversity
cp -r $dirs->{alpha}/5.diff_alpha_diversity/* 04.alpha_diversity/05.diff_alpha_diversity
rm 04.alpha_diversity/05.diff_alpha_diversity/*/file.txt

CMD
	}

	my $index = "05";
	if ($sample_num >= 3)
	{
		$cmd .= <<CMD;
## beta diversity
mkdir $index.beta_diversity
mkdir $index.beta_diversity/01.unifrac
#mkdir $index.beta_diversity/02.PCA
mkdir $index.beta_diversity/03.PCoA
mkdir $index.beta_diversity/04.NMDS
mkdir $index.beta_diversity/05.UPGMA

# unifrac distance file
ln -sf $dirs->{beta}/1.unifrac/weighted_unifrac.tsv $index.beta_diversity/01.unifrac/weighted_unifrac.xls
ln -sf $dirs->{beta}/1.unifrac/unweighted_unifrac.tsv $index.beta_diversity/01.unifrac/unweighted_unifrac.xls
files=(unweighted_unifrac.heatmap.pdf unweighted_unifrac.heatmap.png weighted_unifrac.heatmap.pdf weighted_unifrac.heatmap.png)
push_files "\${files[*]}" $dirs->{beta}/1.unifrac $index.beta_diversity/01.unifrac

# PCA 
ln -sf $dirs->{beta}/2.PCA $index.beta_diversity/02.PCA

# PCoA
files=(unweighted_unifrac.PCoA.nonames.pdf unweighted_unifrac.PCoA.pdf weighted_unifrac.PCoA.nonames.pdf weighted_unifrac.PCoA.pdf unweighted_unifrac.PCoA.nonames.png unweighted_unifrac.PCoA.png weighted_unifrac.PCoA.nonames.png weighted_unifrac.PCoA.png)
push_files "\${files[*]}" $dirs->{beta}/1.unifrac $index.beta_diversity/03.PCoA

# NMDS
files=(unweighted_unifrac.NMDS.nonames.pdf unweighted_unifrac.NMDS.pdf weighted_unifrac.NMDS.nonames.pdf weighted_unifrac.NMDS.pdf unweighted_unifrac.NMDS.nonmaes.png unweighted_unifrac.NMDS.png weighted_unifrac.NMDS.nonmaes.png weighted_unifrac.NMDS.png)
push_files "\${files[*]}" $dirs->{beta}/1.unifrac $index.beta_diversity/04.NMDS

# UPGMA
files=(unweighted_unifrac.UPGMA.tre weighted_unifrac.UPGMA.tre)
push_files "\${files[*]}" $dirs->{beta}/1.unifrac $index.beta_diversity/05.UPGMA
for i in `ls $dirs->{beta}/1.unifrac/*.*.UPGMA_stack.p??`
do
	ln -sf \$i $index.beta_diversity/05.UPGMA
done
CMD
		
		if ($conf->{two_pass_diffs} || $conf->{multi_pass_diffs} || $conf->{adonis_groups_diff})
		{
			$cmd .= <<CMD;
# for diff beta diversity
mkdir $index.beta_diversity/06.diff_beta_diversity
mkdir $index.beta_diversity/06.diff_beta_diversity/weighted_unifrac
mkdir $index.beta_diversity/06.diff_beta_diversity/unweighted_unifrac

cp -rf $dirs->{beta}/3.diff_beta_diversity/weighted_unifrac $index.beta_diversity/06.diff_beta_diversity
cp -rf $dirs->{beta}/3.diff_beta_diversity/unweighted_unifrac $index.beta_diversity/06.diff_beta_diversity
rm $index.beta_diversity/06.diff_beta_diversity/*/file.txt
ln -sf $dirs->{beta}/3.diff_beta_diversity/all.beta_diversity.groups_diff.summary.xls $index.beta_diversity/06.diff_beta_diversity

# for diff beta diversity (Adonis and Anosim)
ln -sf $dirs->{beta}/3.diff_beta_diversity/unweighted_unifrac.tsv.ssim.diff.xls $index.beta_diversity/06.diff_beta_diversity/unweighted_unifrac
ln -sf $dirs->{beta}/3.diff_beta_diversity/weighted_unifrac.tsv.ssim.diff.xls $index.beta_diversity/06.diff_beta_diversity/weighted_unifrac
CMD
		}
		$index ++;
	}
	
	# push function annot result
	if ($conf->{anno_soft} ne "none")
	{
		$cmd .= <<CMD;
# function annot 
mkdir -p $index.function_annot
mkdir -p $index.function_annot/1.annot_result

CMD
		
		if ($conf->{anno_soft} eq "PICRUSt")
		{
			$cmd .= <<CMD;
# for PICRUSt
ln -sf $dirs->{function}/all.ko.NSTI.tab $index.function_annot/1.annot_result/all.ko.NSTI.xls
ln -sf $dirs->{function}/all.pathway.annot.xls $index.function_annot/1.annot_result/all.pathway.annot.xls
ln -sf $dirs->{function}/all.ko.annot.xls $index.function_annot/1.annot_result/all.ko.annot.xls
ln -sf $dirs->{function}/all.pathway.L3.spf $index.function_annot/1.annot_result/all.pathway.L3.spf
CMD
		}
		elsif ($conf->{anno_soft} eq "Tax4Fun")
		{
			
		}
		elsif ($conf->{anno_soft} eq "FUNGuild")
		{
			$cmd .= <<CMD;
# FUNGuild
ln -sf $dirs->{function}/*.guilds.xls $index.function_annot/1.annot_result
CMD
		}

		$index ++;
	}
	
	# push diff analysis result
	if ($conf->{two_groups_diff} || $conf->{multi_groups_diff})
	{
		$cmd .= <<CMD;
# for taxonomy diff analysis
mkdir -p $index.taxa_diff

CMD
		my $metastat_pack_cmd = <<CMD;
# for Metastats
mkdir -p $index.taxa_diff/metastat
mkdir -p $index.taxa_diff/metastat/2.Phylum $index.taxa_diff/metastat/3.Class $index.taxa_diff/metastat/4.Order $index.taxa_diff/metastat/5.Family $index.taxa_diff/metastat/6.Genus $index.taxa_diff/metastat/7.Species

taxones=(Root Domain Phylum Class Order Family Genus Species)
for i in {2..7}
do 
	for j in $dirs->{compare}/metastat/\${taxones[\$i]}/*.metastat.result.xls
	do
		ln -sf \$j $index.taxa_diff/metastat/\$i.\${taxones[\$i]}/
	done
done
CMD
		my $lefse_pack_cmd    = <<CMD;

# for LefSe
mkdir -p $index.taxa_diff/Lefse
for i in $dirs->{compare}/LefSe/*.Lefse.svg
do
	name=`basename \$i .Lefse.svg`
	files=(\$name.Lefse.cladogram.png \$name.Lefse.png \$name.Lefse.cladogram.svg \$name.Lefse.svg \$name)
	push_files "\${files[*]}" $dirs->{compare}/LefSe $index.taxa_diff/Lefse
done
CMD
		
		$cmd .= $metastat_pack_cmd if ($conf->{metastats_groups_diff});
		$cmd .= $lefse_pack_cmd if ($conf->{lefse_groups_diff});

		$index ++;
	}
	
	# push analysis related with env factor
	if ($conf->{env_factor})
	{
		$index ++;
	}
	
	# create html report 
	my $cmd2 = "$Bin/meta16S_report.v2.pl $conf_file $pack_dir";

	create_shell("SH08_1.pack.sh",$cmd);
	run_shell("$dirs->{shell}/SH08_1.pack.sh",1) if ($OPTS{steps} =~ /8/);
	
	create_shell("SH08_2.report.sh",$cmd2);
	run_shell("$dirs->{shell}/SH08_2.report.sh",1) if ($OPTS{steps} =~ /8/);
}


#-------------------------------------------------------------------------------
#  create shell 
#-------------------------------------------------------------------------------
sub create_shell
{
	my ($file,$cmd) = @_;

	open my $fh , ">" , "$dirs->{shell}/$file" or die "can't open file $!";
	print $fh $cmd;
	close $fh;

	return 1;
}

#-------------------------------------------------------------------------------
# run the shell 
#-------------------------------------------------------------------------------
sub run_shell
{
	my ($shell,$cpu,%opts) = @_;
	my $queue = $opts{'-queue'} || $conf->{'queue'};
	my $mem = $opts{'-mem'} || 2;

	$cpu = $#samples+1 if ($cpu > @samples);
	my $cmd;
	if ($OPTS{run} eq "multi" && $cpu > 1)
	{
		$cmd = "perl $src/multi-process.pl -cpu $cpu $shell";
	}
	elsif ($OPTS{run} eq "multi" && $cpu == 1)
	{
		$cmd = "sh $shell 1>$shell.log";
	}
	elsif ($OPTS{run} eq "qsub" && $cpu == 1)
	{
		$cmd = "qsub -cwd -S /bin/sh -sync y -q $queue $shell";
	}
	elsif ($OPTS{run} eq "qsub")
	{
		$cmd = "perl $conf->{soft}->{'qsub-sge'} --queue=$queue --convert no --resource vf=${mem}G --maxjob $cpu $shell";
	}
	elsif ($OPTS{run} eq "no")
	{
		return 1;
	}
	else 
	{
		ERROR("Your -run set is error, either be 'multi' or 'qsub'");
	}

	timeLOG("create shell file, [$shell], runing ... ");
	timeLOG("CMD: $cmd");
	system($cmd);
}

sub run_shell_cmd
{
	my ($shell,$cpu,%opts) = @_;
	my $queue = $opts{'-queue'} || $conf->{'queue'};

	$cpu = $#samples+1 if ($cpu > @samples && $shell !~ /blast/);
	my $cmd;
	
	if ($OPTS{run} eq "multi" && $cpu > 1)
	{
		$cmd = "perl $src/multi-process.pl -cpu $cpu $shell";
	}
	elsif ($OPTS{run} eq "multi" && $cpu == 1)
	{
		$cmd = "sh $shell 1>$shell.log";
	}
	elsif ($OPTS{run} eq "qsub" && $cpu == 1)
	{
		$cmd = "qsub -cwd -S /bin/sh -sync y $shell";
	}
	elsif ($OPTS{run} eq "qsub")
	{
		$cmd = "perl $conf->{soft}->{'qsub-sge'} --maxjob $cpu --convert no --queue $queue $shell";
	}
	elsif ($OPTS{run} eq "no")
	{
		return 1;
	}
	else 
	{
		ERROR("Your -run set is error, either be 'multi' or 'qsub'");
	}
	
	return $cmd;
}

sub read_tags_num
{
	my $file = shift ;
	my %tags_num;
	
	open my $fh , $file or die "can't open stat_tag.xls file $!";
	
	<$fh>;
	while(<$fh>)
	{
		chomp;
		my ($sample,$num) = (split /\t/)[0,1];
		$tags_num{$sample} = $num;
	}

	close $fh;

	return \%tags_num;
}

sub usage
{
	print <<HELP;

Usage:   perl $0 --conf <*.conf> [--steps] [--run]

Options: --conf  STR    the config file, needed
         --steps STR    the steps you want to run, default all, see the README.md for detail
         --run          how to run the shell created by the program, no/qsub/multi [no]
         --help         print this infomation

HELP
	exit;
}
