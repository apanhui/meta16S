#!usr/bin/perl -w
#-------------------------------------------------
#  Program name:
#  Release date: 2017-6
#  Contact     : ffhuang@genedenovo.com
#-------------------------------------------------
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename dirname);
use Data::Dumper;

use lib "/Bio/User/aipeng/bin/meta16S/lib";
use FindBin;
$FindBin::RealBin = "/Bio/User/aipeng/bin/meta16S";
use DEBUG;
use CONF;

my %OPTS;
GetOptions(\%OPTS,
	"dist:s",
	"alpha:s",
	"outdir:s",
	"conf:s",
	"help",
	);

$OPTS{outdir} ||= ".";
$OPTS{outdir} = abs_path($OPTS{outdir});
`mkdir $OPTS{outdir}` unless -e $OPTS{outdir};

## read conf
my $conf = load_conf($OPTS{conf});

## extract group information
my $group = $conf->{groups};
for my $order (keys %{$group}){
	my ($g,$s)=split /:/,$group->{$order};
	my @sam = split /,/,$s;
	@{ $conf->{Group}{$g} } = @sam;
}

## filter
## Set the number of repeat must be greater than 3
my @twodiff = split /,/,$conf->{two_groups_diff};
my @mutildiff = split /,/,$conf->{multi_groups_diff};
$conf->{two_groups_diff} = Check_diff(\@twodiff);
$conf->{multi_groups_diff} = Check_diff(\@mutildiff);

my $ref = $conf->{two_groups_diff};
if (ref($ref) eq "ARRAY"){
	my $tmp = "two_groups_diff : " . (join ",",@{ $conf->{two_groups_diff} });
	timeLOG($tmp);
}
else{
	my $tmp = "warning :: two_groups_diff $conf->{two_groups_diff} !!";
	timeLOG($tmp);
}

$ref = $conf->{multi_groups_diff};
if (ref($ref) eq "ARRAY"){
	my $tmp =  "multi_groups_diff: " . (join ",",@{ $conf->{multi_groups_diff} }) ."\n";
	timeLOG($tmp);
}
else{
	my $tmp = "warning :: multi_groups_diff $conf->{multi_groups_diff} !!";
	timeLOG($tmp);
}


#print Dumper($conf);


if ($OPTS{dist}){
	Read_Dist($OPTS{dist});
}
if ($OPTS{alpha}){
	Read_alpha($OPTS{alpha});

}

## sub program
#------------------------------------------
#  unweighted_unifrac or weighted_unifrac
#------------------------------------------

sub Read_Dist{
	my ($dist)=@_;

	my %Group;
	my %tmp = %{ $conf->{Group} };
	for my $g (keys %tmp){
		my @sam = @{ $tmp{$g} };
		for my $s (@sam){
			$Group{$s} = $g;
		}
	}

	open my $fh,"<",$dist or die $!;
	my %DistInfo;
	my %index;
	timeLOG("Start read file : $dist");
	while (<$fh>){
		chomp;
		my @tmp = split /\t/;
		if ($.==1){
			for my $i ( 1 .. $#tmp ) {
				my $sam = $tmp[$i];
				$index{$i} = $sam;
			}
		}
		else{
			for my $i ( 1 .. $#tmp ) {
				my $head = $index{$i};
				next if $head eq $tmp[0];
				my $group1 = $Group{ $tmp[0] };
				my $group2 = $Group{ $head };
				next if $group1 ne $group2;
				next if $DistInfo{$group1}{"$tmp[0]_$head"};
				next if $DistInfo{$group1}{"${head}_$tmp[0]"};

				$DistInfo{$group1}{"$tmp[0]_$head"} = $tmp[$i];
			}
		}
	}
	close $fh;

	my $file;
	 # two group
	my $ref = $conf->{two_groups_diff};
	if (ref($ref) eq "ARRAY"){
		my @twodiff = @{ $conf->{two_groups_diff} };
		$file = extract_diff_dist(\@twodiff,\%DistInfo);
		timeLOG("Running : Rscript $FindBin::Bin/Ttest-wilcox.R $file $OPTS{outdir}");
		system("Rscript $FindBin::Bin/Ttest-wilcox.R $file $OPTS{outdir}");
	}

	 # mutil group
	$ref = $conf->{multi_groups_diff};
	if (ref($ref) eq "ARRAY"){
		my @multidiff = @{ $conf->{multi_groups_diff} };
		$file = extract_diff_dist(\@multidiff,\%DistInfo);
		timeLOG("Running : Rscript $FindBin::Bin/Kruskal-TukeyHSD.R $file $OPTS{outdir}");
		system("Rscript $FindBin::Bin/Kruskal-TukeyHSD.R $file $OPTS{outdir}");
	}

}

sub extract_diff_dist{
	my ($aref,$DistInfo) = @_;
	my $file = "$OPTS{outdir}/file.txt";
	open my $ofh_file, ">", $file or die $!;
	for my $diff (@{$aref}){
		my @diff = split /&/, $diff;
		$diff =~ s/&/-VS-/g;
		open my $ofh, ">", "$OPTS{outdir}/${diff}.xls" or die $!;
		print $ofh_file "$OPTS{outdir}/${diff}.xls\n";
		print $ofh "Group\tSample\tIndex\n";

		for my $g (@diff) {
			my %sam = %{ $DistInfo->{$g} };
			for (keys %sam){
				print $ofh "$g\t$_\t$DistInfo->{$g}{$_}\n";
			}
		}
	}
	return $file;
}


#---------------------------------------------
#  alpha diverisity
#--------------------------------------------
sub Read_alpha{
	my ($alpha)=@_;
	open my $fh,"<",$alpha or die $!;
	my %index;
	my %AlphaInfo;
	timeLOG("Start read file :$alpha");
	while(<$fh>){
		chomp;
		my @tmp = split /\t/,$_;
		if ($.==1){
			for my $i (1..$#tmp){
				my $alpha = $tmp[$i];
				`mkdir $OPTS{outdir}/$alpha` unless -e "$OPTS{outdir}/$alpha";
				$index{$i}=$alpha;
			}
		}
		else{
			for my $i (1..$#tmp){
				my $head = $index{$i};
				$AlphaInfo{$head}{$tmp[0]}=$tmp[$i];
			}
		}
	}
	timeLOG("Competed read file :$alpha");
###  get diff info
	for my $alpha_index (keys %AlphaInfo){

		 # two group
		my $ref = $conf->{two_groups_diff};
		my $file;
		if (ref($ref) eq "ARRAY"){
			my @twodiff = @{ $ref };
			$file = extract_diff(\@twodiff,$alpha_index,\%AlphaInfo);
			timeLOG("Running : Rscript $FindBin::Bin/Ttest-wilcox.R $file $OPTS{outdir}/$alpha_index");
			system("Rscript $FindBin::Bin/Ttest-wilcox.R $file $OPTS{outdir}/$alpha_index");
		}

		 # mutil group
		$ref = $conf->{multi_groups_diff};
		if (ref($ref) eq "ARRAY" ){
			my @multidiff = @{ $conf->{multi_groups_diff} };
			$file = extract_diff(\@multidiff,$alpha_index,\%AlphaInfo);
			timeLOG("Running : Rscript $FindBin::Bin/Kruskal-TukeyHSD.R $file $OPTS{outdir}/$alpha_index");
			system("Rscript $FindBin::Bin/Kruskal-TukeyHSD.R $file $OPTS{outdir}/$alpha_index");
		}

	}
}

sub extract_diff{
	my ($aref,$type,$AlphaInfo)=@_;
	my $file = "$OPTS{outdir}/$type/file.txt";
#	print "$file\n";
	open my $ofh_file,">",$file or die $!;
	for my $diff (@{$aref}){
		my @diff = split /&/,$diff;
		$diff=~s/&/-VS-/g;
		open my $ofh,">","$OPTS{outdir}/${type}/${diff}.xls" or die $!;
		print $ofh_file "$OPTS{outdir}/${type}/${diff}.xls\n";
		print $ofh "Group\tSample\tIndex\n";

		for my $g (@diff){
			my @sam = @{ $conf->{Group}{$g} };
			for my $s (@sam){
				print $ofh "$g\t$s\t$AlphaInfo->{$type}{$s}\n";
			}
		}
		close $ofh;
	}
	return $file;
}

#-------------------------------
#  filter
#-------------------------------

sub Check_diff{
	my ($aref)=@_;
	my @diff;
	for (@{$aref}){
		my @tmp = split /&/,$_;
		my $flag = Check(\@tmp);
		if ($flag==0){
			print STDERR "Warning::$_,repeat number less than 3, no difference analysis!\n";
		}
		else{
			push @diff,$_;
		}
	}

	if (scalar @diff ==0){
		my $tag = "No comparison";
		return $tag;
	}
	else{
		return \@diff;
	}
}

sub Check{
	my ($aref)=@_;
	my $flag=1;
	for my $g (@{$aref}){
		if (scalar @{ $conf->{Group}{$g} }<3){
			$flag=0;
			last;
		}
	}
	return $flag;
}
