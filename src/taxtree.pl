#!/Bio/bin/perl
#-----------------------------------------------+
#    [APM] This script was created by amp.pl    |
#    [APM] Created time: 2017-05-12 14:06:10    |
#-----------------------------------------------+
# name: taxtree.pl
# func: draw taxtree with taxonomy profiling v2 file by SBV
# version: 0.1

use strict;
use warnings;

use List::Util qw/sum/;
use Getopt::Std;

my %opts = (k=>"temp",m=>0.5,r=>20);
getopts('g:k:m:r:h',\%opts);

&usage unless @ARGV == 1;
&usage if ($opts{h});

my $sbv = "/home/aipeng/work/develepment/SBV/bin/sbv.pl";

my $profiling_file = shift @ARGV;
my @groups = read_groups($opts{g});
my $leaves_num = 0;

&profiling2nhx($profiling_file);

my $height = $opts{r} * 2.5 * $leaves_num + 40;

&taxtree("$opts{k}.samples_taxtree.conf","$opts{k}.nhx","$opts{k}.samples.percent","$opts{k}.samples_taxtree");
&taxtree("$opts{k}.groups_taxtree.conf","$opts{k}.nhx","$opts{k}.groups.percent","$opts{k}.groups_taxtree") if ($opts{g});

sub profiling2nhx
{
	my $file = shift;
	my %relations;

	open my $fh , $file or die "can't open profiling file $file, $!";
	open my $ofh_nhx   , ">" , "$opts{k}.nhx" or die $!;
	open my $ofh_val   , ">" , "$opts{k}.samples.percent" or die $!;

	# create group percnet file handle and print the header info
	open my $ofh_group , ">" , "$opts{k}.groups.percent" or die $! if ($opts{g});
	if ($opts{g})
	{
		my @group_names = map { $_->[0] } @groups;
		print $ofh_group ${[ join "\t" , ("Class",@group_names,"total_percent") ]}[0] , "\n";
	}
	
	# read the first line of profiling file 
	my $header = <$fh>;
	chomp $header;
	my ($taxon_flag,$total_flag,@fields) = split /\t/,$header;
	my $sample_num = ($#fields+1)/2;
	my @samples = map { $fields[$_] =~ s/_tags//; $fields[$_] } 0 .. $sample_num-1;
	my %indexs = map { $samples[$_] => $_ } 0 .. $#samples;
	print $ofh_val ${[ join "\t" , ("Class",@samples,"total_percent") ]}[0] , "\n";
	
	# read the second line of profiling file (which defined the root Mirobe values)
	my $root_line = <$fh>;
	my ($root,$total_tags) = (split /\t/ , $root_line)[0,1];
	$relations{$root}{father} = undef;

	while(<$fh>)
	{
		chomp;
		my ($taxon,$total,@values) = (split /\t/)[0,1,2..$sample_num+1];
		my @array = split /;/ , $taxon;
		my $percent = sprintf "%.2f" , $total*100/$total_tags;
		next if ($percent < $opts{m});

		# save the nodes of 1-7 levels taxonomy
		foreach my$i(1 .. $#array)
		{
			my $child  = $array[$i];
			my $father = $array[$i-1];
			$relations{$child}{father} = $father;
			$relations{$father}{children}{$child} = 1;
		}

		# create samples' perentage file
		my $last = $#array == 1 ? $array[1] : "$array[-2];$array[-1]";
		@values = map { sprintf "%.2f" , $_*100/$total } @values;
		my $line = join "\t" , ($last,@values,$percent);
		print $ofh_val $line , "\n";
		
		# create groups' perentage file
		if ($opts{g})
		{
			my @group_vals;
			foreach my $group (@groups)
			{
				my ($group_name,@samples) = @$group;
				my @percents = map { $values[$indexs{$_}] } @samples;
				push @group_vals , sum(@percents);
			}
			print $ofh_group ${[ join "\t" , ($last,@group_vals,$percent) ]}[0] , "\n";
		}
	}
	
	print_nhx($root,\%relations,$ofh_nhx,1);
	
	close $fh;
	close $ofh_nhx;
	close $ofh_val;
	close $ofh_group if ($opts{g});
}

sub print_nhx
{
	my ($node,$hash,$ofh_nhx,$isLast) = @_;
	
	if (defined $$hash{$node}{father})
	{
		if (! $$hash{$node}{children}) # leaf 
		{
			$leaves_num ++;
			if ($isLast)
			{
				print $ofh_nhx $node;
			}
			else 
			{
				print $ofh_nhx "$node,";
			}
		}
		else # non-leaf nodes 
		{
			my @children = keys %{$$hash{$node}{children}};
			print $ofh_nhx "(";
			foreach my$i(0 .. $#children)
			{
				print_nhx($children[$i],$hash,$ofh_nhx,$i==$#children);
			}
			if ($isLast)
			{
				print $ofh_nhx ")$node";
			}
			else 
			{
				print $ofh_nhx ")$node,";
			}
		}
	}
	else 
	{
		my @children = keys %{$$hash{$node}{children}};

		print $ofh_nhx "(";
		foreach my$i(0 .. $#children)
		{
			print_nhx($children[$i],$hash,$ofh_nhx,$i==$#children);
		}
		print $ofh_nhx ")$node;";
	}
}

sub read_groups
{
	my $group_file = shift ;

	return () unless $group_file;
	
	my @groups;

	open my $fh_group , $group_file or die $!;
	while(<$fh_group>)
	{
		chomp;
		my ($group_name,$samples) = $_ =~ /(\S+)\s*=\s*(\S+)/;
		my @samples = split /,/ , $samples;
		push @groups , [$group_name,@samples];
	}
	close $fh_group;
	
	return @groups;
}

sub taxtree
{
	my ($conf_file,$nhx_file,$percent_file,$outname) = @_;

	my $conf = <<CONF;
width  = 1200
height = $height
margin = 20

<taxtree>

file = $nhx_file
percent = $percent_file

radius = $opts{r}

<<include legend_taxtree.conf>>
</taxtree>
<<include colors.conf>>
<<include styles/styles.taxtree.conf>>
CONF

	open my $ofh_conf , ">" , $conf_file or die $!;
	print $ofh_conf $conf;
	close $ofh_conf;
	
	system("perl $sbv taxtree --strict --conf $conf_file --out $outname");
	system("/Bio/usr/bin/convert $outname.svg $outname.png");
}

sub usage
{
	print <<HELP;
Usage:   perl $0 [options] <profiling.v2.xls>

Options: -g STR    set the groups define file, optional
         -k STR    set the output file keyname, default is 'temp'.
         -m FLOAT  set the minmum taxonomy percent to be displayed, default is 1.
         -r INT    set the radius size of the pie, default is 20
         -h        print the help info.
HELP
	exit;
}
