#!/Bio/bin/perl
#-----------------------------------------------+
#    [APM] This script was created by amp.pl    |
#    [APM] Created time: 2017-02-22 15:45:42    |
#-----------------------------------------------+
# name: uparse2mothur.pl
# func: turn uparse result to Muthor format 
# version:: 0.1

use strict;
use warnings;

my $otus = shift @ARGV;
my $tab = shift @ARGV;
my $newotus = shift @ARGV;
my $newtab = shift @ARGV;

# read the shared file 
open TAB,$tab or die "can't open file, $tab $!";
open OUT,">$newtab" or die "can't open file, $newtab $!";

my %names;
my $num;
chomp( my $header = <TAB> );

while(<TAB>)
{
	chomp;
	my @tmp = split /\t/;
	$tmp[0] = 0.03;
	
	unless ($num)
	{
		$num = length ($tmp[2]);
		my $prefix = '0' x ($num+1);
		my @head = split /\t/,$header;
		@head = map {
			if ($_ > 2)
			{
				my ($order) = $head[$_] =~ /(\d+)$/;
				my $len = length ($order);
				my $newid = "Otu" . substr($prefix,$len) . $order;
				$names{$head[$_]} = $newid;
				$newid;
			}
			else 
			{
				$head[$_]
			}
		} 0 .. $#head;
		
		print OUT ${[ join "\t",@head ]}[0] , "\n";
	}

	print OUT ${[ join "\t",@tmp ]}[0] , "\n";
}

close TAB;
close OUT;

# re write the OTU represstive sequence with new id 
open SEQ,$otus or die "can't open file, $otus $!";
open OUT,">$newotus" or die "can't open file, $newotus $!";

while(<SEQ>)
{
	if (/>(\S+)\n/)
	{
		print OUT ">$names{$1}\n";
	}
	else 
	{
		print OUT $_;
	}
}

close SEQ;
close OUT;
