#!usr/bin/perl -w
use strict;
my $list_file=shift;
my $OtuLabel=shift;
open(IN,$list_file)||die "cannot open:$!";
my %otu_name;    my %otu;
my $label;     my $otu_count;
while(<IN>)
{
	chomp;   my @a=split(/\t/,$_);
	if($a[0] eq "label")
	{
		foreach my $out(2 .. $#a)
		{
			$otu_name{$out}=$a[$out];
		}
	}
	else
	{
		if($a[0] eq "$OtuLabel")
		{
			foreach my $out(2 .. $#a)
			{
				$otu{$otu_name{$out}}=$a[$out];
			}
			$label=$a[0];
			$otu_count=$a[1];
		}
	}
}
close IN;

my $map_files=shift;
my %pre_cluster;
open(IN,$map_files)||die"cannot open:$!";
my $mark=0;   my $key;
while(<IN>)
{
	chomp;   my @a=split(/\t/,$_);
	if($_=~ /^ideal_seq/ || $_ eq "")
	{
		$mark=1;
	}
	else
	{
		if($mark==1)
		{
			$key=$a[0];
			push @{$pre_cluster{$key}},$a[0];
			$mark=0
		}
		else
		{
			push @{$pre_cluster{$key}},$a[0];
		}
	}
}
close IN;

my $all_names=shift;
my %unique;
open(IN,$all_names)||die "cannot open:$!";
while(<IN>)
{
	chomp;    my @a=split(/\t/);
	$unique{$a[0]}=$a[1];
}
close IN;

my $all_count_file=shift;
my %unique_count;
open(IN,$all_count_file)||die "cannot open: $all_count_file $!";
while(<IN>)
{
	chomp;    my @a=split(/\t/);   $unique_count{$a[0]}=$a[1];
}
close IN;

my $all_fasta_file=shift;
open(IN,$all_fasta_file)||die "cannot open:$!";
my %fa;   my $key1;
while(<IN>)
{
	chomp;
	if(/^>/)
	{
		s/>//;    my @a=split(/\s+/,$_);   $key1=$a[0];
	}
	else
	{
		$fa{$key1}.=$_;
	}
}
close IN;

my $all_sample_reformat_file=shift;
open(OUT,">$all_sample_reformat_file")||die"cannot open:$!";
my %otu_unique_count;
print OUT "$label\n$otu_count\n\n";
foreach my $out(sort keys %otu)
{
	my @a=split(/,/,$otu{$out});
	print OUT "$out\n";
	foreach my $in(@a)
	{
		#print "-$in\n";
		foreach my $in1(sort keys %{ {map {$_ => 1} @{$pre_cluster{$in}}} })
		{
			#print "--$in1\n";
			$otu_unique_count{$out}{$in1}++;
			my @b=split(/,/,$unique{$in1});
			foreach my $in2(@b)
			{
		#		print "---$in2\n";
				print OUT "$in2\n";
			}
		}
	}
	print OUT "\n";
}
close OUT;

my $otu_representative_fasta=shift;
open(OUT,">$otu_representative_fasta")||die"cannot open:$!";
foreach my $out(sort keys %otu_unique_count)
{
	print OUT ">$out\n";
	foreach my $in(sort {$unique_count{$b}<=>$unique_count{$a}} keys %{$otu_unique_count{$out}})
	{
		print OUT "$fa{$in}\n";
		last;
	}
}
close OUT;
