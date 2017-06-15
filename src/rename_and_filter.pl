#!usr/bin/perl -w
use strict;
open IN,"$ARGV[0]" or die $!;
my $outdir=$ARGV[1];
my $sample=$ARGV[2];
my $len_min=$ARGV[3];
my $len_max=$ARGV[4];
my $fix_num=$ARGV[5];
my $cut_num=$ARGV[6];

chomp $sample;
open OUT,">$outdir/$sample.fasta" or die $!;
open OUT2,">$outdir/$sample.list" or die$!;
open OUT3,">$outdir/$sample.groups" or die$!;

my $random=1+rand(0.15);

my $i="t0000001";    my $co=1;
while(<IN>)
{
	chomp(my $id=$_);
	chomp(my $seq=<IN>);
	if(length($seq)>$len_min && length($seq)<$len_max)
	{
		next if($fix_num>0 && $co>$fix_num);
		next if($cut_num>0 && $co>$cut_num*$random);
		print OUT ">$i\_$sample\n$seq\n";
		print OUT2 "$id\t$i\n";
		print OUT3 "$i\_$sample\t$sample\n";
		$i++;
		$co++;
	}
}
close IN;
close OUT;
close OUT2;
close OUT3;
