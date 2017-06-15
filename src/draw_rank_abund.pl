#!/usr/bin/perl -w
#use strict;
my $real=shift;
my $outdir=shift;
open IN,"<$real" or die "can't open $real";
<IN>;
my %hash2;
my %hash1;
while (<IN>)
{
	chomp;
	$_=~s/\t$//;
	my @inf=split;
	my $i=0;
	my $num=@inf;
	$num=$num-1;
	for (3..$num)
	{
		$i++;
		if($inf[$_] != 0)
		{
			$hash2{$inf[1]}{$i}=$inf[$_];
			$hash1{$inf[1]}+=$inf[$_];
		}
	}
}
close IN;
my %sum;
open OUT,">$outdir/RankAbdDistr";
print OUT "sample\txl\tvalue\n";
for (keys %hash2)
{
        my $aa=$_;
        my $cc=$hash2{$_};
	my $num=0;
        for (sort{$hash2{$aa}{$b}<=>$hash2{$aa}{$a}} keys %$cc)
        {
		$num++;
		my $log=log($hash2{$aa}{$_}/$hash1{$aa})/log(10);
                print OUT "$aa\t$num\t$log\n";
        }
}
close OUT;
my $samplenum=keys %hash2;
my $col=int($samplenum/21)+1;
$Rscript = <<R;
library(ggplot2)
palette(colors())
m <- read.table("$outdir/RankAbdDistr", header=T, sep="\\t")
mycolor=round(seq(50,657,length.out=$samplenum),0)
max=round(max(m\$value))
min=round(min(m\$value))
lim=c(-5:0)
lim1=10^lim
ggplot(m, aes(x = xl, y= value, group = sample, color = sample )) + geom_line()+
ggtitle("Rank Abundance Distribution Curves")+
labs(x="Rank", y="Relative Abundance")+
theme(plot.title = element_text(size = 20))+
theme(axis.title.x = element_text(size = 15))+
theme(axis.title.y = element_text(size = 15))+
guides(col=guide_legend(ncol=$col))+
scale_y_continuous(breaks=c(-5:0),label=lim1)+
scale_color_manual(values = mycolor)+
theme_bw()
ggsave(filename="Rank_abd_Distr.png",path="$outdir",width=9, height=6)
ggsave(filename="Rank_abd_Distr.pdf",path="$outdir",width=9, height=6)
quit("no")

R
open OUT,"> $outdir/RankAbdDistr.R" || die $!;
print OUT "$Rscript\n";
close OUT;
system "R CMD BATCH $outdir/RankAbdDistr.R $outdir/RankAbdDistr.Rout";
print STDERR "Program End!\n";

