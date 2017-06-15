#!usr/bin/perl -w
use strict;
use Getopt::Long;

#-------------------------------------------------------------------+
#   Program name: Trim_seq.pl                                       |
#   Author      :                                                   |
#   Release date: 2017/3/2                                          |
#   Contact     : ffhuang@genedenovo.com                            |
#-------------------------------------------------------------------+

=pod

=head1 USAGE
  
     Usage : perl Trim_seq.pl [options] <fastq>

        necessary:
            -name     <string>   Sample name
        optional :
            -lqual     <int>     The low quality value filtering,default=3
            -clen     <int>      The max length of Continuous low quality ,default=3
            -minp     <int>      The min percentage  of the length of Continuous high quality  than the length of the sequence ,                                  default = 75
            -hqual     <int>      The high quality value,default=20
            -maxl     <int>      Maximum length of the tag,default=500
            -minl     <int>      Minimum length of the tag,default=300
            -phread   <int>      <33|64>,default=33
            -outfile  <string>   default "./out.fasta"
            -h                   output help information to screen

=head1 Example
    
    perl Trim_seq.pl -name D7 -outfile ./out.fasta fastq

=cut

my ( $low_qual, $high_qual, $continuous_length_cutoff, $max_length_of_tag,
    $min_percent, $min_length_of_tag, $phread, $outfile, $name, $help );

GetOptions(
    "lqual:i"   => \$low_qual,
    "hqual:i"   => \$high_qual,
    "conl:i"    => \$continuous_length_cutoff,
    "minp:i"    => \$min_percent,
    "maxl:i"    => \$max_length_of_tag,
    "minl:i"    => \$min_length_of_tag,
    "phread:i"  => \$phread,
    "name:s"    => \$name,
    "outfile:s" => \$outfile,
    "h"         => \$help,
);

#### default settings
$low_qual                 ||= 3;
$high_qual                ||= 30;
$continuous_length_cutoff ||= 3;
$max_length_of_tag        ||= 500;
$min_length_of_tag        ||= 300;
$phread                   ||= 33;
$outfile                  ||= "./out.fasta";
$min_percent              ||= 75;

my ( $fastq ) = @ARGV;

die `pod2text $0` if $help;
die `pod2text $0` if ( !$fastq || !$name );

############   main

open my $fastq_h, "<:gzip", $fastq or die "Can't open $fastq\n";
open my $out_h,  ">", $outfile   or die "Can't open $outfile\n";

my $tag_name = "t0000001";
FIRST:
while ( my $name1 = <$fastq_h> ) {
    chomp( my $seq = <$fastq_h> );
    <$fastq_h>;    #ignore
    chomp( my $qual = <$fastq_h> );

    if ( $name1 !~ /^@/ ) {
        die "ERROR:File format error\n";
    }

    my @qual = split //, $qual;
    my @ascii = map { ord($_) - $phread } @qual;
    my @encode = map { $_ < $low_qual ? 0 : $_ > $high_qual ? 2 : 1 } @ascii;
    my $encode = join "", @encode;

    my $star = 0;             #Starting from 0
    my %effective;
    my $flag = 0;
SECOND:
    while ( $encode =~ /(0{$continuous_length_cutoff,})/g ) {   #0 is low qual
        $flag = 1;
        my $low_seq_len = length $1;           #the length of low quality
        my $low_end     = pos($encode) - 1;    #Starting from 0
        my $low_star = $low_end - $low_seq_len + 1;    # Starting from 0
        my $len = ( $low_star - 1 ) - $star + 1; ##the length of trim sequence
        my $sequence = substr( $seq, $star, $len );
        my $sub_code = substr( $encode, $star, $len );
        next SECOND if $sequence eq "";

        if ( $len < $max_length_of_tag && $len > $min_length_of_tag )
        {                                        # filter length
            my $boolean = Check( $len, $sub_code, $min_percent );
            if ($boolean) {
                $flag = 2;
                $effective{$sequence} = $len;
            }
        }
        $star = $low_end + 1;
    }
    if ( $flag == 0 ) {
        my $len = length $seq;
        my $boolean = Check( $len, $encode, $min_percent );
        if ($boolean) {
            print $out_h ">$tag_name\_$name\n$seq\n";
            $tag_name++;
        }
    }
    elsif ( $flag == 2 ) {
        my @max_len =
            sort { $effective{$b} <=> $effective{$a} } keys %effective;
        my $max_seq = $max_len[0];
        print $out_h ">$tag_name\_$name\n$max_seq\n";
    }
}

close $fastq;
close $out_h;

#########  sub

sub Check {
    my ( $len, $code, $min_percent ) = @_;
    my $qual_len = int( $min_percent * $len / 100 );
    if ( $code =~ m/2{$qual_len,}/g ) {
        return 1;
    }
    else {
        return 0;
    }
}

