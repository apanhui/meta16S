package General;

use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(timeLOG LOG printn WARN ERROR read_list);

sub timeLOG
{
	my $line = join "" , @_;
	my $time = `date +"%F %T"`;
	chomp $time;
	print STDERR "[$time] $line\n";
}

sub LOG
{
	my $line = join "" , @_;
	print STDERR $line , "\n";
}

sub WARN
{
	my $line = join "" , @_;
	print STDERR "WARN: $line\n";
}

sub ERROR
{
	my $line = join "" , @_;
	print STDERR "ERROR: $line\n";
	exit 1;
}

sub printn
{
	my $line = join "" , @_;
	print $line , "\n";
}

# open file error
sub OFE
{
	my $file = shift;
	die "OFE: can't open file '$file', $!";
}

# read list file
sub  read_list
{
	my $file = shift;
	open FH,$file or die;
	my @arr = <FH>;
	close FH;
	chomp @arr;
	return @arr;
}

