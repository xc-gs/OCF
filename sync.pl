#!/usr/bin/perl

#
# Author: Ahfyth
# Date  : Dec 12, 2018
#
# This program is designed to extract the cfDNA coverage information
# for tissue-specific open chromatin regions
#

use strict;
use warnings;

if( $#ARGV < 1 )
{
	print STDERR "\nUsage: $0 <in.ol> <out.prefix>\n";
	exit 2;
}

my %value;
my %total;
open IN, "$ARGV[0]" or die( "$!" );
while( <IN> )
{
	chomp;
	my @l = split /\t/;	##chr1	1014459	1016459	LABEL	chr1	1015296	1015412	116
	next if $l[4] eq '.' || $l[0]!~/^chr\d/;
	$l[1] ++;
	$l[5] ++;
	my $s = ( $l[5]>=$l[1] ) ? ($l[5]-$l[1]) : 0;
	my $e = ( $l[6]<=$l[2] ) ? ($l[6]-$l[1]+1) : ($l[2]-$l[1]+1);
	$value{ $l[3] }->[ $s ] ++;
	$value{ $l[3] }->[ $e ] --;
	$total{ $l[3] } += $e - $s;
}
close IN;

foreach my $k ( keys %total )
{
	open OUT, ">$ARGV[1].$k.sync.value" or die( "$!" );

	my $v = $value{$k};
	my $t = $total{$k}/10000;

	print OUT "#Locus\tValue\tFraction%%\n";
	$v->[0] = 0 unless defined $v->[0];
	my $j = $v->[0];
	my $c = $#$v>>1;
	print OUT join("\t", -$c, $j, $j/$t), "\n";
	foreach my $k ( 1..($#$v-1) )
	{
		$v->[ $k ] += $v->[ $k-1 ];
		print OUT join("\t", $k-$c, $v->[$k], $v->[$k]/$t), "\n";
	}
	close OUT;
}

