#!/usr/bin/perl
#______________________________________________________________________________
# Title     : sortdccp.pl
# Usage     :
# Function  :
# Example   :
# Keywords  :
# Options   :
# Author    : jong@bio.cc,
# Category  :
# Returns   :
# Version   : 1.0
#------------------------------------------------------------------------------
#
# sortdccp_z.perl: read a DCCP file from STDIN, sort by Z-score, write to STDOUT
#

use strict;

my($nkeep)=@ARGV; # optional to filter nkeep best per cd1-cd2
if($nkeep<1) { $nkeep=999999; } # keep all

my $TRUE=1;
my $FALSE=0;
my @array;
my $first=$TRUE;
my $block="";
while(<STDIN>) {
#DCCP   1   2149.4 0.0 157    34.4         100      1                1mup  1mup
	if(/^ DCCP\s+(\d+)/) {
		next if($1!=1); # remove redundant DCCP   2 ... lines
		# save old block
		if(!$first) { push(@array,$block); } else { $first=$FALSE; }
		# initialize new block
		$block = $_;
	} else {
		$block .= $_;	# append to block
	}
}
# last line
push(@array,$block);

my $i;
my @datakeys;
foreach $i ($[..$#array) {
	$_=$array[$i]; 
	# sort on Dali-score
	#push(@datakeys,(/^ DCCP\s+\d+\s+([\-\.\d]+)/));
	# sort on Z-score
	push(@datakeys, (/^.{26}\s+([\-\.\d]+)/));
}

sub bydatakeys { $datakeys[$b] <=> $datakeys[$a]; }
my @sortarray=@array[sort bydatakeys $[..$#array]; 	# sort array on value

# output
my %seen;
foreach (@sortarray) { 
	my $block=$_;
	my(@x)=split(/\n/);
	my(@y)=split(/\s+/,$x[0]);
	my $cda=pop(@y);
	my $cdb=pop(@y);
	my $key="$cda\_$cdb";
	#warn "key=$key seen=$seen{$key}\n";
	$seen{$key}++;
	next if($seen{$key}>$nkeep);
	print $block; 
}



