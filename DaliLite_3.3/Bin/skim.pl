#!/usr/bin/perl
#______________________________________________________________________________
# Title     : skim.perl
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
# skim.perl
# puts first line in file1, rest in file2
#
($file1,$file2)=@ARGV;
$i=0;
while (<STDIN>) {
	$i++;
	if($i == 1) { open(OUT,"> $file1"); }
	elsif($i == 2) { close(OUT); open(OUT,"> $file2"); }
	print OUT $_;
}
close(OUT);
