# input: list, subsetlist
# output: list filtered to those in subsetlist

use strict;

my($subsetlist)=@ARGV;

# memorize subset 
my %keep;
open(IN,"<$subsetlist");
while(<IN>) {
	my($cd1)=/^(\w+)/;
	$keep{$cd1}=1;
}
close(IN);

while (<STDIN>) {
	my($cd1)=/^(\w+)/;
	if($keep{$cd1}) { print $_; }
}

exit(0);

