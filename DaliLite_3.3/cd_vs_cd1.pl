# get chains of newpdb-entries

# usage: cd_vs_cd1.pl cdlist < cd1list > subset_of_cd1list

use strict;

my($cdlist)=@ARGV;
my %x;

# memorize cd
open(IN,"<$cdlist") || die "Can't open $cdlist\n";
while(<IN>) {
	my($cd)=/^(\w+)/;
	$x{$cd}=1;
}
close(IN);

# read cd1list from STDIN, pass through matches to cdlist
while(<STDIN>) {
	my($cd)=/^(\w{4})/;
	if($x{$cd} > 0) { print $_; }
}

exit();
