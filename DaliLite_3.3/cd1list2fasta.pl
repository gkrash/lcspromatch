# grep sequence from DALIDATDIR

use strict;
my($DALIDATDIR)=@ARGV;

my $cd1;
while(<STDIN>) {
	my($cd1)=/^(\w+)/;
	print ">$cd1\n";
	my $long=$cd1;
	if(length($long)<5) { $long.='_'; }
	my(@lines)=`grep '^\-sequence' $DALIDATDIR\/$long\.dat`;
	$_=$lines[0];
        s/^\S+\s+\"//;
        $_=~s/[a-z]/C/; # disuplhides in DSSP sequence
        $_=~s/U/X/; # selenomet U triggers warning from blast
        print $_;
}

exit;
