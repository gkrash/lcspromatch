# divide STDIN into 8 files new_[1-8].list

use strict;
my $i=0;

open(OUT1,">new_1.list");
open(OUT2,">new_2.list");
open(OUT3,">new_3.list");
open(OUT4,">new_4.list");
open(OUT5,">new_5.list");
open(OUT6,">new_6.list");
open(OUT7,">new_7.list");
open(OUT8,">new_8.list");

my $old_pdbid='?';

while(<STDIN>) {
	my($pdbid,$chainid)=/^(\w{4})(\w*)/;
	if($pdbid ne $old_pdbid) { $i++; }
	my $cd1=$pdbid.$chainid;
	$old_pdbid=$pdbid;
	if($i>8) { $i=1; }
	if($i==1) { print OUT1 "$cd1\n"; }
	elsif($i==2) { print OUT2 "$cd1\n"; }
        elsif($i==3) { print OUT3 "$cd1\n"; }
        elsif($i==4) { print OUT4 "$cd1\n"; }
        elsif($i==5) { print OUT5 "$cd1\n"; }
        elsif($i==6) { print OUT6 "$cd1\n"; }
        elsif($i==7) { print OUT7 "$cd1\n"; }
        elsif($i==8) { print OUT8 "$cd1\n"; }
	else { warn "#ERROR $cd1\n"; }
}

close(OUT1);
close(OUT2);
close(OUT3);
close(OUT4);
close(OUT5);
close(OUT6);
close(OUT7);
close(OUT8);

exit();

