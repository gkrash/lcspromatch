# convert column X to Z-scores
# optional dbsize overrides #points

use strict;

my($col,$dbsize)=@ARGV; # cols 0,1,2,.. which is normalized?

my(@lines)=<STDIN>;

my $n=0;
my $sum=0.0;
my $sum2=0.0;
foreach (@lines) {
	my(@data)=split(/\s+/);
	my $x=$data[$col];
	$sum+=$x;
	$sum2+=$x*$x;
	$n++;
}

# output original data and z-score
if($n<1) { $n=1; }
if($dbsize>$n) { $n=$dbsize; }
my $mean=$sum/$n;
my $stdev=sqrt(($sum2-$n*$mean*$mean)/$n);
if($stdev<1e-60) { $stdev=1e-60; }
foreach (@lines) {
	my(@data)=split(/\s+/);
	my $x=$data[$col];
	my $z=($x-$mean)/$stdev;
	print join("\t",$z,@data),"\n";
}

