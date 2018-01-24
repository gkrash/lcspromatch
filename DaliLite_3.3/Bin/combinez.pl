# read ranklist_sse ranklist_gtg ranklist_wolf, reorder by max(z)
#>> filter to pdb90.list

use strict;

my %best;

&readlist('ranklist_gtg',1,3,1);
&readlist('ranklist_wolf',4,0,0);

# output combined max(z)-list
my $cd2;
foreach $cd2 (keys %best) {
	print join("\t",$cd2,$best{$cd2}),"\n";
}

sub readlist {
	my($filename,$cd2col,$zcol,$addbias)=@_;
	open(IN,"<$filename") || warn "Can't open $filename\n";
	my $i=0;
	while(<IN>) {
		chomp;
		my(@data)=split(/\t/);
		my $cd2=$data[$cd2col];
		next if ($cd2=~/\?/);
		my $z=$data[$zcol];
		if($addbias>0) { $z+=&gtg_bias($i); }
		if($best{$cd2}<$z) { $best{$cd2}=$z; }
		$i++;
	}
	close(IN);
}

sub gtg_bias {
	# add bias to gtg z-scores, because they are good hits!
	my($rank)=@_;
#	my $bias=100-$rank;
	my $bias=4; # add constant
	if($bias<0) { $bias=0; }
	return($bias);
}

