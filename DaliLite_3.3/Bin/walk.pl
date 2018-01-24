use strict;

# global hashes; key is cd2/3
my %lali_induced;
my %lali_observed;
my %inducer;
my %maxtried_lali;
my %map21;
my %map23;
my %z12;
my %z12_domain; # domain mapping of match
my %z23;
my %iteration; # iteration of improved z12

use lib '/data/backup/zope_data/dali_server/DaliLite_3.0/Bin/'; # HARD-CODED!
use FSSP;

if($#ARGV<6) { die "Usage: $0 cd1 dccpfile zcutoff DALIDATDIR_1 DALIDATDIR_2 DALI_SERVER_HOME SQL_DATABASE\n"; }

my($cd1,$dccpfile,$zcutoff,$DALIDATDIR_1,$DALIDATDIR_2,$DALI_SERVER_HOME,$SQL_DATABASE)=@ARGV;
my $MAX_HITS=1000; # max size of first shell; walk to twice outputted amount (500) so as to drop fewer hits with medium Z-scores
my $MAX_DALICON=$MAX_HITS; # max number of comparisons performed
my $MIN_LALI=25;
my $BATCH_SIZE=100; # small batch-size => largest domain dominates!
my $tmp_dccpfile='tmp.dccp';
my $TMP_DALICON='tmp.dalicon';
my $MAX_ITER=50;
my $LFITZ='T';
my $MAX_RADIUS=1000000; # max size of second shell; stop search when number of induced structures exceeds this 

my $table="$SQL_DATABASE\.dali_dccp";
my $fragment_table="$SQL_DATABASE\.dali_segments"; # use dali_segments.from1,from2,blocklen

# dalicon requires dali.default
system "$DALI_SERVER_HOME\/dalidefault.csh";
if (!-e 'dali.default') { die "dali.default missing!\n"; }

FSSP::Connect();

my $iteration=0;
my $ndalicon=0;

# domain structure
my($ndom,@domain_mapping)=&domain_mapping();
print "# $ndom domains: @domain_mapping\n"; 

# load existing cd1.dccp
my($nhit)=&read_dccp($dccpfile); # returns map21{cd2}, lali_observed{cd2}, maxtried_lali{cd2}, $z12{cd2}
print "# $nhit hits returned\n";

# MAIN LOOP
# terminate if more than MAXHITS
# select cd3 which has largest lali_induced{cd3}-lali_observed{cd3}
# 	generate cd1-cd3 alignment 
# 	add cd3 as cd2
# 	induce cd3

my $n_cd3=$nhit;
while($n_cd3>0) {
	$iteration++;

	# pool domain lists
	my $idom;
	my %xx;
###	foreach $idom (1..$ndom) { 
		my(@y)=&sortlist_z($idom);
		my $n=$#y;
		if($n>$BATCH_SIZE) { $n=$BATCH_SIZE; }
		foreach(@y[0..$n]) { $xx{$_}=1; }
###	}
	my(@y)=keys %xx;
	print "# keys: @y\n";
	my $m=$#y;
	last if($m<0);
	my $m0=$m;
	if($m0>$BATCH_SIZE) { $m0=$BATCH_SIZE; }
	my(@keys)=@y[0..$m0];

	print "# calling run_dalicon with list: @keys\n";
	($n_cd3)=&run_dalicon(@keys);
	last if($n_cd3<1);
	print "# $n_cd3 seed alignments sent to dalicon\n";
	$ndalicon+=$n_cd3;
	my($n)=&read_dccp($tmp_dccpfile);
	# nhit := number of unique neighbours
	my(@x)=keys %lali_observed;
	$nhit=$#x+1;
	print "# $n hits returned from dalicon (total $nhit hits / $m radius) iteration $iteration \n";

	system "cat $tmp_dccpfile >> $dccpfile";

	if($nhit>=$MAX_HITS) { 
		print "# Terminating walk at $nhit hits > MAX_HITS=$MAX_HITS\n";
		$n_cd3=0; # to exit loop
	}

	if($ndalicon>=$MAX_DALICON) {
		print "# Terminating walk at $ndalicon dalicon runs > MAX_DALICON=$MAX_DALICON\n";
		$n_cd3=0; # to exit loop
	}

	if ($iteration>$MAX_ITER) {
	        print "# Terminating walk at $iteration iterations > MAX_ITER=$MAX_ITER\n";
		$n_cd3=0; # to exit loop
        }	       

	if($m>$MAX_RADIUS) {
		print "# Terminating walk at $m induced neighbours > MAX_RADIUS=$MAX_RADIUS\n";
		$n_cd3=0; # to exit loop
	}
}

FSSP::DisConnect();

&print_stack($iteration);

# clean up
system "rm -f $TMP_DALICON $tmp_dccpfile";


########################## SUBROUTINES ########################################

sub sortlist_z { #>> select best per each z12_domain
	my($idom)=@_;
	my %x;
	foreach (keys %lali_induced) { 
		#warn "# sortlist_z $idom key: $_ domain: $z12_domain{$inducer{$_}}\n"; 
###		next if($z12_domain{$inducer{$_}} != $idom);
		next if($lali_induced{$_}<$MIN_LALI);
		my $x=$z12{$inducer{$_}};
		my $y=$z23{$_};
		if ($y<$x) { $x=$y; } # min z-score
		next if($x<$zcutoff); # poor hit
		if( ($lali_induced{$_} > $maxtried_lali{$_}) && ($lali_induced{$_} > $lali_observed{$_}) ) { $x{$_}=$x; }
	}
	my(@keys)=sort { $x{$b} <=> $x{$a} } keys %x;
	return(@keys);	
}

sub print_stack {
	my($iteration)=@_;
	my %x;
	foreach (keys %z12,keys %z23) { $x{$_}=1; }
	my $cd3;
	my $rank=0;
	foreach $cd3 (keys %z23) {
		$rank++;
		my $cd2=$inducer{$cd3};
		print join("\t",$iteration,$rank,$cd3,$iteration{$cd3},$z23{$cd3},$lali_induced{$cd3},$maxtried_lali{$cd3},$lali_observed{$cd3},$cd2,$z12{$cd3}),"\n";
	}
}

sub domain_mapping { 
	# read cd1.dat, return domain_no array [1..nres]
	my $long=$cd1;
	if(length($cd1)<5) { $long.='_'; }
	my $filename="$DALIDATDIR_1\/$long\.dat";
	open(IN,"<$filename") || die "Can't open $filename\n";
	my $n=0;
	my $nres=0;
	my $ndom=0;
	my $idom=0;
	my @mapping;
	while(<IN>) { 
		if(/^>>>>\s+\w+\s+(\d+)/) {
			$n++;
			if($n==1) { 
				$nres=$1; 
				foreach (0..$nres) { push(@mapping,0); } 
			} elsif($n==3) { 
				$ndom=$1; 
			}
		} 
		if( ($n==3) && /\*.{17}\s+(\d+.*)$/) {
			$idom++;
			my(@x)=split(/\s+/,$1);
			while($#x>0) {
				my $from=shift(@x);
				my $x;
				while ($from>999) {
					$_=$from;
					($from,$x)=/(\d+)(\d{4})/;
					unshift(@x,$x);
				}
				my $to=shift(@x);
				while($to>999) {
					$_=$to;
					($to,$x)=/(\d+)(\d{4})/;
					unshift(@x,$x);
				}
				foreach ($from..$to) { $mapping[$_]=$idom; }
			}
		}
		last if(/^-/);
	}
	close(IN);
	return($idom,@mapping);
}

sub run_dalicon {
	my(@cd3list)=@_;
        # do batch of 100 best
        my $cd3;
        open(OUT,">$TMP_DALICON");
        print OUT "$cd1\n";
        my $n_cd3=0;
        foreach $cd3 (@cd3list) {
		#last if($n_cd3>$BATCH_SIZE); 
		my $x=$lali_induced{$cd3}-$maxtried_lali{$cd3};

		# output cd3,z13, preali to 'tmp.dalicon'
		my($l13,$m13)=&generate_map13($map21{$inducer{$cd3}},$map23{$cd3});
		#warn "# generate_map13 returned $cd3 $l13 $m13\n";
		printf OUT "%-5.5s\*\n", $cd3;
		print OUT "$l13\n$m13";
		$n_cd3++;
                print "# dalicon iteration=$iteration cd3=$cd3 z12=$z12{$inducer{$cd3}} z23=$z23{$cd3} x=$x n_cd3=$n_cd3\n";
		if($l13>$maxtried_lali{$cd3}) { $maxtried_lali{$cd3}=$l13; }
        }
        print OUT "END\nEND\n";
        close(OUT);

        # run dalicon
        if($n_cd3>0) {
                        my $cmd="rm -f fort.91 ; cat $TMP_DALICON | perl $DALI_SERVER_HOME\/Bin/pipeforwarder.pl $DALIDATDIR_1 $DALIDATDIR_2 $LFITZ | $DALI_SERVER_HOME\/Bin/dalicon > /dev/null ";
                        #warn "$cmd\n";
                        system $cmd;
                        if(-e 'fort.91') {
                                my $cmd="perl $DALI_SERVER_HOME\/Bin/sortdccp.pl < fort.91 | $DALI_SERVER_HOME\/Bin/pipeforwarder.pl $DALIDATDIR_1 $DALIDATDIR_2 DCCP 2.0 10 | $DALI_SERVER_HOME\/Bin/dp > $tmp_dccpfile";
                                #warn "$cmd\n";
                                system $cmd;
                        }
        }
	return($n_cd3);
}

sub generate_map13 {
	my($map21,$map23)=@_;
	# generate induced alignment map13 from map21 x map23
	my %m21;
	my ($cd2starts,$cd1starts,$lengths)=split(/ /,$map21);
	my(@cd2starts)=split(/,/,$cd2starts);
	my(@cd1starts)=split(/,/,$cd1starts);
	my(@lengths)=split(/,/,$lengths);
	my $from2;
	my $i;
	foreach $from2 (@cd2starts) {
		my $from1=shift(@cd1starts);
		my $l=shift(@lengths);
		foreach $i (0..$l-1) { $m21{$from2+$i}=$from1+$i; }
	}
	my %m23;
	my ($cd2starts,$cd3starts,$lengths)=split(/ /,$map23);
	my(@cd2starts)=split(/,/,$cd2starts);
	my(@cd3starts)=split(/,/,$cd3starts);
	my(@lengths)=split(/,/,$lengths);
	my $from2;
	my $i;
	foreach $from2 (@cd2starts) {
		my $from3=shift(@cd3starts);
		my $l=shift(@lengths);
		foreach $i (0..$l-1) { $m23{$from2+$i}=$from3+$i; }
	}
	my $i2;
	my $l13=0;
	my @x;
	my @y;
	foreach $i2 (sort { $a <=> $b } keys %m21) {
		my $i1=$m21{$i2};
		my $i3=$m23{$i2};
		if($i3>0) { push(@x,$i1); push(@y,$i3); push(@x,$i1); push(@y,$i3); $l13++; } # ranges!
	}
	my $m13="@x\n@y\n"; # tmp.dalicon input!
	return($l13,$m13);
}

sub induce { 
	my($cd2)=@_;
	if($z12{$cd2}<$zcutoff) { return; } # too poor quality
	#print "# this is induce: $cd2\n";
	my $cmd="SELECT 1 AS db, 1 AS orderflag, cd2 AS cd3, zscore,align_id FROM $table WHERE zscore\>$zcutoff AND cd1=\'$cd2\' ORDER BY zscore DESC";
        my $cmd1="SELECT 1 AS db, 2 AS orderflag, cd1 AS cd3, zscore, align_id FROM $table WHERE zscore\>$zcutoff AND cd2=\'$cd2\' ORDER BY zscore DESC";
        FSSP::Execute($cmd);
        my $result=FSSP::BuildResultStructure();
        FSSP::Execute($cmd1);
        my $result1=FSSP::BuildResultStructure();
        my $row;
	my $n=$#$result+$#$result1;
        foreach $row (@$result,@$result1) {
                        my $cd3=$row->{cd3};
                        my $z23=$row->{zscore};
			my $alignid=$row->{align_id};
			my $orderflag=$row->{orderflag};

			# get map23
			my $cmd2="SELECT from1,from2,blocklen FROM $fragment_table WHERE align_id=$alignid";
			FSSP::Execute($cmd2);
			my $result2=FSSP::BuildResultStructure();
			my $row2;
			my @xstarts;
			my @ystarts;
			my @lengths;
			my $lali=0;
			foreach $row2 (@$result2) {
				push(@xstarts,$row2->{from1});
				push(@ystarts,$row2->{from2});
				push(@lengths,$row2->{blocklen});
				$lali+=$row2->{blocklen};
			}
			my $map23=join(' ',join(',',@xstarts), join(',',@ystarts), join(',',@lengths));
			my $l23=$lali;
			if($orderflag==2) { $map23=join(' ',join(',',@ystarts), join(',',@xstarts), join(',',@lengths)); }

			# generate induced alignment map13 from map21 x map23
			my($l13,$m13)=&generate_map13($map21{$cd2},$map23);

			# case: z12{$cd3} is defined => only keep if induced alignment is longer than maxtried
			# case: z12{$cd3} is undefined and z23{$cd3} is defined => only keep if induced alignment is longer than maxtried and current induced
			# case: z12{$cd3} is undefined and z23{$cd3} is undefined => keep induced alignment
			my $keep=0;
			if(defined($z12{$cd3})) {
				if($l13 > $maxtried_lali{$cd3}) { $keep=1; }
			} else {
				if(!defined($z23{$cd3})) { $keep=1; }
				else {
					if(($l13 > $maxtried_lali{$cd3}) && ($l13 > $lali_induced{$cd3})) { $keep=1; }
				}
			}
			# overwrite cd3 alignment
			if($keep) {
				$inducer{$cd3}=$cd2;
				$lali_induced{$cd3}=$l13;
				$map23{$cd3}=$map23;
				$z23{$cd3}=$z23;
			}
        }
}

sub read_dccp {
        my($dccpfile)=@_;
        open(IN,"perl $DALI_SERVER_HOME\/Bin/sortdccp.pl 1 < $dccpfile | ") || die "can't open perl $DALI_SERVER_HOME\/Bin/sortdccp.pl 1 < $dccpfile\n";
        my $n=0;
        my $block='';
        my $cd2;
	my %order;
	my %alignment12;
	my %z;
        while(<IN>) {
                if (/DCCP/) {
                        if($block ne '') { $alignment12{$cd2}=$block; }
                        my(@x)=split(/\s+/);
                        my $cdb=pop(@x);
                        my $cda=pop(@x);
                        $cd2=$cdb;
                        $order{$cd2}=1;
                        if($cdb eq $cd1) { $cd2=$cda; $order{$cd2}=2; }
                        my($nfrag)=pop(@x); pop(@x);
                        my($z)=pop(@x);
                        $z{$cd2}=$z;
                        if($z>=2) { $n++; }
                        $block="$nfrag\n";
			#warn "# read_dccp got cda=$cda cdb=$cdb z=$z cd2=$cd2 z12=$z12{$cd2} nhit=$n\n";
                } elsif(/alignment/) {
                        next;
                 } else {
                        # fix nres>1000 here, later split on whitespace
# dccp   1    384.7 6.1 139     0.6          78     21                mol1  1yqwr
# alignment
# 338   341 345   350 384   388 898   902 904   908 915   9191104  11081123  1128
#1129  11451150  11541155  11631166  11691170  11751206  12091210  12131304  1308
#1310  13131314  13201338  13471357  13751376  1379
#   8    11  34    39  42    46  54    58  59    63  64    68  82    86  87    92
#  95   111 112   116 118   126 163   166 168   173 198   201 265   268 290   294
# 300   303 309   315 416   425 426   444 446   449
                        s/^\s+//;
                        my $line=$_;
                        foreach (split(/\s+/,$line)) {
                                if($_>9999) {
                                        my($a,$b)=/^(\d+)(\d{4})$/;
                                        $block.="$a $b ";
                                } else {
                                        $block.="$_ ";
                                }
                        }
                }
        }
        close(IN);
        # last alignment
        if($block ne '') { $alignment12{$cd2}=$block; }

	# convert alignment12{cd2} and order{cd2} to map21{cd2}, lali_observed{cd2}, maxtried_lali{cd2}
	foreach $cd2 (keys %alignment12) {
		#warn "# read_dccp cd2=$cd2 z=$z{$cd2} z12=$z12{$cd2}\n";
		if($z{$cd2}>$z12{$cd2}) { # if z-score is better, then overwrite & induce
			$iteration{$cd2}=$iteration;
			# domain mapping
			my %vote;
			# build alignment		
			my $order=$order{$cd2};
			my @xstarts;
			my @ystarts;
			my @lengths;
			my $lali=0;
	                $_=$alignment12{$cd2};
	                my(@x)=split(/\s+/);
	                my($nfrag)=$x[0];
	                my $i=1;
	                while($i<=$nfrag*2) {
                                my $from1=$x[$i];
                                my $to1=$x[$i+1];
                                my $jres=$x[$i+$nfrag*2];
				my $length=$to1-$from1+1;
				if ($order==2) { 
					my $x=$from1;
					$from1=$jres;
					$jres=$from1;
				}
				push(@xstarts,$from1);
				push(@ystarts,$jres);
				push(@lengths,$length);
				$lali+=$length;
                                $i+=2;	
				my $j=$from1;
				foreach (0..$length) { 
					$vote{$domain_mapping[$j]}++;
					$j++;
				}
	                }
			# which domain got most votes?
			#my $maxdom=1;
			#foreach (2..$ndom) {
			#	if($vote{$_}>$vote{$maxdom}) { $maxdom=$_; }
			#}
			# export values in global variables 
			#warn "# votes: $vote{1}, $vote{2}, $vote{3}\n";
			$z12{$cd2}=$z{$cd2};
			#$z12_domain{$cd2}=$maxdom;
			$map21{$cd2}=join(' ',join(',', @ystarts), join(',', @xstarts), join(',', @lengths)); # @cd2starts @cd1starts @lengths
			$lali_observed{$cd2}=$lali;
			$lali_induced{$cd2}=$lali;
			if(!defined($maxtried_lali{$cd2})) { $maxtried_lali{$cd2}=0; }
			$inducer{$cd2}=$cd1; # because cd2 now has direct alignment to cd1
			#print "# read_dccp: cd2=$cd2 lali=$lali z12=$z12{$cd2}\n";

			# induce new neighbours using improved alignments to cd2 
			&induce($cd2);
		}
	}

        return($n);
}


