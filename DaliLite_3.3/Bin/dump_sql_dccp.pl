# dump_sql_dccp.pl 1ppt 
# STDOUT: @cd1list (cd1 are representatives!)
# 1ppt.dccp, 1pptA.dccp

use lib '/data/backup/zope_data/dali_server/DaliLite_3.0/Bin/';
use FSSP;
use strict;

my($query)=@ARGV;

my $dccp_table='dalidb_dccp';
my $segments_table='dalidb_segments';

# connect to SQL database
FSSP::Connect();

my @queries;
# find representative of query == best match
# case: query is pdbid, unspecified chain(s)
if(length($query)==4) {
	my $cmd="SELECT DISTINCT(cd2) FROM $dccp_table WHERE SUBSTRING(cd2,1,4)=\"$query\" ";
	FSSP::Execute($cmd);
	my $result=FSSP::BuildResultStructure();
	my $row;
	foreach $row (@$result) { push(@queries,$row->{cd2}); }
} else { # chain specified in query
	push(@queries,$query);
}
foreach $query (@queries) {
  my $cmd="SELECT cd1 FROM $dccp_table WHERE cd2=\"$query\" ORDER BY zscore DESC LIMIT 1";
  #warn "$cmd\n";
  FSSP::Execute($cmd);
  my $result=FSSP::BuildResultStructure();
  my $row=shift(@$result);
  my $cd1=$row->{cd1};
  # error: query not found in dalidb
  if($cd1 eq '') { 
	print "ERROR: query $query not found in database\n"; }
  else { # write output to query.dccp
	my $dccpfile="$query\.dccp";
	print "# cd1=$cd1 query=$query\n";
	open(OUT,">$dccpfile");
	# get neighbours of cd1; they include query
	my $cmd="SELECT cd1,cd2,zscore,align_id FROM $dccp_table WHERE cd1=\"$cd1\" ORDER BY zscore DESC";
	#warn "$cmd\n";
	FSSP::Execute($cmd);
	my $result=FSSP::BuildResultStructure();
	my $row;
	foreach $row (@$result) {
		my $align_id=$row->{align_id};
		my $cd2=$row->{cd2};
		my $zscore=$row->{zscore};
		my $lali=$row->{lali};
		my $cmd1="SELECT from1,from2,blocklen FROM $segments_table WHERE align_id\=$align_id";
		#warn "$cmd1\n";
		FSSP::Execute($cmd1);
		my $result1=FSSP::BuildResultStructure();
		my @from1; my @to1; my @from2; my @to2;
		my $row1;
		my $lali=0;
		foreach $row1 (@$result1) { 
			my $l=$row1->{blocklen}; 
			my $from1=$row1->{from1};
			my $from2=$row1->{from2};
			$lali+=$l;
			push(@from1,$from1); 
			push(@to1,$from1+$l-1); 
			push(@from2,$from2); 
			push(@to2,$from2+$l-1);
		}
		my $nfrag=$#from1+1-$[;
		my $ali1='';
		my $ali2='';
		my $ifrag=0;
		while($#from1>-1) { 
			$ali1.=sprintf("%4d  %4d",shift(@from1),shift(@to1));
			$ali2.=sprintf("%4d  %4d",shift(@from2),shift(@to2));
			$ifrag++;
			if($ifrag==8) { $ali1.="\n"; $ali2.="\n"; $ifrag=0; }
		}
		if($ifrag>0) { $ali1.="\n"; $ali2.="\n"; }
		# DCCP format
		my $header=sprintf(" DCCP   1     99.9 0.0%4d  %6.1f           0   %4d                %-5.5s %-5.5s\n alignment\n",$lali,$zscore,$nfrag,$query,$cd2);
		print OUT $header;
		print OUT $ali1;
		print OUT $ali2;
	}
	close(OUT);
  }
}

FSSP::DisConnect();

exit();
