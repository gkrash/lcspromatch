# wait in loop until all processes are gone ...

use strict;
use lib '/data/backup/zope_data/dali_server/DaliLite_3.0/Bin/'; # HARD-CODED!
use FSSP;



my $cmd="ps -fu wwwrun";

my $i=0;

#goto L1;


# wait for DaliLite -u to finish
while(1) {
	my(@x)=`$cmd`;
	my $n=0;
	foreach (@x) { if(/DaliLite -u/) { $n++; } }
	$i++;
	print "i=$i: $n DaliLite -u processes running\n";
	last if($n<=0); # ps is one process
	sleep(1000);
}

L1:

# load into mysql
FSSP::Connect();

my $cmd="SELECT MAX(align_id) AS align_id FROM dali.dali_dccp";
warn "cmd=$cmd\n";
FSSP::Execute($cmd);
my $result=FSSP::BuildResultStructure();
my $row=shift(@$result);
my $max_align_id=$row->{align_id};
if($max_align_id<1) { die "Fatal error accessing dali database\n"; }

my $cmd="cat /data/backup/zope_data/dali_server/DaliLite_3.0/new_[1-8]/*.dccp | perl /data/backup/zope_data/dali_server/DaliLite_3.0/load_dccp.pl $max_align_id 1 > /data/backup/zope_data/dali_server/DaliLite_3.0/dali_dccp 2> /data/backup/zope_data/dali_server/DaliLite_3.0/dali_segments";
warn "# cmd=$cmd\n";
system $cmd;

my $cmd="LOAD DATA LOCAL INFILE \'/data/backup/zope_data/dali_server/DaliLite_3.0/dali_dccp\' INTO TABLE dali.dali_dccp";
FSSP::Execute($cmd);
warn "# cmd=$cmd\n";

my $cmd="LOAD DATA LOCAL INFILE \'/data/backup/zope_data/dali_server/DaliLite_3.0/dali_segments\' INTO TABLE dali.dali_segments";
FSSP::Execute($cmd);
warn "# cmd=$cmd\n";

FSSP::DisConnect();

exit();

# update pdb90.list
system "cp pdb90.list ./DB/";

