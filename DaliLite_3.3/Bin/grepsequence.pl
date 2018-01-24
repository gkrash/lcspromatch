# write pdb.fasta 

use strict;

#my $DATADIR='/dataprod/RealData/sjkaaria/bin/DaliLite_2.4.4/DAT/';
my($DATADIR)=@ARGV;

opendir(DIR,$DATADIR) || die "Can't open $DATADIR\n";
my @filelist=(grep(!/^\./,readdir(DIR)));
closedir(DIR);

my $file;
foreach $file (@filelist) {
	#warn "$file\n";
	$_=$file;
	s/\.dat//;
	next if(/_/); # skip blank chain identifier
	s/_//;
	my $id=$_;
	my(@lines)=`grep '^\-sequence' $DATADIR\/$file`;
	$_=$lines[0];
	s/^\S+\s+\"//;
	if(length($_)>1) {
		print ">$id\n";
		print $_;
	}
}

exit;
