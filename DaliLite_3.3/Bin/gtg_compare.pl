# compare gtg attributes to database
# query from STDIN: perl pairsdb_mapper.pl --sequence=string --gtg --name=1abcD | grep gtg | cut -f 3,6-7,13-16 | python gtg_attributes.py 1
# database = pdb90.attributes_3 [simple attributes clusid.aaix] or pdb90.attributes_1 [CAA1s]

use strict;

if($#ARGV<0) { die "Usage: $0 /afs/bi/bioinfo/user/luholm/pdb90.attributes_1\n"; }

my($db)=@ARGV; # '/afs/bi/bioinfo/user/luholm/pdb90.attributes_1';
my %hash;
my $query_cd;
my $old_query='?';

while(<STDIN>) {
	next if(/^#/);
	my($query_cd,$ires,$aa,$nid,$jres,$aaix,$clusid)=split(/\s+/);
        if($query_cd ne $old_query) {
                if($old_query ne '?') { &dbscan($old_query); }
                undef(%hash);
        }
	my $key=join('_',$clusid,$aaix);
	$hash{$key}=$ires;
	$old_query=$query_cd;
}
if($old_query ne '?') {
        &dbscan($old_query);
}

exit();


sub dbscan {
	my($query_cd)=@_;
	open(IN,"<$db") || die "Can't open gtg attributes database $db\n";
	my $old_cd='?';
	my %seen;
	my %match;
	my $n=0;
	while(<IN>) {
		my($cd,$ires,$aa,$nid,$jres,$aaix,$clusid)=split(/\s+/);
		if($cd ne $old_cd) { $n++; undef(%seen); }
		$old_cd=$cd;
	        my $key=join('_',$clusid,$aaix);
		my $query_ires=$hash{$key};
		if($query_ires>0) { 
			next if(defined($seen{$query_ires})); # query_ires counted once only
			$match{$cd}++; 
		}
		$seen{$query_ires}=1;
	}
	close(IN);

	# get mean,stdev
	my $cd;
	my $sum2=0.0;
	my $sum=0.0;
	foreach $cd (keys %match) {
		my $score=$match{$cd};
		$sum+=$score;
		$sum2+=$score*$score;
	}
	if($n==0) { $n=1; }
	my $mean=$sum/$n;
	my $stdev=sqrt(($sum2-$n*$mean*$mean)/$n);
	if($stdev<1e-6) { $stdev=1e-6; }

	# output matches, bogus Z-score
	my $cd;
	foreach $cd (keys %match) {
		my $score=$match{$cd};
		substr($cd,1,3)=~tr/[A-Z]/[a-z]/;
		my $z=($score-$mean)/$stdev;
		print join("\t",$query_cd,$cd,$score,$z),"\n";
	}
}

