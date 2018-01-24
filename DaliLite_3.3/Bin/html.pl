# HTML output to STDOUT: summary block, pairwise alignment block
# text output to STDERR: summary block, range block
#
# >> query_pdbresnos field!! 
# 
use strict;

if($#ARGV<0) { die "Usage $0 DALIDATDIR_1 DALIDATDIR_2 jobid \n"; }
my($DALIDATDIR_1, $DALIDATDIR_2, $jobid)=@ARGV;
# query (cd1) in DALIDATDIR_1, targets (cd2) in DALIDATDIR_2

my $gap='-';
my $identity_symbol='|';
my $lw=60;

my $summary='';
my $summaryhtml='';
my $alignments='';
my $ranges='';
my %sequence;
my %dssp;
my %compnd;
my(@u);
my(@t);
my $pdbresnos; # zeroth element is nres
my $hasmol2=0;

# read one alignment from STDIN
# append summary line to $summary
# append pairwise alignment to $alignments
my $cd2;
my $nfrag;
my $ires;
my $jres;
my @seq1;
my @seq2;
my @ss1;
my @ss2;
my @ranges;
my $cd1;
my $qstarts='';
my $sstarts='';
my $lengths='';
my $nhit=0;
my $query_sequence='';
my $query_dssp='';
while(<STDIN>) {
	if(/^-summar\s+\"(.*)\"/) { # this is last line in dbhit-block
		#print join("\n","<PRE>",join('',@ss1),join('',@seq1),join('',@ss2),join('',@seq2),"</PRE>"),"\n";

		my $x=$1;
		$x=~s/     0      0\s+\d+ S   //;
		my $y=sprintf("%4d: ",$nhit);
		$x=~s/^....../$y/;

                $_=$x;
		my($z)=/\:\s+\S+\s+(\S+)/;

                # rotations
                my $u=join(',',@u);
                my $t=join(',',@t);

		# link to PDB download
		my $xhtml=$x;
		my($cd2x)=substr($cd2,0,4);
		my($chainid)=substr($cd2,4,1);
		if($cd2x ne 'mol2') { $xhtml=~s/^(.{40})(.*)$/$1\<A HREF=\"http:\/\/ekhidna.biocenter.helsinki.fi\/dali_server\/qz-test\?jobid=pdb&pdbid\=$cd2x&u\=$u&t\=$t\"\>PDB\<\/A\>$2/; }
		else { $xhtml=~s/^(.{40})(.*)$/$1\<A HREF=\"http:\/\/ekhidna.biocenter.helsinki.fi\/dali_server\/qz-test\?jobid\=$jobid\&pdbid\=$cd2x&u\=$u&t\=$t\"\>PDB\<\/A\>$2/; } 
                # link cd2 to DaliDB
		if($cd2x ne 'mol2') { $xhtml=~s/^(.{7})(\S+)(.*)$/$1<A HREF=\"..\/..\/..\/\/dali\/daliquery?pdbid\=$cd2x\&chainid\=$chainid\"\>$2\<\/A\>$3/; }
		else { $hasmol2=1; } # exclude db-link info from header

		# jump to pairwise alignment
		$xhtml=~s/^(\s+)(\d+)/$1\<A HREF=#alignment\-$2\>$2\<\/A>/;

		# checkbox
		my $seq2=join('',@seq2); $seq2=~s/\W//g; # remove leading '?'
		my $dssp2=join('',@ss2); $dssp2=~s/\W//g; # remove leading '?'
		$xhtml=~s/^/\<INPUT TYPE=\"CHECKBOX\" NAME=\"cd2list\" VALUE=\"cd2=$cd2 sequence=$seq2 dssp=$dssp2 qstarts=$qstarts sstarts=$sstarts lengths=$lengths u=$u t=$t\">/;

		$summary.=$x."\n";
		$summaryhtml.=$xhtml."\n";

		# unaligned C-terminus
                my $i;
                foreach $i ($ires+1 .. length($sequence{$cd1})) {
                        ($seq1[$i])=~tr/[A-Z]/[a-z]/; # unaligned
                        ($ss1[$i])=~tr/HEL/hel/; # unaligned
                }
                foreach $i ($jres+1 .. length($seq2)) {
                        ($seq2[$i])=~tr/[A-Z]/[a-z]/; # unaligned
                        ($ss2[$i])=~tr/HEL/hel/; # unaligned
                }
		# pairwise alignment block
		$alignments.="<a name=alignment\-$nhit\><h3> No $nhit\: Query=$cd1 Sbjct=$cd2 Z-score=$z</h3>\n";
		$alignments.="\n<A HREF=\"\#$cd1\"\>back to top\<\/A>\n<PRE>\n";
		my @lines;
		my $seq1; foreach (@seq1) { if(/[\w]/) { $seq1.=$_; } } # cygwin adds ^M 
		my $ss1; foreach (@ss1) { if(/[\w]/) { $ss1.=$_; } } 
		my $seq2; foreach (@seq2) { if(/[\w]/) { $seq2.=$_; } } 
		my $ss2; foreach (@ss2) { if(/[\w]/) { $ss2.=$_; } }
		#$alignments.="<PRE>\nss1  $ss1\nseq1 $seq1\nseq2 $seq2\nss2  $ss2\n</PRE>\n";

		push(@lines,'?'.$ss1."\n",'?'.$seq1."\n",'?'.$seq2."\n",'?'.$ss2."\n");
		my(@longali)=&poppush(@lines);
		#$alignments.="\n<PRE>\n @longali\n";
		my @mali;
		foreach (@longali) { s/\?//; s/\n//g; push(@mali,$_); }
		# mark identical amino acids
		my(@x)=split(//,$mali[1]);
		my(@y)=split(//,$mali[2]);
		my $star='';
		my $x;
		foreach $x (@x) {
			my $y=shift(@y);
			if($x eq $y && $x =~ /[A-Z]/) { $star.=$identity_symbol; } else { $star.=' '; }
		}
		$mali[4]=$mali[3]; $mali[3]=$mali[2]; $mali[2]=$star;		
		my(@prolog)=('DSSP  ','Query ','ident ','Sbjct ','DSSP  ');
		my $i=0;
		my $qres=0;
		my $sres=0;
		my $l=length($mali[$[]);
		while($i<$l) {
			my $row;
			foreach $row (0..4) {
				my $resno='';
				my $nc=0;
				if($prolog[$row]=~/^[QS]/) {
					$_=substr($mali[$row],$i,$lw);
					$nc=s/[a-zA-Z]//g;
				}
				if($prolog[$row]=~/^Q/) { $qres+=$nc; $resno=sprintf("%5d",$qres); } 
				elsif($prolog[$row]=~/^S/) { $sres+=$nc; $resno=sprintf("%5d",$sres); }
				$alignments.=$prolog[$row].substr($mali[$row],$i,$lw).$resno."\n";
			}
			$alignments.="\n\n";
			$i+=$lw;
		}
		$alignments.="</PRE>";

		# ranges block
		# print delimiter 'cos pairwise alignments can have multiple identical cd1-cd2 pairs!
		#$ranges.="# No $nhit\: z-score=$z\n";
		foreach(@ranges) {
			my $x=sprintf("%4d: ",$nhit);
			s/^-ranges\s+/$x/;
			s/\"//g;
			$ranges.=$_;
		}
	} elsif(/^-matrix\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+([\d\/\-]+)/) {
		push(@u,$1,$2,$3);
		push(@t,$4);
	} elsif(/^-protein\s+(\w+)/) {
		$cd2=$1;
		$nhit++;
		$ires=0;
		$jres=0;
		$qstarts='';
		$sstarts='';
		$lengths='';
		undef(@u);
		undef(@t);
		undef(@ranges);
		# initialize Query for poppush
		@seq1=split(//,'?'.$sequence{$cd1}); 
		@ss1=split(//,'?'.$dssp{$cd1}); 		
        } elsif(/^-sequen\s+\"(\S+)\"/) { # follows immediately after protein, precedes ranges
                my $x=$1;
		$x=~s/[a-z]/C/g; # disulphides in DSSP
		@seq2=split(//,'?'.$x);
        } elsif(/^-struct\s+\"(.+)\"/) { # follows immediately after protein, precedes ranges
		my $x=$1;
                @ss2=split(//,'?'.$x);
	} elsif(/^-query/) {
		s/^-query\s+//;
		s/[\-]//;
		s/\s+//;
		$cd1=$_;
		&get_data($cd1,$DALIDATDIR_1);
		$query_sequence=$sequence{$cd1}; # query must be in uppercase for poppush
		$query_dssp=$dssp{$cd1};
		@seq1=split(//,'?'.$sequence{$cd1}); 
		@ss1=split(//,'?'.$dssp{$cd1}); 		
		$summary="# Query: $cd1\n# No:  Chain   Z    rmsd lali nres  \%id PDB  Description\n";
	} elsif(/^-ranges/) {
		push(@ranges,$_);
		my($from1,$to1,$from2,$to2)=/(\d+) -\s*(\d+) <=>\s+(\d+) -\s*(\d+)/;
		my $i;
                foreach $i ($ires+1..$from1-1) {
                        ($seq1[$i])=~tr/[A-Z]/[a-z]/; # aligned
                        ($ss1[$i])=~tr/[HEL]/[hel]/; # aligned
                }
                foreach $i ($jres+1..$from2-1) {
                        ($seq2[$i])=~tr/[A-Z]/[a-z]/; # aligned
                        ($ss2[$i])=~tr/[HEL]/[hel]/; # aligned
                }
		$ires=$to1;
		$jres=$to2;
		#print "mark range ($from1,$to1,$from2,$to2)\n";
		my $q=$from1-1; # segment numbering starts from zero!
		my $s=$from2-1; # segment numbering starts from zero!
		my $l=$to1-$from1+1;
		if($qstarts eq '') { $qstarts="$q"; $sstarts="$s"; $lengths="$l"; } else { $qstarts.=",$q"; $sstarts.=",$s"; $lengths.=",$l"; }
	} elsif(/^-nres/) {
		s/^\S+\s+//;
		$pdbresnos=$_;
	}
}

# always: HTML output goes to STDOUT, text output goes to STDERR
warn "$summary\n";
warn "# Structural equivalences\n";
warn "$ranges\n";

# HTML header
my $x='to pre-computed structural neighbours in the Dali Database, ';
if($hasmol2) { $x=''; }
print<<EOB;
<TITLE>Dali: $cd1, $compnd{$cd1}</TITLE><A name=$cd1>
<h1>Query: <a name=\"$cd1\">$cd1\</a></h1>
$compnd{$cd1}
<P><A NAME="summary">
Select neighbours (check boxes) for viewing as multiple structural alignment or 3D superimposition. 
The list of neighbours is sorted by Z-score. Similarities with a Z-score lower than 2 are spurious.
Each neighbour has links to pairwise structural alignment with the query structure, $x
 and to the PDB format coordinate file where the neighbour
is superimposed onto the query structure. 

<form METHOD=\"POST\" ACTION=\"http://ekhidna.biocenter.helsinki.fi/dali_server/qz-test\" TARGET="_blank"/>
<input TYPE=\"hidden\" NAME=\"jobid\" VALUE\=$jobid\>
<input TYPE=\"hidden\" NAME=\"query_sequence\" VALUE=\"$query_sequence\"> 
<input TYPE=\"hidden\" NAME=\"query_dssp\" VALUE=\"$query_dssp\"> 
<input TYPE=\"hidden\" NAME=\"cd1\" VALUE=$cd1\>
<input TYPE=\"hidden\" NAME=\"pdbresnos\" VALUE=\"$pdbresnos\"\>
<input TYPE=\"hidden\" NAME=\"cd2list\" VALUE=\"dummy\">
<input TYPE=\"submit\" NAME=\"D1\" VALUE=\"Structural Alignment\"/> 
<INPUT TYPE=\"CHECKBOX\" NAME=\"expand_gaps\" VALUE=1 CHECKED\"/> Expand gaps 
<input TYPE=\"submit\" NAME=\"D3\" VALUE=\"3D Superimposition (Jmol Applet)\"/> 
<input TYPE=\"reset\" VALUE=\"Reset Selection\"/>
<h2>Summary</h2>
<PRE>
    No:  Chain   Z    rmsd lali nres  \%id PDB  Description
$summaryhtml
</form>
</PRE>
<A NAME="alignments">
<h2>Pairwise Structural Alignments </h2>
Notation: three-state secondary structure definitions by DSSP (reduced to H=helix, E=sheet, L=coil) are shown above the amino acid sequence. Structurally equivalent residues are in uppercase, structurally non-equivalent residues (e.g. in loops) are in lowercase. Amino acid identities are marked by vertical bars.
$alignments
</PRE>
<HR>
EOB

exit();

###############################################################################

sub get_data {
	my($cdx,$DALIDATDIR)=@_;
	my($cd)=substr($cdx,0,4);
	my($chain)=substr($cdx,4,1);
	#warn "get_data: $cdx gave $cd and $chain\n"; 
	if($chain !~ /\w/) { $chain='_'; }
	my $datafile="$DALIDATDIR\/$cd$chain\.dat";
	open(IN,"<$datafile") || warn "Can't open $datafile\n";
	while(<IN>) {
		if(/^\-dssp\s+\"(.*)\"*$/) { $dssp{$cdx}=$1; }
		elsif(/^\-sequence\s+\"(.*)\"*$/) { $sequence{$cdx}=$1; }
		elsif(/^\-compnd\s+\"(.*)\"/) { $compnd{$cdx}=$1; }
	}
	close(IN);
	# disulphides in DSSP
	$sequence{$cdx}=~s/[a-z]/C/g;
	#warn "get_data: $cdx > $cd $chain $dssp{$cdx} $sequence{$cdx} $compnd{$cdx}\n";
}

###############################################################################

sub poppush {
	my(@multali)=@_;
	my($iseq,$maxl,$nseq,@where,$finished,$tuck,@tuck,$x,@current,@l,$l);
	my(@pending,$y,$i,$j,$c);
	my @longali;

	$nseq=$#multali; 
	foreach($[..$nseq) { 
		$longali[$_]=substr($multali[$_],$[,1); 
		$where[$_]=$[;
	}
	$finished=0;
	while($finished==0) {
		$finished=1;
		$tuck=0;
		foreach $iseq ($[..$nseq) {
		  $x=$multali[$iseq]; $current[$iseq]='';
		  $where[$iseq]++; $tuck[$iseq]=0;
		  if($where[$iseq]<=length($x)-1+$[) {
			$finished=0;
			$i=$where[$iseq]; 
			$c=substr($x,$i,1); $_=$c;
			while(/[a-z]/) {
				$pending[$iseq].=$c;
				$i++; $where[$iseq]++; 
				$c=substr($x,$i,1); $_=$c;
			} 
			if(/\\/) { # C-term insert
				$i++; $where[$iseq]++; 
				$c=substr($x,$i,1);
				$pending[$iseq].=$c;
			} elsif(/[$gap]/) { $current[$iseq].=$c; 
			} else {
				# tuck longest insert first
				$tuck=1; $tuck[$iseq]=1;
				$current[$iseq].=$c;
			}
		  }
		}
		if($tuck) {
			$i=0; $maxl=0; 
			foreach $iseq ($[..$nseq) { # free space=insert+gap-row
			  if($tuck[$iseq]) {
				$y=$longali[$iseq];
				$j=0; $_=substr($y,$j-1,1); 
				while(/[$gap]/) { $j--; $_=substr($y,$j-1,1); }
				$x=length($pending[$iseq])+$j;
				if($x>$maxl) { 
					$maxl=$x; $i=$iseq; 
				}
			  }
			}
			$iseq=$i;
			&tuck($iseq,$nseq,\@longali,$pending[$iseq]); 
			$pending[$iseq]=''; $tuck[$iseq]=0;
			foreach  $iseq ($[..$nseq) {
			  if($tuck[$iseq]) {
				&tuck($iseq,$nseq,\@longali,$pending[$iseq]);
				$pending[$iseq]='';
			  }
			}
		}
		foreach $iseq ($[..$nseq) { 
			$longali[$iseq].=$current[$iseq]; 
		}
	}
	$maxl=0;
	foreach $iseq ($[..$nseq) { 
		$x=$longali[$iseq]; 
		$i=0; $_=substr($x,$i-1,1); 
		while(/[$gap]/) { $i--; $_=substr($x,$i-1,1); }
		$l=length($x); $y=substr($x,$[,$l+$i);
		$x=$y.$pending[$iseq]; $longali[$iseq]=$x;
		$l=length($x); $l[$iseq]=$l;
		if ($l>$maxl) { $maxl=$l; }
	}
	foreach $iseq ($[..$nseq) { foreach($l[$iseq]+1..$maxl) { $longali[$iseq].=$gap; }}

	return(@longali);
}

#############################################################################

sub tuck {
	my($iseq,$nseq,$longali,$insert)=@_;
	my ($i,$j,$k,$in,$y,$x,$xseq);
	my($leni,$leng,$q,$z,@pending);

	$leni=length($insert); return if($leni==0);

	$x=@$longali[$iseq];
	$i=0; $_=substr($x,$i-1,1); while(/[$gap]/) { $i--; $_=substr($x,$i,1);}
	if($i==0) { $leng=0; } else { $leng=-1-$i; }

	$in=''; $k=$leni-$leng; foreach(1..$k) { $in.=$gap; } 

	$i=0; $x=@$longali[$[]; $k=length($x); $y=0;
	if($in ne '' && $y==0) {
		$_=substr($x,$i-1,1); 
		while(-$i<$leng&&/[A-Z]/) {
			$i--;$_=substr($x,$i,1);
			if(/[$gap]/) { $y=1; last; } 
		}
	}
	if($y) { # right-justify in-gaps if an insert starts right there
	  foreach $xseq ($[..$nseq) {
		# print "check $xseq of $nseq: $pending[$xseq]\n";
		next if($pending[$xseq] eq '');
		$q=@$longali[$xseq]; $z=length($q); 
		$_=substr($q,$z-1,1); # print "$z / $xseq gave $_\n";
		if(/[A-Z]/) { $y=0; } last if(!$y);
	  }
	}
	if($y) {
		foreach($[..$nseq) { 
			my $pre=substr(@$longali[$_],$[,$k+$i);
			my $suf=substr(@$longali[$_],$i);
			@$longali[$_]=$pre.$in.$suf; 
		}
	} else { 
		foreach $iseq ($[..$nseq) { 
			@$longali[$iseq].=$in; 
		} 
	}

	if($leng<$leni) {
		$x=@$longali[$iseq];
		$i=0; $_=substr($x,$i-1,1); while(/[$gap]/) { $i--; $_=substr($x,$i,1);}
		$leng=-1-$i;
	}

	my $pre=substr(@$longali[$iseq],$[,length(@$longali[$iseq])-$leni);
	$x=$pre.$insert;
	@$longali[$iseq]=$x;
}

##############################################################################

