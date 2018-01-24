#!/usr/bin/perl
#
# use as: ... | pipeforwarder.perl arguments | ...
#
# forwards all lines preceded by command-line arguments
#
#
# print argument list
#
#print "\n***** Arguments given to $0 were: @ARGV\n\n";

foreach $i (0 .. $#ARGV) {
	print $ARGV[$i],"\n";
}

#
# pass through lines from standard input
#

while(<STDIN>){
	print ($_);
}
#
# what happens if STDIN is dry ?
#
