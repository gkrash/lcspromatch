#!/usr/bin/perl
#
# use as: forwarder.perl arguments | ...
#
# forwards all command-line arguments on separate lines
#

#
# print argument list
#

foreach $i (0 .. $#ARGV) {
	print $ARGV[$i],"\n";
}


