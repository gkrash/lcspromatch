# Makefile for "DaliLite", the protein structure comparison program.
#
#  Copyright (C) 2007 Liisa Holm.
#
# You must set HOMEDIR and ESCAPED_HOMEDIR (lines 17-18) !!
#
# make install			# sets home directory in DaliLite script, makes all 
# make all 			# compile and link all programs used by DaliLite
# make test1			# testing DaliLite -pairwise
# make test2			# testing advanced DaliLite options
# make clean			# remove all .o and binary files

SHELL = /bin/sh

# Directories <--- modify these!!
# WARNING: the maximum acceptable length of HOMEDIR is 70 characters!
HOMEDIR=/home/krash/git/DaliLite_3.3/
ESCAPED_HOMEDIR=\/home\/krash\/git\/DaliLite_3.3\/
 
# Compiler options
FC = g77                   # Fortran compiler
OPT = # -O		   # uncomment this for optimization
DEB = # -g		   # uncomment this for debugging
DEFS =			   # place other -D types of things here
CFLAGS = $(OPT) $(DEB) $(DEFS)
LFLAGS =

# assuming that Fortran compiler does not know ran function
RANF = ran.f
RANO = ran.o

# objects
OBJS = clean.o dccpto90.o dimple.o dp.o filter95fitz.o fssp.o uagadccp.o \
       gagagene.o gagahash.o gagarubber.o gagatool.o lean.o main.o \
       parsi-1.o parsi-admin.o parsi-align.o parsi-score.o parsi-stack.o \
       pipe.o pipedccp.o puu.o puutos.o $(RANO) replace.o selfdccp.o \
       seqfitz.o soap4.o ssap.o subfitzfast.o testi.o triplet.o u3b-8.o \
       util.o wolf.o AccSurf.o CalcAccSurf.o Contacts.o Date.o \
       DsspCMBI.o Vector.o p2clib.o

EXES = soap4 fssp wolf seqfitz dp pipedccp replace selfdccp \
	filter95fitz pipe puutos dalicon parsi puu dssp wolf_grid

# Top-level targets
all: soap4 fssp wolf90 wolf seqfitz dp pipedccp replace selfdccp \
	filter95fitz pipe puutos dalicon parsi puu dssp \
	wolf_grid

install: all
	# create customised copy of DaliLite script in homedir
	perl -pe 's/\$$DALI_SERVER_HOME=\"(\S+)\"/\$$DALI_SERVER_HOME=\"$(ESCAPED_HOMEDIR)\"/' < DaliLite > $(HOMEDIR)/DaliLite 
	chmod 755 $(HOMEDIR)/DaliLite

dssp: DsspCMBI.o p2clib.o Date.o Vector.o AccSurf.o CalcAccSurf.o Contacts.o
	echo "Running script to compile the CMBI version of DSSP, please wait..."
	gcc -c -O -Wall -DGCC AccSurf.c
	gcc -c -O -Wall -DGCC CalcAccSurf.c
	gcc -c -O -Wall -DGCC Contacts.c
	gcc -c -O -Wall -DGCC Date.c
	gcc -c -O -Wall -DGCC DsspCMBI.c
	gcc -c -O -Wall -DGCC Vector.c
	gcc -c -O -Wall -DGCC p2clib.c
	gcc -o dssp DsspCMBI.o p2clib.o Date.o Vector.o AccSurf.o CalcAccSurf.o Contacts.o -lm
	echo "Type dssp PDBSourcefile DSSPDestinationfile to run the program..."

# Cygwin appends .exe to all executables, but DaliLite expects bare -o name

soap4: parsizes.for soap4.f u3b-8.f subfitzfast.f util.f
	$(FC) $(CFLAGS) soap4.f u3b-8.f subfitzfast.f util.f -o soap4

fssp: fssp.f util.f u3b-8.f
	$(FC) $(CFLAGS) fssp.f util.f u3b-8.f -o fssp

wolf90: parsizes.for wolf.f 
	$(FC) $(CFLAGS) wolf.f -o wolf90

wolf: parsizes.for wolf_original.f subfitzfast.f u3b-8.f util.f
	$(FC) $(CFLAGS) wolf_original.f subfitzfast.f u3b-8.f util.f -o wolf


seqfitz: parsizes.for seqfitz.f subfitzfast.f u3b-8.f util.f
	$(FC) $(CFLAGS) seqfitz.f subfitzfast.f u3b-8.f util.f -o seqfitz

dp: parsizes.for dp.f util.f 
	$(FC) $(CFLAGS) dp.f util.f -o dp

pipedccp: parsizes.for pipedccp.f
	$(FC) $(CFLAGS) pipedccp.f -o pipedccp

selfdccp: parsizes.for selfdccp.f util.f
	$(FC) $(CFLAGS) selfdccp.f util.f -o selfdccp

filter95fitz: parsizes.for filter95fitz.f subfitzfast.f u3b-8.f util.f
	$(FC) $(CFLAGS) filter95fitz.f subfitzfast.f u3b-8.f util.f -o filter95fitz

pipe: parsizes.for pipe.f util.f
	$(FC) $(CFLAGS) pipe.f util.f -o pipe

puutos: parsizes-puutos.for puutos.f util.f 
	$(FC) $(CFLAGS) puutos.f util.f -o puutos

dalicon: gagasizes.for gagarubber.f clean.f dimple.f gagadccp.f \
		gagagene.f gagahash.f gagatool.f util.f ssap.f triplet.f \
		u3b-8.f lean.f testi.f $(RANF) subfitzfast.f
	$(FC) $(CFLAGS) gagarubber.f clean.f dimple.f gagadccp.f \
		gagagene.f gagahash.f gagatool.f util.f ssap.f triplet.f \
		u3b-8.f lean.f testi.f $(RANF) subfitzfast.f -o dalicon

parsi: parsizes.for parsi-admin.f parsi-align.f parsi-score.f \
		parsi-stack.f parsi-1.f $(RANF) util.f
	$(FC) $(CFLAGS) parsi-admin.f parsi-align.f parsi-score.f \
		parsi-stack.f parsi-1.f $(RANF) util.f -o parsi

puu: sizes-puu.for main.f puu.f util.f $(RANF)
	$(FC) $(CFLAGS) main.f puu.f util.f $(RANF) -o puu

# Miscellaneous
clean:
	-rm -f $(OBJS) $(EXES)

# Simple test for pairwise runmode
test1:
	echo "#####" ; \
	echo "# Testing -pairwise mode" ; \
	echo "#####" ; \
	echo ""; \
	if [ ! -d "$(HOMEDIR)/test1" ]; then \
                mkdir $(HOMEDIR)/test1; \
		cp -f $(HOMEDIR)/8abp.brk $(HOMEDIR)/test1/; \
		cp -f $(HOMEDIR)/2lbp.brk $(HOMEDIR)/test1/; \
        fi; \
	cd $(HOMEDIR)/test1 ; \
	rm -f dali.lock ; \
	../DaliLite -pairwise 8abp.brk 2lbp.brk >log 2>err; \
	echo "#####" ; \
	echo "# Use a web browser to view the result page $(HOMEDIR)/test1/index.html" ; \
	echo "#####" ; \

# Clean up after test
clean_test1:
	echo ""; \
	if [ -d "$(HOMEDIR)/test1" ]; then \
		cd $(HOMEDIR)/test1 ; \
		rm -f dali.lock CA* *.html replist list* log err matrix.txt mol* ranges.txt; \
		echo "Done."; \
	fi;

# More tests
test2:
	echo "" ; \
	echo "#####" ; \
	echo "# Testing -align and -format modes" ; \
	echo "#####" ; \
	echo ""; \
	if [ ! -d "$(HOMEDIR)/test2" ]; then \
		mkdir $(HOMEDIR)/test2 ; \
		cp -f $(HOMEDIR)/1fkx.brk $(HOMEDIR)/test2/ ; \
                cp -f $(HOMEDIR)/2kau.brk $(HOMEDIR)/test2/ ; \
	fi; \
	cd $(HOMEDIR)/test2 ; \
	rm -f dali.lock ; \
	echo "# (1) Reading PDB files: " ; \
	echo "$(HOMEDIR)/DaliLite -readbrk 1fkx.brk 1fkx >log 2>err" ; \
	$(HOMEDIR)/DaliLite -readbrk 1fkx.brk 1fkx >log 2>err ; \
	echo "$(HOMEDIR)/DaliLite -readbrk 2kau.brk 2kau  >>log 2>>err" ; \
	$(HOMEDIR)/DaliLite -readbrk 2kau.brk 2kau  >>log 2>>err ; \
	echo "#  For the 2kau test, chains A, B and C are generated." ; \
	echo "#  So, 1fkx_.dat, 2kauA.dat, 2kauB.dat and 2kauC.dat will be" ; \
	echo "#  generated in $(HOMEDIR)/DAT." ; \
	echo "" ; \
	echo "# (2) Testing structural alignment by DaliLite:" ; \
	echo "#  Here, we choose 2 structures to align: 1fkx and 2kauC which" ; \
	echo "#  have been generated in step (1)." ; \
	echo "$(HOMEDIR)/DaliLite -align 1fkx 2kauC >>log 2>>err" ; \
	$(HOMEDIR)/DaliLite -align 1fkx 2kauC >>log 2>>err ; \
	echo "#  The result will be DCCP file 1fkx.dccp ." ; \
	echo "$(HOMEDIR)/DaliLite -align 2kauC 1fkx >>log 2>>err"; \
	$(HOMEDIR)/DaliLite -align 2kauC 1fkx >>log 2>>err ; \
	echo "#  The result will be DCCP file 2kauC.dccp ." ; \
	echo "" ; \
	echo "# (3) Testing HTML file generation:" ; \
	echo "#  There must exist a file 'replist' which contains the list of" ; \
	echo "#  representatives. No other structures will occur in the output," ; \
	echo "#  even if they occur in the dccp file. A DSSP file for each" ; \
	echo "#  structure must also be available." ; \
	if [ ! -e "$(HOMEDIR)/test2/replist" ]; then \
                cp -f $(HOMEDIR)/replist $(HOMEDIR)/test2/ ; \
        fi; \
	echo "$(HOMEDIR)/DaliLite -format 1fkx 1fkx.dccp replist 1fkx.html >>log 2>>err" ; \
	echo "#  $(HOMEDIR)/test2/1fkx.html is the output file name." ; \
	$(HOMEDIR)/DaliLite -format 1fkx 1fkx.dccp replist 1fkx.html >>log 2>>err ; \
	echo "$(HOMEDIR)/DaliLite -format 2kauC 2kauC.dccp replist 2kauC.html >>log 2>>err" ; \
	echo "#  $(HOMEDIR)/test2/2kauC.html is the output file name." ; \
	$(HOMEDIR)/DaliLite -format 2kauC 2kauC.dccp replist 2kauC.html >>log 2>>err ;
	#####
	# Use a web browser to view the results in $(HOMEDIR)/test2
	#
	# You can compare the results with our results stored in
	# $(HOMEDIR)/test2/example_output.  There is a log file test.log which is
	# the screen output from test.sh.  Note that you may get different numbers
	# on different systems because optimization algorithms used are stochastic.
	# Apart from the above mentioned outputs, there are DAT files generated
	# and stored in $(HOMEDIR)/DAT directory. If there exist 4 files:
	# 
	#  1fkx_.dat  2kauA.dat  2kauB.dat  2kauC.dat
	#
	# Your testing was good. If not, something went wrong which you would have
	# noticed already in the testing steps.
	#####

clean_test2:
	if [ -d "$(HOMEDIR)/test2" ]; then \
		cd $(HOMEDIR)/test2 ; \
		rm -f dali.lock dali.default *.dssp list1 list2 *.dccp *.html log err; \
		echo "Done."; \
	fi;

# database search filtering
wolf_grid: parsizes.for wolf_grid.f util.f
	$(FC) $(CFLAGS) wolf_grid.f util.f -o wolf_grid

