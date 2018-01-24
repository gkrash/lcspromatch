# work directory
cd ~/zdata/dali_server/DaliLite_3.3

# import.csh new_cd -- creates DAT entries -- pdbnew.list is list of new pdb entries
## pdbnew.list must be created manually: 
## ls -ls /data/databases/pdb | grep ' Jan  2 ' | perl -pe 's/^.*pdb//' | perl -pe 's/\.ent\.gz.*$//' | sort > pdbnew.list
perl -pe 's/^(\w+)/\.\/import.csh $1/' < pdbnew.list > x.csh
chmod +x x.csh
./x.csh

# pdb.fasta
perl Bin/grepsequence.pl ./DAT/ > pdb.fasta

# pdb.list; exclude blank chainid
ls ./DAT | perl -pe 's/\.dat//' | grep -v '_' > pdb.list
#grep '^>' pdb.fasta | perl -pe 's/^>//' > pdb.list

# pdb90.fasta
Bin/cd-hit -i pdb.fasta -o pdb90.fasta

# formatdb
formatdb -i pdb90.fasta
formatdb -i pdb.fasta 

# pdb90.list
grep '^>' pdb90.fasta | perl -pe 's/^>//' > pdb90.list

# WOLF90.dat
grep '^>' pdb90.fasta | perl -pe 's/^>//' | grep '^[1-9]' > list1
perl Bin/forwarder.pl ./DAT/ 12.0 | Bin/wolf_grid
mv fort.1 ./DB/WOLF90.dat

# pdbnew.list is list of new pdb entries
# make list of new chains
cp pdb.list ./DB/pdb.list
perl cd_vs_cd1.pl pdbnew.list < ./DB/pdb.list > new.list
# new pdb90-representatives
perl cd_vs_cd1.pl pdbnew.list < ./pdb90.list > do.list
# as fastafile
perl cd1list2fasta.pl ./DAT/ < do.list > do.fasta

cp pdb.fasta* pdb90.fasta* ./DB/ 

# pdb90.map incremental update
# exclude selenocysteine warnings from blastall, they truncate sequence
python Blast_GTG_multi.py DB/nrdb40_v2.fasta do.fasta | grep -v blastall | grep gtg > tmp.map
cat tmp.map >> ./DB/pdb90.map

# incremental pdb90.attributes_1
python Bin/gtg_attributes.py 1 ./DB/ < tmp.map > tmp
cat tmp >> ./DB/pdb90.attributes_1

# DaliLite -u new chains (does not use gtg)
# -u arguments are do-list new-list (for parallel processing)
# e.g. perl DaliLite -u new.list new.list 
# parallel update:
rm -f new_*/*
perl split8.pl < new.list # creates new_[1-8].list
./u.csh                   # 8 parallel jobs run in ./new_[1-8]/

# wait for all jobs to finish
perl ./perpetuum.pl

# usage statistics
#perl /data/backup/zope_data/dali_server/log_html.pl < /data/backup/zope_data/dali_server/log > /data/backup/zope_results/dali/stat.html

# PDB statistics
#perl /data/backup/zope_data/dali_server/pdbstat.pl > /data/backup/zope_results/dali/pdbstat.html
