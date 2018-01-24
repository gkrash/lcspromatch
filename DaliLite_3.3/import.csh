# run DaliLite -r on /data/databases/pdb/pdb<cd>.ent.gz

zcat /data/databases/pdb/pdb$1\.ent.gz > $1\.brk
rm -f dali.lock
perl DaliLite -r $1\.brk $1

rm -f $1\.brk $1\.dssp
