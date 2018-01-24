# map query sequence to pairsdb_v2 nid
# output input line for gtg_attributes.py
# sequence database must be in format nid|accessionnumber|identifier

## HARD-CODED: /data/bin/blastall

import os, sys, commands, re, datetime, random
import string


# Run blast on given sequence
def RunBlast(seq,qid=''):
        QID = qid
        MAX_HITS = 1
        # Seq should already be as single uppercase sequence!
        statement = "echo %(seq)s | /data/bin/blastall -pblastp -d%(db)s -v%(hits)d -b%(hits)d -FF -f1000 -W3" %{'seq':seq, 'hits':MAX_HITS, 'db':RSDB}
        header = []
        ph = os.popen( statement )
        # capture data line by line
        data = ""
        while 1:
                data = ph.readline()
                if not data: break
                if data[0:9] == "Sequences":
                        data = ph.readline()
                        while data == '\n': data = ph.readline() #advance to first hit summary
                        break
                elif data[0:20] == " ***** No hits found":
                        data = ""
                        break
        # In result-summary lines, chars 1-66 = header, 67-72 = score, 73-80 = evalue
        if not data: return 0
        count = 1
        while 1:
                line = data.split('|',1)
                try: hit_nid = int(line[0])
                except(ValueError): print "Error! Blastall hit-summary data were corrupt, no nid at start of line."; return 0
                if not QID: QID = line[1].split('|')[0]
                # Leave a * marker for later insertion of alignment string (first radiobutton will be checked by default)
                if count == 1: radio = '<input TYPE="radio" NAME="qnidqidqali" VALUE="%d|%s|*" CHECKED />' %(hit_nid,QID)
                else: radio = '<input TYPE="radio" NAME="qnidqidqali" VALUE="%d|%s|*"/>' %(hit_nid,QID)
                anchor = '<a href="#%s">%d</a>' %(count,hit_nid)
                description = line[1].split()
                evalue = description.pop(-1)
                score = description.pop(-1)
                description = ' '.join(description)
                header.append( [radio, anchor, description, score, evalue] )
                data = ph.readline()
                if data == '\n': break
                count += 1
        # Grab rest of blast data and split on '>' char
        data = ph.read()
        ph.close()
        data = string.replace(data,'->','--') # hack
        data = data.split('  Database: nrdb40_v2.fasta')[0].split('>')[1:]
        # Grab compressed alignment (as underscore-separated string) and add to input tag in summary data
        for x in range(0,len(header),1):
                line = data[x].split('|',1)
                hit_nid=line[0]
                qfrom,qto,hit_from,hit_to,ali = _GetCompressedAli(data=data[x],asarray=1,fromto=1)  # ali is list of tuples [ (qstart,sstart,blocklen) ]
                qstarts=''
                sstarts=''
                blocks=''
                for tuple in ali:
                        qstarts+="%i " %tuple[0]
                        sstarts+="%i " %tuple[1]
                        blocks+="%i " %tuple[2]
                # print gtg_attributes input line
                print "%s\tgtg\t%s\t%s\t%s\t%s\t%s" %(qid,hit_nid,qstarts,sstarts,blocks,seq)



# INTERNAL FUNCTIONS
def _GetCompressedAli(data,asarray=1,fromto=1):
        GAPCHAR = '-'
        qseq,sseq = "",""
        qfrom,qto,sfrom,sto = 0,0,0,0
        mydata = data.split("Score")[1].splitlines() # only want the first alignment to the sbjct for each hit !!!
        while '' in mydata: mydata.remove('')
        # 1: delineate ali segments and get from/to values
        for x in range(0,len(mydata),1):
                if mydata[x][0:5] == "Query":
                        line = mydata[x].split()[1:]
                        if qfrom == 0: qfrom = int(line[0])
                        qseq += line[1]
                        qto = int(line[2])
                elif mydata[x][0:5] == "Sbjct":
                        line = mydata[x].split()[1:]
                        if sfrom == 0: sfrom = int(line[0])
                        sseq += line[1]
                        sto = int(line[2])
        qali,sali = "",""
        mylen = len(qseq)
        # Quick sanity check
        if mylen != len(sseq): print "Input data corrupt!"; return 0

        # get aligned blocks
        blocklen=0
        ali=[]
        i=0
        qi=qfrom-1
        si=sfrom-1
        qstart=qi
        sstart=si
        while i<mylen:
		#print "#i,ali ",i,si,qi,blocklen,qseq[i],sseq[i],ali
                if qseq[i]==GAPCHAR or sseq[i]==GAPCHAR:
                        ali.append( (qstart,sstart,blocklen) )
                        blocklen=0
			if qseq[i]==GAPCHAR:
				si+=1
	                        while qseq[i+1]==GAPCHAR:
	                                i+=1
	                                si+=1
			else:
				qi+=1
				while sseq[i+1]==GAPCHAR:
	                                i+=1
        	                        qi+=1
                        qstart=qi
                        sstart=si
                else:
                        blocklen+=1
                        si+=1
                        qi+=1
 	        i+=1
        if blocklen>0: ali.append( (qstart,sstart,blocklen) )

        return(qfrom,qto,sfrom,sto,ali)


def RunFromFastafile(fastafile):
        file=open(fastafile,'r')
        qid=''
        QSEQ=''
        for line in file:
                if line[0] == '>':
                        # run previous sequence
                        if len(QSEQ)>0:
                                QSEQ=string.strip(QSEQ)
                                RunBlast(QSEQ,qid)
                        qid=line[1:].split()[0]
                        QSEQ=''
                else:
                        QSEQ+=line
        if len(QSEQ)>0:
                QSEQ=string.strip(QSEQ)
                RunBlast(QSEQ,qid)

# Test GetFold() Method:
# input fasta file, call RunBlast on each sequence

if __name__ == '__main__':
        RSDB=sys.argv[1]
        fastafile=sys.argv[2]
        RunFromFastafile(fastafile)

        sys.exit(0)

