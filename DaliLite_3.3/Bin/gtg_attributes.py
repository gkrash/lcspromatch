import sys, os, re
import array

###############################################################################
# Classes:
###############################################################################

class MyReadTable:
 """Retrieve topology data for a cluster. 
 
 [nid] -> front, back, seqstart-pointer, nres (type 'L')
	byte address = nid*16 

 sequence array (type 'c')
 	byte address = seqstart-pointer[nid]

 [leaf] -> cluster, arrayElement-pointer, nid, residue (type 'L')
 	byte address = vertex*16

 [vertex] -> nrecords in amino acid frequency vector array 

 amino acid frequency vector array (type 'B')
 	byte address = nrecords*21
	[0] = index of most frequent aa
	[1-20] = ACDEFGHIKLMNPQRSTVWY frequency as percentage

 [cluster] -> clusterstart-pointer (type 'L') -> cluster array tree
 	byte address = cluster*4
 
 Cluster tree array (type 'l')
	Each cluster is in one block array[0:nt*24]. 
	The array is conceptually organised as records of 6 elements (of 4 bytes). 

	Row 0 = header:
	===============
 	byte address = clusterstart-pointer*24
	[0] = record number (can be used for checking)
	[1] = cluster index (can be used for checking)
	[2] = root row number*6 = array element
	[3] = number of rows (nt)
	[4] = number of inter-cluster-edges from this cluster
	[5] = number of leafs in this cluster

	Rows 1..ns = leaf nodes:
	========================
	byte address = clusterstart-pointer*24 + arrayElement-pointer*4
 	[0] = vertex				# to link s..i, i-j..t pathlets
	[1] = parent_row*6 = array element
 	[2] = -nid 				# redundant with vertex array 
 	[3] = -(ires + 1,000,000*aatypeindex)	# redundant with vertex array 
 	[4] = pointer to inter-cluster-edge block in edge-database
	[5] = number of inter-cluster-edges in subtree
 
	Rows ns+1..nt-1 = internal nodes:
	===============================
 	[0] = vertex 				# internal vertex labels are not unique! 0, 21, 42, 63, ...
	[1] = parent_row*6 = array element
	[2] = leftChild_row*6 = array element
	[3] = rightChild_row*6 = array element
	[4] = pathWeight
	[5] = number of inter-cluster-edges in subtree
  
 """
 def __init__(self):
 	pass

 def getClusterHeader(self,cluster_index):
	# get pointer from fh_cluster_ptrin
	ptr=4*cluster_index
	fh_cluster_ptrin.seek(ptr)
	temp=array.array('I') ##64-bit?##temp=array.array('L')
	temp.fromfile(fh_cluster_ptrin,1)
	nrecords=temp[0]
	if nrecords<1: 
		print "# Warning from getClusterHeader: cluster %i not in database" %cluster_index
		return(0,0,0,0,0,0) 
	# read header record from fh_cluster_datain
	ptr=24*nrecords
	fh_cluster_datain.seek(ptr)
	temp=array.array('i') ##64-bit?## temp=array.array('l')
	temp.fromfile(fh_cluster_datain,6)
	return(temp)

###############################################################################

class MyCluster:
 """
 x=MyCluster(clusterIndex) -> x.nt, x.ni, x.ns
 x.getTable() -> cluster topology table
 """
 def __init__(self,clusterindex):
	 nrecords,x,root,nt,ni,ns=c.getClusterHeader(clusterindex)
	 self.clusterindex=clusterindex
	 self.nrecords=nrecords
	 self.nt=nt
	 self.ni=ni
	 self.ns=ns
	 self.root=root
 
 def getTable(self): 
	 nt_bytes=self.nt*6
	 ptr=self.nrecords*24
	 fh_cluster_datain.seek(ptr)
	 temp=array.array('i') ##64-bit?## temp=array.array('l')
	 temp.fromfile(fh_cluster_datain,nt_bytes)
	 return(temp)

###############################################################################

class MyAacomp:
 def __init__(self):
 	pass

 def get_anc40_aafreqnrecords(self,irow,table):
 	"""
 A leaf's anc40 is the last ancestral node that has pathWeight > 40. 
 	"""
 	while (table[irow+4]>param_anc_threshold or table[irow+2]<0) and table[irow+1]>0:
		if table[irow+1]==0: break # root
		prow=table[irow+1]
		if table[prow+4]<=param_anc_threshold: break
		irow=prow
	aafreqnrecords=(table[0]-1)+irow/6 # nrecords in aafreq.store, byteaddress=nrecords*21
	return(aafreqnrecords)

 def get_aaix_of_leaf(self,irow,table):
	 aafreqnrecords=(table[0]-1)+irow/6
	 temp=self.get_aacomp_from_aafreqptr(aafreqnrecords)
	 return(temp[0]) # most frequent aaix

 def get_aacomp_from_aafreqptr(self,aafreqnrecords):
	fh_aafreq_datain.seek(aafreqnrecords*21) # nrecords in aafreq.store, byteaddress=nrecords*21 
	temp=array.array('B')
	temp.fromfile(fh_aafreq_datain,21)
	#print "get_aacomp_from_aafreqptr",aafreqnrecords,temp
	return(temp)
	
###############################################################################


###############################################################################

def getCAA(leaf_vertex):
        # get CAA_id and aaix for leaf_vertex
        indexArray=array.array('L')
        address=leaf_vertex*8*2 # 64-bit
        try: fh_CAA_vertex_ptrin.seek(address)
        except: return((0,0)) # vertex not in database
        indexArray.fromfile(fh_CAA_vertex_ptrin,2)
        data_address=indexArray[0]*4*2 # 64-bit
        try: fh_CAA_vertex_datain.seek(data_address)
        except: return((0,0)) # address not in database
        dataArray=array.array('L')
        dataArray.fromfile(fh_CAA_vertex_datain,indexArray[1]) # list of CAAs containing query vertex
	#print "# leaf_vertex, dataArray",leaf_vertex,dataArray
        # return list of CAAs with aaix
        tuplelist=[]
        for CAA_id in dataArray:
                # lookup aaix of CAA_id
                aaix=get_aaix_for_CAAid(CAA_id)
                tuplelist.append( (CAA_id, aaix) )
        return(tuplelist) # CAA_id, aaix

def get_aaix_for_CAAid(CAA_id):
                # get all members of CAA
                indexArray=array.array('L')
                address=(CAA_id*12-12)*2 # 64-bit
                #print "# CAA_id=%i address=%i" %(CAA_id,address)
                try: fh_CAA_member_ptrin.seek(address)
                except:
                        #print "CAA_id not in database:",CAA_id
                        return(0)
                indexArray.fromfile(fh_CAA_member_ptrin,3)
                #print "query CAA_id retrieved", CAA_id,indexArray
                return(indexArray[2])

###############################################################################

class MyGetSequence: # global fh_nid_ptrin/datain
	"""emulates gtglib's g.getSequence object

	x=MyGetSequence(nid) -> x.front, x.back, x.aasequence()
	
	"""
	def __init__(self,nid):
		address=(nid-1)*16*2 # 64-bit
		fh_nid_ptrin.seek(address)
		indexArray=array.array('L') 
		indexArray.fromfile(fh_nid_ptrin,4)
		self.front=indexArray[0]
		self.back=indexArray[1]
		self.seqstart=indexArray[2]
		self.nres=indexArray[3]
		#print "MyGetSequence(%i): " %nid,indexArray
		
	def aasequence(self):
		fh_nid_datain.seek(self.seqstart)
		seqArray=array.array('c')
		seqArray.fromfile(fh_nid_datain,self.nres)
		aasequence=seqArray.tostring()
		return(aasequence)

class MyVertex: 
	"""
	x=MyVertex(leaf vertex) -> x.asClusterIndex, x.asNid, x.asResidue, x.arrayElement (points to cluster topology table) 
	"""
	def __init__(self,vertex):
		address=vertex*16*2 # 64-bit
		fh_vertex_ptrin.seek(address)
		indexArray=array.array('L')
		try: indexArray.fromfile(fh_vertex_ptrin,4)
		except: 
#			print "#MyVertex error (MyVertex)",address,vertex	# 137727259, 2203636144
			indexArray.fromlist([0,0,0,0]) # HACK
		self.asClusterIndex=indexArray[0]
		self.arrayElement=indexArray[1]
		self.asNid=indexArray[2]
		self.asResidue=indexArray[3]

def input_one_line():
#	line=raw_input("Query,v4_nid,query_starts,sbjct_starts,blocklengths,v2_rep40,query_starts,sbjct_starts,blocklengths,query_sequence:")
#	print
	line=raw_input()
	return(line)

def get_alignment(query_starts,sbjct_starts,block_length):
	# ali1[query_residue]=sbjct_residue
	# ali2[sbjct_residue]=query_residue
	ali1={}
	ali2={}
	n=len(block_length)
	for i in range(0,n):
		q=query_starts[i]
		s=sbjct_starts[i]
		for k in range(0,block_length[i]):
			q+=1
			s+=1
			ali1[q]=s
			ali2[s]=q
	return(ali1,ali2)

def get_attributes_intracluster(query_name,gtg_nid,ali2,query_sequence):
	gtg_sequence=MyGetSequence(gtg_nid)
        front=gtg_sequence.front
        back=gtg_sequence.back
        for gtg_vertex in range (front,back):
		x=MyVertex(gtg_vertex)
		sbjct_residue=x.asResidue
		if not ali2.has_key(sbjct_residue): continue # not aligned
		cluster_id=x.asClusterIndex
		array_element=x.arrayElement
		query_residue=ali2[sbjct_residue]
		print "%s\t%i\t%s\t%i\t%i\t%i\t%i\t%i" %(query_name,query_residue,query_sequence[query_residue-1],gtg_nid,sbjct_residue,gtg_vertex,cluster_id,array_element)

def get_attributes(query_name,gtg_nid,ali2,query_sequence):
	gtg_sequence=MyGetSequence(gtg_nid)
        front=gtg_sequence.front
        back=gtg_sequence.back
	#print "get_attributes: gtg_nid=%i front=%i back=%i nres=%i" %(gtg_nid,front,back,len(query_sequence))
	nres=len(query_sequence)
        for gtg_vertex in range (front,back):
		x=MyVertex(gtg_vertex)
		sbjct_residue=x.asResidue
		if not ali2.has_key(sbjct_residue): continue # not aligned
		query_residue=ali2[sbjct_residue]
		if query_residue>nres: break # security hack
		query_aaix=aa2ix[query_sequence[query_residue-1]]
		tuplelist=getCAA(gtg_vertex) # match CAA-aaix
		if len(tuplelist)>param_exclude_CAAs_larger: continue # exclude oversized CAAs
		for CAA_vertex,CAA_aaix in tuplelist:
			#print "# ",CAA_vertex,CAA_aaix,gtg_vertex,query_aaix,sbjct_residue,query_residue
			if CAA_aaix<>query_aaix: continue
			print "%s\t%i\t%s\t%i\t%i\t%i\t%i" %(query_name,query_residue,query_sequence[query_residue-1],gtg_nid,sbjct_residue,CAA_aaix,CAA_vertex)

def get_simple_attributes(query_name,gtg_nid,ali2,query_sequence):
	# attribute == clusid.aaix
        gtg_sequence=MyGetSequence(gtg_nid)
        front=gtg_sequence.front
        back=gtg_sequence.back
        for gtg_vertex in range (front,back):
                x=MyVertex(gtg_vertex)
                sbjct_residue=x.asResidue
                if not ali2.has_key(sbjct_residue): continue # not aligned
                query_residue=ali2[sbjct_residue]
                query_aaix=aa2ix[query_sequence[query_residue-1]]
		clusid=x.asClusterIndex
		print "%s\t%i\t%s\t%i\t%i\t%i\t%i" %(query_name,query_residue,query_sequence[query_residue-1],gtg_nid,sbjct_residue,query_aaix,clusid)

def get_consensus_anc40(query_name,gtg_nid,ali2,query_sequence):
	# get anc40 profile of each vertex in query
        gtg_sequence=MyGetSequence(gtg_nid)
        front=gtg_sequence.front
        back=gtg_sequence.back
        for gtg_vertex in range (front,back):
                x=MyVertex(gtg_vertex)
                sbjct_residue=x.asResidue
                if not ali2.has_key(sbjct_residue): continue # not aligned
                query_residue=ali2[sbjct_residue]
                query_aaix=aa2ix[query_sequence[query_residue-1]]
                clusid=x.asClusterIndex
		irow=x.arrayElement
		x=MyCluster(clusid)
		if x.nrecords<1: continue
		table=x.getTable()
		s_aafreqptr=e.get_anc40_aafreqnrecords(irow,table)
		s_aacomp=e.get_aacomp_from_aafreqptr(s_aafreqptr)
                c='X'
		if s_aacomp[s_aacomp[0]]>50: c=alfabet[s_aacomp[0]]
		print "%s\t%i\t%s\t%i\t%i\t%i\t%i\t%s\t" %(query_name,query_residue,query_sequence[query_residue-1],gtg_nid,sbjct_residue,query_aaix,clusid,c),
		for x in s_aacomp[1:]: print "%i " %x,
		print

if __name__ == "__main__":
	USAGE="""
USAGE: python gtg_attributes.py <CAA_level> <DATA_DIR>

where
	CAA_level=1 outputs CAA1-attributes (multiple values per residue)
        CAA_level=3 outputs clusid (single value per residue)
        CAA_level=4 outputs both clusid and anc30 aa-profile

"""

	try:
		DATADIR=sys.argv[2] # '/data/liisa/'
		DATADIR+='/' # just in case
		param_CAA_level=int(sys.argv[1]) # 0 or 1 or 2 or 3 (3=simple clusid.aaix attributes)
	except:
                print USAGE
                sys.exit(1)

	param_exclude_CAAs_larger=10000
	param_anc_threshold=30

	aa2ix={'A':1, 'C':2,'D':3,'E':4,'F':5,'G':6,'H':7,'I':8,'K':9,'L':10,'M':11,'N':12,'P':13,'Q':14,'R':15,'S':16,'T':17,'V':18,'W':19,'Y':20,'X':0, 'Z': 0, 'B': 0,
		'a':2,'b':2,'c':2,'d':2,'e':2,'f':2,'g':2,'h':2,'i':2,'j':2,'k':2,'l':2,'m':2,'n':2,'o':2,'p':2,'q':2,'r':2,'s':2,'t':2,'u':2,'v':2,'w':2,'x':2,'y':2,'z':2}
	alfabet='XACDEFGHIKLMNPQRSTVWY'
	
	# connect databases
	#print "DATADIR=%s" %DATADIR
	fh_nid_ptrin=open(DATADIR+'nidindex.ptr','rb')
	fh_nid_datain=open(DATADIR+'nidindex.store','rb')
	fh_vertex_ptrin=open(DATADIR+'vertex.ptr','rb') 

	# CAA databases
        if param_CAA_level==0:
                fh_CAA_vertex_datain=open(DATADIR+"CAA0_leafs.store","rb")
                fh_CAA_vertex_ptrin=open(DATADIR+"CAA0_leafs.ptr","rb")
                fh_CAA_member_datain=open(DATADIR+"CAA0_members.store","rb")
                fh_CAA_member_ptrin=open(DATADIR+"CAA0_members.ptr","rb")
        elif param_CAA_level==1:
                fh_CAA_vertex_datain=open(DATADIR+"CAA1_leafs.store","rb")
                fh_CAA_vertex_ptrin=open(DATADIR+"CAA1_leafs.ptr","rb")
                fh_CAA_member_datain=open(DATADIR+"CAA1_members.store","rb")
                fh_CAA_member_ptrin=open(DATADIR+"CAA1_members.ptr","rb")
        elif param_CAA_level==2:
                fh_CAA_vertex_datain=open(DATADIR+"CAA2_leafs.store","rb")
                fh_CAA_vertex_ptrin=open(DATADIR+"CAA2_leafs.ptr","rb")
                fh_CAA_member_datain=open(DATADIR+"CAA2_members.store","rb")
                fh_CAA_member_ptrin=open(DATADIR+"CAA2_members.ptr","rb")
	elif param_CAA_level>3:
		e=MyAacomp()
		c=MyReadTable()
		fh_aafreq_datain=open(DATADIR+"aafreq.store","rb")
		fh_cluster_ptrin=open(DATADIR+'cluster.ptr',"rb")
        	fh_cluster_datain=open(DATADIR+'cluster.store',"rb")

	while(1):
		try: line=input_one_line()
		except: break

		data=line.split("\t")
		query_name=data[0]
		gtg_nid=int(data[2])
		#except: continue
		#print "gtg_nid: ",gtg_nid
		if gtg_nid<1: continue # no mapping
		query_starts=map(int,data[3].split(None))
		sbjct_starts=map(int,data[4].split(None))
		block_lengths=map(int,data[5].split(None))
		query_sequence=data[6]
		#print "query_starts: ",query_starts
		#print "sbjct_starts: ",sbjct_starts
		#print "block_lengths: ",block_lengths

		ali1,ali2=get_alignment(query_starts,sbjct_starts,block_lengths)
	        #print "ali1: ",ali1
		#print "ali2: ",ali2
		if param_CAA_level<3: get_attributes(query_name,gtg_nid,ali2,query_sequence)
		elif param_CAA_level==3: get_simple_attributes(query_name,gtg_nid,ali2,query_sequence)
		else: get_consensus_anc40(query_name,gtg_nid,ali2,query_sequence)
	
	if param_CAA_level<3:
		fh_nid_ptrin.close()
		fh_nid_datain.close()
		fh_vertex_ptrin.close()
	        fh_CAA_vertex_datain.close()
	        fh_CAA_vertex_ptrin.close()
	        fh_CAA_member_datain.close()
	        fh_CAA_member_ptrin.close()
	elif param_CAA_level>3:
		fh_cluster_datain.close()
		fh_cluster_ptrin.close()
		fh_aafreq_datain.close()
