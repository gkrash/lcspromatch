c This module/program is part of DaliLite (c) L. Holm 1999
c
c
c----------------------------------------------------------------------
c
	function findfile(code,listfile,iunit)
	implicit none
	character*80 listfile,findfile
	character*4 code
	integer iunit
c
	integer i
	character*4 c
	character*80 f
	logical lfound
c
c	generic function to scan listfile for code and read complete filename
c	format of listfile is a4,5x,a80
c	use e.g. ~/lib/dssp.list or ~/lib/pdb.list
c
	lfound=.false.
	do i=1,80
		findfile(i:i)=' '
	end do
	write(*,500) code, listfile
	open(iunit,file=listfile,status='old')
	i=0
10	read(iunit,510,end=19) c,f
	i=i+1
	if(c.eq.code) then
		findfile(1:60)=f
		write(*,520) i,f
		lfound=.true.
		goto 19
	end if
	goto 10
19	if(.not.lfound) write(*,530) i
	close(iunit)

500	format(/'Looking for ',a4,' in ',a60)
510	format(a4,5x,a80)
520	format('position:',i4,2x,a60)
530	format('No match found in ',i4,' entries')

	return
	end
c
c-----------------------------------------------------------------------------
c
	subroutine setup(code1,chainid1,nres1,seq1,struc1,relacc1,resno1,
     $		ca1,d1,hdr1,ierr,infile1,lpreali,pr1,pr2,nfrag,nres0,tsil,
     $		dalidatpath)
	implicit none
	include 'gagasizes.for'
	character*4 code1,resno1(maxres0)
	character*80 infile1,hdr1,findfile,dalidatpath
	character chainid1,seq1(maxres0),struc1(maxres0)
	integer acc1(maxres0),nres1,relacc1(maxres),nres0,tsil(maxres)
	real ca1(3,maxres0)
	integer*2 d1(maxres,maxres)
	integer ierr
c
	character c1,c2
	integer iseg,nseg,range(2,1000),pr1(2,1000),pr2(2,1000)
	integer nfrag,i,ires,j
	logical lpreali
	character*5 cd1
c
	ierr=0
	nfrag=0
	write(*,*)
     $	' enter code+chain identifier, preali/range flags (a4,a1,a1,a1) '
	read(*,500) code1,chainid1,c1,c2
	write(*,*) 'setup: ',code1,chainid1,c1,c2
	lpreali=(c1.eq.'*')
	if(code1(1:3).eq.'END') then
		ierr=-99
		return
	end if
!	infile1=findfile(code1,findfiledefault,iounit1)
c
c	prealignment (only makes sense with second protein input)
c	example:
c		1mli
c		1ndkA*
c		2 1 10 12 50 11 20 42 60
c		1aps *
c		1 1 50 3 52
c		END
c
	if(lpreali) then
		write(*,*) 'enter prealignment'
		read(*,*) nfrag,(pr1(1,i),pr1(2,i),i=1,nfrag),
     $			(pr2(1,i),pr2(2,i),i=1,min(1000,nfrag))
		if(nfrag.gt.1000) then
			write(*,*) 'WARNING: truncated to 1000 fragments'
			nfrag=1000
		end if
		write(*,*) 'using prealignment',nfrag
		write(*,510) (pr1(1,i),pr1(2,i),i=1,nfrag)
		write(*,510) (pr2(1,i),pr2(2,i),i=1,nfrag)
	end if
c
c	domain ranges
c
	if(c2.eq.'*') then
		write(*,*) 'enter ranges'
		read(*,*) nseg,(range(1,iseg),range(2,iseg),
     $			iseg=1,min(100,nseg))
		if(nseg.gt.100) then
			write(*,*) 'WARNING: truncated to 100 segments'
			nseg=100
		end if
		write(*,*) 'using range(s)',nseg
		write(*,510) (range(1,iseg),range(2,iseg),iseg=1,nseg)
	end if
c
	cd1(1:4)=code1
	cd1(5:5)=chainid1
	write(*,*) 'call setup_new ',cd1,nres0,' from ',code1,chainid1
	call setup_new(cd1,nres0,ca1,99,dalidatpath,seq1)
	write(*,*) 'back from setup_new ',cd1,nres0
	if(nres0.gt.maxres) nres0=maxres

	if(ierr.eq.1) goto 99
	do i=1,80
		hdr1(i:i)=' '
	end do
	do i=1,nres0
		struc1(i)=' '
		resno1(i)='    '
	end do
	if(c2.ne.'*') then
		nseg=1
		range(1,1)=1
		range(2,1)=nres0
	end if
	nres1=0
	do iseg=1,nseg
		do ires=range(1,iseg),range(2,iseg)
			nres1=nres1+1
			if(nres1.gt.maxres) then
				write(*,*) ' too many residues -- skip at',nres1
				nres1=maxres
				goto 19
			end if
			tsil(nres1)=ires
			do j=1,3
				ca1(j,nres1)=ca1(j,ires)
			end do
		end do
	end do
	write(*,*) 'in setup nres0,nres1 ',nres0,nres1
19	call getdist(ca1,nres1,d1,maxres)
c
500	format(a4,a1,a1,a1)
510	format(20i4)
c
99	return
	end
c
c----------------------------------------------------------------------
c
	subroutine getrelacc(nres,acc,relacc,seq)
c
c	input DSSP acc, return relative accessibility*100 % in relacc
c
	implicit none
	integer nres,acc(nres),relacc(nres)
	character seq(nres)
c
	integer i
	character aa
	real x,reference

	do i=1,nres
		aa=seq(i)
		if(aa.ge.'a') aa='C'
		x=reference(aa)
		if(x.gt.0.0) then
			relacc(i)=nint(100.0*acc(i)/x)
		else
			relacc(i)=100
		end if
	end do

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine getdist(ca,nres,d,maxres)
	implicit none
	integer nres,maxres
	real ca(3,nres)
	integer*2 d(maxres,maxres),distance
c
	integer i,j

	do i=1,min(maxres,nres)
		d(i,i)=0.0
		do j=i+1,nres
			d(i,j)=distance(ca(1,i),ca(2,i),ca(3,i),
     $				ca(1,j),ca(2,j),ca(3,j))
			d(j,i)=d(i,j)
		end do
	end do

	return
	end
c
c----------------------------------------------------------------------
c
	function hashfilnam(code,chainid)
	implicit none
	character chainid,fchain
	character*4 code
	character*80 hashfilnam
c
	integer i

	do i=1,80
		hashfilnam(i:i)=' '
	end do
	hashfilnam(1:4)=code
	hashfilnam(5:5)=fchain(chainid)
	hashfilnam(6:10)='.hash'

	return
	end
c
c----------------------------------------------------------------------
c
	function pairfilnam(code1,chainid1,code2,chainid2,extension)
	implicit none
	character chainid1,chainid2,fchain
	character*4 code1,code2,extension
	character*80 pairfilnam
c
	integer i

	do i=1,80
		pairfilnam(i:i)=' '
	end do
	pairfilnam(1:4)=code1
	pairfilnam(5:5)=fchain(chainid1)
	pairfilnam(6:9)=code2
	pairfilnam(10:10)=fchain(chainid2)
	pairfilnam(11:11)='.'
	pairfilnam(12:15)=extension

	return
	end
c
c----------------------------------------------------------------------
c
	function fchain(chainid)
	implicit none
	character fchain,chainid

	fchain=chainid
	if(chainid.eq.' ') fchain='_'
	if(chainid.eq.'*') fchain='_'

	return
	end
c
c----------------------------------------------------------------------
c
	function getp(x)
	implicit none
	real getp,x

	if(x.gt.0.0) then
		getp=1.0
	else if(x.lt.-5) then
		getp=0.0
	else
		getp=exp(x)
	end if

	return
	end
c
c----------------------------------------------------------------------------=
c
	function file_exists(infile)
	implicit none
	include 'gagasizes.for'
	logical file_exists
	character*80 infile
c
	file_exists=.false.
	if(machine.eq.'V') then
		open(iounit1,file=infile,err=19,status='old')
	else
		open(iounit1,file=infile,err=19,status='old')
	end if
	file_exists=.true.
19	close(iounit1)

	return
	end
c
c	---------------------------------------------------------------------
c
	function code_parser(infile, machine)
	implicit none
	character*80 infile,code
	character*4 code_parser
	character machine
c
	integer i,j,k,l

c	parse infile to get code: j is start of code, k is end of code
c
c	VAX syntax {dir:}{[.subdir{.subdir}]}code{.brk}
c                      j   k       k       j      k  l
c	Unix syntax {.,/dir/subdir/}code{.suffixes}
c                                 j      k       l
c
	if (machine.eq.'V') then
		do i=1,80
			code(i:i)=' '
		end do
		j=1
		i=1
		do while(infile(i:i).eq.' ')
			i=i+1
		end do
		l=i-1
			do i=1,len(infile)
			if(infile(i:i).eq.':') j=i+1
			if(infile(i:i).eq.']') j=i+1
			if(infile(i:i).eq.'.') k=i-1
			if(infile(i:i).ne.' ') l=l+1
		end do
		if(k.eq.0) k=l
		l=0
		do i=j,k
			l=l+1
			code(l:l)=infile(i:i)
		end do
	else
		do i=1,80
			code(i:i)=' '
		end do
		j=1
		i=1
		do while(infile(i:i).eq.' ')
			i=i+1
		end do
		l=i-1
		k=0
		do i=1,len(infile)
			if(infile(i:i).eq.'/') j=i+1
			if(infile(i:i).eq.'.') k=i-1
			if(infile(i:i).ne.' ') l=l+1
		end do
		if(k.eq.0) k=l
		l=0
		do i=j,k
			l=l+1
			code(l:l)=infile(i:i)
		end do
	end if
c
c	keep 4 characters
c
	do i=1,min(4,l)
		code_parser(i:i)=code(i:i)
	end do

	return
	end
c
c----------------------------------------------------------------------
c
	function distance(a1,a2,a3,b1,b2,b3)
	implicit none
	real a1,a2,a3,b1,b2,b3
	integer*2 distance

	distance=nint(10.0*sqrt((a1-b1)*(a1-b1)+(a2-b2)*(a2-b2)+
     $		(a3-b3)*(a3-b3)))

	return
	end
c
c-----------------------------------------------------------------------
	subroutine xstring(hdr,lhdr,line)
c
c	extracts useful information from PDB compnd & source records
c
	implicit none
	character*80 hdr,line
	integer lhdr
c
	integer i

	i=6
	do while(i.lt.80)
		i=i+1
		if(line(i:i).eq.'(') then
			do while ((line(i:i).ne.')').and.(i.lt.80))
				i=i+1
			end do
			if(line(i:i).eq.')') i=i+1
		end if
		if((line(i-1:i).ne.'  ').and.(lhdr.lt.80)) then
			lhdr=lhdr+1
			hdr(lhdr:lhdr)=line(i:i)
		end if
	end do

	return
	end
c
c-----------------------------------------------------------------------
c
	function reference(aa)
	implicit none
	character aa
	real reference

	reference=-1.0
	if(aa.eq.'A') then
		reference=106.3
	else if(aa.eq.'P') then
		reference=135.9
	else if(aa.eq.'S') then
		reference=123.1
	else if(aa.eq.'C') then
		reference=138.5
	else if(aa.eq.'T') then
		reference=141.7
	else if(aa.eq.'V') then
		reference=148.4
	else if(aa.eq.'I') then
		reference=171.5
	else if(aa.eq.'L') then
		reference=163.6
	else if(aa.eq.'D') then
		reference=149.5
	else if(aa.eq.'N') then
		reference=149.8
	else if(aa.eq.'H') then
		reference=182.2
	else if(aa.eq.'F') then
		reference=200.9
	else if(aa.eq.'Y') then
		reference=212.4
	else if(aa.eq.'W') then
		reference=245.4
	else if(aa.eq.'M') then
		reference=193.7
	else if(aa.eq.'E') then
		reference=182.8
	else if(aa.eq.'Q') then
		reference=186.6
	else if(aa.eq.'K') then
		reference=200.8
	else if(aa.eq.'R') then
		reference=239.5
	else if(aa.eq.'G') then
		reference=83.6
	end if

	return
	end
c
c----------------------------------------------------------------------
c


C======================================================================
C this library contains subroutines which are calling system specific
C things, like get the actual date, time, open a file etc.
C ===> have one system-lib.for for the VMS, UNIX is nix.... machines
C      and link them.
C======================================================================

C======================================================================
c SUBROUTINE GETDATE RS89
c returns date in a string of implied length
c UNIX version
      subroutine getdate(date)
      character date*(*)
      character ctemp*24
      character day*2, month*3, year*2

	     call fdate(ctemp)
 	     month = ctemp(5:7)
 	     day = ctemp(9:10)
 	     year = ctemp(23:24)
  	    date = (((day // '-') // month) // '-') // year

      return
      end
c..END GETDATE................................................
C======================================================================
C this library contains subroutines which are calling system specific
C things, like get the actual date, time, open a file etc.
C ===> have one system-lib.for for the VMS, UNIX is nix.... machines
C      and link them.
C======================================================================


C======================================================================
c	SUBROUTINE GETDATE(CDATE)
c returns date in a string of implied length -- VAX version
c	CHARACTER CDATE*(*)
c	CHARACTER*9 CTEMP
c	CTEMP=' '
c	CALL DATE(CTEMP)
c	CDATE(1:9)=CTEMP(1:9)
c	RETURN
c	END
C======================================================================

	function val(instring,s)
c
c	converts first digit-segment of instring to an integer, s is sign
c
	implicit none
	integer val,ic,s
	character*(*) instring
c
	integer k,l,r,v
	real m
c
c	ignore leading non-digits
c
	l=1
	do while(((instring(l:l).lt.'0').or.(instring(l:l).gt.'9'))
     $		.and.(l.le.len(instring)))
		l=l+1
	end do
c
c	check sign
c
	s=1
	if(l.gt.1) then
		if(instring(l-1:l-1).eq.'-') s=-1
	end if
c
c	ignore trailing non-digits
c
	r=l
	do while((instring(r+1:r+1).ge.'0').and.(instring(r+1:r+1).le.'9')
     $		.and.(r.lt.len(instring)))
		r=r+1
	end do
c	write(*,*) l,r
	m=0.1
	v=0
	do k=r,l,-1
		m=m*10
		v=v+ic(instring(k:k))*m
	end do
	val=v

	return
	end
c
c-----------------------------------------------------------------------------
c
	function ic(a)
	implicit none
	character a,t(0:9)
	integer ic,i
	data t/'0','1','2','3','4','5','6','7','8','9'/

	ic=0
	do i=0,9
		if(a.eq.t(i)) ic=i
	end do

	return
	end
c
c-----------------------------------------------------------------------
c

	subroutine readdssp(filnam,chainid,iunit,
     $		compnd,nres,resno,seq,struc,ca,acc,ierr)
c
c	frontend to getdssp to transport wanted subset of DSSP information
c
c	parameter:
c		maxres
c
c	input:	filnam
c		chainid
c		iunit
c
c	output:	header, compnd, source, author -- echo PDB cards
c		nres
c		nchain
c		nss, nssintra, nssinter -- SS-bridges
c		area -- total accessible surface area
c		nhb, nhbp -- # backbone H-bonds, same per-cent
c		resno -- PDB residue number string with appended insertion
c			 indicator, e.g. 1ACX 63A
c		seq -- one-letter sequence
c		struc -- one-letter secondary structure
c		acc -- accessible surface area per residue
c		dsspinfo -- structure string with renumbered beta-sheet resnos
c				* marks H-bonds to outside given chain
c		dsspresno, fsspresno -- local work arrays
c		tco,kappa,alfa,phi,psi -- backbone conformational angles
c		ca -- CA coordinates
c		ierr = 0 -- normal
c		     = -1 -- file not found
c		     = -2 -- too many residues, chain truncated to maxres
c
	implicit none
	include 'gagasizes.for'
	integer iunit,nres,nchain,nss,nssintra,nssinter,nhb,acc(maxres0)
	character*80 filnam,header,compnd,source,author
	character chainid,seq(maxres0),struc(maxres0)
	real nhbp,tco(maxres0),kappa(maxres0),alfa(maxres0),phi(maxres0),
     $		psi(maxres0),ca(3,maxres0)
	character*4 resno(maxres0)
	character*33 dsspinfo(maxres0)
	integer area,ierr,dsspresno(maxres0)
c
	call getdssp(filnam,chainid,iunit,maxres0,header,compnd,source,
     $		author,nres,nchain,nss,nssintra,nssinter,area,
     $		resno,seq,struc,acc,dsspinfo,dsspresno,
     $		tco,kappa,alfa,phi,psi,ca,ierr)
	if(ierr.eq.0) write(*,*) nres,' residues read ',filnam(1:60)
	if(ierr.eq.-1) write(*,*) ' DSSP file not found ',filnam(1:60)
	if(ierr.eq.-2) write(*,*) ' DSSP file truncated at ',maxres,'residues'

	return
	end
c
c-----------------------------------------------------------------------
c
c
c----------------------------------------------------------------------
c
        subroutine setup_new(cd1,nres,ca,iunit,dalidatpath,seq)
        implicit none
        include 'gagasizes.for'
        integer iunit,nres
        real ca(3,maxres)
        character*5 cd1
        character seq(maxres)
	character*80 dalidatpath
c
        character*80 filnam,line,constructfilnam
        character node_type(maxdom),tag1
        integer i,j,k,idom,nseg,ndom
	real x,y,z
	character*10 tag
c
        filnam=constructfilnam(cd1,dalidatpath,'.dat')
	write(*,*) 'cd1,dalidatpath,filnam ',cd1,dalidatpath,filnam
        open(iunit,file=filnam,status='old',err=219)
	read(iunit,530,end=219) line
        read(line(11:20),*) nres,nseg
	write(*,*) 'in setup_new nres,nseg ',nres,nseg
	if(nres.gt.maxres) then 
		write(*,*) 'WARNING: too many residues ',nres,
     $			' > ',maxres,' in ',cd1	
	end if
	do j=1,nseg
 		read(iunit,530) line
	end do
	read(iunit,*) ((ca(k,j),k=1,3),j=1,min(maxres,nres))
300	read(iunit,530,end=309) line
	if(line(1:1).eq.'>') then
		read(line,500) tag1,ndom
	else
		goto 300
	end if
309	write(*,*) 'tag1 = ',tag1
	do i=1,ndom
		read(iunit,530) line
	end do
310	read(iunit,530,end=319) line
	if(line(1:1).eq.'>') then
		read(line,500) tag1,ndom
	else
		goto 310
	end if
319	write(*,*) 'tag1 = ',tag1
	do i=1,ndom
		read(iunit,530) line
	end do
c       -sequence
218	read(iunit,530,end=219) line
	if(line(1:9).eq.'-sequence') then
		read(line,540) tag,(seq(i),i=1,min(maxres,nres))
	else
		goto 218
	end if
219	write(*,*) 'tag = ',tag
	close(iunit)

500     format(a1,9x,i5)
530     format(a80)
540     format(a10,1x,<maxres>a1)

        return
        end
c
c------------------------------------------------------------------------------
c
