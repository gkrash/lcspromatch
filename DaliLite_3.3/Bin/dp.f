c This module/program is part of DaliLite (c) L. Holm 1999
c
c
c------------------------------------------------------------------------------
c
c	dp.f
c
c	derivative of ~/daliserv/domainparser-DCCP-WOLF.f
c	reads wolf file from standard input, passes on to DCCP file if Zmax>1.0
c
c	cat 'old' DCCP alignment, calculate Z-scores per domain-domain pairs
c
c	-> parameters: zcut, DCCP/WOLFITZ, wolfpath, Nbest [assume sorted input]
c
c------------------------------------------------------------------------------
c
	program domainparser
	implicit none
	include 'parsizes.for'
c
	integer i,nbest
	real wght(0:100),zcut
	character*4 mode
	character*80 filnam,dalidatpath_1,dalidatpath_2
c
c	parameters
c
c	type *,'enter dalidatpath_1, dalidatpath_2, mode [DCCP/WOLF], zcut, nbest'
	read(*,500) dalidatpath_1
	read(*,500) dalidatpath_2
	read(*,510) mode
	read(*,*) zcut
	read(*,*) nbest
c
c	scan through database
c
	call weights(wght)
c
	if(mode.eq.'DCCP') then
		call read_DCCP(wght,zcut,dalidatpath_1,dalidatpath_2,nbest)
	else if(mode.eq.'WOLF') then
		call read_wolf(5,wght,zcut,dalidatpath_1,dalidatpath_2,nbest)
	else
		write(*,*) 'domainparser in unknown mode'
	end if

500	format(a80)
510	format(a4)

	end
c
c------------------------------------------------------------------------------
c
	subroutine read_wolf(iunit,wght,zcut,dalidatpath_1,dalidatpath_2,nbest)
	implicit none
	include 'parsizes.for'
c
	character*5 cd1,cd2,oldcd1,oldcd2,cd
	real score
	integer nblock,l1(maxres),l2(maxres),r1(maxres),r2(maxres)
	integer i,ndom1,ndom2,domns1(maxdom),domns2(maxdom)
	integer domseglist1(maxseg,2,maxdom),domseglist2(maxseg,2,maxdom)
	real ca1(3,maxres),ca2(3,maxres),wght(0:100)
	integer*2 d1(maxres,maxres),d2(maxres,maxres)
	integer nres1,nres2,ide,lali
	real rmsd,zcut
	integer iunit,ali(maxres),nbest,ibest
	character*80 dalidatpath_1,dalidatpath_2
c
	oldcd1='?????'
	oldcd2='?????'
	ibest=0
100	read(iunit,500,end=999) cd1,cd2,nres1, (ali(i),i=1,nres1)
!		write(*,*) 'WOLFITZ INPUT ',cd1,nres1
!		write(*,*) cd1,cd2,nres1,(ali(i),i=1,nres1)
		nblock=0
		do i=1,nres1
		  if(ali(i).ne.0) then
			nblock=nblock+1
			l1(nblock)=i
			r1(nblock)=i
			l2(nblock)=ali(i)
			r2(nblock)=ali(i)
		  end if
		end do
		cd=cd1
		if(cd(5:5).eq.'_') cd(5:5)=' '
		if(cd2(5:5).eq.'_') cd2(5:5)=' '
		rmsd=-1.0
		lali=nblock
		ide=-1
		score=-infinit
		ibest=ibest+1
		call innerloop(cd,cd2,l1,l2,r1,r2,nblock,oldcd1,ndom1,domns1,
     $			domseglist1,nres1,ca1,d1,ndom2,domns2,domseglist2,
     $			nres2,ca2,d2,wght,oldcd2,score,rmsd,lali,ide,zcut,
     $			dalidatpath_1,dalidatpath_2,nbest,ibest)
!		write(*,*) 'innerloop done ',cd1,cd2,score
	goto 100
999	continue

500	format(8x,a5,a5,45x,i5,<maxres>i4)
cOLFITZ 2hsp_1nscA    1    1    2    1   25       2.5   23    5   71   0   0   0
	return
	end
c
c------------------------------------------------------------------------------
c
	subroutine read_DCCP(wght,zcut,dalidatpath_1,dalidatpath_2,nbest)
	implicit none
	include 'parsizes.for'
c
	character*5 cd1,cd2,oldcd1,oldcd2
	character*10 line
	real score
	integer nblock,l1(maxres),l2(maxres),r1(maxres),r2(maxres)
	integer i,ndom1,ndom2,domns1(maxdom),domns2(maxdom)
	integer domseglist1(maxseg,2,maxdom),domseglist2(maxseg,2,maxdom)
	real ca1(3,maxres),ca2(3,maxres),wght(0:100)
	integer*2 d1(maxres,maxres),d2(maxres,maxres)
	integer nres1,nres2,ide,lali,nbest,ibest
	real rmsd,zcut
	character*80 dalidatpath_1,dalidatpath_2
c
	oldcd1='?????'
	oldcd2='?????'
	ibest=0
100	read(*,600,err=100,end=999) line,score,rmsd,lali,ide,nblock,cd1,cd2
	if(score.lt.0) score=-infinit
	if(line(1:9).eq.'DCCP   1 '.or.line(2:10).eq.'DCCP   1 ') then
10		read(*,605,end=999) line
		if(line(1:9).ne.'alignment'.and.line(2:10).ne.'alignment')
     $			goto 10
		read(*,610,err=100) (l1(i),r1(i),i=1,nblock)
		read(*,610,err=100) (l2(i),r2(i),i=1,nblock)
		ibest=ibest+1
		call innerloop(cd1,cd2,l1,l2,r1,r2,nblock,oldcd1,ndom1,domns1,
     $			domseglist1,nres1,ca1,d1,ndom2,domns2,domseglist2,
     $			nres2,ca2,d2,wght,oldcd2,score,rmsd,lali,ide,zcut,
     $			dalidatpath_1,dalidatpath_2,nbest,ibest)
	end if
	goto 100

600	format(1x,a9,f8.1,f4.1,i4,16x,i4,3x,i4,16x,a5,1x,a5)
601	format(1x,a9,f8.1,f4.1,i4,f8.1,8x,i4,3x,i4,16x,a5,1x,a5)
605	format(1x,a9)
610	format(8(i4,2x,i4))

999	return
	end
c
c------------------------------------------------------------------------------
c
	subroutine innerloop(cd1,cd2,l1,l2,r1,r2,nblock,oldcd1,ndom1,domns1,
     $		domseglist1,nres1,ca1,d1,ndom2,domns2,domseglist2,nres2,ca2,
     $		d2,wght,oldcd2,score,rmsd,lali,ide,zcut,dalidatpath_1,
     $		dalidatpath_2,nbest,ibest)
	implicit none
	include 'parsizes.for'
	character*5 cd1,cd2,oldcd1,oldcd2
	real score
	integer nblock,l1(maxres),l2(maxres),r1(maxres),r2(maxres)
	integer i,ndom1,ndom2,domns1(maxdom),domns2(maxdom)
	integer domseglist1(maxseg,2,maxdom),domseglist2(maxseg,2,maxdom)
	real ca1(3,maxres),ca2(3,maxres),wght(0:100)
	integer*2 d1(maxres,maxres),d2(maxres,maxres)
	integer nres1,nres2,idom1,idom2,ide,lali,nbest,ibest
	real zmax,rmsd,zcut,z,x,x1
	character*80 dalidatpath_1,dalidatpath_2
	character seq1(maxres),seq2(maxres)
c
!	write(*,*) 'innerloop ',cd1,cd2,nres1,nres2,score,zcut
	if(cd1.ne.oldcd1) then
	  oldcd1=cd1
	  call setup(cd1,ndom1,domns1,domseglist1,nres1,ca1,d1,99,dalidatpath_1,
     $		seq1)
c	  file-not-found handle
	  if(ndom1.eq.0) return
	end if
	if(cd2.ne.oldcd2) then
	  oldcd2=cd2
	  call setup(cd2,ndom2,domns2,domseglist2,nres2,ca2,d2,99,dalidatpath_2,
     $		seq2)
c	  file-not-found handle
	  if(ndom2.eq.0) return
	end if
c
c	idom1 * idom2
c
	zmax=0.0
	do idom1=1,ndom1
		do idom2=1,ndom2
!			write(*,*) 'dopair ',idom1,idom2
			call dopair(idom1,idom2,domns1,domns2,
     $	domseglist1,domseglist2,l1,r1,l2,r2,nblock,d1,d2,ca1,ca2,nres1,
     $	nres2,wght,cd1,cd2,zmax,z,x)
			if(idom1.eq.1.and.idom2.eq.1) x1=x
!			write(*,*) 'dopair done ',idom1,idom2
		end do
	end do
c
c	print protein pair if any domain pair had z>zcut
c
	if(zmax.ge.zcut.or.ibest.lt.nbest) then
		if(score.eq.-infinit) score=x1	! return Daliscore in WOLF-mode
		if(lali.lt.0.or.ide.lt.0) call getide(cd1,cd2,nblock,
     $			l1,r1,l2,r2,lali,ide,seq1,seq2,nres1)
		write(*,601) 'DCCP   1 ',score,rmsd,abs(lali),
     $			zmax,ide,nblock,cd1,cd2
		write(*,605) 'alignment'
		write(*,610) (l1(i),r1(i),i=1,nblock)
		write(*,610) (l2(i),r2(i),i=1,nblock)
	end if

600	format(1x,a9,f8.1,f4.1,i4,16x,i4,3x,i4,16x,a5,1x,a5)
601	format(1x,a9,f8.1,f4.1,i4,f8.1,8x,i4,3x,i4,16x,a5,1x,a5)
605	format(1x,a9)
610	format(8(i4,2x,i4))

!	write(*,*) 'innerloop done'

	return
	end
c
c------------------------------------------------------------------------------
c
	subroutine getide(cd1,cd2,nblock,l1,r1,l2,r2,lali,intide,seq1,seq2,
     $		nres1)
	implicit none
	include 'parsizes.for'
	character*5 cd1,cd2
	integer nblock,l1(maxres),l2(maxres),r1(maxres),r2(maxres),lali
	integer intide
c
	character seq1(maxres),seq2(maxres),aligned(maxres)
	integer nres1,ali(maxres),i,j,k
	real xide
c
	write(*,*) 'getide: ',cd1,cd2!,nblock,lali,intide,nres1

	do i=1,maxres
		ali(i)=0
	end do
	do i=1,nblock
		do j=l1(i),r1(i)
			k=l2(i)+j-l1(i)
			ali(j)=abs(k)
		end do
	end do
	! aligned sequence string
	do i=1,maxres
 		k=ali(i)
		if(k.eq.0) then
			aligned(i)='.'
		else
			aligned(i)=seq2(k)
		end if
	end do
	! sequence identity
	xide=0.0
	lali=0
	do i=1,nres1
		if(ali(i).ne.0) then
			if(aligned(i).ge.'a'.and.seq1(i).ge.'a') then
                                xide=xide+100.0
			else if(aligned(i).ge.'a'.and.seq1(i).eq.'C') then
                                xide=xide+100.0
			else if(aligned(i).eq.'C'.and.seq1(i).ge.'a') then
                                xide=xide+100.0
			else if(aligned(i).eq.seq1(i)) then
                                xide=xide+100.0
                        end if
			lali=lali+1
                end if
	end do

	intide=xide/max(lali,1)

	return
	end
c
c------------------------------------------------------------------------------
c
	subroutine dopair(idom1,idom2,domns1,domns2,domseglist1,
     $		domseglist2,l1,r1,l2,r2,nblock,d1,d2,ca1,ca2,nres1,nres2,
     $		wght,cd1,cd2,zmax,z,x1)
c
c	extract idom1, idom2 from raw alignment
c	trim alignment
c	write out zscore,confined alignment
c
	implicit none
	include 'parsizes.for'
	integer nblock,l1(maxres),l2(maxres),r1(maxres),r2(maxres)
	integer idom1,idom2
	integer domns1(maxdom),domns2(maxdom),nres1,nres2
	integer domseglist1(maxseg,2,maxdom),domseglist2(maxseg,2,maxdom)
	real ca1(3,maxres),ca2(3,maxres),wght(0:100),zmax,x1
	character*5 cd1,cd2
	integer*2 d1(maxres,maxres),d2(maxres,maxres)
c
	real totscore,zscore
c
	integer i,j,ali1(maxres),lali,len1,len2
	logical lactive1(maxres),lactive2(maxres)
	real z,x
c
c	mark domains
c
!	write(*,*) 'dopair',idom1,idom2,nres1,nres2,cd1,cd2
	do i=1,nres1
		lactive1(i)=.false.
		ali1(i)=0
	end do
	len1=0
	do i=1,domns1(idom1)
		do j=domseglist1(i,1,idom1),domseglist1(i,2,idom1)
			lactive1(j)=.true.
			len1=len1+1
		end do
	end do
	do i=1,nres2
		lactive2(i)=.false.
	end do
	len2=0
	do i=1,domns2(idom2)
		do j=domseglist2(i,1,idom2),domseglist2(i,2,idom2)
			lactive2(j)=.true.
			len2=len2+1
		end do
	end do
c
c	construct ali1
c
	lali=0
	do i=1,nblock
		do j=l1(i),r1(i)
		  if(lactive1(j).and.lactive2(l2(i)+j-l1(i))) then
			ali1(j)=l2(i)+j-l1(i)
			lali=lali+1
		  end if
		end do
	end do
!	write(*,*) 'lali',lali
c
	if(lali.eq.0) then
!		write(*,*) idom1,idom2,lali,' lali'
		return
	end if
c
c	calculate score
c
	x=totscore(ali1,d1,d2,nres1,wght)
!	write(*,*) 'totscore =',x
	z=zscore(len1,len2,x)
!	write(*,*) 'zscore = ',z
!	write(*,500) z,x,lali,len1,len2,cd1,idom1,cd2,idom2
	if(idom1.eq.1.and.idom2.eq.1) x1=x

!
	if(z.gt.0.90) then	!! write out all domain-domain pairs
!				!! 	for calibration !
c		! find blocks !
 		write(12,630) x,0.0,lali,l1(1),r1(nblock),l2(1),r2(nblock),
     $			-9,0,0,nblock,idom1,idom2,nres1,nres2,cd1,cd2
!		write(12,640) 1
!		write(12,650) (l1(i),r1(i),i=1,nblock)
!		write(12,650) (l2(i),r2(i),i=1,nblock)
!
	end if
	zmax=max(z,zmax)
c
500	format('ZZZZ',f10.2,f10.1,3i5,2(2x,a5,i3))
630	format(' DCCP   1 ',f8.1,f4.1,6i4,2i2,3i3,2i4,2x,a5,1x,a5)
640	format(' alignment (resnos)',i10,':')
650	format(8(i4,' -',i4))
c
	return
	end
c
c------------------------------------------------------------------------------
c
	function totscore(ali1,d1,d2,nres1,wght)
        implicit none
        include 'parsizes.for'
	real totscore,wght(0:100)
        integer nres1
        integer*2 d1(maxres,maxres),d2(maxres,maxres)
        integer ali1(maxres)
c
        real scorefun
c
        integer i,j,k,l,q,r,a(maxres),n
        real x
c
        totscore=0.0
	n=0
	do i=1,nres1
		if(ali1(i).ne.0) then
			n=n+1
			a(n)=i
		end if
	end do
	do i=1,n
		k=a(i)
		do j=1,n
			l=a(j)
			q=abs(ali1(k))
			r=abs(ali1(l))
			x=scorefun(d1(k,l),d2(q,r),wght)
			totscore=totscore+x
		end do
	end do
	totscore=totscore
c
500	format(8(i4,' -',i4))
c
        return
        end
c
c----------------------------------------------------------------------
c
	function scorefun(a,b,wght)
	implicit none
	include 'parsizes.for'
	real scorefun,wght(0:100)
	integer*2 a,b
c
	real x,y,d0
	logical lela
	parameter(lela=.true.)
	parameter(d0=0.20)
c !!! 	elastic uses weights !!!
	x=float(abs(a-b))/10
	if(lela) then
		y=float(a+b)/20
		if(y.gt.100) then
			scorefun=0.0
		else
			if(y.gt.0) then
				scorefun=wght(nint(y))*(d0-x/y)
			else
				scorefun=wght(nint(y))*d0
			end if
		end if
	end if

	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine weights(wght)
	implicit none
	include 'parsizes.for'
	integer i
	real wght(0:100)
c
	real x,enveloperadius
c
	enveloperadius=20.0
	x=1/(enveloperadius*enveloperadius)
	do i=0,100
		wght(i)=exp(-x*i*i)
	end do

	return
	end
c
c----------------------------------------------------------------------
c
	function distanceint2(a1,a2,a3,b1,b2,b3)
	implicit none
	real a1,a2,a3,b1,b2,b3
	integer*2 distanceint2

	distanceint2=nint(10.0*sqrt((a1-b1)*(a1-b1)+(a2-b2)*(a2-b2)+
     $		(a3-b3)*(a3-b3)))

	return
	end
c
c-----------------------------------------------------------------------
c
	subroutine getdist(ca,d1,nres1)
	implicit none
	include 'parsizes.for'
	real ca(3,maxres)
	integer*2 d1(maxres,maxres),distanceint2,x
	integer nres1
c
	integer i,j
c
	do i=1,nres1
		d1(i,i)=0
		do j=1,i-1
			x=distanceint2(ca(1,i),ca(2,i),ca(3,i),
     $				ca(1,j),ca(2,j),ca(3,j))
			d1(i,j)=x
			d1(j,i)=x
		end do
	end do
c
	return
	end
c
c------------------------------------------------------------------------------
c
	function zscore(l1,l2,score)
	implicit none
	real zscore,score
	integer l1,l2
c
	real n12,mean,sigma,x
c
        n12=sqrt(float(l1*l2))
        x=min(n12,400.0)
        mean=7.9494+0.70852*x+2.5895e-4*x*x-1.9156e-6*x*x*x
	if(n12.gt.400.0) mean=mean+(n12-400.0)*1.0              ! hack !
        sigma=0.50*mean
        zscore=(score-mean)/max(1.0,sigma)


	return
	end
c
c------------------------------------------------------------------------------
c
	subroutine setup(cd1,ndom,domns,domseglist,nres,ca,d1,iunit,
     $		dalidatpath,seq)
	implicit none
	include 'parsizes.for'
	integer ndom,domns(maxdom),domseglist(maxseg,2,maxdom),iunit,nres
	real ca(3,maxres)
	character*5 cd1
	integer*2 d1(maxres,maxres)
	character seq(maxres)
c
	character*80 filnam,line,constructfilnam,dalidatpath
	character node_type(maxdom)
	integer i,j,k,idom,nseg
c
!	write(*,*) 'setup ',cd1
	do i=1,nres
		seq(i)='?'
	end do
	filnam=constructfilnam(cd1,dalidatpath,'.dat')
	open(iunit,file=filnam,status='old',err=19)
		i=0
200		read(iunit,530) line
		if(line(1:4).ne.'>>>>') goto 200
		read(line,550) ndom
		i=i+1
		if(i.eq.1) then
			read(line(11:20),*) nres,nseg
			do j=1,nseg
				read(iunit,530) line
			end do
			read(iunit,*) ((ca(k,j),k=1,3),j=1,nres)
		end if
		if(i.lt.3) goto 200
		do idom=1,ndom
			read(iunit,710,end=219) j,node_type(j),domns(j),
     $				((domseglist(i,k,j),k=1,2),i=1,domns(j))
		end do
c		-sequence
		read(iunit,540,end=219) (seq(i),i=1,nres)
219	close(iunit)
c
c	calculate distance matrix
c
	call getdist(ca,d1,nres)
c
c	keep only (+,*)-units
c
	ndom=0
	do i=1,j
		if(i.eq.1.or.node_type(i).eq.'+'.or.node_type(i).eq.'*') then
			ndom=ndom+1
			domns(ndom)=domns(i)
			do k=1,domns(i)
				domseglist(k,1,ndom)=domseglist(k,1,i)
				domseglist(k,2,ndom)=domseglist(k,2,i)
			end do
		end if
	end do
c
c	normal exit
c
	return
c
c	error exit
c
19	write(*,*) 'ERROR in fetchdomseglist: could not open file',filnam
	ndom=0
c
530	format(a80)
540	format(10x,1x,<maxres>a1)
550	format(10x,i5)
710	format(i4,1x,a1,13x,400i4)

	return
	end
c
c------------------------------------------------------------------------------
c
        subroutine getlist(filnam,iunit,list,nprot)
        implicit none
        include 'parsizes.for'
        character*5 list(maxprot)
        integer nprot,iunit
        character*(*) filnam
c
        character*5 cd
c
        nprot=0
        open(iunit,file=filnam,status='old')
10      read(iunit,500,end=19) cd
        if(nprot.eq.maxprot) then
c               write(*,*) 'WARNING: skip reading list after maxprot',maxprot
                goto 19
        end if
        nprot=nprot+1
        list(nprot)=cd
        goto 10
19      close(iunit)
!        write(*,*) nprot,' proteins in list from ',filnam
c
500     format(a5)
c

        return
        end
c
c----------------------------------------------------------------------
c
