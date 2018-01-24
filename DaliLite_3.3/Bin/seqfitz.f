c This module/program is part of DaliLite (c) L. Holm 1999
c
	program seqfitz
c
c	f77 seqfitz.f subfitz.f u3b-8.f
c
	implicit none
	include 'parsizes.for'
	character*5 cd1,cd2
	character*80 dalidatpath_1,dalidatpath_2,findfile,filnam
	character seq1(maxres),seq2(maxres)
	integer nres1,nres2
	integer i,j,score,ali(maxres),lali,iter,maxiter
	integer*2 trace(maxres,maxres,2)
	integer*2 table0(maxres,maxres)
	real rms,rcut,u(3,3),t(3),ca1(3,maxres),ca2(3,maxres)
c
	write(*,*),'enter cd1, cd2, dalidatpath_1, _2, fitzrcut, fitzmaxiter'
	read(*,500) cd1
	read(*,500) cd2
	read(*,510) dalidatpath_1
	read(*,510) dalidatpath_2
	read(*,*) rcut
	read(*,*) maxiter
c
	call setup(cd1,nres1,ca1,99,dalidatpath_1,seq1)
	call setup(cd2,nres2,ca2,99,dalidatpath_2,seq2)
c
	call fillseq(nres1,nres2,seq1,seq2,table0)
	call nw_maxsim(nres1,nres2,table0,score,ali,trace)
	lali=0
	do i=1,nres1
		if(ali(i).gt.0) lali=lali+1
	end do
c
	if(lali.gt.2) then
		call getut(nres1,ali,ca1,ca2,u,t,lali,rms)
		call transrotate(ca1,nres1,u,t)
		call fitz(ca1,nres1,ca2,nres2,ali,rcut,maxiter,rms,lali,iter)
	end if
	write(*,520),cd1,cd2,(0,j=1,5),rms,lali,iter,nres1,(ali(j),j=1,nres1)

500	format(a5)
510	format(a80)
520	format('WOLFITZ ',2a5,5i5,f10.1,3i5,<maxres>i4)

	end
c
c-----------------------------------------------------------------------
c
 	subroutine fillseq(nres1,nres2,seq1,seq2,table0)
	implicit none
	integer nres1,nres2
	character seq1(nres1),seq2(nres2)
	integer*2 table0(nres1,nres2)
c
	integer i,j
	character c
c
	do i=1,nres1
		c=seq1(i)
		do j=1,nres2
			if(seq2(j).eq.c) then
				table0(i,j)=1
			else
				table0(i,j)=0
			end if
		end do
	end do

	return
	end
c
c-----------------------------------------------------------------------
c
      function distance(a1,a2,a3,b1,b2,b3)
        implicit none
        real distance,a1,a2,a3,b1,b2,b3
c
        distance=sqrt((a1-b1)*(a1-b1)+(a2-b2)*(a2-b2)+(a3-b3)*(a3-b3))
c
        return
        end
c
c----------------------------------------------------------------------
c
        subroutine setup(cd1,nres,ca,iunit,dalidatpath,seq)
        implicit none
        include 'parsizes.for'
        integer iunit,nres
        real ca(3,maxres)
        character*5 cd1
        character seq(maxres)
c
        character*80 filnam,line,constructfilnam,dalidatpath
        character node_type(maxdom)
        integer i,j,k,idom,nseg
c
        filnam=constructfilnam(cd1,dalidatpath,'.dat')
        open(iunit,file=filnam,status='old',err=219)
                i=0
200             read(iunit,530,end=219) line
		if(line(1:4).ne.'>>>>'.and.line(1:9).ne.'-sequence') goto 200
                i=i+1
                if(i.eq.1) then
                        read(line(11:20),*) nres,nseg
                        do j=1,nseg
                                read(iunit,530) line
                        end do
                        read(iunit,*) ((ca(k,j),k=1,3),j=1,nres)
			goto 200
                end if
		if(line(1:9).ne.'-sequence') goto 200
c               -sequence
                read(iunit,540,end=219) (seq(i),i=1,nres)
219     close(iunit)

500	format(10x,i5)
530	format(a80)
540     format(10x,1x,<maxres>a1)

	return
	end
c
c------------------------------------------------------------------------------
c
