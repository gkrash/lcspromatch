c This module/program is part of DaliLite (c) L. Holm 1999
c
	program replace
c
c	input:	 cd1,cd2,dccpfile1,dccpfile2, list0
c	output:	 prealignments for dalicon (pipe format)
c
	implicit none
	include 'parsizes.for'
	character*5 cd1,cd2,cd,list0(maxprot) ! list0 has space for all-of-PDB
        integer nblock,l1(maxres),l2(maxres),r1(maxres),r2(maxres)
        integer i,j,k,ierr,ali1(maxres),ali(maxres),getix,nprot
	character*80 filnam1,filnam2,listfile
	real zscore,zmax
	logical lfirst
	!
	! new=cd1 old=cd2
	!
c	type *,'enter code+chainid of new, old'
	read(*,500) cd1
	read(*,500) cd2
c	type *,'enter DCCP filenames of new,old'
	read(*,510) filnam1
	read(*,510) filnam2
c	type *,'enter list of valid codes'
	read(*,510) listfile
	!
	! read list0
	!
	nprot=0
	open(90,file=listfile,status='old')
10	read(90,500,end=19) cd
	nprot=nprot+1
	if(nprot.gt.maxprot) stop 'replace: pdb size exceeded'
	list0(nprot)=cd
	goto 10
19	close(90)
	!
	! step 1: find cd1-cd2 in cd1.dccp
	!
	open(90,file=filnam1,status='old',err=999)
c	type *,'read ',filnam1
	ierr=0
	zmax=0.0
	do while(ierr.eq.0)
	  call read_nextDCCP(cd2,cd,zscore,nblock,l1,r1,l2,r2,ierr,90)
c	  type *,cd2,cd,nblock,ierr,zscore,zmax,cd1,getix(cd2,list0,nprot)
	  if(ierr.eq.0) then
C	    if(getix(cd2,list0,nprot).gt.0) then
		if(cd.eq.cd1) then		! select cd2-cd1 pair
		  if(zscore.gt.zmax) then	! keep best-z ali only
			zmax=zscore
			do i=1,maxres
				ali1(i)=0
			end do
			do i=1,nblock		! cd2-cd1
c				type *,'ali1',i,l1(i),r1(i),l2(i),r2(i)
				k=l2(i)-l1(i)
				do j=l1(i),r1(i)
					ali1(j)=k+j
				end do
			end do
		  end if
		end if
c	    end if
	  end if
	end do
	close(90)
	!
	! step 2: loop through each ali in cd2.dccp
	! write directly in DCCP format 
	!
	lfirst=.true.
	open(90,file=filnam2,status='old',err=999)
c	type *,'read ',filnam2
	ierr=0
	do while(ierr.eq.0)
	  call read_nextDCCP(cd2,cd,zscore,nblock,l1,r1,l2,r2,ierr,90)
c	  type *,cd1,cd,nblock,ierr
	  if(ierr.eq.0) then
		do i=1,maxres
			ali(i)=0
		end do
c		do i=1,nblock
c			type *,'ali',i,l1(i),r1(i),l2(i),r2(i)
c		end do
		do i=1,nblock			! cd2-x
			k=l2(i)-l1(i)
			do j=l1(i),r1(i)
				ali(j)=k+j
			end do
		end do
		!
		! generate cd1-x from cd2-cd1 & cd2-x
		!
		nblock=0
		do i=1,maxres
		  if(ali(i).ne.0.and.ali1(i).ne.0) then
			nblock=nblock+1
			l1(nblock)=ali1(i)
			r1(nblock)=ali1(i)
			l2(nblock)=ali(i)
			r2(nblock)=ali(i)
		  end if
		end do
c		do i=1,nblock
c			type *,'new',i,l1(i),r1(i),l2(i),r2(i)
c		end do
c		type *,'mergeblocks: ',cd1,cd2
		call mergeblocks(nblock,l1,r1,l2,r2)
		call output(cd1,cd,nblock,l1,r1,l2,r2,zscore,91) 
	  end if
	end do
	close(90)

500	format(a5)
510	format(a80)

999	end
c	
c-------------------------------------------------------------------------------c
        function getix(cd,list,nprot)
        implicit none
        integer getix,nprot
        character*5 list(nprot),cd
c
        integer i
c
        getix=0
        do i=1,nprot
                if(list(i).eq.cd) then
                        getix=i
                        return
                end if
        end do
  
        return
        end
c
c----------------------------------------------------------------------------
c
	subroutine mergeblocks(nblock,l1,r1,l2,r2)
	implicit none
	include 'parsizes.for'
	integer nblock,l1(maxres),r1(maxres),l2(maxres),r2(maxres)
	integer i,j,k

	! nothing to do
	if(nblock.eq.1) return
	! merge blocks
	do i=2,nblock
		if(l1(i).eq.r1(i-1)+1.and.l2(i).eq.r2(i-1)+1) then
			j=i-1
c			type 500,i,j,l1(j),r1(j),l2(j),r2(j),
c     $					l1(i),r1(i),l2(i),r2(i)
			l1(i)=l1(i-1)
			l2(i)=l2(i-1)
			l1(i-1)=0
			l2(i-1)=0
		end if
	end do
	k=nblock
	nblock=0
	do i=1,k
		if(l1(i).ne.0) then
			nblock=nblock+1
			l1(nblock)=l1(i)
			r1(nblock)=r1(i)
			l2(nblock)=l2(i)
			r2(nblock)=r2(i)
		end if
	end do

500	format('merge',10i6)

	return
	end
c	
c-------------------------------------------------------------------------------c
 	subroutine output(cd1,cd2,nblock,l1,r1,l2,r2,zscore,iunit)
	implicit none
	include 'parsizes.for'
	character*5 cd1,cd2
        integer nblock,l1(maxres),l2(maxres),r1(maxres),r2(maxres)
	real zscore
	integer i,iunit
c
	write(iunit,601) 'DCCP   1 ',-9.9,-9.9,-9,
     $   	zscore,-9,nblock,cd1,cd2
	write(iunit,605) 'alignment'
	write(iunit,610) (l1(i),r1(i),i=1,nblock)
	write(iunit,610) (l2(i),r2(i),i=1,nblock)
 
600     format(1x,a9,f8.1,f4.1,i4,16x,i4,4x,i3,16x,a5,1x,a5)
601     format(1x,a9,f8.1,f4.1,i4,f8.1,8x,i4,4x,i3,16x,a5,1x,a5)
605     format(1x,a9)
610     format(8(i4,2x,i4))
 
	return
	end
c
c------------------------------------------------------------------------------
c
        subroutine read_nextDCCP(cd1,cd2,zscore,
     $          nblock,l1,r1,l2,r2,ierr,iunit)
        implicit none
        include 'parsizes.for'
c
c       reads dali-scores & Z-scores from domainparser output !
c
        character*5 cd1,cd2,c1,c2
        character*9 line
        real score,zscore
        integer nblock,l1(maxres),l2(maxres),r1(maxres),r2(maxres)
        integer i,ierr,iunit
        integer ide,lali
        real rmsd
c
cDCCP   1   1651.8 0.0 175    24.6         100      1                1hgeB 1hgeB
c
100     read(iunit,601,err=100,end=999) line,score,rmsd,lali,zscore,
     $          ide,nblock,c1,c2
        if(line(1:9).eq.'DCCP   1 ') then
                if(c1.eq.cd1) then      ! expected order
10                      read(iunit,605,err=100,end=999) line
                        if(line(1:9).ne.'alignment') goto 10
                        read(iunit,610,err=100) (l1(i),r1(i),i=1,nblock)
                        read(iunit,610,err=100) (l2(i),r2(i),i=1,nblock)
                        cd2=c2
                else if(c2.eq.cd1) then ! reverse pair !
20                      read(iunit,605,err=100,end=999) line
                        if(line(1:9).ne.'alignment') goto 20
                        read(iunit,610,err=100) (l2(i),r2(i),i=1,nblock)
                        read(iunit,610,err=100) (l1(i),r1(i),i=1,nblock)
                        cd2=c1
                else 
                        goto 100
                end if
        else
                goto 100
        end if
 
600     format(1x,a9,f8.1,f4.1,i4,16x,i4,4x,i3,6x,2i4,2x,a5,1x,a5)
601     format(1x,a9,f8.1,f4.1,i4,f8.1,8x,i4,4x,i3,16x,a5,1x,a5)
605     format(1x,a9)
610     format(8(i4,2x,i4))
 
!       normal exit
        ierr=0
        return
 
!       error exit
999     ierr=-1
        return
        end
c	
c-------------------------------------------------------------------------------c
