c This module/program is part of DaliLite (c) L. Holm 1999
c
	program dccpto90
c
c	use as: grep DCCP wolfall.dccp | dccpto90
c
c	write raw scores to fort.90, Z-scores to fort.91 !
c
	implicit none
	include 'parsizes.for'
	integer i,j,nprot2,s,score(maxprot,maxprot),getix,n,ns(maxprot)
	character*5 list2(maxprot),cd1,cd2
	integer sout(maxprot,maxprot),k,l
	real x,z

	call getlist('list2',99,list2,nprot2)
c	write(*,*),(list2(i),i=1,nprot2)
	do i=1,nprot2
		ns(i)=0
		do j=1,nprot2
			score(j,i)=0
		end do
	end do

100	i=0
	j=0
        read(*,600,err=99,end=19) x,z,cd1,cd2
	i=getix(cd1,list2,nprot2)
	j=getix(cd2,list2,nprot2)
	if(i.eq.0) goto 100
	if(j.eq.0) goto 100
!	s=max(score(i,j),nint(x))	! use raw scores
	s=max(score(i,j),nint(z*10))	! use Z-scores
	score(i,j)=s
	score(j,i)=s
99	write(*,*),cd1,cd2,x,i,j
	goto 100

!DCCP   1    439.0-1.0  58    15.7          -1     58                1fxd  1fxd
600     format(1x,9x,f8.1,4x,4x,f8.1,8x,4x,4x,3x,16x,a5,1x,a5)


19	n=0
	do i=1,nprot2
		write(*,*),'selfscore: ',list2(i),score(i,i)
		if(score(i,i).gt.0) n=n+1
	end do
c
29	continue
c	write out
	k=0
	do i=1,nprot2
	  if(score(i,i).gt.0) then
		k=k+1
		l=0
		do j=1,nprot2
		  if(score(j,j).gt.0) then
			l=l+1
			sout(k,l)=score(i,j)
		  end if
		end do
	  end if
	end do
	write(*,*),k,l,n
	do i=1,2
		write(90,*) n
		do j=1,nprot2
			if(score(j,j).gt.0) write(90,510) list2(j)
		end do
	end do
	write(90,*) ((sout(j,i),j=1,n),i=1,n)

500	format(13x,a5,35x,i5,5x,i5,<maxres>i4)
510	format(2a5)

	end
c
c----------------------------------------------------------------------
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
        write(*,*) nprot,' proteins in list from ',filnam
c
500     format(a5)
c

        return
        end
c
c----------------------------------------------------------------------
c
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
c----------------------------------------------------------------------
c
