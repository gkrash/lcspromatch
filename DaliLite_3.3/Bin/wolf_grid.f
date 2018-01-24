c This module/program is part of DaliLite (c) L. Holm 1999,2008
c
	program wolf_grid
c
c	f77 wolf_grid.f
c
c       reads list1
c       writes out SSE-vectors in all frames
c
c       use in pipe to create wolf-db, to compare query versus wolf-db
c
c
	implicit none
	include 'parsizes.for'
	character secstr(maxseg)
	integer nseg
c
	integer nprot2
	character*5 cd,list2(maxprot)
c
	logical lverb,lgrid
c
	character*80 dalidatpath
c
	real yca(3,maxres),x(3,maxres)
	integer ny
	real midpoint(3,maxseg+3),direction(3,maxseg+3),neidist(maxseg,maxseg)
	real neiborcutoff, fitzrcut
c
        real r,midx(3),dirx(3)
        integer gx,gy,gz,iseg,jseg,kseg,iprot,fung,i
        character atype,ctype
c
c	write(*,*),'# enter dalidatpath, neiborcutoff'
	read(*,520) dalidatpath
	read(*,*) neiborcutoff
c	write(*,*),'# imported parameters are:',dalidatpath,neiborcutoff
c
	lverb=.false.
	nprot2=0
	call getlist('list1',90,list2,nprot2)
c
c       loop over list2: write all boxed segments
c        nocc=0
        do iprot=1,nprot2
                cd=list2(iprot)
                call setupprotein(cd,nseg,midpoint,direction,
     $                  secstr,90,
     $                  yca,ny,neidist,dalidatpath)
                if(nseg.gt.2) then
                 do iseg=1,nseg
                  atype=secstr(iseg)
                  do jseg=1,nseg
                    if(iseg.ne.jseg) then
                      r=neidist(iseg,jseg)
                      if(r.lt.neiborcutoff) then
                 call preparex(x,iseg,jseg,nseg,midpoint,direction)
                 call twist(x,3+nseg+nseg)
                        do kseg=1,nseg
                          if(kseg.ne.iseg) then
                                do i=1,3
        midx(i)=(x(i,3+kseg)+x(i,3+nseg+kseg))/2
        dirx(i)=x(i,3+nseg+kseg)-midx(i)
                                end do
                                ctype=secstr(kseg)
                                r=(x(1,3+kseg)+x(1,3+nseg+kseg))/2.0
                                gx=fung(r)
                                r=(x(2,3+kseg)+x(2,3+nseg+kseg))/2.0
                                gy=fung(r)
                                r=(x(3,3+kseg)+x(3,3+nseg+kseg))/2.0
                                gz=fung(r)
                                if(lgrid(gx,gy,gz,20)) then
c output text to STDOUT
c        write(*,530) cd,iseg,jseg,kseg,
c     $  atype,ctype,gx,gy,gz,(midx(i),i=1,3),(dirx(i),i=1,3)
c output binary to fort.1
        write(1) cd,iseg,jseg,kseg,
     $  atype,ctype,gx,gy,gz,(midx(i),i=1,3),(dirx(i),i=1,3)
                                end if

                          end if
                        end do
                      end if
                    end if
                  end do
                 end do
                else
                        write(*,*) '# too few SSEs: ',cd,nseg
                end if
        end do

520	format(a80)
530     format('#boxing ',a5,3i5,2(1x,a1),3i5,6f10.1)

	end
c
c----------------------------------------------------------------------
c
	subroutine preparex(x,iseg,jseg,nx,midpoint,direction)
	implicit none
	include 'parsizes.for'
	integer nx,iseg,jseg
	real x(3,*),midpoint(3,*),direction(3,*)
c
	integer i,kseg
c
	do i=1,3
		x(i,1)=midpoint(i,iseg)
		x(i,2)=midpoint(i,iseg)+direction(i,iseg)
		x(i,3)=midpoint(i,jseg)
	end do
	do kseg=1,nx
		do i=1,3
		  x(i,3+kseg)=midpoint(i,kseg)-direction(i,kseg)
		  x(i,3+nx+kseg)=midpoint(i,kseg)+direction(i,kseg)
		end do
	end do

	return
	end
c
c----------------------------------------------------------------------
c

c
c----------------------------------------------------------------------
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

c
c----------------------------------------------------------------------
c
	subroutine vec(ca,nres,nseg,segmentrange,midpoint,direction)
	implicit none
	include 'parsizes.for'
	integer nseg,segmentrange(2,maxseg),nres
	real ca(3,maxres),midpoint(3,maxseg+3),direction(3,maxseg+3)
c
	integer iseg,ires,left,rite,mid,i,l,r
	real nmid(3),cmid(3)
c
	do iseg=1,nseg
		left=segmentrange(1,iseg)
		rite=segmentrange(2,iseg)
		mid=(left+rite)/2
		l=mid-left+1
		r=rite-mid+1
		do i=1,3
			nmid(i)=0
			cmid(i)=0
		end do
		do ires=left,mid
			do i=1,3
				nmid(i)=nmid(i)+ca(i,ires)/l
			end do
		end do
		do ires=mid,rite
			do i=1,3
				cmid(i)=cmid(i)+ca(i,ires)/r
			end do
		end do
		do i=1,3
			midpoint(i,iseg)=(nmid(i)+cmid(i))/2
			direction(i,iseg)=cmid(i)-nmid(i)
		end do
	end do

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine readproteindata(nseg,secstr,segmentrange,ca,nres,iunit)
	implicit none
	include 'parsizes.for'
        integer ndom,node_child(2,maxdom),nres,nseg,iunit
        character node_type(maxdom),secstr(maxseg)
	integer domns(maxdom),domseglist(maxseg,maxdom),na,nb
	real ca(3,maxres)
	integer segmentrange(2,maxseg),checkrange(2,maxseg),checkx(maxseg)
c
	integer i,j,idom,iseg
c
	read(iunit,500) nres,nseg,na,nb,(secstr(i),i=1,nseg)
	do iseg=1,nseg
	  read(iunit,510) i,(segmentrange(j,i),j=1,2),(checkrange(j,i),j=1,2),
     $		checkx(i)
	end do
	read(iunit,520) ((ca(j,i),j=1,3),i=1,nres)
	read(iunit,500) ndom
	do idom=1,ndom
	  read(iunit,530) i,node_type(i),(node_child(j,i),j=1,2),domns(i),
     $		(domseglist(j,i),j=1,domns(i))
	end do

500	format(10x,4i5,2x,<maxseg>a1)
510	format(6i10)
520	format(10f8.1)
530	format(i4,1x,a1,1x,3i4,<maxseg>i4)

	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine setupprotein(cd,nseg,midpoint,direction,secstr,iunit,
     $		ca,nres,neidist,dalidatpath)
	implicit none
	include 'parsizes.for'
	character*5 cd
	integer nseg,iunit
	character secstr(maxseg)
	real midpoint(3,maxseg+3),direction(3,maxseg+3),neidist(maxseg,maxseg)
c
	character*80 filnam,constructfilnam,dalidatpath
	integer segmentrange(2,maxseg),nres,jseg,iseg
	real ca(3,maxres),distance,r
c
	filnam=constructfilnam(cd,dalidatpath,'.dat')
c	write(*,*),'setup: ',filnam
	open(iunit,file=filnam,status='old')
	call readproteindata(nseg,secstr,segmentrange,ca,nres,iunit)
	close(iunit)
c	write(*,*) cd,nseg,nres,' ',(secstr(i),i=1,nseg)
	call vec(ca,nres,nseg,segmentrange,midpoint,direction)
!	remember all distances in neidist
	do iseg=1,nseg
		do jseg=1,nseg
			if(iseg.ne.jseg) then
	r=distance(midpoint(1,iseg),midpoint(2,iseg),midpoint(3,iseg),
     $	midpoint(1,jseg),midpoint(2,jseg),midpoint(3,jseg))
				neidist(iseg,jseg)=r
			else
				neidist(iseg,jseg)=0.0
			end if
		end do
	end do

500	format(a4,a1)

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
10	read(iunit,500,end=19) cd
	if(nprot.eq.maxprot) then
		write(*,*) '# WARNING: skip reading list after maxprot',maxprot
		goto 19
	end if
	nprot=nprot+1
	list(nprot)=cd
	goto 10
19	close(iunit)
c	write(*,*) nprot,' proteins in list from ',filnam
c
500	format(a5)
c

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine twist(x,looplen)
c
c	puts x(*,1) in origo, x(*,2) along y axis, x(*,3) in positive yz plane
c	output: transrotated coordinates x
c
	implicit none
	integer looplen
	real x(3,looplen)
c
	real pii
	parameter(pii=3.141592653589793)
c
	integer i,j
	real u,sinu,cosu,y(3)
c
c	set origo
c
	do i=looplen,1,-1
		do j=1,3
			x(j,i)=x(j,i)-x(j,1)
		end do
	end do
c
c	rotate around x
c	-- if y=0 do first 90 deg around z
c
	if(abs(x(2,2)).lt.1e-6) then
		do i=2,looplen
			y(1)=-x(2,i)
			y(2)=x(1,i)
			y(3)=x(3,i)
			do j=1,3
				x(j,i)=y(j)
			end do
		end do
	end if
c
c	-- if y=0 & x=0 do first 90 deg around x
c
	if(abs(x(2,2)).lt.1e-6.and.abs(x(1,2)).lt.1e-6) then
		do i=2,looplen
			y(2)=-x(3,i)
			y(3)=x(2,i)
			y(1)=x(1,i)
			do j=1,3
				x(j,i)=y(j)
			end do
		end do
	end if
	if(abs(x(2,2)).gt.1e-6) then
		u=atan(x(3,2)/x(2,2))
	else
		u=0.0
	end if
	sinu=sin(u)
	cosu=cos(u)
	do i=2,looplen
		y(2)=cosu*x(2,i)+sinu*x(3,i)
		y(3)=-sinu*x(2,i)+cosu*x(3,i)
		y(1)=x(1,i)
		do j=1,3
			x(j,i)=y(j)
		end do
	end do
c	if(x(2,2).eq.0.0) write(*,*) ' y=0 !!!',u,sinu,cosu
c
c	rotate around z
c
	if(abs(x(2,2)).gt.1e-6) then
		u=-atan(x(1,2)/x(2,2))
	else
		u=0.0
	end if
	if(x(2,2).lt.0.0) u=pii+u
	sinu=sin(u)
	cosu=cos(u)
	do i=2,looplen
		y(1)=cosu*x(1,i)+sinu*x(2,i)
		y(2)=-sinu*x(1,i)+cosu*x(2,i)
		y(3)=x(3,i)
		do j=1,3
			x(j,i)=y(j)
		end do
	end do
c
c	rotate around y
c	-- if z=0 do first 90 deg around y
c
	if(abs(x(3,3)).lt.1e-6) then
		do i=2,looplen
			y(1)=x(3,i)
			y(2)=x(2,i)
			y(3)=-x(1,i)
			do j=1,3
				x(j,i)=y(j)
			end do
		end do
	end if
	if(x(3,3).ne.0.0) then
		u=atan(x(1,3)/x(3,3))
	else
		u=0.0
	end if
	sinu=sin(u)
	cosu=cos(u)
	do i=3,looplen
		y(3)=cosu*x(3,i)+sinu*x(1,i)
		y(1)=-sinu*x(3,i)+cosu*x(1,i)
		y(2)=x(2,i)
		do j=1,3
			x(j,i)=y(j)
		end do
	end do
c
c	-- if z<0 then rotate 180 deg around y
c
	if(x(3,3).lt.0.0) then
		do i=2,looplen
			y(1)=-x(1,i)
			y(2)=x(2,i)
			y(3)=-x(3,i)
			do j=1,3
				x(j,i)=y(j)
			end do
		end do
	end if
c	if(x(3,3).lt.0) write(*,*) ' z<0 !!!',u,sinu,cosu

	return
	end
c
c---------------------------------------------------------------------------
c
	function lgrid(gx,gy,gz,maxgrid)
c
c       check that g* are inside -maxgrid..+maxgrid
c       .true. if ok, .false. if outside
c
        implicit none
        logical lgrid
        integer gx,gy,gz,maxgrid
c
        lgrid=.true.
        if(gx.gt.maxgrid) lgrid=.false.
        if(gx.lt.-maxgrid)lgrid=.false.
        if(gy.gt.maxgrid) lgrid=.false.
        if(gy.lt.-maxgrid)lgrid=.false.
        if(gz.gt.maxgrid) lgrid=.false.
        if(gz.lt.-maxgrid)lgrid=.false.

        return
        end
c
c------------------------------------------------------------------------------
c
        function fung(x)
        implicit none
        integer fung
        real x

        fung=nint(x/2)

        return
        end
c
c------------------------------------------------------------------------------
c
