c This module/program is part of DaliLite (c) L. Holm 1999
c
	program wolf
c
c	f77 wolf.f subfitz.f u3b-8.f
c
c	use as:
c		cp list.sep list2
c		echo "codeX" > list1
c		wolf-one | sort -nr +1 > codeX.wolf
c
	implicit none
	include 'parsizes.for'
	character secstr(maxseg),secstrx(maxseg)
	integer nseg
	real midy(3,maxseg+3),diry(3,maxseg+3)
c
	integer nprot1,nprot2,ns
	character*5 cd,list1(maxprot),list2(maxprot),cd1
c
	integer*4 box(-20:20,-20:20,-20:20),link_a(10000),link_c(10000),
     $		link_b(10000)
	real link_from(3,10000),link_to(3,10000)
	integer lastlink,link_next(10000)
	integer i,j,protcount,bestpair(4)
	logical lverb
c
	character*80 dalidatpath_1,dalidatpath_2
c
	real xca(3,maxres),yca(3,maxres),x(3,maxres),y(3,maxres),rms
	integer nx,ny,ali(maxres),k,lali,niter,nsegx
	real midx(3,maxseg+3),dirx(3,maxseg+3),neidist(maxseg,maxseg)
	integer nocc,occ(3,60000), fitzmaxiter
	real neiborcutoff, fitzrcut
c
	write(*,*),'enter dalidatpath_1, _2, fitzrcut, fitzmaxiter, neiborcutoff'
	read(*,520) dalidatpath_1
	read(*,520) dalidatpath_2
	read(*,*) fitzrcut
	read(*,*) fitzmaxiter
	read(*,*) neiborcutoff
	write(*,*),'imported parameters are:',dalidatpath_1,
     $		dalidatpath_2,fitzrcut,fitzmaxiter,
     $		neiborcutoff
c
	lverb=.false.
	nprot1=0
	nprot2=0
	call getlist('list1',90,list1,nprot1)
	call getlist('list2',90,list2,nprot2)
c
	call initgrid(box,20,20,20)
c
c	search-loop
c
	nocc=0
	cd=list1(1)
	cd1=cd
	if(cd1(5:5).eq.' ') cd1(5:5)='_'
	call setupprotein(cd,nsegx,midx,dirx,secstrx,90,xca,nx,neidist,
     $		dalidatpath_1)
	call loadprotein(box,20,20,20,nsegx,midx,dirx,ns,neiborcutoff,
     $		link_a,link_b,
     $		link_c,link_from,link_to,lastlink,link_next,neidist,nocc,occ)
	if (nsegx.gt.2) then
	  do i=1,nprot2
		cd=list2(i)
		call setupprotein(cd,nseg,midy,diry,secstr,90,
     $			yca,ny,neidist,dalidatpath_2)
		if(nseg.gt.2) then
		  if(lverb) write(*,*),ns,' AB pairs boxed',lastlink
		  call compare(nseg,midy,diry,secstr,cd,nseg,secstrx,
     $			link_next,link_from,link_to,link_a,link_b,link_c,
     $			box,20,20,20,neiborcutoff,protcount,bestpair,neidist)
		  write(*,500),cd,protcount,(bestpair(j),j=1,4)
c
c		  fitz C(alphas) in best-per-protein trial superimpositions
c
		  if(bestpair(1).eq.0.or.bestpair(2).eq.0) goto 39
		  if(lverb) write(*,*),'fitzing ',cd,i

		  call preparex(x,bestpair(3),bestpair(4),0,midx,dirx)
		  do j=1,nx
			do k=1,3
				x(k,j+3)=xca(k,j)
			end do
		  end do
		  call twist(x,3+nx)	! prepare sets twist frame, transform CA
		  do j=1,nx
			do k=1,3
				x(k,j)=x(k,j+3)
			end do
		  end do
		  call preparex(y,bestpair(1),bestpair(2),0,midy,diry)
		  do j=1,ny
			do k=1,3
				y(k,j+3)=yca(k,j)
			end do
		   end do
		  call twist(y,3+ny)	! prepare sets twist frame, transform CA
		  do j=1,ny
			do k=1,3
				y(k,j)=y(k,j+3)
			end do
		  end do

		  call fitz(x,nx,y,ny,ali,fitzrcut,fitzmaxiter,rms,lali,niter)
		  write(*,510),cd1,cd,protcount,(bestpair(j),j=3,4),
     $		    (bestpair(j),j=1,2),rms,lali,niter,nx,(ali(j),j=1,nx)

		end if
39	  end do
	end if

500	format(a5,i20,4i10)
510	format('WOLFITZ ',2a5,5i5,f10.1,3i5,<maxres>i4)
520	format(a80)

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
	subroutine compare(nseg,midpoint,direction,secstr,protein_code,
     $		protein_nseg,protein_secstr,link_next,link_from,link_to,
     $		link_a,link_b,link_c,box,boxdim1,boxdim2,boxdim3,
     $		rcut,protcount,bestpair,neidist)
	implicit none
	include 'parsizes.for'
	integer nseg,boxdim1,boxdim2,boxdim3,bestpair(4)
	real midpoint(3,maxseg),direction(3,maxseg),rcut
	integer*4 box(-boxdim1:boxdim1,-boxdim2:boxdim2,-boxdim3:boxdim3),
     $		link_a(10000),link_b(10000),link_c(10000)
 	real link_from(3,10000),link_to(3,10000),neidist(maxseg,maxseg)
	integer link_next(10000),protcount
	character secstr(maxseg)
c
	integer protein_nseg,count(maxseg,maxseg)
	character*5 protein_code
	character protein_secstr(maxseg)
c
	integer cur_a,cur_c,cur_link,iseg,jseg,kseg,gx,gy,gz,fung,i,j
	real mid(3),r,distance,x(3,maxseg+maxseg+3),cosi,dir(3),midx(3),dirx(3)
	character atype,ctype
	logical lgrid,less,lest
	integer dgx,dgy,dgz,gx0,gy0,gz0,cur_b,useg,vseg,lseg
c
	dgx=2
	dgy=2
	dgz=2
c	write(*,*) 'this is compare',nseg,rcut
	protcount=0
	do j=1,4
		bestpair(j)=0
	end do
	do iseg=1,nseg
	  atype=secstr(iseg)
	  do jseg=1,nseg
	    lest=(iseg.lt.jseg)
	    if(iseg.ne.jseg) then
	      r=neidist(iseg,jseg)
	      if(r.lt.rcut) then
		call initcomparison(protein_nseg,count)
		call preparex(x,iseg,jseg,nseg,midpoint,direction)
		call twist(x,3+nseg+nseg)
		! require midpoint distance < 3 A and cosi(c,c')>0.5
		do kseg=1,nseg
		  if(kseg.ne.iseg) then
			less=(iseg.lt.kseg)
			do i=1,3
				midx(i)=(x(i,3+kseg)+x(i,3+nseg+kseg))/2
				dirx(i)=x(i,3+nseg+kseg)-midx(i)
			end do
			ctype=secstr(kseg)
			r=(x(1,3+kseg)+x(1,3+nseg+kseg))/2.0
			gx0=fung(r)
			r=(x(2,3+kseg)+x(2,3+nseg+kseg))/2.0
			gy0=fung(r)
			r=(x(3,3+kseg)+x(3,3+nseg+kseg))/2.0
			gz0=fung(r)
			do gx=gx0-dgx,gx0+dgx
			  do gy=gy0-dgy,gy0+dgy
			    do gz=gz0-dgz,gz0+dgz
			      if(lgrid(gx,gy,gz,20)) then
				cur_link=box(gx,gy,gz)
c				write(*,*) cur_link,gx,gy,gz
				do while(cur_link.ne.0)
					cur_a=link_a(cur_link)
					cur_b=link_b(cur_link)
					cur_c=link_c(cur_link)
c	write(*,*) 'test:',iseg,jseg,kseg,atype,ctype,cur_link,
c     $		link_a(cur_link),link_c(cur_link),
c     $		protein_secstr(cur_a),protein_secstr(cur_c)
	if(protein_secstr(cur_a).ne.atype) goto 19
	if(protein_secstr(cur_c).ne.ctype) goto 19
	if(less) then	! simple topology filters
		if(cur_a.ge.cur_c) goto 19
	else
		if(cur_a.le.cur_c) goto 19
	end if
	if(lest) then
		if(cur_a.ge.cur_b) goto 19
	else
		if(cur_a.le.cur_b) goto 19
	end if
	do i=1,3
		mid(i)=(link_from(i,cur_link)+link_to(i,cur_link))/2
		dir(i)=link_to(i,cur_link)-mid(i)
	end do
	if(distance(mid(1),mid(2),mid(3),midx(1),midx(2),midx(3))
     $		.gt.4.0) goto 19
	if(cosi(dir(1),dir(2),dir(3),dirx(1),dirx(2),dirx(3)).lt.0.5) goto 19
				  count(cur_a,cur_b)=count(cur_a,cur_b)+1
19				  cur_link=link_next(cur_link)
				end do
			      end if
			    end do
			  end do
			end do
		  end if
		end do
c
c		report counts
c
		i=0
		useg=0
		vseg=0
		do lseg=1,protein_nseg
		  do kseg=1,protein_nseg
			if(count(lseg,kseg).gt.i) then
				i=count(lseg,kseg)
				useg=lseg
				vseg=kseg
			end if
		  end do
		end do
c		write(*,*) 'counts: ',protein_code,iseg,jseg,useg,vseg,i
		if(i.gt.protcount) then
			protcount=i
			bestpair(1)=iseg
			bestpair(2)=jseg
			bestpair(3)=useg
			bestpair(4)=vseg
		end if
	      end if
	    end if
	  end do
	end do

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine loadprotein(box,boxdim1,boxdim2,boxdim3,
     $		nseg,midpoint,direction,ns,rcut,
     $		link_a,link_b,link_c,link_from,link_to,
     $		lastlink,link_next,neidist,nocc,occ)
	implicit none
	include 'parsizes.for'
	integer nocc,occ(3,60000)
	integer boxdim1,boxdim2,boxdim3
	integer*4 box(-boxdim1:boxdim1,-boxdim2:boxdim2,-boxdim3:boxdim3)
	integer nseg,ns
	real midpoint(3,maxseg),direction(3,maxseg),rcut
	integer*4 link_a(10000),link_c(10000)
	integer*4 link_b(10000)
	real link_from(3,10000),link_to(3,10000),neidist(maxseg,maxseg)
	integer lastlink,link_next(10000)
c
	integer iseg,jseg,i
	real r,x(3,3+maxseg+maxseg)
c
c	clean up box
c
c	write(*,*),'unload',nocc
	do i=1,nocc
		box(occ(1,i),occ(2,i),occ(3,i))=0
	end do
	nocc=0
	lastlink=0
c
c
c
	ns=0
	do iseg=1,nseg
	  do jseg=1,nseg
	    if(jseg.ne.iseg.and.neidist(iseg,jseg).lt.rcut) then
	      r=neidist(iseg,jseg)
	      if(r.lt.rcut) then
		call preparex(x,iseg,jseg,nseg,midpoint,direction)
		call twist(x,3+nseg+nseg)
		call boxit(nseg,x,iseg,jseg,lastlink,link_a,link_b,
     $			link_c,link_from,link_to,link_next,
     $			box,boxdim1,boxdim2,boxdim3,nocc,occ)
		ns=ns+1
	      end if
	    end if
	  end do
	end do

500	format('before',12f5.1)
505	format('after ',12f5.1)
510	format('line',3f10.1,' to ',3f10.1,';')

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine boxit(nseg,x,aseg,bseg,lastlink,link_a,
     $		link_b,link_c,link_from,link_to,
     $		link_next,box,boxdim1,boxdim2,boxdim3,nocc,occ)
	implicit none
	include 'parsizes.for'
	integer nseg,aseg,bseg
	real x(3,maxseg+maxseg+3)
	integer boxdim1,boxdim2,boxdim3
	integer*4 box(-boxdim1:boxdim1,-boxdim2:boxdim2,-boxdim3:boxdim3)
	integer*4 link_a(10000),link_c(10000)
	integer*4 link_b(10000)
	real link_from(3,10000),link_to(3,10000)
	integer lastlink,link_next(10000),nocc,occ(3,60000)
c
	integer fung
	logical lgrid
c
	integer gx,gy,gz,iseg,i,nx
	real mid(3)
c
	nx=0
	do iseg=1,nseg
	  if(iseg.ne.aseg.and.lastlink.lt.10000) then
		lastlink=lastlink+1
		if(lastlink.gt.10000) stop 'link overflow'
		link_a(lastlink)=aseg
		link_b(lastlink)=bseg
		link_c(lastlink)=iseg
		do i=1,3
			link_from(i,lastlink)=x(i,3+iseg)
			link_to(i,lastlink)=x(i,3+nseg+iseg)
			mid(i)=(x(i,3+iseg)+x(i,3+nseg+iseg))/2.0
		end do
		gx=fung(mid(1))
		gy=fung(mid(2))
		gz=fung(mid(3))
		if(lgrid(gx,gy,gz,20)) then
			link_next(lastlink)=box(gx,gy,gz)
			if(box(gx,gy,gz).eq.0) then
				nocc=nocc+1
				occ(1,nocc)=gx
				occ(2,nocc)=gy
				occ(3,nocc)=gz
			end if
			box(gx,gy,gz)=lastlink
c			write(*,*) 'boxed',lastlink,aseg,iseg,gx,gy,gz
		else
			lastlink=lastlink-1
			nx=nx+1
		end if
	  end if
	end do
c	if(nx.gt.0) write(*,*) nx,' elements ignored by boxit'

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine initcomparison(protein_nseg,count)
	implicit none
	include 'parsizes.for'
	integer protein_nseg,count(maxseg,maxseg)
c
	integer j,k
c
	do j=1,protein_nseg
		do k=1,protein_nseg
			count(j,k)=0
		end do
	end do

	return
	end
c
c----------------------------------------------------------------------
c
	function cosi(z1,z2,z3,x1,x2,x3)
	implicit none
	real cosi,z1,z2,z3,x1,x2,x3
c
	real lz,lx,norm
c
	lz=sqrt(z1*z1+z2*z2+z3*z3)
	lx=sqrt(x1*x1+x2*x2+x3*x3)
	norm=max(1e-6,lz*lx)
	cosi=(z1*x1+z2*x2+z3*x3)/norm
c
	return
	end
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
	subroutine initgrid(box,boxdim1,boxdim2,boxdim3)
	implicit none
	integer boxdim1,boxdim2,boxdim3
	integer*4 box(-boxdim1:boxdim1,-boxdim2:boxdim2,-boxdim3:boxdim3)
c
	integer i,j,k
c
	do k=-boxdim3,boxdim3
		do j=-boxdim2,boxdim2
			do i=-boxdim1,boxdim1
				box(i,j,k)=0
			end do
		end do
	end do

	return
	end
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
	write(*,*),'setup: ',filnam
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
c		write(*,*) 'WARNING: skip reading list after maxprot',maxprot
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
