c This module/program is part of DaliLite (c) L. Holm 1999,2008
c
	program wolf
c
c	f77 wolf.f 
c
c	use as: wolf_grid | wolf
c
c       reads WOLFDB
c
	implicit none
	include 'parsizes.for'
        integer maxlines
        parameter (maxlines=1000000)
c       largest protein: 1y1vA 1426 residues, 69 SSEs, 4041 frame-lines
	integer*2 box(-20:20,-20:20,-20:20),link_list(maxlines)
c       parsizes.for has maxseg=200
        integer q_aseg(maxlines),q_bseg(maxlines),q_cseg(maxlines)
        character*1 q_atype(maxlines),q_ctype(maxlines)
        real q_mid(3,maxlines),q_dir(3,maxlines)
        character*5 query_cd
c
        character*5 cd,old_cd
        real mid(3),dir(3)
        character*1 atype,ctype
        integer aseg,bseg,cseg,nseg
c
        integer i,n,old_aseg,old_bseg
        integer gx,gy,gz
        real scores(maxseg,maxseg),best_seen
        logical less
        integer dx,dy,dz,gx0,gy0,gz0
        real score,cosi,r2,distance_squared,cosine
c
        real bestscore,bestperframe_score(maxseg,maxseg)
        integer best_iseg,best_jseg
        integer bestperframe_ijseg(2,maxseg,maxseg)
	character*80 filename
c
c	binary database
c
	write(*,*) 'enter WOLF90 database'
	read(*,500) filename
500	format(a80)

c       initialise box
c
        call initgrid(box,20,20,20)
c
c       read query from fort.1 (binary)
c
        n=1
10      read(1, end=19) query_cd,q_aseg(n),q_bseg(n),q_cseg(n),
     $          q_atype(n),q_ctype(n),gx,gy,gz,
     $          (q_mid(i,n),i=1,3),(q_dir(i,n),i=1,3)

c       load linked list of line-numbers in box
        i=box(gx,gy,gz)
        link_list(n)=i 
        box(gx,gy,gz)=n
        n=n+1
        goto 10 ! query input loop
19      continue
c
c       read targets stream from WOLFDB (binary)
c
        old_cd='?????'
        old_aseg=0
        old_bseg=0
        nseg=0
        open(90,file=filename,form='unformatted',
     $          status='old')
20      read(90,end=29) cd,aseg,bseg,cseg,atype,ctype,
     $          gx0,gy0,gz0,(mid(i),i=1,3),(dir(i),i=1,3)

	if(aseg.gt.nseg) nseg=aseg
        if(bseg.gt.nseg) nseg=bseg
        if(cseg.gt.nseg) nseg=cseg
c       new target        
        if(cd.ne.old_cd) then
                if(old_cd.ne.'?????') 
     $		call output(query_cd,old_cd,nseg,
     $          bestperframe_ijseg,
     $          bestperframe_score)
                call init_scores(bestperframe_score,nseg)
                call init_counts(bestperframe_ijseg,nseg)
                call init_scores(scores,nseg)
                old_aseg=0
                old_bseg=0
        end if
        old_cd=cd
c       new target.frame
        if(aseg.ne.old_aseg.and.bseg.ne.old_bseg) then
          if(old_aseg.gt.0.and.old_bseg.gt.0) then
            call save_frame(scores,nseg,best_iseg,best_jseg,
     $          bestscore) 
            bestperframe_ijseg(1,old_aseg,old_bseg)=best_iseg
            bestperframe_ijseg(2,old_aseg,old_bseg)=best_jseg
            bestperframe_score(old_aseg,old_bseg)=bestscore
          end if
        end if
        old_aseg=aseg
        old_bseg=bseg
c       cumulate comparison scores per query.frame
c       find best one-to-one segment match of target.segment
        less=aseg.lt.bseg
        best_seen=0.0
        do dx=-2,2
                gx=gx0+dx
		if(gx.lt.-20.or.gx.gt.20) goto 31
                do dy=-2,2
                        gy=gy0+dy
			if(gy.lt.-20.or.gy.gt.20) goto 32
                        do dz=-2,2
                              gz=gz0+dz
	if(gz.lt.-20.or.gz.gt.20) goto 33
                              i=box(gx,gy,gz)
                      do while(i.gt.0)
        if(atype.ne.q_atype(i)) goto 39
        if(ctype.ne.q_ctype(i)) goto 39
        if(less.and.q_aseg(i).ge.q_cseg(i)) goto 39
        if(.not.less.and.q_aseg(i).le.q_cseg(i)) goto 39
        r2=distance_squared(mid(1),mid(2),mid(3),
     $          q_mid(1,i),q_mid(2,i),q_mid(3,i))
        if(r2.gt.16.0) goto 39
        cosi=cosine(dir(1),dir(2),dir(3),
     $          q_dir(1,i),q_dir(2,i),q_dir(3,i))
        if(cosi.lt.0.5) goto 39
        score=1.0-r2/16.0+(cosi-0.5)*2
        if(score.gt.best_seen) best_seen=score
39                                      i=link_list(i)
                              end do
33                      end do
32              end do
31      end do
c       only allow one-to-one segment matching of target segment
        scores(aseg,bseg)=scores(aseg,bseg)+best_seen
        goto 20 ! WOLFDB input loop
29      close(90)
        if(old_aseg.ne.0.and.old_bseg.ne.0) then
          if(old_aseg.gt.0.and.old_bseg.gt.0) then
            call save_frame(scores,nseg,best_iseg,best_jseg,
     $          bestscore) 
            bestperframe_ijseg(1,old_aseg,old_bseg)=best_iseg
            bestperframe_ijseg(2,old_aseg,old_bseg)=best_jseg
            bestperframe_score(old_aseg,old_bseg)=bestscore
          end if
        end if
        if(old_cd.ne.'?????') call output(query_cd,old_cd,nseg,
     $          bestperframe_ijseg,
     $          bestperframe_score)

550     format(8x,a5,3i5,2(1x,a1),3i5,6f10.1)

        end
c
c----------------------------------------------------------------------
c
        subroutine save_frame(scores,nseg,best_iseg,best_jseg,bestscore)
        implicit none
        include 'parsizes.for'
        real bestscore
        integer best_iseg,best_jseg
        real scores(maxseg,maxseg)
        integer nseg
c
        integer iseg,jseg

        bestscore=0.0
        best_iseg=0
        best_jseg=0
        do iseg=1,nseg
                do jseg=1,nseg
                        if(scores(iseg,jseg).gt.bestscore) then
                               bestscore=scores(iseg,jseg)
                               best_iseg=iseg
                               best_jseg=jseg
                        end if
                end do
        end do

        return
        end
c
c----------------------------------------------------------------------
c
        subroutine output(query_cd,target_cd,nseg,bestperframe_ijseg,
     $          bestperframe_score)
        implicit none
        include 'parsizes.for'
        character*5 query_cd,target_cd
        integer nseg,iseg,jseg
        integer bestperframe_ijseg(2,maxseg,maxseg)
        real bestperframe_score(maxseg,maxseg)
c
        real x,protcount
        integer lseg,kseg,useg,vseg,bestpair(4)
c       
c       report best match per target_cd
c       
	x=0.0
	useg=0
	vseg=0
        iseg=0
        jseg=0
	do lseg=1,nseg
	  do kseg=1,nseg
		if(bestperframe_score(lseg,kseg).gt.x) then
			x=bestperframe_score(lseg,kseg)
			useg=lseg
			vseg=kseg
                        iseg=bestperframe_ijseg(1,lseg,kseg)
                        jseg=bestperframe_ijseg(2,lseg,kseg)
		end if
	  end do
	end do
        write(*,500) query_cd,iseg,jseg,
     $          target_cd,useg,vseg,x

500     format(a5,2i5,1x,a5,2i5,f10.1)

        return
        end
c
c----------------------------------------------------------------------
c
	function cosine(z1,z2,z3,x1,x2,x3)
	implicit none
	real cosine,z1,z2,z3,x1,x2,x3
c
	real lz,lx,norm
c
	lz=sqrt(z1*z1+z2*z2+z3*z3)
	lx=sqrt(x1*x1+x2*x2+x3*x3)
	norm=max(1e-6,lz*lx)
	cosine=(z1*x1+z2*x2+z3*x3)/norm
c
	return
	end
c
c----------------------------------------------------------------------
c
	function distance_squared(a1,a2,a3,b1,b2,b3)
	implicit none
	real distance_squared,a1,a2,a3,b1,b2,b3
c
	distance_squared=(a1-b1)*(a1-b1)+(a2-b2)*(a2-b2)+(a3-b3)*(a3-b3)
c
	return
	end
c
c----------------------------------------------------------------------
c
	subroutine initgrid(box,boxdim1,boxdim2,boxdim3)
	implicit none
	integer boxdim1,boxdim2,boxdim3
	integer*2 box(-boxdim1:boxdim1,-boxdim2:boxdim2,-boxdim3:boxdim3)
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
c------------------------------------------------------------------------------
c
        subroutine init_scores(scores,nseg)
        implicit none
        include 'parsizes.for'
        integer nseg
        real scores(maxseg,maxseg)
c
        integer i,j

        do i=1,nseg
                do j=1,nseg
                        scores(i,j)=0.0
                end do
        end do

        return
        end
c
c------------------------------------------------------------------------------
c
        subroutine init_counts(scores,nseg)
        implicit none
        include 'parsizes.for'
        integer nseg
        integer scores(2,maxseg,maxseg)
c
        integer i,j,k

        do i=1,nseg
                do j=1,nseg
                        do k=1,2
                                scores(k,i,j)=0
                        end do
                end do
        end do

        return
        end
c
c------------------------------------------------------------------------------
c

