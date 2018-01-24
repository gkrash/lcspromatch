	program run_u3b

	implicit none

	integer maxatom
	parameter(maxatom=100000)
	integer n
	real w(maxatom),x(3,maxatom),y(3,maxatom)
	real u(3,3), t(3), rms
	integer i,j,mode, ier
	parameter(mode=1)

	! input = n, (x,y,w,i=1,n)
	read(*,*) n
        rms=0.0
        ier=-2
        if(n.ge.3) then
                do i=1,n
 	                read(*,*) (x(j,i),j=1,3),(y(j,i),j=1,3),w(i)
                end do
                call u3b(w,x,y,n,mode,rms,u,t,ier)
        end if 
        write(*,500) rms,u,t,ier

500     format('rms ',f10.1/'u ',9f10.3/'t ',3f10.3/'ier ',i10)

        end
	


