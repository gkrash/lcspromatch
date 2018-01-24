c This module/program is part of DaliLite (c) L. Holm 1999
c
	program main
	call puu
	end

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

