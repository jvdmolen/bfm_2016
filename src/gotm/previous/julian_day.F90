!$Id: time.F90,v 1.8 2005-06-27 13:44:07 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Convert a calendar date to true Julian day
!
! !INTERFACE:
   subroutine julian_day(yyyy,mm,dd,julian)
!
! !DESCRIPTION:
!  Converts a calendar date to a Julian day.
!  Based on a similar routine in \emph{Numerical Recipes}.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer                             :: yyyy,mm,dd
!
! !OUTPUT PARAMETERS:
   integer(8)                          :: julian
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!
! !LOCAL VARIABLES:
   integer(8), PARAMETER        :: IGREG=15+31*(10+12*1582)
   integer(8)                   :: ja,jy,jm
!
!-----------------------------------------------------------------------
!BOC
   jy = yyyy
   if(jy .lt. 0) jy = jy+1
   if (mm .gt. 2) then
      jm = mm+1
   else
      jy = jy-1
      jm = mm+13
   end if
   julian = int(floor(365.25*real(jy))+floor(30.6001*real(jm))+real(dd+1720995))
   if (dd+31*(mm+12*yyyy) .ge. IGREG) then
      ja = int(0.01*real(jy))
      julian = julian+2-ja+int(0.25*real(ja))
   end if

   return
   end subroutine julian_day
!EOC


