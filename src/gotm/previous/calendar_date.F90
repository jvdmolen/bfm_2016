!$Id: calendar_date.F90,v 1.8 2005-06-27 13:44:07 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Convert true Julian day to calendar date
!
! !INTERFACE:
   subroutine calendar_date(julian,yyyy,mm,dd)
!
! !DESCRIPTION:
!  Converts a Julian day to a calendar date --- year, month and day.
!  Based on a similar routine in \emph{Numerical Recipes}.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer(8),intent(in)                          :: julian
!
! !OUTPUT PARAMETERS:
   integer,intent(out)                            :: yyyy,mm,dd
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!
! !LOCAL VARIABLES:
   integer(8), parameter        :: IGREG=2299161
   integer(8)                   :: ja,jb,jc,jd,je
   REALTYPE                     :: x
!
!-----------------------------------------------------------------------
!BOC
   if(julian .ge. IGREG ) then
      x = (real(julian-1867216)-0.25)/36524.25
      ja = julian+1+int(x)-int(0.25*x)
   else
      ja = julian
   end if

   jb = ja+1524
   jc = int(6680 + (real(jb-2439870)-122.1)/365.25)
   jd = 365*jc+int((0.25*real(jc)))
   je = int(real(jb-jd)/30.6001)

   dd = int(jb-jd-int(30.6001*real(je)))
   mm = int(je-1)
   if (mm .gt. 12) mm = mm-12
   yyyy = int(jc - 4715)
   if (mm .gt. 2) yyyy = yyyy-1
   if (yyyy .le. 0) yyyy = yyyy-1

   return
   end subroutine calendar_date
!EOC



