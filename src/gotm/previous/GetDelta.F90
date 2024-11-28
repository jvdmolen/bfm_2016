#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the numeric timestep from gotm
!
! !INTERFACE:
   function GetDelta() 
!
! !DESCRIPTION:
!  Transfer the integration time step from GOTM to the BFM
!  Unit conversion from seconds to days
!
! !USES:
   use constants, only: RLEN
   use bio_var, only: dt,secs_pr_day
   use time, only: timestep
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi 
!
!EOP
!-----------------------------------------------------------------------
!BOC
   real(RLEN) :: GetDelta

   if ( dt< 0.0) then
      GetDelta = timestep/secs_pr_day
   else
      GetDelta = dt/secs_pr_day
   endif
!  STDERR copy_of_dt,timestep,GetDelta

   return
   end function GetDelta
!EOC
!-----------------------------------------------------------------------


