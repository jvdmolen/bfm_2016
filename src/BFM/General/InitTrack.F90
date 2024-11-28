#include "DEBUG.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: InitTrack
!
! DESCRIPTION
!   Initialization of tracking
!   Allocation of memory for variables, reading of data files 

! !INTERFACE
  SUBROUTINE InitTrack
!
! USES:
  use track, only:calc_all_track_rates

  call calc_all_track_rates(0)

end
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: fill_sfl_if_tracking
!
! DESCRIPTION
!   Initialization of tracking
!   Allocation of memory for variables, reading of data files 

! !INTERFACE
  SUBROUTINE fill_sfl( parallel,nr,numc,sfl)

!
! USES:
  use constants, only: RLEN
  use track, only:fill_sfl_if_tracking

   logical                    :: parallel
   integer,intent(IN)         :: nr
   integer,intent(IN)         :: numc
   real(RLEN),intent(INOUT)   :: sfl(0:numc)

   call fill_sfl_if_tracking( parallel,nr,numc,sfl)


end
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model! 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
