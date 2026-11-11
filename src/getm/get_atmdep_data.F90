#include "cppdefs.h"

#ifdef BFM_GOTM

!-----------------------------------------------------------------------
! Update atmospheric-deposition forcing for model step n.
!-----------------------------------------------------------------------
subroutine get_atmdep_data(n)

   use atmdep, only: use_atmdep, atmdep_flux
   use ncdf_atmdep, only: get_atmdep_data_ncdf
   implicit none

   integer, intent(in) :: n

   if (.not. use_atmdep) return
   call get_atmdep_data_ncdf(n,atmdep_flux)

end subroutine get_atmdep_data

#endif
