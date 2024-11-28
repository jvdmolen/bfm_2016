#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: cnv_3d() - convert 3D-scalar fields to real*4.
!
! !INTERFACE:
   subroutine cnv_3d(imin,jmin,imax,jmax,kmin,kmax,mask,var,missing, &
                     il,ih,jl,jh,kl,kh,ws)
   IMPLICIT NONE
!
! !DESCRIPTION:
! 
! !BFM
!   code added to overcome warning from netcdf :
!   "NetCDF: Numeric conversion not representable"
!   This happen in case that double real (REAL*8) value which 
!   are saved as "float" (REAL*4) in a netcdf-results file are
!   larger than 1.0E125 or small then 1.0E-125 
!
! !USES:
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: imin,jmin,imax,jmax
   integer, intent(in)                 :: kmin(I2DFIELD)
   integer, intent(in)                 :: kmax
   integer, intent(in)                 :: mask(E2DFIELD)
   REALTYPE, intent(in)                :: var(I2DFIELD,kl:kh)
   REALTYPE, intent(in)                :: missing
   integer, intent(in)                 :: il,jl,ih,jh,kl,kh
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: ws(I2DFIELD,kl:kh)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k
   REAL(kind=4)              :: short
   REALTYPE                  :: r
   REALTYPE,parameter        :: max_value=2.0D+00**(maxexponent(short)-1)
   REALTYPE,parameter        :: min_value=2.0D+00**(minexponent(short)+1)
!EOP
!-----------------------------------------------------------------------
!BOC
   do k=kl,kh
      do j=jl,jh
         do i=il,ih
            if (mask(i,j) .gt. 0 .and. k .ge. kmin(i,j) ) then
               r=var(i,j,k)
               if (abs(r).lt.min_value) r=_ZERO_
               if (abs(r).gt.max_value) r=sign(max_value,r)
               ws(i,j,k) =r
            else
               ws(i,j,k) = missing
            end if
         end do
      end do
   end do
   return
   end subroutine cnv_3d
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
