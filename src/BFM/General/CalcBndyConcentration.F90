#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: CalcBndyConentration.f90
!
! DESCRIPTION
!   !
    ! Calculation of Column values....
    ! This routine is meant to use in 3D-model to calculate the river concentration.

!
! !INTERFACE
  SUBROUTINE CalcBndyConcentration(area,side,numc,kmax,mode,depend, &
                           multi,rho,salt,temp, cc_old,cc_inout,error)

!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,LOGUNIT,NZERO,ZERO,DONE
  use mem, ONLY: ppO3h,ppO3c,ppP1c,ppP1s,ppR3c,ppP6c,iiP6,ppPcc,OCDepth,ppB1c,ppBac
  use string_functions, ONLY:getseq_number
  use mem_phaeo,only: LIMIT_BOUND
  use global_interface,ONLY:PhaeocystisCalc
  use mem_globalfun,   ONLY: insw_vector

!  
!
! !AUTHORS
!   P. Ruardij
!
! !REVISION_HISTORY
!
! COPYING
!   
!   Copyright (C) 2006 P. Ruardij and M. Vichi
!   (rua@nioz.nl, vichi@bo.ingv.it)
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation;
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!EOP
!-------------------------------------------------------------------------!
!BOC
!
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE
  character(len=*),intent(IN)            :: area
  character(len=*),intent(IN)            :: side
  integer,intent(IN)                     :: numc,kmax,mode,depend
  real(RLEN),intent(IN)                  :: multi
  real(RLEN),intent(IN)                  :: temp(0:kmax)
  real(RLEN),intent(IN)                  :: salt(0:kmax)
  real(RLEN),intent(IN)                  :: rho(0:kmax)
  real(RLEN),intent(IN)                  :: cc_old(0:kmax,1:numc)
  real(RLEN),intent(INOUT)               :: cc_inout(0:kmax,1:numc)
  integer,intent(OUT)                    :: error
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! local variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN)                             :: salt_corr(0:kmax)
  real(RLEN)                             :: temp_corr(0:kmax)
  real(RLEN)                             :: old_DIC(0:kmax)
  real(RLEN)                             :: p(0:kmax)
  real(RLEN)                             :: mDepth
  character(len=4),dimension(4)          :: bdy_side 
  integer                                :: m,n(1),j
  logical                                :: llNorthSea

  data bdy_side/'west','north','east','south'/ 

  error=0
  
  if ( mode.eq.0) return
  llNorthSea=.true.
  if (trim(area) .eq. "WaddenSea") llNorthSea=.false.
  if(mode.eq.depend.and.multi.gt.0.0) then
        cc_inout(0:kmax,mode)=cc_inout(0:kmax,mode) *multi
  elseif (mode.gt.0.and.depend>0.and.depend<=numc) then
        cc_inout(0:kmax,mode)=cc_inout(0:kmax,depend) &
            *cc_old(0:kmax,mode)/(NZERO+cc_old(0:kmax,depend))
        if (multi.gt.0.0) then
          cc_inout(0:kmax,mode)=max(cc_inout(0:kmax,mode),cc_inout(0:kmax,depend)*multi)
        endif
     return
  elseif (mode.eq.ppbac) then
        cc_inout(0:kmax,mode)=min(cc_inout(0:kmax,mode),cc_inout(0:kmax,ppB1c))
  elseif (mode.eq.ppO3h.or.mode.eq.ppO3c) then
       j=getseq_number(side,bdy_side,4,.true.)
       if (j.eq.3.and.llNorthSea ) then
          ! From Lesly Salt
          cc_inout(:,ppO3h)= salt(:)*28.1D+00+ 1391.0D+00
       else
       ! Calculated according 
       !Lee et al. GeoPhysical Research Letters, Vol. 33 L19605
       ! Equation is only valid for range salt range 31-37 and temp-range 0-20
       salt_corr=min(max(31.0D+00,salt(:)),37.0D+00)
       temp_corr=min(max(0.0D+00,temp(:)),20.0D+00)
       mDepth=OCDepth(1)
       p=OCDepth/mDepth
       cc_inout(:,ppO3h)=(DONE-p)*cc_inout(:,ppO3h)+ &
         p*(2305.0D+00+53.97D+00*(salt_corr(:)-35.0D+00) &
          +2.74D+00*(salt_corr(:)-35.0D+00)**2 - 1.16D+00*(temp_corr(:)-20.0D+00)- &
           0.04D+00*(temp_corr(:)-20.0D+00)**2 )
       endif
  elseif (mode.eq.ppR3c) then
    p=DONE
    where(cc_inout(:,ppP6c) >1.0D-10)p=cc_inout(:,ppP6c)/(NZERO+ cc_inout(:,ppPcc))
    p=insw_vector(1.99D+00-p)*cc_inout(:,ppP6c)
    call PhaeocystisCalc(LIMIT_BOUND,iiP6,cc_inout(:,ppR3c),p,ZERO)
  endif
  select case (mode)
#ifdef INCLUDE_PELCO2
    case (ppO3h)
       !See above
    case (ppO3c)
       old_DIC(1:kmax)=cc_inout(1:kmax,ppO3c)
       mDepth=OCDepth(1)
       p=OCDepth/mDepth
       call CalcCO2SatInField(kmax,numc,rho(1:kmax),temp(1:kmax), &
                                        salt(1:kmax),cc_inout(1:kmax,1:numc))
       cc_inout(1:kmax,ppO3c)=cc_inout(1:kmax,ppO3c)*(DONE-p)+p*old_DIC(1:kmax)
       temp_corr(1:kmax)=abs(old_DIC(1:kmax)/cc_inout(1:kmax,ppO3c)-1.0)
       n=maxloc(temp_corr);m=n(1)
       if (abs(temp_corr(m))> 0.10) error=1;
       if (error.eq.1) then
         write(LOGUNIT,*) &
          'Boundary for O3c reset to equilibrium value in at least in one layer'
         write(LOGUNIT,*) 'old:',old_DIC(m),'new:',cc_inout(m,ppO3c);
         forall (j=1:kmax,temp_corr(j)> 0.10) &
         cc_inout(j,ppO3c)= max(0.95* cc_inout(j,ppO3c), &
              min(1.05*cc_inout(j,ppO3c),0.5*(cc_inout(j,ppO3c)+old_DIC(j))))
       endif
#endif
  end select

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
