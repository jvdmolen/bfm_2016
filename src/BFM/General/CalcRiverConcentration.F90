!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: CalcRiverConentration.f90
!
! DESCRIPTION
!
! Calculation of Column values....
! This routine is meant to use in 3D-model to calculate river concentrations.

!
! !INTERFACE
  SUBROUTINE CalcRiverConcentration(mode,kmax,numc,ETW,cc_old_river,cc_river)
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,NZERO
  use mem, ONLY: ppO3c,ppO3h,ppN4n,ppN3n,ppR6c,ppR6p,ppR6n,ppN1p, &
                 ppB1c,ppBac,ppR9x
  use mem_Param,ONLY:p_qR6cQ9x
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
  integer,intent(IN)         :: mode
  integer,intent(IN)         :: kmax
  integer,intent(IN)         :: numc
  real(RLEN),intent(IN)      :: ETW(1:kmax)
  real(RLEN),intent(IN)      :: cc_old_river(1:kmax,1:numc)
  real(RLEN),intent(INOUT)   :: cc_river(1:kmax,1:numc)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! external functions
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

   real(RLEN)               :: ERHO(1:kmax)
   real(RLEN)               :: ESW(1:kmax)

#ifdef INCLUDE_PELCO2

   if ( mode.eq.0) return
   select case (mode)
     case ( ppO3h)
     case ( ppO3c)
       ERHO=1027.0;
       ESW=0.3;
       call CalcCO2SatInField(kmax,numc,ERHO,ETW,ESW,cc_river)
     case ( ppBac)
       cc_river(:,mode)=(cc_river(:,ppB1c)) &
                /(NZERO+cc_old_river(:,ppB1c))*cc_old_river(:,mode)
     case ( ppR6c)
       cc_river(:,mode)=(cc_river(:,ppN1p))*(1.0/0.70-1.0) /0.000391
!      cc_river(ppR9x,:)=cc_river(:,mode)/p_qR6cQ9x
     case ( ppR6n)
       cc_river(:,mode)=(cc_river(:,ppN3n)+cc_river(:,ppN4n))*(1.0/0.93-1.0)
     case ( ppR6p)
       cc_river(:,mode)=cc_river(:,ppN1p)*(1.0/0.70-1.0)
   end select
#endif

  end

!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
