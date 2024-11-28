#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: CO2
!
! DESCRIPTION
!   !   This is the top level of the benthic submodel.
!       All the biological processes affecting the benthic dynamics are called
!       in a specified sequence according to the calculation flags.
!
!

!
! !INTERFACE
  subroutine CO2Dynamics(control)
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,LOGUNIT,ZERO
  use mem, ONLY:ppO3c,ppO3h,ppG3c,ppG3h
  use mem_Param, ONLY: CalcBenthicFlag,p_dry_ben,p_clDxm


!
!
! !AUTHORS
!   Piet Ruardij ERSEM group
!       
!
!
! !REVISION_HISTORY
!   !
!
! COPYING
!
!   Copyright (C) 2006 P. Ruardij, the mfstep group, the ERSEM team
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

  integer,intent(IN)         ::control
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  integer            ::i,j

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! CO2 ppO3c > 0 : all routines and variables are present for CO2processes
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#ifdef INCLUDE_PELCO2

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! CO2 ppO3c > 0 : all routines and variables are present for CO2processes
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  select case (control)
    case(1)
      ! This is called before the call in PelDynamics in the routine Ecology 
      ! Aim: calculation of the partial CO2's and pH
      if ( ppO3c > 0 )  call PelCO2Dynamics()
      ! Calculation of the Ph in the 3 benthic layers.
#ifdef INCLUDE_BENCO2
      if ( CalcBenthicFlag==3.and. ppG3c>0 .and. ppG3h> 0 ) call BenpHDynamics
#endif
    case(2)
      !This called after the last routine where biochemical rates are defined.
      !in pelagic
      if ( ppO3c > 0 )  call AlkalinityDynamics(1)
    case(3)
      ! ppG3c > 0 : BenCO2 dynamics is active....
#ifdef INCLUDE_BENCO2
      if ( CalcBenthicFlag==3 ) then
        if  (ppG3c > 0 ) call BenCO2TransportDynamics
        if ( ppG3h > 0 ) call BenAlkalinityDynamics

        call BenProfiles
        if ( ppG3c > 0 .and. ppG3h> 0 ) call BenCO2Profiles
      endif
#endif
    end select
#endif


  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
