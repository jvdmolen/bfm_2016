#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: SurfaceCO2Processes
!
! DESCRIPTION
!   Model describes reaeration between air and water column.
!       as forced by temperature, wind and chemical enhancement.
!	second routine to be applied according to Schneider et al., 1999
!
!       The equation and correlation used in this routine
!       are found in the
!
!		R. Wanninkhof (1992), Relationship between windspeed and gas
!		exchange over the oecean
!               J. GeoPhys. Res. 97, 7373-7382
!
!	notes: K0 = co2/pco2 => dimension = 1.e-6mol/(l*1.e-6atm)
!	exchange coefficient: deltapCO2 * k660 * K0
!	=> 1.e-6atm * cm/hr * 1.e-6mol/(l * 1.e-6atm) = cm/hr * 1.e-6mol / l
!	Temp in degrees C
!	test parameter see CalculateCO2system.f (DIC=2133,
!	AC=2260, pco2=341), O7.c = AC-2210,
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine SurfaceCO2Processes(BoxNumber,BoxNumberXY)
!
! !USES:

  ! For the following Pelagic-states fluxes are defined: O3c
  ! The following global scalar vars are used: Wind
  ! The following Pelagic 1-d global boxvars are used: ETW, ERHO, pCO2, &
  ! Depth
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,LOGUNIT
#ifdef INCLUDE_PELCO2
  use mem,  ONLY: D2STATE
  use mem, ONLY: ppO3c, Wind, ETW, ESW, ERHO, pCO2, jsurO3c, &
    Depth, iiPel, flux
  use constants,only:MW_C
! use gotm_error_msg,only:set_warning_for_getm
  use global_interface,   ONLY: SurfaceCO2Diffusion
  use mem_param,   ONLY: AssignAirPelFluxesInBFMFlag
  use mem_CO2,only:pCO2_air

#endif
!
! !AUTHORS
!   16 March 1999 Original version by H. Thomas
!
!
!
! !REVISION_HISTORY
!   
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
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

  integer, INTENT(IN)    :: BoxNumber
  integer, INTENT(IN)    :: BoxNumberXY

#ifdef INCLUDE_PELCO2
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN)  :: jsurDIC
  real(RLEN)  :: r

       call SurfaceCO2Diffusion(2,ETW(BoxNumber),ESW(BoxNumber),& 
       ERHO(BoxNumber),Wind,Depth(BoxNumber),pCO2_air-pCO2(BoxNumber), &
                                                          jsurDIC)

       jsurO3c(BoxNumberXY)= jsurDIC*MW_C
       r=jsurO3c(BoxNumberXY)/Depth(BoxNumber)
       if (AssignAirPelFluxesInBFMFlag) call flux(BoxNumber,iiPel,ppO3c,ppO3c,r)

#endif

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
