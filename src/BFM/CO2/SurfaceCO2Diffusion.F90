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
  subroutine SurfaceCO2Diffusion(mode,t,s,rho,Wind,depth,rdiff,flux)
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
  use CO2System, ONLY:CalcK0Ac
  use gotm_error_msg,only:set_warning_for_getm

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following global functions are used:CalcSchmidtNumberCO2
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_interface,   ONLY: CalcSchmidtNumberCO2

#endif

!  
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
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

  integer,INTENT(IN)        :: mode
  real(RLEN), INTENT(IN)    :: t
  real(RLEN), INTENT(IN)    :: s
  real(RLEN), INTENT(IN)    :: rho
  real(RLEN), INTENT(IN)    :: Wind
  real(RLEN), INTENT(IN)    :: depth
  real(RLEN), INTENT(IN)    :: rdiff
  real(RLEN), INTENT(OUT)   :: flux

#ifdef INCLUDE_PELCO2
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN)  :: p_schmidt
  real(RLEN)  :: kex
  real(RLEN)  :: mkg_to_mmm3   ! mol/kg --> mmol/m3  (density to volume measure) 
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! air-sea exchange
      !  ( In lower boxes co2-input take place via vertical transport)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


       !
       ! The authors assumed a Schimdt number of CO2 (=reference) of 660.0
       !

        p_schmidt  =   CalcSchmidtNumberCO2(t)/ 660.0D+00


        !
        ! Calculate wind dependency:
        !	including conversion cm/hr => m/dag :
        !

        kex  =  ( 0.31D+00* (Wind)**(2.0D+00))/ sqrt(  p_schmidt)* 0.24D+00  
        ! kex is in cm / hr
        ! kex *.24 is then m/day

        !units of CO2 flux: kex * delta pCO2 * K0 = 10-3 mol m-2 d-1
        ! unit: (m/day) * (1.e-6atm) * (1.e-6mol / (l *1.e-6atm)) = 10-3 mol m-2
        ! d-1
        !
    
        mkg_to_mmm3=rho*1000.0;

!       CO2AirSea(1) = -(( kex*( pCO2(BoxNumber)- pCO2_air)* &
!         K0Ac* mkg_to_mmm3)/ 1000.0D+00)  ! flux co2 in mol/m2/d

         select case (mode)
           case (0) ; flux=kex*depth                              !m2/d
           case (1) ; flux=kex*rdiff                              !C/m2/d
           case (2) ; flux=kex*rdiff * CalcK0Ac(s,t)* mkg_to_mmm3 !C/m2/d:w
         end select
         if (isnan(flux)) then 
            write(LOGUNIT,*)'SurfaceCO2Diffusion kex,rdiff,wind,p_schmidt,t:',&
                                                kex,rdiff,Wind,p_schmidt,t
            write(LOGUNIT,*)'SurfaceCO2Diffusion s,t,CalcK0Ac(s,t):',&
               s,t , CalcK0Ac(s,t)
            write(LOGUNIT,*)'SurfaceCO2Diffusion mkg_to_mmm3,p_schmidt:', &
                          mkg_to_mmm3,p_schmidt
             call set_warning_for_getm 
         endif
#endif

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
