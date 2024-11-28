#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: SurfaceO2Diffusion
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
  subroutine SurfaceO2Diffusion(mode,t,Wind,depth,rdiff,flux)
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

  use global_mem, ONLY:RLEN


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following global functions are used:CalcSchmidtNumberCO2
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_interface,   ONLY: CalcSchmidtNumberOx
  use constants,ONLY: HOURS_PER_DAY
  use mem_WindOxReaeration_3

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
  real(RLEN), INTENT(IN)    :: Wind
  real(RLEN), INTENT(IN)    :: depth
  real(RLEN), INTENT(IN)    :: rdiff
  real(RLEN), INTENT(OUT)   :: flux

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN)  :: p_schmidt
  real(RLEN)  :: kex
  logical     :: check

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! air-sea exchange
        !  ( In lower boxes co2-input take place via vertical transport)
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        !
        ! The authors assumed a Schimdt number of CO2 (=reference) of 660.0
        !

        p_schmidt  =   CalcSchmidtNumberOx(t)/ 660.0D+00

        if (mode ==1 ) then

           !
           ! Calculate wind dependency:
           !	including conversion cm/hr => m/dag :
           !

           kex  = p_k* (Wind)**(2.0D+00)/ sqrt(  p_schmidt) 
        else
 
          if (Wind < 3.6) then
!JM             kex = 0.17  * Wind / p_schmidt**-0.6667;
             kex = 0.17  * Wind / p_schmidt**(-0.6667);
          elseif ( Wind < 13.0 ) then
             kex = (2.85  * Wind- 9.65) / sqrt(p_schmidt);
          else
             kex = (5.90  * Wind- 49.3) / sqrt(p_schmidt);
          endif
     
          ! Correct Exchange Coefficient from cm/h to m/d
          !
             kex =  kex * HOURS_PER_DAY * 0.01;

        endif


        !
        check=(mode.gt.0)
        select case (check)
           case (.false.) ; flux=kex*depth                              !m2/d
           case (.true. ) ; flux=kex*rdiff                              !C/m2/d
        end select


  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
