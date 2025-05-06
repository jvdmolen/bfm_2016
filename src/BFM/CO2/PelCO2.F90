#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: PelCO2Dynamics
!
! DESCRIPTION
!   !

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  SUBROUTINE PelCO2Dynamics()

#ifdef INCLUDE_PELCO2
!
! !USES:
  ! The following Pelagic-states are used (NOT in fluxes): O3h, O3c
  ! The following global scalar vars are used: &
  ! BoxNumberZ, NO_BOXES_Z, BoxNumberY, NO_BOXES_Y, BoxNumberX, NO_BOXES_X, &
  ! BoxNumber, BoxNumberXY
  ! The following Pelagic 1-d global boxvars are modified:Ac,DIC, K0Ac, K1Ac, &
  ! K2Ac, KsAc, CO2, HCO3, CO3, pCO2, pH
  ! The following Pelagic 1-d global boxvars got a value: KbAc, KsiAc, KwAc, &
  ! KfAc
  ! The following Pelagic 1-d global boxvars  are used: ETW, ESW, ERHO
  ! The following Pelagic 2-d global boxvars got a value: KpAc
  ! The following 0-d global parameters are used: MethodCalcCO2
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,LOGUNIT
  use mem,  ONLY: O3h, O3c, D3STATE
  use mem, ONLY: BoxNumberZ, NO_BOXES_Z, BoxNumberY, &
    NO_BOXES_Y, BoxNumberX, NO_BOXES_X, BoxNumber, BoxNumberXY, Ac, DIC, &
    N3n,N4n,N1p,N5s,CO2, HCO3, CO3,CAc, pCO2, pH, ETW, ESW, ERHO
  use CO2System,ONLY: CalcCO2System
  use mem_CO2    
  USE BFM_ERROR_MSG, ONLY: BFM_ERROR,set_warning_for_getm
  use global_interface,only:SurfaceCO2Processes
  use constants,only:MW_C
!  
!
! !AUTHORS
!   H. Thomas   and  P. Ruardij
!
! !REVISION_HISTORY
!   calculates CO2 speciation major aim CO2 aequeous for CO2air-sea
!  exchange
! nutrient, carbon units 10^-6mol per liter = mmol per m^3
! sN1p_in = PO4 = N1.p, sN5s_in= SIO4=N5.s, O3.c/12.011=DIC, O3.h=Ac, ta=TA
! st = sulfate
! conversion to kgs required!!!!!!!!!!!!!!!
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
  integer            ::j,k=0
  integer            ::error=0
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! external functions
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer, external  :: D3toD1
  integer, external  :: D2toD1

  j=0;k=0;error=0
  DO BoxNumberZ=NO_BOXES_Z,1,-1
    DO BoxNumberY=1,NO_BOXES_Y
      DO BoxNumberX=1,NO_BOXES_X
        BoxNumber=D3toD1(BoxNumberX,BoxNumberY,BoxNumberZ)
        BoxNumberXY=D2toD1(BoxNumberX,BoxNumberY)

        Ac(BoxNumber)   =   O3h(BoxNumber)
        DIC(BoxNumber)  =   O3c(BoxNumber)/ MW_C
        error= CalcCO2System(MethodCalcCO2, &
                ESW(BoxNumber),ETW(BoxNumber),ERHO(BoxNumber),&
                N1p(BoxNumber),N5s(BoxNumber),Ac(BoxNumber),&
                CO2(BoxNumber),HCO3(BoxNumber),CO3(BoxNumber),CAc(BoxNumber),&
                DIC_in=DIC(BoxNumber), &
                pCO2_out=pCO2(BoxNumber),pH_out=pH(BoxNumber))
        if ( error > 0 ) then
           j=j+1
           if (k ==0) then
             write(LOGUNIT,*) 'BoxNumber',BoxNumber
             write(LOGUNIT,*) 'DIC',DIC(BoxNumber),O3c(BoxNumber)
             write(LOGUNIT,*) 'Ac',Ac(BoxNumber)
             write(LOGUNIT,*) 'N1p',N1p(BoxNumber)
             write(LOGUNIT,*) 'N5s',N5s(BoxNumber)
             write(LOGUNIT,*) 'N3n',N3n(BoxNumber)
             write(LOGUNIT,*) 'N4n',N4n(BoxNumber)
             write(LOGUNIT,*) 'rho',ERHO(BoxNumber)
             write(LOGUNIT,*) 'pCO2',pCO2(BoxNumber)
             write(LOGUNIT,*) 'pH',pH(BoxNumber)
!             k=1
           endif
           error= CalcCO2System(MethodCalcCO2, &
             ESW(BoxNumber),ETW(BoxNumber),ERHO(BoxNumber),&
             N1p(BoxNumber),N5s(BoxNumber),Ac(BoxNumber),&
             CO2(BoxNumber),HCO3(BoxNumber),CO3(BoxNumber),CAc(BoxNumber),&
             pH_in=9.D+00,DIC_in=DIC(BoxNumber), &
             pCO2_out=pCO2(BoxNumber))
write(LOGUNIT,*)'BoxNumber',BoxNumber
write(LOGUNIT,*)'Ac',Ac(BoxNumber)
write(LOGUNIT,*)'O3h',O3h(BoxNumber)
           O3h(BoxNumber)=Ac(BoxNumber)
           write(LOGUNIT,*) 'O3h=AC reset to',O3h(BoxNumber)
        endif

        if ( BoxNumber==NO_BOXES_Z ) then
          call SurfaceCO2Processes(BoxNumber,BoxNumberXY)
        endif

      end DO
      ! of compute
    end DO
  ENDDO
 if (j>0) then
   write(LOGUNIT,*) 'For this column ',j,' times pH was outside bounds (4-11)'
!stop
   call set_warning_for_getm()
 endif

#endif

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
