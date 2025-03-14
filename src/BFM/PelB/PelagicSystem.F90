#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: PelagicSystem
!
! DESCRIPTION
!   This is the Pelagic Submodel. 
!   All the pelagic biogeochemical modules are called in sequence
!   according to the logical switches
!        
!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine PelagicSystemDynamics(mode)
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,LOGUNIT,ZERO,NZERO,DONE
! use BFM_ERROR_MSG,only:set_warning_for_getm
use mem, ONLY:R6c,P6c,Pcc
  use mem, ONLY: iiP2,iiP6, iiC, iiN,iiP,iiS,iiL,D3STATE, &
      ppPhytoPlankton,iiPhytoPlankton, ppMesoZooPlankton,iiMesoZooPlankton, &
      ppMicroZooPlankton,iiMicroZooPlankton,sunPI, NO_BOXES,eO2Mo2
  use mem_Param, ONLY: CalcPhytoPlankton,CalcMicroZooPlankton, &
    CalcMesoZooPlankton, CalcBacteria, CalcPelChemistry
  use mem_Phaeo,ONLY:NEW_COLONIES
  use mem_phyto,ONLY:CalcPhytoCopy
#ifdef INCLUDE_MACROPHYT
 use mem,ONLY:iiMacroContent,ppMacroContent,ppMacroStructure
 use mem_Param, ONLY:CalcMacroPhyto 
 use mem_MacroPhyto,ONLY: test_MacroPhyto_status
 use global_interface,only:MacroPhytoDynamics
#endif

! use mem,  ONLY: Output2d_1,Output2d_2,Output2d_3, Output2d_4,Source_D3_vector
 use mem,only:ppR6c,iiPel,ppN1p
! use SourceFunctions,only: Source_D3_withstate
! use mem_PelBac, ONLY:p_qpBc=>p_qpc
use mem,only:ppR2c

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following global functions are used:CalcChlorophylla, &
  ! CalcOxygenSaturation
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_interface,   ONLY: PhaeocystisCalc,PhytoDynamics,MesoZooDynamics, &
                          MicroZooDynamics,CalcVerticalExtinction,CalcLight  

 
  !-=-=-=-= -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following group processes are &
  ! used: PhytoDynamics, MesoZooDynamics, MicroZooDynamics
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!  
!
! !AUTHORS
!   ERSEM team
!
! !REVISION_HISTORY

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
  integer,intent(IN)::mode
  integer           ::i,ic,ip,in,il,is,j   ,iout
  real(RLEN)        ::rzero=0.0D+00,hold,r0
  character(len=10)  ::msg=""
  real(RLEN),dimension(NO_BOXES)    :: rx_any,sx_any,px_any !,c_xprevious_c,c_xprevious_p

!write(LOGUNIT,*)'start pelagicsystemdynamics',mode

  select case (mode)
  case (0,1) !------------------------------------------------------------------
!write(LOGUNIT,*)'CalcVertExt'
    call CalcVerticalExtinction(1)
!write(LOGUNIT,*)'Calclight'
    Call CalcLight
!write(LOGUNIT,*)'CalcChl'
    call  CalcChlorophylla( )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Compute oxygen variables: cmO2o eO2mO2
  ! calculate oxygen OxygenReaeration
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!write(LOGUNIT,*)'Calcoxsat'
    call  CalcOxygenSaturation
!write(LOGUNIT,*)'Calcoxreair'
    call OxygenReaerationDynamics
!write(LOGUNIT,*)'after OxygenReaearation'
    if (mode .eq.1 ) then
!write(LOGUNIT,*)'call testmacrophytostatus 1'
#ifdef INCLUDE_MACROPHYT
!write(LOGUNIT,*)'call testmacrophytostatus 2'
     if (CalcMacroPhyto) call test_MacroPhyto_status
#endif
!write(LOGUNIT,*)'call pelglobaldynamics'
    call PelGlobalDynamics(1)
!write(LOGUNIT,*)'pelglobaldynamics'
!stop
   endif
 

  case(2) !-------------------------------------------------------------------
!write(LOGUNIT,*)'pelagicsystem case2'
!stop

    !predictor-corrector
    call LimitNutrientUptake


!   call FindNaNInRates(iiPel,ppN1p,'After PelagicSystems:Phyto case 2')

      
#ifdef INCLUDE_MACROPHYT
!write(LOGUNIT,*)'calcmacrophyto'
   if (CalcMacroPhyto ) then
     do i=1,iiMacroContent
        is=ppMacroStructure(i,iiC)
        ic=ppMacroContent(i,iiC)
        in=ppMacroContent(i,iiN)
        ip=ppMacroContent(i,iiP)
        il=ppMacroContent(i,iiL)
        call MacroPhytoDynamics( i, is, ic, in, ip, il)
!       write(msg,'(I4)') i
     enddo
   endif
#endif
!write(LOGUNIT,*)'findnaninrates'
!stop
   call FindNaNInRates(iiPel,ppN1p,'Before PelagicSystems:MicroZoo')
!write(LOGUNIT,*)'microzooloop'
!stop
    do i=1,iiMicroZooPlankton
      if ( CalcMicroZooPlankton(i)) then
        ic=ppMicroZooPlankton(i,iiC)
        in=ppMicroZooPlankton(i,iiN)
        ip=ppMicroZooPlankton(i,iiP)
        call MicroZooDynamics( i, ic, in, ip)
       write(msg,'(I4)') i
!write(LOGUNIT,*)'inside microzooloop'
!stop
   call FindNaNInRates(iiPel,ppN1p,'After PelagicSystems:MicroZoo'//msg)
!write(LOGUNIT,*)'betweenfindnans'
!stop
   call FindNaNInRates(iiPel,ppR2c,'After PelagicSystems:MicroZoo'//msg)
!write(LOGUNIT,*)'afterbetweenfindnans'
!stop
      end if
    enddo

!write(LOGUNIT,*)'mesozooloop'
!stop

    do i=1,iiMesoZooPlankton
      if ( CalcMesoZooPlankton(i)) then
        ic=ppMesoZooPlankton(i,iiC)
        in=ppMesoZooPlankton(i,iiN)
        ip=ppMesoZooPlankton(i,iiP)
        call MesoZooDynamics( i, ic, in, ip)
      end if
    enddo
!write(LOGUNIT,*)'after mesozooloop'
!stop
    call FindNaNInRates(iiPel,ppN1p,'After PelagicSystems:MesoZoo')
!write(LOGUNIT,*)'after findnan'
!stop

    if ( CalcBacteria) call PelBacDynamics
!   call FindNaNInRates(iiPel,ppN1p,'After PelagicSystems:Bact')
!write(LOGUNIT,*)'after calcbacteria'
!stop
   call FindNaNInRates(iiPel,ppR2c,'After PelagicSystems:MicroZoo')
!write(LOGUNIT,*)'after findnan'
!stop

  case(3) !-----------------------------------------------------------------

!write(LOGUNIT,*)'pelagicsystemdynamics phytoplankton'

    do i=1,iiPhytoPlankton
      if ( CalcPhytoPlankton(i)) then
        ic=ppPhytoPlankton(i,iiC)
        in=ppPhytoPlankton(i,iiN)
        ip=ppPhytoPlankton(i,iiP)
        is=ppPhytoPlankton(i,iiS)
        il=ppPhytoPlankton(i,iiL)
        call PhytoDynamics( i, ic, in, ip, is, il)
       write(msg,'(I4)') i
!      call FindNaNInRates(iiPel,ppR6c,'After PelagicSystems:Phyto'//msg)
      end if
    enddo

!write(LOGUNIT,*)'pelagicsystemdynamics, pelchem'

    if ( CalcPelChemistry)  call PelChemDynamics
!   call FindNaNInRates(iiPel,ppN1p,'After PelagicSystems:PelChem')
  
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Calculate the biomass  of singleflagellates which are transferred into 
    ! colonies. For the calculation the total net nutrient uptake is needed
    ! A large nutrient uptake induces singel flagellates to make colonies.
    ! r is a dummy variable
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    if ( CalcPhytoCopy(iiP6) .and. CalcPhytoPlankton(iiP2)) then
!write(LOGUNIT,*)'pelagicsystemdynamicx,pheo'
       sx_any=sunPI(:,iiP2)
       call PhaeocystisCalc(NEW_COLONIES,iiP2,rx_any,sx_any,rzero)
    endif
!write(LOGUNIT,*)'pelagicsystemdynamics, call pelglobaldynamics'
    call PelGlobalDynamics(2)
  end select 

!write(LOGUNIT,*),'end pelagicsystemdynamics'

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

