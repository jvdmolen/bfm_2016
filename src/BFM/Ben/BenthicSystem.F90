#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenthicSystem
!
! DESCRIPTION
!   !   This is the top level of the benthic submodel.
!       All the biological processes affecting the benthic dynamics are called
!       in a specified sequence according to the calculation flags.
!
!

!   This file is generated directly from OpenSesame model code, using a code
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine BenthicSystemDynamics
!
! !USES:
  ! The following Benthic-states are used (NOT in fluxes): Y1c, Y1n, Y1p, Y2c, &
  ! Y2n, Y2p, Y4c, Y4n, Y4p, Y5c, Y5n, Y5p, H1c, H1n, H1p, H2c 
  ! The following groupmember vars  are used: iiY1, iiY2, iiY4, iiY5, iiH1, iiH2
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,LOGUNIT,DONE,NZERO,ZERO
  use mem, ONLY:  ppY1c, ppY1n, ppY1p, ppY2c, ppY2n, ppY2p, &
    ppY4c, ppY4n, ppY4p, ppY5c, ppY5n, ppY5p,  &
    iiY1, iiY2, iiY4,iiY5, iiH1, iiH2, iiH3, iiHN,iiBen, &
    iiY3,iiYy3,  ppY3c, ppY3n, ppY3p, ppYy3c, ppYy3n, ppYy3p, &
    ppBP1c,ppBP1n,ppBP1p,ppBP1l,ppBP1s,iiBP1,ppK1p,ppK14n, ppDfm, &
    ppH1c,ppH2c,ppH3c,ppHNc, ppH1n,ppH2n,ppH3n,ppHNn, ppH1p,ppH2p,ppH3p,ppHNp 
  use mem_Param, ONLY: CalcBenOrganisms, CalcSuspensionFeeders,  &
                       CalcBenBacteria,CalcBenPhyto,combine_anabac
     use BFM_ERROR_MSG,only:set_warning_for_getm
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_interface, ONLY: BenOrganismDynamics, BenBacDynamics, &
     FilterFeederDynamics,BenPhytoDynamics,FindNanInRates
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


! use mem,ONLY: Output2d_1,Output2d_2,Output2d_3, Output2d_4
  use mem,ONLY: Dfm,NO_BOXES_XY,ppDfm,Source_D2_vector,ppBP1c,BP1c,ppQ6s,jPIQ6s, &
                   iiProduction,iiConsumption,ppQ1c,ppG3c,ppQ6c,ppG2o,D2SINK
  use SourceFunctions,only:Source_D2_withstate

!
!
! !AUTHORS
!   ERSEM group
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

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   REAL(RLEN),dimension(NO_BOXES_XY)::r1,p1
   INTEGER                          :: iout
   real(RLEN),external  ::GetDelta;


  ! mode ==1 : Bioturbation Dynamics is already called  in the 
  ! call BenthicNutrient3Dynamics(1) where some global variables calculated in 
  ! BenthicNutrient3Dynamics get already an value (pH, M5s,....)


  call findLarge(Dfm,No_BOXES_XY,1.0D+08,iout)
  if ( iout.gt.0) write(LOGUNIT,*) 'Dfm is large',Dfm

  call BenGlobalDynamics(1)
! call FindNaNInRates(iiBen,ppK14n,'BenthicSystem:after start')
! call FindNaNInRates(iiBen,ppDfm,'BenthicSystem:after start')


  if ( CalcBenOrganisms(iiY1)) &
    call BenOrganismDynamics( iiY1, ppY1c, ppY1n, ppY1p)
! call FindNaNInRates(iiBen,ppK14n,'BenthicSystem:after iiY1')
  

  if ( CalcBenOrganisms(iiY5)) &
    call BenOrganismDynamics( iiY5, ppY5c, ppY5n, ppY5p)
! call FindNaNInRates(iiBen,ppK14n,'BenthicSystem:after iiY5')


  if ( CalcBenOrganisms(iiY4)) &
    call BenOrganismDynamics( iiY4, ppY4c, ppY4n, ppY4p)
! call FindNaNInRates(iiBen,ppK14n,'BenthicSystem:after iiY4')

   r1=Source_D2_vector(ppDfm,0)
    if (r1(1)>DONE) write(LOGUNIT,*) 'After CalcB>enOrganism:Rate Dfm>1',r1(1)



  if ( CalcBenOrganisms(iiY2)) &
    call BenOrganismDynamics( iiY2, ppY2c, ppY2n, ppY2p)


! call FindNaNInRates(iiBen,ppK14n,'BenthicSystem:after iiY2')

  if ( CalcSuspensionFeeders(iiY3) .and. CalcSuspensionFeeders(iiYy3)) &
    call FilterFeederDynamics(iiYy3,ppYy3c,ppYy3n,ppYy3p)
  if ( CalcSuspensionFeeders(iiY3)) &
    call FilterFeederDynamics(iiY3,ppY3c,ppY3n,ppY3p)

! call FindNaNInRates(iiBen,ppK14n,'BenthicSystem:after FilterFeeder')


  call BioturbationDynamics

  !mode =1 : preprocessing
  call BenLimitNutrientDynamics

  if ( CalcBenBacteria(iiHN)) call BenNBacDynamics


! call FindNaNInRates(iiBen,ppK14n,'BenthicSystem:after iiHN')

  if ( CalcBenPhyto(iiBP1)) &
  call BenPhytoDynamics( iiBP1, ppBP1c, ppBP1n, ppBP1p, ppBP1s, ppBP1l)
!  write(LOGUNIT,*) 'BenPhytoDyna'
   r1=Source_D2_vector(ppDfm,0)
    if (r1(1)>DONE) write(LOGUNIT,*) 'After BenPhyto:Rate Dfm>1',r1(1)
! OUtput2d_2=jPIQ6s(1,:)

! call FindNaNInRates(iiBen,ppK14n,'BenthicSystem:after iiBP1')
  if ( CalcBenBacteria(iiH1)) &
    call BenBacDynamics( iiH1, ppH1c, ppH1n, ppH1p)

! D2SINK(ppH2c,:,:)=ZERO

  if ( CalcBenBacteria(iiH2))  &
    call BenBacDynamics( iiH2, ppH2c, ppH2n, ppH2p)
! call FindNaNInRates(iiBen,ppK14n,'BenthicSystem:after iiH2')


  if ( CalcBenBacteria(iiH3).and.(.not.combine_anabac))  &
    call BenBacDynamics( iiH3, ppH3c, ppH3n, ppH3p)


  !mode =2 : postprocessing
  call BenGlobalDynamics(2)

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
