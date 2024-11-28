#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenthicNutrient3
!
! DESCRIPTION
!  Description of the diagenetic processes in the sediment
!
!   This file is generated directly from OpenSesame model code, using a code
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine BenthicNutrient3Dynamics
!
! !USES:
  ! The following global scalar vars got a value: LocalDelta
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,LOGUNIT
  use mem_BenthicNutrient3
  use mem_Param,only:p_qro
  use constants,only:p_qnUc
  use mem_BenPhyto,only:CalculateAverageNutrientConcInLightLayer
  use mem,only: ruBPn3,ruBTc,ruBTn,reBTc,reBTn,reATc,reATn, ruBTp,reBTp,reATp, &
    reK6o,ruBTs,reBTs,rrBTo,reBTo,rrDTo,rrATo,Ru_Benn,jbotR1c,jbotR1n, &
    NO_BOXES_XY, KQ1c,KQun,iiC,iiN,iiP,iiS,&
    ppK3n,ppK4n,ppK14n,ppK24n,ppK1p,ppK1p,ppK11p,ppK21p,ppK5s,ppG2o, &
    ppG3c,ppG13c,ppK6r,ppK16r,ppK26r, ppR1c,ppR1n,ppHNc,ppHNn,ppH2n,ppH2p, &
    iiProduction,iiConsumption,iiTotal, &
    ppBenPhyto,iiBenPhyto,ppBenbacteria,iiBenbacteria,ppBenOrganisms, &
    iiBenOrganisms, ppSusPensionFeeders,iiSusPensionFeeders
  use global_interface,only:BenQ1TransportCalculate,BenQ1TransportDynamics
  use SourceFunctions,only: Source_D2_withstate,Source_D2_withgroup

!   use mem,  ONLY: Output2d_1,Output2d_2,Output2d_3, Output2d_4, &
!    Source_D2_vector,ppQ1c,Q1c
   use mem,ONLY:ppQ1c,iiBen,Source_D2_vector,ppK14n,K14n
!
!
! !AUTHORS
!   P. Ruardij
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
!   the Free Software Foundation
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
  integer iout
  real(RLEN),dimension(NO_BOXES_XY)    ::r,s


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN), external  :: GetDelta
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Get actual time step for the calculation of the transient profile of &
  ! the nutrients
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rrBTo=Source_D2_withstate(ppG2o,ppG2o,iiConsumption)
  reK6o=Source_D2_withstate(ppK6r,ppK6r,iiProduction)/p_qro
  reBTo=Source_D2_withstate(ppG2o,ppG2o,iiProduction)
  rrDTo=Source_D2_withstate(ppK16r,ppK16r,iiProduction)/p_qro
  rrATo=rrDto+Source_D2_withstate(ppK26r,ppK26r,iiProduction)/p_qro
  !uptake by Benthic Diatoms and bacteria
  ruBTp=Source_D2_withgroup(ppK1p,ppBenPhyto,iiBenPhyto,iiP,iiConsumption) &
   +Source_D2_withgroup(ppK1p,ppBenbacteria,iiBenbacteria,iiP,iiConsumption)
  ! production(=mineralization) by
  ! benphyto (only in winter),bacteria,benthic organisms en filterfeeders
  reBTp=Source_D2_withgroup(ppK1p,ppBenPhyto,iiBenPhyto,iiP,iiProduction) &
   +Source_D2_withgroup(ppK1p,ppBenbacteria,iiBenbacteria,iiP,iiProduction) &
   +Source_D2_withgroup(ppK1p,ppBenOrganisms,iiBenOrganisms,iiP,iiProduction) &
   +Source_D2_withgroup(ppK1p,ppSusPensionFeeders,iiSusPensionFeeders, &
                                                           iiP,iiProduction)
   ! anoxic mineralization by anoxic bacteria ans benthic organisms
  reATp=Source_D2_withstate(ppK11p,ppH2p,iiTotal) &
       +Source_D2_withstate(ppK21p,ppH2p,iiTotal) &
   +Source_D2_withgroup(ppK11p,ppBenOrganisms,iiBenOrganisms,iiP,iiTotal)
  ruBPn3=Source_D2_withgroup(ppK3n,ppBenPhyto,iiBenPhyto,iiN,iiConsumption) &
   +Source_D2_withgroup(ppK3n,ppBenbacteria,iiBenbacteria,iiN,iiConsumption)
  ! Carbon fixzation by benthic diatoms and benthic nitrifiers.
  ruBTc=Source_D2_withgroup(ppG3c,ppBenPhyto,iiBenPhyto,iiC,iiConsumption) &
       +Source_D2_withstate(ppG3c,ppHNc,iiConsumption)
  !nitrate uptake by benphyto and aeorobic bateria
  ruBTn=ruBPn3 &
     +Source_D2_withgroup(ppK4n,ppBenPhyto,iiBenPhyto,iiN,iiConsumption) &
     +Source_D2_withgroup(ppK4n,ppBenbacteria,iiBenbacteria,iiN,iiConsumption)
  ! G3c produced by mineralization and decomposition of urea
  reBTc=Source_D2_vector(ppG3c,iiProduction)
  ! mineralization of N by bacteria (higher org produce urea)
  reBTn=Source_D2_withgroup(ppK4n,ppBenPhyto,iiBenPhyto,iiN,iiProduction)  &
      +Source_D2_withgroup(ppK4n,ppBenbacteria,iiBenbacteria,iiN,iiProduction)
  ! G13c produced by total anoxic mineralization and decomposition of urea
  reATc=Source_D2_vector(ppG13c,iiTotal)
  ! anoxic mineralization
  reATn=Source_D2_withstate(ppK14n,ppH2n,iiTotal) &
      +Source_D2_withstate(ppK24n,ppH2n,iiTotal)
  ! production and consumption of silica
  ruBTs=Source_D2_withgroup(ppK5s,ppBenPhyto,iiBenPhyto,iiS,iiConsumption)
  reBTs=Source_D2_withgroup(ppK5s,ppBenPhyto,iiBenPhyto,iiS,iiProduction)


! write(LOGUNIT,*) ' BenAmmoniumDynamics'
  r=Source_D2_vector(ppK14n,0)
  call findlarge(r,NO_BOXES_XY,50.0D+00,iout)
  if ( iout>0)  write(LOGUNIT,*) 'BN3 At start:rateK14n=',r
  call BenAmmoniumDynamics
  r=Source_D2_vector(ppK14n,0)
  call findlarge(r,NO_BOXES_XY,50.0D+00,iout)
  if ( iout>0)  write(LOGUNIT,*) 'BN3 After BenAmmonium:rateK14n=',r
! write(LOGUNIT,*) ' BenNitrateDynamics'
  call BenNitrateDynamics
  r=Source_D2_vector(ppK14n,0)
  call findlarge(r,NO_BOXES_XY,50.0D+00,iout)
  if ( iout>0)  write(LOGUNIT,*) 'BN3 After BenNitrate:rateK14n=',r
! write(LOGUNIT,*) ' BenAnoxicDynamics'
  call BenAnoxicDynamics
! write(LOGUNIT,*) ' BenOxygenDynamics'
  call BenOxygenDynamics(1)
! write(LOGUNIT,*) ' BenDenitriDepthDynamics'
  call BenDenitriDepthDynamics
  r=Source_D2_vector(ppK14n,0)
  call findlarge(r,NO_BOXES_XY,50.0D+00,iout)
  if ( iout>0)  write(LOGUNIT,*) 'BN3 After BenDenitriDepth:rateK14n=',r
! write(LOGUNIT,*) ' BenNitrogenShiftingDynamics'
  call BenNitrogenShiftingDynamics
! write(LOGUNIT,*) ' BenPhosphateDynamics'
  call BenPhosphateDynamics(1)
! write(LOGUNIT,*) ' BenSilicaDynamics(3)'
  call BenSilicaDynamics(3)
! write(LOGUNIT,*) ' BenQ1TransportCalculate'
  if (jbotR1n(1).gt.1000.0) then
     write(LOGUNIT,*)' BenthicNutrient3: jbotR1n is >1000',jbotR1n
  endif
  if (jbotR1c(1).gt.1000.0) then
     write(LOGUNIT,*)' BenthicNutrient3: jbotR1c is >1000',jbotR1c
  endif
  !Make gradient for a fraction of LOC: Urea in sediment for 2 layers:Q1un,Q11un
  call BenQ1TransportCalculate(KQun,iiN,R1_Benx=Ru_Benn)
  ! Calculates rates on basis of gradient for N fraction in urea
  call BenQ1TransportDynamics(KQun,iiN,ppR1x=ppR1n,R1_Benx=Ru_Benn)
  ! Calculates rates on basis of gradient for C fraction in urea
  call BenQ1TransportDynamics(KQun,iiC,ppR1x=ppR1c,R1_Benx=Ru_Benn, &
                                                                p_qnUc=p_qnUc)
  ! Make gradient for other LOC in sediment for 2 layers: (Q1c-Qun/p_qnUc,
  ! Q11c-Q1un/p_qnUc).It is assumed that LOC has such a small diffustion,
  ! that exchange with waer column is neglectable.
  ! The main aim of determining these rates is the determination of exchange
  ! due to changes in the thickness of the aerobic layer
  call BenQ1TransportCalculate(KQ1c,iiC,p_qnUc=p_qnUc)
  ! Calculates rates on basis of gradient for for the 3 constituents in LOC:
  !Q1,Q1n,Q1p
  call BenQ1TransportDynamics(KQ1c,iiC)
! write(LOGUNIT,*) ' BenQ1TransportDynamics'
  call FindNaNInRates(iiBen,ppQ1c,'Benthic Nutrient:at the end ')

  call CalculateAverageNutrientConcInLightLayer

! write(LOGUNIT,*) ' CalculateAverageNutrientConcInLightLayer'
  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
