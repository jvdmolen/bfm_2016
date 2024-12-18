#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Ecology
!
! DESCRIPTION
!   !	This submodel calls all other submodels
!
!

!   This file is generated directly from OpenSesame model code, using a code
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine EcologyDynamics
!
! !USES:
  ! The following 0-d global parameters are used: CalcPelagicFlag, &
  ! CalcBenthicFlag
  ! The following global constants are used: RLEN
  ! The following constants are used: BENTHIC_RETURN, BENTHIC_BIO, BENTHIC_FULL

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use bio,only: ActualN
  use bio_var,Only:julianday,var_names,cc_before_transport
  use global_mem, ONLY:RLEN,LOGUNIT,ZERO,DONE,NZERO
  use mem,  ONLY: iiPel,iiBen, iiReset,flux,Depth,NO_BOXES,iiConsumption, &
                   LocalDelta, max_change_per_step,max_rate_per_step
  use constants,  ONLY: BENTHIC_RETURN, BENTHIC_BIO, BENTHIC_FULL
  use mem_Param,  ONLY: CalcPelagicFlag, CalcBenthicFlag,&
      mass_conservation_check,p_max_state_change,p_max_state_rate
  use Track, ONLY:calc_all_track_rates,check_track_states
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following global functions are used:ResetTotMassVar
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_interface,   ONLY: ResetTotMassVar
  use mem,ONLY:Source_D3_vector,Source_D2_vector,NO_BOXES_XY,ppN1p,N1p, &
    D3SOURCE,D3SINK,K14n,ResetSource_D3_vector, ppK14n,ResetSource_D2_vector, &
    N1p, iiTotal,NO_D3_BOX_STATES,iiPhytoPlankton,ppPhytoPlankton, Source, &
    iiL,iiC,D3STATE,D6m,D7m,D8m,D9m,ppN3n,ppN4n,ppN5s,ppR6c,ppR6n,ppR1n,ppB1n, &
    R6c,R6n,R2c,ppY2c,jbotR1c,jbotR1n,jbotR1p,K4n,ppK24n,ppK3n,K3n, &
    PELBOTTOM,jbotR6c,jbotQ6c,ppP1c,ppR3c,R3c,P6c,Pcc,ppP6c,ppP1n,ppH2c,ppH2p, &
    ppK5s,ppH2n,ppH1c,ppH1p,ppH1n,ppR2c,Dfm,Dcm,ppDcm,ppDfm,ppBP1l,sediPI, &
    sediR6,ppQ21p,ppQ21n,ppHnc,D2Sink,ppK4n,ppQ6s
  use mem_phyto,only:p_qchlc
  use SourceFunctions,only: Source_D3_withstate

  use mem,ONLY: H1c,H2c,Hnc,H1n,H2n

  use BFM_ERROR_MSG,only:set_warning_for_getm
  use mem_BenPhyto,ONLY:test_BenPhyto_status
  use mem_Phyto,ONLY:test_Phyto_status
  use constants,ONLY:INITIALIZE,ADD
  use botflux,ONLY:setbotflux,checkbotflux
  use mem,ONLY:Bac,B1c
  use mem,ONLY:P1s,P5s,R6s,N5s
#ifdef INCLUDE_MACROPHYT
 use mem_Param, ONLY:CalcMacroPhyto
 use mem_MacroPhyto,ONLY:AddMacroPhytAddToMassConservationNP
 use bio_bfm,only:print_structure
#endif


 use mem,  ONLY: Output2d_1,Output2d_2,Output2d_3, Output2d_4
!
!
! !AUTHORS
!   ERSEM team	
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
  real(RLEN),external  ::GetDelta
  real(RLEN),dimension(NO_BOXES)   :: r,s
  real(RLEN)                       ::scalar,Dcm_old
  integer,save                     ::follow=0
  integer                          :: iout,i,j,k,l,mm

write(LOGUNIT,*) 'subroutine ecologydynamics'
!stop

  call ResetSource_D2_vector(ppQ6s)
  Output2d_1(1)=sum((P1s+P5s+R6s+N5s)*Depth)
  Output2d_2(1)=sum(R6s*Depth)

  call findnan(R2c,NO_BOXES,iout)
  if ( iout>0) then
         write(LOGUNIT,*) 'Ecology at start(1):NAN in R2c layer',iout
         call set_warning_for_getm
  endif
  call findnega(H1n,NO_BOXES_XY,iout)
  if ( iout>0) then
         write(LOGUNIT,*) 'Ecology at start(1):nega in H1n layer',iout
         call set_warning_for_getm
  endif
  call findnega(H2n,NO_BOXES_XY,iout)
  if ( iout>0) then
         write(LOGUNIT,*) 'Ecology at start(1):nega in Hnc layer',iout
         call set_warning_for_getm
  endif
  call findnan(Hnc,NO_BOXES,iout)
  if ( iout>0) then
         write(LOGUNIT,*) 'Ecology after Pelagic(1):NAN in Hnc layer',iout
         call set_warning_for_getm
  endif
  call findnan(K3n,NO_BOXES_XY,iout)
  if ( iout>0) then
         write(LOGUNIT,*) 'Ecology after Pelagic(1):NAN in K3n'
         call set_warning_for_getm
  endif

#ifdef DEBUG
  call  flux(1,iiReset,1,1,0.00D+00)
#endif
  if ( follow.eq.0) then
     follow=julianday
     write(LOGUNIT,*) 'Start of run..........follow=',follow
  endif

!write(LOGUNIT,*) 'start ecology'

  LocalDelta  =   GetDelta( )
  max_change_per_step=p_max_state_change/LocalDelta
  max_rate_per_step=p_max_state_rate/LocalDelta

  call check_track_states(1,'at start of the BFM calculations')

  call ResetTotMassVar( )


! if ( iout>0) then
!        write(logunit,*) 'Ecology at start:large in Dcm layer',Dcm(iout)
!        call set_warning_for_getm
! endif
  if ( follow ==julianday ) then
!   Dcm=Dfm
!   Pcc=P6c
!   if ( D6m(1).lt.0.0) D6m(1)=0.3
!   if ( D7m(1).lt.0.0) D7m(1)=0.3
!   if ( D8m(1).lt.0.0) D8m(1)=0.3
!   if ( D9m(1).lt.0.0) D9m(1)=0.3
!   if ( Dcm(1).gt.0.1) Dcm(1)=0.001
!    Bac=min(Bac,B1c)
  endif
  call findlarge(Dcm,NO_BOXES_XY,0.1D+00,iout)
  Dcm_old=Dcm(1)
  if ( iout>0) then
         write(logunit,*) 'Ecology at:large in Dcm layer',Dcm(iout)
         call set_warning_for_getm
  endif
  call findlarge(N1p,NO_BOXES,50.0D+00,iout)
  if ( iout>0) then
         write(logunit,*) 'Ecology at start:large in N1p layer',iout,N1p(iout)
         call set_warning_for_getm
  endif
  call findnan(K14n,NO_BOXES_XY,iout)
  if ( iout>0) then
         write(logunit,*) 'Ecology at start:Nan in K14n layer',iout,K14n(iout)
         call set_warning_for_getm
  endif

!write(LOGUNIT,*)'ecology before reset flux interface'
!stop

  ! Reset plagic-benthic fluxes-interface..............
  if (CalcBenthicFlag==BENTHIC_BIO .or. CalcBenthicFlag==BENTHIC_FULL)  &
                                  call setbotflux(INITIALIZE,NO_BOXES)
  !Calculate
  ! -initilialize (cellect) variables
  ! - calculate light climate in pelagic
  ! - calculate Chlorophyll concentrations
  ! - caluclate internal quaota's
  ! -oxygen saturation
  select case (CalcPelagicFlag)
    case (.true.)
!write(LOGUNIT,*)'ecology, pelagic flag true'
      if (CalcBenthicFlag==BENTHIC_BIO .or. CalcBenthicFlag==BENTHIC_FULL)  &
                         call test_BenPhyto_status(1)
!write(LOGUNIT,*)'call calcthermovars'
      call CalculateThermoVarsBFM
!write(LOGUNIT,*)'call pelagicsystemdynamics'
      call PelagicSystemDynamics(1)
!write(LOGUNIT,*)'after pelagicsystemdynamics'
!stop
      ! -pH in water column and in sediment
    case (.false.)
!write(logunit,*)'ecology, pelagic flag false'
      ! Calculate the minimum such that pelagic serve as 
      !dynamic boundary condition,
      call PelagicSystemDynamics(0)
  end select
!write(LOGUNIT,*)'findlarge',iiPel,ppN1p
!stop
!JM a pointer error occurs in this call. Don't understand what's going wrong. Skip for now!!!!!!!!!
  call FindLargeInRates(iiPel,ppN1p,'Ecology:after Pelagic_1')
!write(LOGUNIT,*)'CO2dynamics'
!stop
  call CO2Dynamics(1)
  call findnan(R6c,NO_BOXES,iout)
  if ( iout>0) then
         write(LOGUNIT,*) 'Ecology after Pelagic(1):NAN in R6c layer',iout
         call set_warning_for_getm
  endif

!write(LOGUNIT,*)'ecology: before benthicdiatomtest'
!stop

  !test if benthic diatoms are present on this grid point
!write(LOGUNIT,*) 'next call pelagic system dynamics'
!stop

  if ( CalcPelagicFlag) then
!write(LOGUNIT,*) 'calcpelagicflag'
!stop
    !test which phyton groups are presenta on this grid point
    call test_Phyto_status
!write(LOGUNIT,*) 'aftertest'
!stop
    ! Calculate Pelagic processes inclusive change in DIC
    call PelagicSystemDynamics(2)
!write(LOGUNIT,*) 'after pelagicsystem'
!stop
    call FindNaNInRates(iiPel,ppR3c,'Ecology:after Pelagic_2')
    call FindNaNInRates(iiPel,ppR2c,'Ecology:after Pelagic_2')
    call FindLargeInRates(iiPel,ppN1p,'Ecology:after Pelagic_2')
  endif
!write(LOGUNIT,*) 'past it'
!stop

  if (CalcBenthicFlag==BENTHIC_BIO .or. CalcBenthicFlag==BENTHIC_FULL) then
    !Calculate for for benthic model: use concentration form the lowest layer
    call PelForcingForBenDynamics(1)
    !check on negative values in layers,if so reset concentration to positive 
    ! values
    call BenPhosphateDynamics(0)
  endif
!write(LOGUNIT,*) 'after benphosphatedynamics',CalcBenthicFlag
!stop
  if (CalcBenthicFlag > 0 ) then

    select case ( CalcBenthicFlag)
      case ( BENTHIC_RETURN )  ! Simple benthic return
        call BenthicReturn1Dynamics

      case ( BENTHIC_BIO )  ! Intermediate benthic return
        call BenthicSystemDynamics

        call BenthicNutrient2Dynamics

      case ( BENTHIC_FULL )  ! Full benthic nutrients
!write(LOGUNIT,*) 'benthic_full'
!stop
        call PelForcingForBenDynamics(2)
!write(LOGUNIT,*) 'after pelforcingforbendynamics'
!stop
        !Calculate Cand coupled nutrient fluxes of functional groups in benthic system
        call BenthicSystemDynamics
!write(LOGUNIT,*) 'after benthicsystemdynamics'
!stop
        call FindNanInRates(iiPel,ppK3n,'Ecology:full benthic system')

        !Diagenetic model
                call BenthicNutrient3Dynamics
                call FindNanInRates(iiPel,ppK3n,'Ecology:full benthic nutrients')

!write(LOGUNIT,*) 'after benthicsnutrient3Ddynamics'
!stop
        
        !For good funtioning of benthic nutrient model it is necessary
        ! that are no empty buffere of the esstial nutrients.
        ! if not ,correct and give an message.
        ! Calculate Benthic DIC and Alkalinity fluxes.
        call CO2Dynamics(3)
!write(LOGUNIT,*) 'after CO2Ddynamics'
!stop
        if ( CalcPelagicFlag) then
          !Resuspension and settling of bnethic filterfeeder larvae
          ! (fluxes between Yy3 and Z2)
          call Y3Z2CoupDynamics
!write(LOGUNIT,*) 'after Y3Z2'
!stop

          call FindNanInRates(iiPel,ppR3c,'Ecology:after Y3Y2Dyn')
          call FindNanInRates(iiPel,ppR2c,'Ecology:after Y3Y2Dyn')

          ! Resuspension and of detritus   (fluxes beteween R6 and Q6)
          call ResuspensionPartDetritus
!write(LOGUNIT,*) 'after resuspensionPart'
!stop
          ! Resuspension of diatoms  (fluxes between P5 and BP1)
          call ResuspensionBenPhyto
!write(LOGUNIT,*) 'after resuspensionBenPhyto'
!stop
          ! Settling of phytoplankton,detritus
        endif
    end select
  endif
!write(LOGUNIT,*) 'after allbenthiccalls'
!stop
  call FindLargeInRates(iiPel,ppN1p,'Ecology:before PSF(3)')
  if ( CalcPelagicFlag) then
    call PelagicSystemDynamics(3)
    call CO2Dynamics(2)
    if (CalcBenthicFlag==BENTHIC_BIO .or. CalcBenthicFlag==BENTHIC_FULL)  & 
                                            call SettlingDynamics

  endif

!write(LOGUNIT,*) 'after calcpelagicflag'
!stop
   !Fill variables to test if no masis lows or created.
   ! ( can only test in a 1D-setup)

   ! Add calculated fluxes at the interface in the benthic system
   ! to the lowest  layer of the pelagic
   ! collected rates at sediment-waterface

   
! call checkbotflux(1,NO_BOXES,0.1D+00,Depth)
  if (CalcBenthicFlag > 0 ) call setbotflux(ADD,NO_BOXES,Depth)
!write(LOGUNIT,*) 'after setbotflux'
!stop

  call FindNanInRates(iiPel,ppR3c,'Ecology:after setbotflux')
!write(LOGUNIT,*) 'after first check'
!stop
  call FindNanInRates(iiBen,ppK14n,'Ecology:after setbotflux')
  call FindLargeInRates(iiPel,ppN1p,'Ecology:after setbotflux')
!          Output2d_4=Source_D2_vector(ppQ6s,0)

!write(LOGUNIT,*) 'after calcbenthicflag'
!stop

  if ( mass_conservation_check ) then
    call CheckMassConservationNPSDynamics
#ifdef INCLUDE_MACROPHYT
!skip     if (CalcMacroPhyto.and.CalcPelagicFlag )  &
!skip                                 call AddMacroPhytAddToMassConservationNP
#endif
     call CheckMassConservationCDynamics
  endif
  call ControlBenPartNutrientBuffersDynamics
  
  if (abs(Dcm(1)-Dcm_old)/Dcm_old> 1.0D-6) then
     write(LOGUNIT,*) "Dcm changed during rate defs"
  endif
! call print_structure(2,-ppK4n,1)

  call FindNanInRates(iiPel,ppK3n,'Ecology: at end')

!write(LOGUNIT,*)'end ecology'
!stop
  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
