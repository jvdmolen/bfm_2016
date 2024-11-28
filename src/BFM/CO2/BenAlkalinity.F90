#include "INCLUDE.h"
#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenAlkalinity
!
! DESCRIPTION
!   Description of the anoxic diagenitic processes in the sediment
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine BenAlkalinityDynamics
!

#ifdef INCLUDE_BENCO2

! !USES:

  ! For the following Benthic-states fluxes are defined: G13h, G3h
  ! The following Benthic-states are used (NOT in fluxes): D1m, D6m, D2m
  ! The following global scalar vars are used: &
  !    &
  !   BoxNumberXY, InitializeModel, LocalDelta
  ! The following Benthic 1-d global boxvars are modified : KHplus 
  ! The following Benthic 1-d global boxvars got a value: Acae, Acan
  ! The following Benthic 1-d global boxvars are used: rrBTo, KQ1, &
  ! irrenh, ETW_Ben, rrATo, O3h_Ben, shiftD1m
  ! The following Benthic 1-d global boxpars  are used: p_poro
  ! The following 0-d global parameters are used: p_d_tot, p_q10diff
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  USE BFM_ERROR_MSG, ONLY: set_warning_for_getm
  use global_mem, ONLY:RLEN,ZERO,NZERO,DONE,dummy,idummy
  use mem,ONLY: G23h,G13h,G3h, D1m,DH2m,DH3m, D6m,D2m, D2STATE
  use mem, ONLY: ppG23h, ppG13h, ppG3h,ppO3h, &
    LocalDelta,sK4K3,jK3G4n,  jK23G4n, dry_z, &
    NO_BOXES_XY,BoxNumber,BoxNumberXY,PelBoxAbove, &
    InitializeModel, KHplus,KNO3,KNH4,KRED, ruBPn3,irrenh,  &
    ETW_Ben,jK4K3n,jG2K7o,rrATo,rrDTo,O3h_Ben,ctO3m2h, &
    shiftD1m,shiftD2m, iiPel,iiBen, flux,jG33G23h,max_change_per_step,Source, &
    CO2Ae,HCO3Ae,CO3ae,DICae,CO2An,HCO3An,CO3an,DICan
  use constants, ONLY: GET, LABDA_1, DIFFUSION, COEFFICIENT,&
    LAYERS, LAYER1, LAYER2, LAYER3,LAYER4,LAYER5,FOR_ALL_LAYERS, POROSITY, &
    ADSORPTION, DOUBLE_DEFINE, DEFINE, QUADRATIC_TERM, &
    LINEAR_TERM, CONSTANT_TERM, SET_CONTINUITY, FLAG, MASS, SET_BOUNDARY, &
    EQUATION, INPUT_TERM, PARAMETER, ZERO_EXPONENTIAL_TERM, DIST_EXPONENTIAL_TERM, &
    SET_LAYER_INTEGRAL_UNTIL, SET_LAYER_INTEGRAL,LABDA_2, &
    ADD, DERIVATIVE, RFLUX, SHIFT,SEC_PER_DAY,ANY,POSITIVE,NEGATIVE
  use mem_Param,  ONLY: p_poro, p_d_tot, p_d_tot_2,p_q10diff, p_dry_ben, &
                 p_qon_dentri,p_clDxm,p_qro,combine_anabac
  use mem_Diffusion,ONLY: p_diff=>p_diff_O3h
  use botflux,onlY:addbotflux
  use mem_BenAlkalinity
  use mem_BenAmmonium,only:p_slK4K3
  use mem_BenNitrate,only:p_slK3G4
  use mem_CO2,ONLY:p_qhK4K3n,p_qhATo,p_qhK3G4n
  use LimitRates, ONLY:LimitChange,LimitShift,LimitShift_m3

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following bennut functions are used:GetInfoFromSet, &
  ! InitializeSet, DefineSet, CompleteSet, CalculateSet, CalculateTau, &
  ! CalculateFromSet
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use bennut_interface, ONLY: GetInfoFromSet, InitializeSet, DefineSet, &
    CompleteSet, CalculateSet, CalculateTau, CalculateFromSet,GetInfoFromSet


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following global functions are used:eTq
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: IntegralExpDist,ExpDist,SetExpDist,IntegralExp,insw
  

  use global_interface,   ONLY: eTq

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following sesame functions are used:IntegralExp
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!  
!
! !AUTHORS
!   Original version by  P. Ruardij
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

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer     :: i,llfo
  real(RLEN)  :: r,s  !general help variables
  real(RLEN)  :: loss
  real(RLEN)  :: jflux
  real(RLEN)  :: lambda
  real(RLEN)  :: gamma,zeta
  real(RLEN)  :: alpha,beta
  real(RLEN)  :: mdiff0,mdiff1
  real(RLEN)  :: diff0,diff1
  real(RLEN)  :: Tau
  real(RLEN)  :: cG3h
  real(RLEN)  :: zuD1o,zuD2o
  real(RLEN)  :: jcG3O3h,jG3O3h,jG13G3h
  real(RLEN)  :: jG23G13h
  real(RLEN)  :: Dnew
  real(RLEN)  :: jG14G13h,jG24G23h
  real(RLEN)  :: sK13G4n,sK13G4h,sK13G4o
  real(RLEN)  :: sK4K3h,sOh
  real(RLEN)  :: jK16G4o,jK26G4o
  real(RLEN)  :: corr_1,corr_2
  real(RLEN)  :: n11,n12,n13,n14,r15,a15,n15
  !Beaware:
  ! rrATo: total anoxic mineralization, rrDTo: mineralization  in dentri layer
  real(RLEN)  :: jrDTo,jrATo ! respiration resp. in denitri and anoxic layer
  real(RLEN)  :: n21,n22,n31
  real(RLEN)  :: Dxm,Dym
  real(RLEN)  :: limit_rate_0,limit_rate_1,limit_rate_2,limit_rate_3
  real(RLEN)  :: limit_rate_2a,limit_rate_0a,limit_rate_3a,limit_rate_4

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  DO BoxNumberXY=1,NO_BOXES_XY
    BoxNumber=PelBoxAbove(BoxNumberXY)

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=
    ! Check if there is a first order process (=an exponential term ) in oxic 
    ! layer defined in the Anoxic Submodel ( profile in KRED). 
    ! The same type of process will be used for the Alklinity Model.
    r=GetInfoFromSet(KRED(BoxNumberXY),-GET,idummy,11)
    llfo=0;if (r>ZERO) llfo=1;

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! NH4-minrealization
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    if ( combine_anabac) then
      alpha  =   SetExpDist(D6m(BoxNumberXY),Dlm=p_clDxm,Dmm=p_d_tot)
      beta=alpha
      zuD1o = max( 1.D-20, rrATo(BoxNumberXY))/ p_poro(BoxNumberXY)/ &
                      IntegralExpDist( -alpha, p_d_tot_2- D1m(BoxNumberXY))
      zuD2o  =   zuD1o* ExpDist( - alpha, D2m(BoxNumberXY)- D1m(BoxNumberXY))
      jrATo= zuD2o * IntegralExpDist( -alpha, p_d_tot_2- D2m(BoxNumberXY)) &
                                                      *p_poro(BoxNumberXY) 
      jrDTo=rrAto(BoxNumberXY)-jrATo
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    else
      alpha  =   SetExpDist(DH2m(BoxNumberXY),p_clDxm)
      beta  =    SetExpDist(DH3m(BoxNumberXY),p_clDxm)
      jrDTo=rrDTo(BoxNumberXY)
      jrATo=rrATo(BoxNumberXY)-rrDTo(BoxNumberXY)
      zuD1o = max( 1.D-20, jrDTo)/ p_poro(BoxNumberXY)/ &
                    IntegralExpDist( -alpha, D2m(BoxNumberXY)- D1m(BoxNumberXY))
      zuD2o = max( 1.D-20, jrATo)/ p_poro(BoxNumberXY)/ &
                    IntegralExpDist( -beta, p_d_tot_2- D2m(BoxNumberXY))
    endif

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! H+-loss  due to nitrification and deoxidation of H--
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    loss= p_qhK4K3n*( jK4K3n(BoxNumberXY)-ruBPn3(BoxNumberXY) ) &
                                  +  p_qhATo* jG2K7o(BoxNumberXY) 

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Temperature Correction (diffusion)
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    if  (p_diff==ZERO) then
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Assume that Alkalnity diffusion is equal to the one of the HCO3-
      ! Calc Diffusion of HCO3- at actual temperature according Zeebe,RE(2011)
      ! Zeebe,RE(2011)
      ! Geochimica et Cosmochimica Acta 75,2483-2948
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      call  CalcDiffusionCO2(ppG3h,ETW_Ben(BoxNumberXY),CO2ae(BoxNumberXY), &
        HCO3ae(BoxNumberXY), CO3ae(BoxNumberXY),DICae(BoxNumberXY),mdiff0)
      call  CalcDiffusionCO2(ppG3h,ETW_Ben(BoxNumberXY),CO2an(BoxNumberXY), &
          HCO3an(BoxNumberXY),CO3an(BoxNumberXY),DICan(BoxNumberXY),mdiff1)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Correction for porosity and bioirrigation
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      diff0 = mdiff0*SEC_PER_DAY* irrenh(BoxNumberXY)* p_poro(BoxNumberXY)
      diff1 = mdiff1*SEC_PER_DAY* irrenh(BoxNumberXY)* p_poro(BoxNumberXY)
    else
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! (Classic) Temperature Correction (diffusion)
      ! Correction for porosity and bioirrigation
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      diff0 = p_diff*SEC_PER_DAY* irrenh(BoxNumberXY)* p_poro(BoxNumberXY) *& 
                                       eTq( ETW_Ben(BoxNumberXY), p_q10diff)
      diff1=diff0
    endif


    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Get coefficients describing ammonium in the oxic layer :
    ! 1. gamma of the exponential curve
    ! 2. parameter of the nitrification term
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    gamma = GetInfoFromSet( KNH4(BoxNumberXY), GET, LABDA_1, 11)
    sK4K3h = p_qhK4K3n*GetInfoFromSet( KNH4(BoxNumberXY), GET, LABDA_2, 11)

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Get coefficients of all terms of equation valid for the
    ! first layer of ammonium (integration constants)
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    n11 = -sK4K3h * GetInfoFromSet( KNH4(BoxNumberXY), GET, COEFFICIENT, 11)
    n12 = -sK4K3h * GetInfoFromSet( KNH4(BoxNumberXY), GET, COEFFICIENT, 12)
    a15 = -sK4K3h * GetInfoFromSet( KNH4(BoxNumberXY), GET, COEFFICIENT, 15)

    if (llfo==1) then
       zeta = GetInfoFromSet( KRED(BoxNumberXY), GET, LABDA_1, 11)
       sOh=GetInfoFromSet( KRED(BoxNumberXY), GET, LABDA_2, 11)*p_qhATo/p_qro
       n13 = -sOh * GetInfoFromSet( KRED(BoxNumberXY), GET, COEFFICIENT, 11)
       n14 = -sOh * GetInfoFromSet( KRED(BoxNumberXY), GET, COEFFICIENT, 12)
       r15 = -sOh * GetInfoFromSet( KRED(BoxNumberXY), GET, COEFFICIENT, 15)
       n15  =  a15 +  r15 &
       + ruBPn3(BoxNumberXY)/D1m(BoxNumberXY)/p_poro(BoxNumberXY)
    else
      n15=sK4K3h* a15-p_qhATo*jG2K7o(BoxNumberXY)/(D1m(BoxNumberXY) &
                                                     *p_poro(BoxNumberXY))
    endif

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Recalculate Mineralization m2 --> m3 porewater
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! anoxidation rate at interface D1m 

    sK13G4n = GetInfoFromSet( KNO3(BoxNumberXY), GET, LABDA_2, 21) 
    sK13G4h =sK13G4n  * p_qhK3G4n
    sK13G4o =sK13G4n / p_qon_dentri
    ! Calculate proportion between corrected flux from the nitrate model and
    ! uncorrected flux
    jK16G4o=jK3G4n(BoxNumberXY) /p_qon_dentri 
    jK26G4o=jK23G4n(BoxNumberXY)/p_qon_dentri
    if ( max(rrATo(BoxNumberXY),jK16G4o) > 1.0D-10 ) then
      r=jrDTo-jK16G4o;corr_1=r/(jrDTo+NZERO)
      r=jrATo-jK26G4o;corr_2=r/(jrATo+NZERO)
    else
      zuD1o=1.0D-6;zuD2o=ZERO;
      corr_1=ZERO;corr_2=ZERO;
      jrATo=ZERO;
    endif

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Get coefficients describing Nitrate in anoxic layer :
    ! 1. lambda of the exponential curve, and the denitrification rate
    ! 2. parameter of the denitrification term (integration constant)
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    lambda = abs(GetInfoFromSet( KNO3(BoxNumberXY), GET, LABDA_1, 21))
    r=max(p_slK3G4* p_qhK3G4n, sK13G4h)
    n21 =  GetInfoFromSet( KNO3(BoxNumberXY), GET, COEFFICIENT, 21)* r
    n22 =  GetInfoFromSet( KNO3(BoxNumberXY), GET, COEFFICIENT, 22)* r
    n31 =  GetInfoFromSet( KNO3(BoxNumberXY), GET, COEFFICIENT, 31)* r

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Initialize and input physical boundaries and forcing:
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! (p_flux_at_deep_end=1) assumption is made for calculation of gradient 
    ! that the boundary condition for this calculation  at 30 cm  is 0
    i=0;if (p_flux_at_deep_end>1)i=1;
    KHplus(BoxNumberXY) = InitializeSet( KHplus(BoxNumberXY), 5, 18+i+llfo*2)

    Dxm=(D1m(BoxNumberXY)+D2m(BoxNumberXY)) *0.5
    Dym=(D2m(BoxNumberXY)+p_d_tot) *0.5
    call DefineSet(KHplus(BoxNumberXY),LAYERS,LAYER1,LAYER2,&
                                                     D1m(BoxNumberXY),Dxm)
    call DefineSet(KHplus(BoxNumberXY),LAYERS,LAYER3,LAYER4,&
                                                     D2m(BoxNumberXY),Dym)
    call DefineSet(KHplus(BoxNumberXY),DIFFUSION, LAYER1, LAYER2, diff0,diff1)
    call DefineSet(KHplus(BoxNumberXY),DIFFUSION, LAYER3, LAYER4, diff1,diff1)
    call DefineSet(KHplus(BoxNumberXY),DIFFUSION, LAYER5, 0, diff1,dummy)

    call DefineSet(KHplus(BoxNumberXY),POROSITY, FOR_ALL_LAYERS,0,&
                                                  p_poro(BoxNumberXY), ZERO)
    call DefineSet(KHplus(BoxNumberXY),ADSORPTION,FOR_ALL_LAYERS,0,p_p, ZERO)

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Give particular solution for all layers:
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    call DefineSet(KHplus(BoxNumberXY), DOUBLE_DEFINE, 11, &
                             ZERO_EXPONENTIAL_TERM, gamma, sK4K3h)

    call DefineSet(KHplus(BoxNumberXY), DOUBLE_DEFINE, 12, &
                             ZERO_EXPONENTIAL_TERM, - gamma, sK4K3h)

    if ( llfo==1) then
      call DefineSet(KHplus(BoxNumberXY), DOUBLE_DEFINE, 13, &
                             ZERO_EXPONENTIAL_TERM, zeta, sOh)

      call DefineSet(KHplus(BoxNumberXY), DOUBLE_DEFINE, 14, &
                             ZERO_EXPONENTIAL_TERM, - zeta, sOh)
    endif

    call DefineSet(KHplus(BoxNumberXY),DEFINE,15,QUADRATIC_TERM, dummy, dummy)
    call DefineSet(KHplus(BoxNumberXY),DEFINE,16,LINEAR_TERM, dummy, dummy)
    call DefineSet(KHplus(BoxNumberXY),DEFINE,17,CONSTANT_TERM, dummy, dummy)

    call DefineSet(KHplus(BoxNumberXY),DEFINE,21,DIST_EXPONENTIAL_TERM, &
                                                                -alpha, dummy)
    call DefineSet(KHplus(BoxNumberXY),DOUBLE_DEFINE,22, &
                                 ZERO_EXPONENTIAL_TERM, -lambda, sK13G4h)
    call DefineSet(KHplus(BoxNumberXY),DOUBLE_DEFINE,23, &
                                  ZERO_EXPONENTIAL_TERM, +lambda, sK13G4h)

    call DefineSet(KHplus(BoxNumberXY), DEFINE,25, LINEAR_TERM, dummy, dummy)
    call DefineSet(KHplus(BoxNumberXY), DEFINE,25, CONSTANT_TERM, dummy,dummy)

    call DefineSet(KHplus(BoxNumberXY), DEFINE,31, DIST_EXPONENTIAL_TERM, &
                                                                 -alpha, dummy)
    call DefineSet(KHplus(BoxNumberXY), DOUBLE_DEFINE, 32, &
                                 ZERO_EXPONENTIAL_TERM, -lambda, sK13G4h)
    call DefineSet(KHplus(BoxNumberXY), DEFINE,34, LINEAR_TERM, dummy, dummy)
    call DefineSet(KHplus(BoxNumberXY), DEFINE,35, CONSTANT_TERM,dummy, dummy)

    call DefineSet( KHplus(BoxNumberXY),DEFINE,41, DIST_EXPONENTIAL_TERM, &
                                                                -beta, dummy)
    call DefineSet(KHplus(BoxNumberXY),DEFINE, 44, LINEAR_TERM, dummy, dummy)
    call DefineSet(KHplus(BoxNumberXY),DEFINE, 45, CONSTANT_TERM,dummy, dummy)
    if (p_flux_at_deep_end>1) &
       call DefineSet(KHplus(BoxNumberXY),DEFINE,51, DIST_EXPONENTIAL_TERM, &
                                                                -beta, dummy)
    call DefineSet(KHplus(BoxNumberXY),DEFINE,55, CONSTANT_TERM, dummy, dummy)

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Insert boundary conditions:
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    !1-8:1-8 boundary conditions:
    call CompleteSet(KHplus(BoxNumberXY), SET_CONTINUITY, FLAG, MASS,&
                                                            dummy,dummy)

    !9:9th boundary condition:
    call CompleteSet(KHplus(BoxNumberXY), SET_BOUNDARY, LAYER1, &
                        EQUATION, ZERO, value=O3h_Ben(BoxNumberXY))

    !10-11:10-11 boundary condition:
    select case  ( InitializeModel ) 
      case (0)
        call CompleteSet(KHplus(BoxNumberXY), SET_LAYER_INTEGRAL, &
                               LAYER2, LAYER3, dummy, value=G13h(BoxNumberXY))
        call CompleteSet(KHplus(BoxNumberXY), SET_LAYER_INTEGRAL_UNTIL, &
                         LAYER4, LAYER5, p_d_tot_2, value=G23h(BoxNumberXY))
      case(1)
        ! intiail coase will not work!!!!!
    end select

   !12-13:12-13:
    call CompleteSet(KHplus(BoxNumberXY),INPUT_TERM,22, &
                                               PARAMETER,ZERO,value=n21)
    call CompleteSet(KHplus(BoxNumberXY),INPUT_TERM,32, &
                                               PARAMETER,ZERO,value=n31)
    if (p_flux_at_deep_end>1) &
   call FixProportionCoeff(KHplus(BoxNumberXY),41,51, &
                                 DONE,ExpDist(-beta , Dym-D2m(BoxNumberXY)))
!   call CompleteSet(KHplus(BoxNumberXY),INPUT_TERM,23, &
!                                                PARAMETER,ZERO,value=n22)
    !14:14
    r=p_qhATo*zuD1o
    call CompleteSet( KHplus(BoxNumberXY),INPUT_TERM,21, &
                                               PARAMETER,ZERO,value=r)
    r=p_qhATo*zuD2o
    call CompleteSet( KHplus(BoxNumberXY),INPUT_TERM,31, &
                                               PARAMETER,ZERO,value=r)

    !15-16:15-16:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! a11 / (labda * labda * diff) = a12 / (labda * labda * diff)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      call CompleteSet(KHplus(BoxNumberXY),INPUT_TERM,12,PARAMETER,ZERO, &
                                                 value=n12)
    if ( llfo==1) then
      call CompleteSet(KHplus(BoxNumberXY),INPUT_TERM,14, &
                                                 PARAMETER,ZERO,value=n14)
      call CompleteSet(KHplus(BoxNumberXY),INPUT_TERM,11, &
                                                 PARAMETER,ZERO,value=n11)
      call CompleteSet(KHplus(BoxNumberXY),INPUT_TERM,13, &
                                                 PARAMETER,ZERO,value=n13)
    else
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! a11 / (labda * labda * diff) = a15 / (2 * diff)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      call CompleteSet(KHplus(BoxNumberXY),INPUT_TERM,11, &
                                                 PARAMETER,ZERO,value=n11)
    endif

    !17:18
    call CompleteSet(KHplus(BoxNumberXY),INPUT_TERM,15, &
                                                 PARAMETER,ZERO,value=n15)

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! 1. Calculate for above defined set of boundary conditions
    !      the gradient of the nutrient
    ! 2. Replace last condition by an alternative new one
    ! 3. Calculate the value belongin by the alternative condition with
    !  the gradient calculated under 1.
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    cG3h = CalculateSet( KHplus(BoxNumberXY), SET_LAYER_INTEGRAL, LAYER1, &
      LAYER1, ZERO, ZERO)

    if ( InitializeModel== 0) then
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate adaptation time absolute and Delta() relative to Tau:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      r=abs(loss)/(NZERO+max(ZERO,G3h(BoxNumberXY)))
      Tau=CalculateTau(r,diff0,p_p,D1m(BoxNumberXY))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Estimate the average value of K4n during the next time step:
      ! Hence this value depend as well as on adaptation time,
      ! ''old'' value, and on ''equilibrium value''
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      cG3h = cG3h+( G3h(BoxNumberXY)- cG3h)* IntegralExp( - LocalDelta/ &
        Tau, DONE)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Recalulate gradient, now using cG3h
      ! Complete first the alternative condition by inputting the value
      ! for G3h (cG3h)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      r = CalculateSet( KHplus(BoxNumberXY), ADD, 0, 0, dummy, cG3h)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! flux affecting G13h and G23h (fixed fluxes)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      jG14G13h= p_qhATo* (corr_1* jrDTo + jK16G4o)
      jG24G23h= p_qhATo* (corr_2* jrATo            + jK26G4o)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! flux at sediment water interface
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      if ( dry_z(BoxNumberXY) <=p_clDxm .and. p_dry_ben) then
        jG3O3h = ZERO 
      else
        jG3O3h = CalculateFromSet(KHplus(BoxNumberXY),DERIVATIVE,RFLUX, &
                                                                ZERO, ZERO) 
      endif

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! flux at D1m
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-

      jG13G3h = CalculateFromSet( KHplus(BoxNumberXY), DERIVATIVE, &
        RFLUX, D1m(BoxNumberXY), ZERO)

      Dnew  =   D1m(BoxNumberXY)+ LocalDelta* shiftD1m(BoxNumberXY)
      jG13G3h = jG13G3h+ CalculateFromSet( KHplus(BoxNumberXY), SHIFT, &
          LAYER1, D1m(BoxNumberXY), Dnew)/ LocalDelta

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-
      ! flux at D2m
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      jG23G13h =  CalculateFromSet( KHplus(BoxNumberXY), DERIVATIVE, RFLUX, &
        D2m(BoxNumberXY), ZERO)
      r=jG23G13h

      Dnew  =   D2m(BoxNumberXY)+ LocalDelta* shiftD2m(BoxNumberXY)
      jG23G13h = jG23G13h+ CalculateFromSet( KHplus(BoxNumberXY), SHIFT, &
        LAYER3, D2m(BoxNumberXY), Dnew)/ LocalDelta

!     if ( abs(jG23G13h) > 1.0D2 ) then
!       write(LOGUNIT,*) 'Rate jG23G13h:' ,jG23G13h
!       write(LOGUNIT,*) 'State G23h,G13h:',G23h(BoxNumberXY),G13h(BoxNumberXY)
!       write(LOGUNIT,*) 'State G23h,G13h:',G3h(BoxNumberXY),O3h_Ben(BoxNumberXY)
!       write(LOGUNIT,*) 'jG23G13h flux:',r
!       write(LOGUNIT,*) 'jG23G13h shift:',jG23G13h-r
!         write(LOGUNIT,*) 'D2m shift:',shiftD2m(BoxNumberXY)
!       write(LOGUNIT,*) 'jG14G13h',jG14G13h
!       write(LOGUNIT,*) 'jG24G13h',jG24G23h
!       call set_warning_for_getm
!     endif
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! flux at underside
      ! 1=no flux, 2= only fluxes_downwards (sink), 3,=full_flux
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       select case (p_flux_at_deep_end)  
         case(1); jflux=ZERO
         case(2); jflux=min(ZERO,CalculateFromSet(KHplus(BoxNumberXY),&
                                      DERIVATIVE,RFLUX, p_d_tot_2, ZERO))
         case(3); jflux=         CalculateFromSet(KHplus(BoxNumberXY),&
                                      DERIVATIVE, RFLUX, p_d_tot_2, ZERO)
      end select
      jG33G23h(BoxNumberXY) = jflux

      ! Damp for too large fluxes

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! limit_rates
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      call LimitChange(ANY,jG33G23h(BoxNumberXY),G23h(BoxNumberXY), &
                                         max_change_per_step,limit_rate_4)

      call LimitShift(jG23G13h,G13h(BoxNumberXY),G23h(BoxNumberXY), &
                                       max_change_per_step,limit_rate_3a)
      limit_rate_3=1
      s=p_d_tot_2-D2m(BoxNumberXY); r=D2m(BoxNumberXY)-D1m(BoxNumberXY)
      call LimitShift_m3(jG23G13h,G13h(BoxNumberXY),G23h(BoxNumberXY), &
             r,s,shiftD2m(BoxNumberXY),max_change_per_step,limit_rate_3)
      limit_rate_3=min(limit_rate_3,limit_rate_3a)
       
      !Calculate first rate at sediment water interface.
      !If volume of the water above sediment limit the shift to the water
      !Consider the first layer (nearly) as a part of the pealgic system:
      call LimitChange(NEGATIVE,jG3O3h,ctO3m2h(BoxNumberXY), &
                                   max_change_per_step,limit_rate_0)
      call LimitChange(ANY,jG3O3h,G3h(BoxNumberXY), &
                                   max_change_per_step,limit_rate_0a)
      limit_rate_0=min(limit_rate_0,limit_rate_0a)

      r=D1m(BoxNumberXY); s=D2m(BoxNumberXY)-D1m(BoxNumberXY)
      call LimitShift_m3(jG13G3h,G3h(BoxNumberXY),G13h(BoxNumberXY), &
                r,s,shiftD1m(BoxNumberXY),max_change_per_step,limit_rate_2)
      call LimitShift(jG13G3h,G3h(BoxNumberXY),G13h(BoxNumberXY), &
                                         max_change_per_step,limit_rate_2a)
      !limit the concentration that input rates cannot  too big!
      limit_rate_2=min(limit_rate_0,limit_rate_2,limit_rate_2a)

      !Forces flux to water column, especially when D1m is low, 
      !to keep G3h at appropiate values....
      s=loss + jG3O3h*insw(jG3O3h)
      r=G3h(BoxNumberXY)+(jG13G3h*limit_rate_2+  &
                                           jG3O3h*insw(-jG3O3h))*LocalDelta
      call LimitChange(NEGATIVE,s,r,    max_change_per_step,limit_rate_1)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! calculate limitation  and correct the fluxes
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      jG13G3h =              jG13G3h   * limit_rate_2
      jG23G13h=              jG23G13h  * limit_rate_3
      jG33G23h(BoxNumberXY)= jG33G23h(BoxNumberXY)* limit_rate_4

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! set rates
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      call flux(BoxNumberXY, iiBen, ppG3h, ppG3h,  -loss )

      call flux(BoxNumberXY, iiBen, ppG13h, ppG13h, jG14G13h )
      call flux(BoxNumberXY, iiBen, ppG23h, ppG23h, jG24G23h)

      !limit flux from sediment is case of high losses 
      jcG3O3h=jG3O3h*min(limit_rate_0,limit_rate_1)*insw(jG3O3h) &
                                      +jG3O3h*insw(-jG3O3h)*limit_rate_0
      !jbot > 0.0 : input to pelagic, loss which cannot be effecuated in 
      ! benthos wil be subtracted in pelagic. Hence the rest of loss 
      ! must be subtractedfrom flux
      ! in case of too low/negative values in bottom force flux from pelagic 
      ! ( keep your fingers crossed!)
      jcG3O3h=jcG3O3h-(DONE-limit_rate_1)*loss
      call addbotflux(ANY,BoxNumberXY,iiBen,ppG3h,iiPel,ppO3h,jcG3O3h)

      call flux(BoxNumberXY,iiBen, ppG13h, ppG3h,  jG13G3h* insw(  jG13G3h) )
      call flux(BoxNumberXY,iiBen, ppG3h, ppG13h,- jG13G3h* insw( -jG13G3h) )

      call flux(BoxNumberXY,iiBen, ppG23h, ppG13h, jG23G13h* insw( jG23G13h))
      call flux(BoxNumberXY,iiBen, ppG13h, ppG23h,-jG23G13h* insw(-jG23G13h))

      call flux(BoxNumberXY, iiBen, ppG23h, ppG23h, jG33G23h(BoxNumberXY) )
   
      s=Source(iiBen,BoxNumberXY,ppG3h)
!     if ( abs(s) > 1.0D2 .or.G3h(BoxNumberXY)> 4000.0 ) then
!        write(LOGUNIT,*) 'Rate +State G3h:' ,s,G3h(BoxNumberXY)
!     endif
!     s=Source(iiBen,BoxNumberXY,ppG13h)
!     if ( abs(s) > 1.0D2 ) then
!        write(LOGUNIT,*) 'Rate+State  G13h:' ,s,G13h(BoxNumberXY)
!     endif
      s=Source(iiBen,BoxNumberXY,ppG23h)
      if ( abs(s)/(NZERO+G23h(BoxNumberXY)) >2.0  ) then
         write(LOGUNIT,*) 'Rate+State  G23h:' ,s,G23h(BoxNumberXY)
         write(LOGUNIT,*) 'jG23G13h:',jG23G13h
         write(LOGUNIT,*) 'limit_rate_3 +3a:',limit_rate_3,limit_rate_3a
         call set_warning_for_getm
      endif
      r=-jG3O3h*limit_rate_1+jG13G3h-loss*limit_rate_1
!     if (-r*LocalDelta> G3h(BoxNumberXY) &
!             .or.G3h(BoxNumberXY)<ZERO.or.s>max_change_per_step) then 
!        if (s>max_change_per_step) write(LOGUNIT,*)'--- TOO large rates-----'
!         write(LOGUNIT,*) 'BA:D1m,shiftD1m=',D1m(BoxNumberXY), &
!                                           shiftD1m(BoxNumberXY)
!         write(LOGUNIT,*) 'BA:D2m,shiftD2m=',D2m(BoxNumberXY), &
!                                           shiftD2m(BoxNumberXY)
!         write(LOGUNIT,*) 'BA:O3h',O3h_Ben(BoxNumberXY)
!         write(LOGUNIT,*) 'BA:ACae,G3h',ACae(BoxNumberXY),G3h(BoxNumberXY)
!         write(LOGUNIT,*) 'BA:sflux of G3h',Source(iiBen,BoxNumberXY,ppG3h),&
!                             Source(iiBen,BoxNumberXY,ppG3h)/G3h(BoxNumberXY)
!         write(LOGUNIT,*) 'BA:diff',ACae(BoxNumberXY)-O3h_Ben(BoxNumberXY)
!         write(LOGUNIT,*) 'BA:netchange from G3h',r,r*LocalDelta 
!         write(LOGUNIT,*) 'BA:loss',loss,limit_rate_1
!         write(LOGUNIT,*) 'BA:jG3O3h',jG3O3h
!         write(LOGUNIT,*) 'BA:jG13G3h',jG13G3h,limit_rate_2
!         write(LOGUNIT,*) 'BA:jG23G13h',jG23G13h
!         write(LOGUNIT,*) 'BA:jG33G23h',jG33G23h(BoxNumberXY)
!         write(LOGUNIT,*) 'BA:jG14G13h',jG14G13h
!         call set_warning_for_getm
!     endif
    end if
  end do
#endif

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
