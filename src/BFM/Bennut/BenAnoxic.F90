#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenAnoxic
!
! DESCRIPTION
!   Description of the anoxic diagenetic processes in the sediment
!       Details on the equations and the method used to calculate
!       the equilibrium and transient profiles can be found in
!       Ruardij et al., 1995. Neth. J. Sea Res. 33(3/4):453-483
!
!

!   This file is generated directly from OpenSesame model code, using a code
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine BenAnoxicDynamics
!
! !USES:

  ! For the following Benthic-states fluxes are defined: K16r, K6r
  ! The following Benthic-states are used (NOT in fluxes): D6m, D1m,D2m
  ! The following global scalar vars are used: &
  !    NO_BOXES_XY,  &
  !   BoxNumberXY, LocalDelta, InitializeModel
  ! The following Benthic 1-d global boxvars are modified:M6r,KRED,
  ! jG2K7o
  ! The following Benthic 1-d global boxvars are used: rrATo, rrBTo, irrenh, &
  ! ETW_Ben, KNO3, N6r_Ben
  ! The following Benthic 1-d global boxpars  are used: p_poro
  ! The following 0-d global parameters are used: p_d_tot, p_clDxm, &
  ! p_qro, p_q10diff, p_qon_dentri

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO,DONE,LOGUNIT,ZERO_KELVIN,dummy
  use mem,  ONLY: K26r, K16r, K6r, D6m,DH2m,DH3m, D1m, D2m, D2STATE
  use mem,ONLY: ppK26r,ppK16r,ppK6r, ppN6r, ppG2o, iiPel,iiBen,NO_BOXES_XY, &
    BoxNumberXY, LocalDelta, InitializeModel, M6r, KRED, KNO3,   &
    rrDTo, rrATo, reBTo,rrBTo,reK6o, jK36K26r, jK3G4n,jK23G4n, jG2K7o, &
    irrenh, ETW_Ben, N6r_Ben, G2_xavail_o, pHae, dry_z, &
    flux,max_change_per_step
  use constants, ONLY: GET, LABDA_1, LABDA_2, COEFFICIENT, SET_LAYER_INTEGRAL,&
    LAYERS, LAYER1, DIFFUSION, FOR_ALL_LAYERS, POROSITY, ADSORPTION, ADD,&
    DEFINE,DOUBLE_DEFINE,ZERO_EXPONENTIAL_TERM,DIST_EXPONENTIAL_TERM, &
    LINEAR_TERM, EXPONENTIAL_TERM, CONSTANT_TERM, SET_CONTINUITY, FLAG, MASS, &
    SET_BOUNDARY, EQUATION, INPUT_TERM, PARAMETER, SET_LAYER_INTEGRAL_UNTIL, &
    ADD, RFLUX, DERIVATIVE, INTEGRAL,LAYER2,  LAYER3, LAYER4,LAYER5, &
    QUADRATIC_TERM,ANY,POSITIVE,NEGATIVE,STANDARD
  use mem_Param,ONLY:p_poro,p_d_tot,p_d_tot_2,p_clDxm,p_qro,p_q10diff, &
                   p_qon_dentri,p_dry_ben,combine_anabac
  use mem_Diffusion,ONLY:p_diff=>p_diff_N6
  use mem_BenAnoxic
  use LimitRates, ONLY:LimitChange
  use mem_BenBac,ONLY: sw_an
  USE BFM_ERROR_MSG, ONLY: set_warning_for_getm
#ifdef INCLUDE_BENCO2
  use mem,ONLY:pHae,G3h
  use mem_CO2,ONLY:p_qhATo
#endif
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following bennut functions are used:GetInfoFromSet, &
  ! InitializeSet, DefineSet, CompleteSet, CalculateSet, CalculateTau, &
  ! CalculateFromSet
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use bennut_interface, ONLY: GetInfoFromSet, InitializeSet, DefineSet, &
    CompleteSet, CalculateSet, CalculateTau, CalculateFromSet
  use botflux,onlY:addbotflux

! use mem,  ONLY: Output2d_1,Output2d_2,Output2d_3, Output2d_4

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following global functions are used:eTq
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_interface,   ONLY: eTq,CalcPelMassInM2

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following sesame functions are used:IntegralExp, insw
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: IntegralExpDist,ExpDist,SetExpDist,IntegralExp,insw


!
!
! !AUTHORS
!   Original version by  P. Ruardij
!
! !REVISION_HISTORY
!   September 1999 by M. Vichi Commented version
!
!
! COPYING
!
!   Copyright (C) 2006 P. Ruardij & M.VIchi
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
  real(RLEN),dimension(NO_BOXES_XY)  :: cN6m2r
  real(RLEN)  :: alpha,beta
  real(RLEN)  :: zuD1,zuBT
  real(RLEN)  :: diff
  real(RLEN)  :: gamma
  real(RLEN)  :: lambda
  real(RLEN)  :: sK13G4
  real(RLEN)  :: n21,n22,n31,a15
  real(RLEN)  :: Tau
  real(RLEN)  :: cK6r
  real(RLEN)  :: jBTK6r,jATK6r,jDTK6r
  real(RLEN)  :: jK6BTr , jK6G4r , jK6N6r , jK26G4r
  real(RLEN)  :: limit_rate_0,limit_rate_1
  real(RLEN)  :: sOS
  real(RLEN)  :: r
  real(RLEN)  :: rx_any,sx_any,cx_any
  real(RLEN)  :: Dxm,Dym,M16r
  real(RLEN)  :: flK6BTr
  real(RLEN)  :: O2o,G2_xcorr_o
  real(RLEN)  :: limit_oxygen,limit_G3h

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  cN6m2r=CalcPelMassInM2(ppN6r)
  do BoxNumberXY=1,NO_BOXES_XY

      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate the pore-water average concentration from the state variables
      ! (Diagnostic variables, not used in calculations)
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      if ( K6r(BoxNumberXY)< ZERO ) then
         K6r(BoxNumberXY)=ZERO
         write(LOGUNIT,*) 'K6r< 0.0, reset on zero'
         call set_warning_for_getm
      endif
      if ( K16r(BoxNumberXY)< ZERO ) then
         K16r(BoxNumberXY)=ZERO
         write(LOGUNIT,*) 'K16r< 0.0, reset on zero'
         call set_warning_for_getm
      endif

      M6r(BoxNumberXY) = K6r(BoxNumberXY)/ p_poro(BoxNumberXY)/&
               ( p_p+ DONE)/(  D1m(BoxNumberXY))
      M16r =             K16r(BoxNumberXY)/ p_poro(BoxNumberXY)/&
               ( p_p+ DONE)/(  D2m(BoxNumberXY)-D1m(BoxNumberXY))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate coefficient for the e-folding distribution of the anoxic
      ! mineralization. D6.m is the average penetration depth for C-detritus
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      if ( combine_anabac) then
        alpha  =   SetExpDist(D6m(BoxNumberXY),Dlm=p_clDxm,Dmm=p_d_tot)
        beta=alpha
      else
        alpha  =   SetExpDist(DH2m(BoxNumberXY),Dlm=p_clDxm,Dmm=p_d_tot)
        beta  =    SetExpDist(DH3m(BoxNumberXY),Dlm=p_clDxm,Dmm=p_d_tot)
      endif


      if ( InitializeModel == 0 ) then
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Convert anoxic mineralization (mmol S/m2/d)
        ! This rate is already assigned to the dynamical equation for K6.r
        ! in BenBacDyanmics for H2:
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

          jATK6r  =   p_qro* rrATo(BoxNumberXY)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Recalculate Mineralization m2 --> m3 porewater
        ! Anoxic mineralization at D1.m
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        select case  (combine_anabac)
          case(.FALSE.)
             jDTK6r  =   p_qro* rrDTo(BoxNumberXY)
             zuD1 = jDTK6r/ p_poro(BoxNumberXY)/ IntegralExpDist( - alpha, &
             D2m(BoxNumberXY)- D1m(BoxNumberXY))
          case(.TRUE.)
             zuD1 = jATK6r/ p_poro(BoxNumberXY)/ IntegralExpDist( - alpha, &
             p_d_tot_2- D1m(BoxNumberXY))
        end select
      else
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! In case of no info about anoxi mineralization at the start
        ! Reconstruct using detritus distribution alpha and oxic mineralization
        ! the anoxic mineralization at the upperside of the denitrification
        ! layer
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        jBTK6r= p_qro * rrBTo(BoxNumberXY)
        zuD1 = 0.01D+00 * jBTK6r / p_poro(BoxNumberXY)/ IntegralExpDist( &
               -alpha, D1m(BoxNumberXY)) *ExpDist(-alpha , D1m(BoxNumberXY))
      endif

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Correction due to environmental regulating factors,
      ! diffusion coefficient: temperature and bioirrigation
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      diff = p_diff* irrenh(BoxNumberXY)* eTq( ETW_Ben(BoxNumberXY), p_q10diff)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate coefficient for the exponential terms of the solution
      ! limitation of reoxidation of first order process at low anoxic
      ! mineralization
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      G2_xcorr_o=G2_xavail_o(BoxNumberXY) &
                     +(reBTo(BoxNumberXY)-rrBTo(BoxNumberXY))*LocalDelta
      if ( p_sOS.gt.ZERO) then
          sOS=max(0.001D+00,p_sOS * eTq( ETW_Ben(BoxNumberXY),p_q10)  &
            *rrATo(BoxNumberXY)/(0.01D+00+rrATo(BoxNumberXY)))
      else
          cx_any=7.0
#ifdef INCLUDE_BENCO2
          if ( pHae(BoxNumberXY).gt.ZERO)  cx_any=pHae(BoxNumberXY)
#endif
          !equation of Millero Environ Sci Technol. 1987,21,439-443
          ! chosen for constant ion strength of 0.7 (seawater) :
          ! 0.3681=0.44*0.7**0.5
          sOS=10.0D+00**(11.78D+00-3.0D3/ &
          (-ZERO_KELVIN+ETW_Ben(BoxNumberXY))+0.16D+00*cx_any+0.3681D+00)*24.D-6
          O2o=G2_xcorr_o/D1m(BoxNumberXY)/p_poro(BoxNumberXY)
          sOs=sOs*O2o
      endif

      limit_G3h=DONE
      flK6BTr=sOS*M6r(BoxNumberXY)         ! unit: mmol S2-/m3 porewater
      jK6BTr=sOS*K6r(BoxNumberXY)     ! unit: mmol S2-/m2
#ifdef INCLUDE_BENCO2
      r=jK6BTr/p_qro*p_qhATo
      call LimitChange(POSITIVE,r,G3h(BoxNumberXY) , &
                                max_change_per_step,limit_G3h)
#endif
      call LimitChange(POSITIVE,jK6BTr,p_qro*G2_xcorr_o , &
                                max_change_per_step,limit_oxygen)
      call LimitChange(POSITIVE,jK6BTr,K6r(BoxNumberXY)       , &
                                           max_change_per_step,r)
      flK6BTr=flK6BTr*min(r,limit_oxygen,limit_G3h)
      sOS=max(0.001D+00,sOS*min(r,limit_oxygen,limit_G3h))
      if  (sw_set>=1) gamma= sqrt(sOS/diff)

      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Get coefficients describing Nitrate in anoxic layer :
      ! 1. lambda of the exponential curve, and the denitrification rate
      ! 2. parameter of the denitrification term (integration constant)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      lambda = abs(GetInfoFromSet( KNO3(BoxNumberXY), GET, LABDA_1,21))
      sK13G4 = GetInfoFromSet( KNO3(BoxNumberXY), GET, LABDA_2, 21) &
                                               * p_qro/ p_qon_dentri
      n21 = - GetInfoFromSet( KNO3(BoxNumberXY), GET, COEFFICIENT, 21)* sK13G4
      n22 = - GetInfoFromSet( KNO3(BoxNumberXY), GET, COEFFICIENT, 22)* sK13G4
      n31 = - GetInfoFromSet( KNO3(BoxNumberXY), GET, COEFFICIENT, 31)* sK13G4

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Initialize the set of differential equations giving:
      ! - n. of layers
      ! - n. of coefficients
      ! - layers depths
      ! - environmental conditions (diffusion, porosity and adsorption coeff.)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      Dxm=(D1m(BoxNumberXY)+D2m(BoxNumberXY)) *0.5
      Dym=(D2m(BoxNumberXY)+p_d_tot) *0.5

      select case (p_flux_at_deep_end)
        case (1)    ; KRED(BoxNumberXY)= InitializeSet( KRED(BoxNumberXY),5,16)
        case default; KRED(BoxNumberXY)= InitializeSet( KRED(BoxNumberXY),5,17)
      end select

     call DefineSet(KRED(BoxNumberXY),LAYERS,LAYER1,LAYER2, &
                                                        D1m(BoxNumberXY),Dxm)
     call DefineSet(KRED(BoxNumberXY),LAYERS,LAYER3,LAYER4, &
                                                        D2m(BoxNumberXY),Dym)

     call DefineSet(KRED(BoxNumberXY),DIFFUSION, FOR_ALL_LAYERS,0, diff, dummy)
     call DefineSet(KRED(BoxNumberXY), POROSITY, FOR_ALL_LAYERS, 0, &
        p_poro(BoxNumberXY), dummy)

     call DefineSet(KRED(BoxNumberXY),ADSORPTION, FOR_ALL_LAYERS,0, p_p, dummy)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Define coefficients for the steady-state solutions in each layer
      ! General solution of the equilibrium profile:
      ! 1st layer:
      ! R(z) = r11*exp(gamma*z) + r12*exp(-gamma*z)
      ! 2nd layer:
      ! R(z) = r21*exp[-alpha*(z-D1.m)] + r22*exp(-lambda*z) +r23*z^2 + r24*z
      !                                                               + r25
      !    r23, r24 = 0 (boundary condition)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      select case  (sw_set )
      case (0)
         ! Assume in oxic layer 0-order process
        call DefineSet(KRED(BoxNumberXY),DEFINE,13, QUADRATIC_TERM,dummy,dummy)
        call DefineSet(KRED(BoxNumberXY),DEFINE,14, LINEAR_TERM, dummy, dummy)
        call DefineSet(KRED(BoxNumberXY),DEFINE,15, CONSTANT_TERM, dummy,dummy)
      case (1,2)
        ! Assume in oxic layer 1-order process
        call DefineSet(KRED(BoxNumberXY),DOUBLE_DEFINE,11,EXPONENTIAL_TERM, &
                                                                  -gamma,sOS)
        call DefineSet(KRED(BoxNumberXY),DOUBLE_DEFINE,12,EXPONENTIAL_TERM, &
                                                                   gamma,sOS)
        call DefineSet(KRED(BoxNumberXY),DOUBLE_DEFINE,15,CONSTANT_TERM, &
                                                                dummy, dummy)
      end select

      call DefineSet(KRED(BoxNumberXY), DEFINE, 21, DIST_EXPONENTIAL_TERM, &
                                                                -alpha, dummy)
      call DefineSet(KRED(BoxNumberXY), DOUBLE_DEFINE, 22, &
                                       ZERO_EXPONENTIAL_TERM, -lambda, sK13G4)
      call DefineSet(KRED(BoxNumberXY),DOUBLE_DEFINE, 23, &
                                        ZERO_EXPONENTIAL_TERM, +lambda, sK13G4)
      call DefineSet(KRED(BoxNumberXY),DEFINE, 24, LINEAR_TERM, dummy, dummy)
      call DefineSet(KRED(BoxNumberXY),DEFINE, 25, CONSTANT_TERM, dummy, dummy)

      call DefineSet(KRED(BoxNumberXY),DEFINE, 31, DIST_EXPONENTIAL_TERM,  &
                                                                - alpha, dummy)
      call DefineSet(KRED(BoxNumberXY), DOUBLE_DEFINE, 32, &
                                        ZERO_EXPONENTIAL_TERM, -lambda, sK13G4)
      call DefineSet(KRED(BoxNumberXY),DEFINE, 34, LINEAR_TERM, dummy, dummy)
      call DefineSet(KRED(BoxNumberXY),DEFINE, 35, CONSTANT_TERM, dummy, dummy)

      call DefineSet(KRED(BoxNumberXY),DEFINE, 41, DIST_EXPONENTIAL_TERM,- &
                                                                beta, dummy)
      call DefineSet(KRED(BoxNumberXY),DEFINE, 44, LINEAR_TERM, dummy, dummy)
      call DefineSet(KRED(BoxNumberXY),DEFINE, 45, CONSTANT_TERM, dummy, dummy)

      if (p_flux_at_deep_end > 1 ) &
      call DefineSet(KRED(BoxNumberXY), DEFINE, 51, DIST_EXPONENTIAL_TERM, &
                                                                -beta, dummy)
      call DefineSet(KRED(BoxNumberXY),DEFINE, 55, CONSTANT_TERM, dummy, dummy)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Insert other boundary conditions and continuity between layers:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      !1-6:
      call CompleteSet( KRED(BoxNumberXY), SET_CONTINUITY, FLAG, MASS, dummy)

      !7:
      call CompleteSet( KRED(BoxNumberXY), SET_BOUNDARY, LAYER1, &
          EQUATION, ZERO, value=N6r_Ben(BoxNumberXY))

      !8-9
      if (p_flux_at_deep_end > 1 ) &
      call FixProportionCoeff(KRED(BoxNumberXY),41,51, &
                                   DONE,ExpDist(-beta , Dym-D2m(BoxNumberXY)))
      call FixProportionCoeff(KRED(BoxNumberXY),22,23,n21,n22)
      call CompleteSet( KRED(BoxNumberXY), INPUT_TERM, 22, &
                             PARAMETER, dummy, value=n21)

      select case (InitializeModel)
        !10 -11
        case(0)
          call CompleteSet(KRED(BoxNumberXY), SET_LAYER_INTEGRAL, &
                            LAYER2, LAYER3,dummy ,     value=K16r(BoxNumberXY))
          call CompleteSet(KRED(BoxNumberXY), SET_LAYER_INTEGRAL_UNTIL, &
                            LAYER4, LAYER5, p_d_tot_2, value=K26r(BoxNumberXY))
        case(1)
          call CompleteSet( KRED(BoxNumberXY), INPUT_TERM, 31, PARAMETER, &
            dummy, value=zuD1* ExpDist(-alpha,Dxm-D1m(BoxNumberXY)))
          call CompleteSet( KRED(BoxNumberXY), INPUT_TERM, 41, PARAMETER, &
            dummy,value=zuD1*ExpDist(-alpha, D2m(BoxNumberXY)-D1m(BoxNumberXY)))
      end select

      !12
      call FixProportionCoeff(KRED(BoxNumberXY),21,31, &
                                   DONE,ExpDist(-alpha , Dxm-D1m(BoxNumberXY)))
      call CompleteSet( KRED(BoxNumberXY), INPUT_TERM, 21, PARAMETER, &
                                   dummy, value=max(1.0D-4,zuD1))

     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     ! Calculate for the above defined set of boundary conditions
     ! the steady-state profiles and return the vertically integrated
     ! concentration.
     !
     ! Technical improvements: in case of utlimate low mineralization rates and
     ! a nearly empty reduction equivalent pool. There is a chance the
     ! estimated equilibrium value is negative. Therfore cK6r is limited to &
     ! values >=0
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

     if ( InitializeModel== 0) then
       !13
       select case (sw_set )
        case(1)
         r=diff*p_poro(BoxNumberXY)
         call CompleteSet( KRED(BoxNumberXY),SET_BOUNDARY, &
                        LAYER1,DERIVATIVE, ZERO,mfac=-r)
         call CompleteSet( KRED(BoxNumberXY),ADD, &
                        LAYER1,DERIVATIVE, D1m(BoxNumberXY), mfac=+r)
         rx_any=flK6BTr*D1m(BoxNumberXY) !unit S2-/(m*m3pw) << S20/M3pw
         call CompleteSet( KRED(BoxNumberXY),ADD+SET_LAYER_INTEGRAL, &
                         LAYER1, LAYER1, dummy, mfac=sOs,value=rx_any)
        case(2)
          zuBT  = p_qro*reK6o(BoxNumberXY)/p_poro(BoxNumberXY)/ D1m(BoxNumberXY)
          zuBT=max(zuBT,1.079D-6)
          ! Calculate coefficient of the zero order term
          a15  =   zuBT/sOs 
          call CompleteSet( KRED(BoxNumberXY), INPUT_TERM, 15, STANDARD, &
            dummy,value=a15)
      end select

      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate the adaptation time to the steady-state profile
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      cK6r = max( ZERO, CalculateSet( KRED(BoxNumberXY), &
                 SET_LAYER_INTEGRAL, LAYER1, LAYER1, ZERO , ZERO))

      Tau  =   CalculateTau(  sOS,  diff,  p_p,  D1m(BoxNumberXY))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Estimate the average value of K6r over the actual time step
      ! (transient value).
      ! This value depends on the adaptation time, the actual time step,
      ! the ''old'' value and the ''equilibrium value''
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      cK6r= cK6r+( K6r(BoxNumberXY)- cK6r)* IntegralExp(-LocalDelta/Tau, DONE)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Derive the equations for the transient profiles, assuming the same
      ! solution as for the steady-state case and using cK6r as new &
      ! constraint.
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      r = CalculateSet( KRED(BoxNumberXY), ADD, 0, 0, dummy, cK6r)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! flux at  sediment/water interface
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      if ( dry_z(BoxNumberXY) <=p_clDxm .and. p_dry_ben) then
        jK6N6r= ZERO
      else
        jK6N6r= CalculateFromSet(KRED(BoxNumberXY),DERIVATIVE,RFLUX,ZERO,dummy)
      endif

      if (sw_set <=1) then
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! processes affecting K6r
      ! To overcome with problems with negative K6r, only the difference between
      ! the actual value and the estimated value on next time step are taken
      ! in account for the calculation of the reoxidation. This means in
      ! practice the value is kept higher than estimated.
      ! Of course negative reoxidation has to be neglected.
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        jK6BTr=max(ZERO,sOS*(K6r(BoxNumberXY) &
           -CalculateFromSet( KRED(BoxNumberXY),  &
                                      INTEGRAL, RFLUX,ZERO, D1m(BoxNumberXY))))
      else
        jK6BTr=max(ZERO,sOS*CalculateFromSet( KRED(BoxNumberXY),  &
                                    INTEGRAL, RFLUX,ZERO, D1m(BoxNumberXY)))
      endif
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! processes affecting K16r
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! loss flux due to denitrification (fixed rate, calculate elsewhere)
      jK6G4r= p_qro/ p_qon_dentri* jK3G4n(BoxNumberXY)

      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      !  processes affecting K26r
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! fixed rate ( calculated elsewhere)
      jK26G4r=jK23G4n(BoxNumberXY) * p_qro/ p_qon_dentri

      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! flux at underside
      ! p_flux_at_deep_end 1=no flux,2= only fluxes_downwards(sink),3=full_flux
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      select case (p_flux_at_deep_end)
        case (1); rx_any = ZERO
        case (2); rx_any = min(ZERO,CalculateFromSet(KRED(BoxNumberXY), &
                                         DERIVATIVE, RFLUX, p_d_tot_2, dummy))
        case (3); rx_any = CalculateFromSet(KRED(BoxNumberXY), &
                                          DERIVATIVE, RFLUX, p_d_tot_2, dummy)
      end select
      jK36K26r(BoxNumberXY)=rx_any

      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! set limits
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! limit for change in K6r
      sx_any=-jK6N6r*insw(jK6N6r)
      cx_any=K6r(BoxNumberXY)+(-jK6BTr)*LocalDelta
      call LimitChange(NEGATIVE,sx_any,cx_any,max_change_per_step,limit_rate_1)

      ! limit for change in N6r
      call LimitChange(NEGATIVE,jK6N6r,cN6m2r(BoxNumberXY), &
                            max_change_per_step,limit_rate_0)

      jK6N6r=min(limit_rate_0,limit_rate_1)*jK6N6r
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! set fluxes
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


      call flux(BoxNumberXY, iiBen, ppK6r, ppK6r,  -jK6BTr )
      jG2K7o(BoxNumberXY)  =   jK6BTr/ p_qro
      call flux(BoxNumberXY, iiBen, ppG2o, ppG2o,-jK6BTr/ p_qro )

      call addbotflux(ANY,BoxNumberXY,iiBen,ppK6r,iiPel,ppN6r, jK6N6r)
      call flux(BoxNumberXY,          iiBen,ppK16r,     ppK16r,-jK6G4r )

      !Denitrification in anxoic layer..
      call flux(BoxNumberXY, iiBen, ppK26r, ppK26r, -jK26G4r )

    else
      call CompleteSet( KRED(BoxNumberXY), INPUT_TERM, 21, PARAMETER, &
                dummy, value=zuD1)

       r = CalculateSet( KRED(BoxNumberXY), 0, 0,  0, dummy, dummy)

       jK6BTr = sOS* CalculateFromSet( KRED(BoxNumberXY), INTEGRAL, &
                                         RFLUX, ZERO, D1m(BoxNumberXY))
      jG2K7o(BoxNumberXY)  =   jK6BTr/ p_qro
    endif


  end do

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
