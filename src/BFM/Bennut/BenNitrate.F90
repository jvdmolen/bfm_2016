#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenNitrate
!
! DESCRIPTION
!   Description of the diagenetic nitrate processes in the sediment
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
  subroutine BenNitrateDynamics
!
! !USES:

  ! For the following Benthic-states fluxes are defined: K3n, G4n
  ! The following Benthic-states are used (NOT in fluxes): D2m, D6m, D1m, K16r
  ! The following global scalar vars are used: &
  !    NO_BOXES_XY,  &
  !   BoxNumberXY, InitializeModel, LocalDelta
  ! The following Benthic 1-d global boxvars are modified : KNO3
  ! The following Benthic 1-d global boxvars got a value: M3n,jK3G4n
  ! The following Benthic 1-d global boxvars are used: rrATo, irrenh,ETW_Ben, &
  ! KNH4, N3n_Ben,Depth_Ben
  ! The following Benthic 1-d global boxpars  are used: p_poro
  ! The following 0-d global parameters are used: p_d_tot, p_q10diff, &
  ! p_qro, p_qon_dentri

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO,NZERO,DONE,dummy,LOGUNIT
  use mem,  ONLY: K3n,K4n,K13n,K23n, D1m, D2m, D6m, K16r,K26r, D2STATE
  use mem, ONLY: ppN3n,ppK3n, ppG4n, ruBPn3,sK4K3,dry_z, &
   NO_BOXES_XY, LocalDelta,max_change_per_step, &
     BoxNumberXY, InitializeModel, KNO3, M3n, jK3G4n,jK23G4n,jK23K13n,&
    jK4K3n, rrATo, irrenh, ETW_Ben, KNH4, N3n_Ben, iiPel,iiBen, flux,Depth_Ben
  use constants, ONLY: GET, LABDA_1,LABDA_2,&
    COEFFICIENT, LAYERS, LAYER1, DIFFUSION, FOR_ALL_LAYERS, POROSITY, &
    ADSORPTION, DOUBLE_DEFINE, ZERO_EXPONENTIAL_TERM, DEFINE, QUADRATIC_TERM, &
    LINEAR_TERM, CONSTANT_TERM, EXPONENTIAL_TERM, SET_CONTINUITY, FLAG, MASS, &
    SET_BOUNDARY, EQUATION, INPUT_TERM, PARAMETER, &
    SET_LAYER_INTEGRAL_UNTIL, LAYER2, ADD, DERIVATIVE, RFLUX,&
    INTEGRAL, MIN_VAL_EXPFUN,LAYER3, ANY,POSITIVE,NEGATIVE
  use mem_Param,ONLY: p_poro, p_q10diff, p_qro, p_qon_dentri,p_d_tot_2, &
                    p_dry_ben, p_clDxm
  use mem_Diffusion,ONLY: p_diff=>p_diff_N3
  use mem_BenNitrate
  use mem_BenAmmonium,only:p_slK4K3
  use LimitRates, ONLY:LimitChange
  use botflux,only:addbotflux

  use gotm_error_msg,only:set_warning_for_getm



  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following bennut functions are used:GetInfoFromSet, &
  ! InitializeSet, DefineSet, CompleteSet, CalculateSet, CalculateTau, &
  ! CalculateFromSet
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use bennut_interface, ONLY: GetInfoFromSet, InitializeSet, DefineSet, &
    CompleteSet, CalculateSet, CalculateTau, CalculateFromSet

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following global functions are used:eTq
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_interface,   ONLY: eTq

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following sesame functions are used:IntegralExp
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: IntegralExp,IntegralExpDist

!
!
! !AUTHORS
!   Original version by  P. Ruardij
!
! !REVISION_HISTORY
!   September 1999 by M. Vichi  Commented version
!
!
! COPYING
!   Copyright (C) 2006 P. Ruardij and M. Vichi
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
  real(RLEN)  :: Dxm
  real(RLEN)  :: sK13G4,sK23G4,sK4K3_l
  real(RLEN)  :: diff
  real(RLEN)  :: gamma_13,gamma_23,lambda,alpha
  real(RLEN)  :: a1,a2,a5
  real(RLEN)  :: cK3n
  real(RLEN)  :: Tau
  real(RLEN)  :: zuD1,jATK16o,jATK26o
  real(RLEN)  :: jK13G14n,jK23G24n,jK3N3n
  real(RLEN)  :: m3pwtm2
  real(RLEN)  :: r
  real(RLEN)  :: sx_any,rx_any,cx_any  ! resp. any specific rate, rate, concentration
  logical     :: isanoxic=.false.,thinoxiclayer
  real(RLEN)  :: limit_rate_0,limit_rate_1,limit_rate_r,limit_rate_n

  real(RLEN), external  :: GetDelta
  external              :: FixProportionCoeff
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  do BoxNumberXY=1,NO_BOXES_XY
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      m3pwtm2=( p_p+ DONE)*p_poro(BoxNumberXY)*D2m(BoxNumberXY)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate the pore-water average concentrations from the state variables
      ! (Diagnostic variables, not used in calculations)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      M3n(BoxNumberXY) = max(ZERO,K3n(BoxNumberXY)/m3pwtm2)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate coefficient for the e-folding distribution of the anoxic
      ! mineralization. D6.m is the average penetration depth for C-detritus
      ! Calculate anaerobic minralizzation in layer D1m-D2m
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      alpha  =   DONE/ D6m(BoxNumberXY)
      zuD1 = rrATo(BoxNumberXY)/ IntegralExpDist( - alpha, &
             p_d_tot_2- D1m(BoxNumberXY))
      jATK16o=zuD1*IntegralExpDist( - alpha, D2m(BoxNumberXY)- D1m(BoxNumberXY))
      jATK26o=rrATo(BoxNumberXY)-jATK16o


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Correction due to environmental regulating factors,
      ! diffusion coefficient: temperature and bioirrigation
      ! denitrification: temperature and coupling with anoxic mineralization
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      diff = p_diff* irrenh(BoxNumberXY)* eTq( ETW_Ben(BoxNumberXY), p_q10diff)

!     ! Do a raw estimation about how much nitrate is potentially inputted
!     ! into the system. This will happen only when nitrate in the sediment
!     ! is lower than in the overlying water.
!     ! This number is used to estimate the denitrifitation rate sK13G3n
!     r=max(ZERO,diff*p_poro(BoxNumberXY)* &
!           (N3n_Ben(BoxNumberXY)-M3n(BoxNumberXY))/ &
!           (Depth_Ben(BoxNumberXY)+D1m(BoxNumberXY))/2.0)

      sK13G4=rrATo(BoxNumberXY)*p_qon_dentri / &
                               (K13n(BoxNumberXY)+K23n(BoxNumberXY))
      sK23G4=sK13G4
      jK13G14n=sK13G4*K13n(BoxNumberXY)
      jK23G24n=sK23G4*K23n(BoxNumberXY)

      !limit denitrification at low concentration of K3n
      call LimitChange(POSITIVE,sK13G4,DONE,max_change_per_step ,limit_rate_n)

      !limit denitrification at low concentration of K16r
      cx_any=max(ZERO,(K16r(BoxNumberXY)/p_qro+jATK16o*LocalDelta)*p_qon_dentri)
      call LimitChange(ANY,jK13G14n,cx_any,max_change_per_step ,limit_rate_r)

      !take the minimum of the two limitations
      sK13G4=min(limit_rate_n,limit_rate_r)*sK13G4
      if (isnan(sK13G4).or.K16r(BoxNumberXY)< ZERO ) then
          write(LOGUNIT,*) ' BenNitrate '
          write(LOGUNIT,*) 'K16r(BoxNumberXY)=',K16r(BoxNumberXY)
          write(LOGUNIT,*) 'rrAto(BoxNumberXY)=',rrAto(BoxNumberXY)
          write(LOGUNIT,*) 'jK4K3n(BoxNumberXY)=',jK4K3n(BoxNumberXY)
          write(LOGUNIT,*) 'K3n(BoxNumberXY)=',K3n(BoxNumberXY)
      endif

      !limit denitrification at low concentration of K23n
      call LimitChange(POSITIVE,sK23G4,DONE,max_change_per_step ,limit_rate_n)

      !limit denitrification at low concentration of K26r
      cx_any=max(ZERO,K26r(BoxNumberXY)/p_qro *p_qon_dentri)
      call LimitChange(ANY,jK23G24n,cx_any,max_change_per_step ,limit_rate_r)

      !take the minimum of the two limitations
      sK23G4=min(limit_rate_n,limit_rate_r)*sK23G4

      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate coefficient for the exponential terms of the solution
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      !Limit first order rate
      sK13G4=min(p_shK3G4,sK13G4)
      sK23G4=min(p_shK3G4,sK23G4)
      !Maximize first order rate only in the calculatin of sqrt(diff/sk3G4)
      gamma_13  = max(MIN_VAL_EXPFUN,  sqrt( max(p_slK3G4, sK13G4)/ diff))
      gamma_23  = max(MIN_VAL_EXPFUN,  sqrt( max(p_slK3G4, sK23G4)/ diff))
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Get coefficients describing ammonium in the oxic layer :
      ! 1. lambda of the exponential curve
      ! 2. parameter of the nitrification term
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Get coefficients of all terms of equation valid for the
      ! first layer of ammonium (integration constants)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      sK4K3_l  =   GetInfoFromSet( KNH4(BoxNumberXY), GET, LABDA_2, 11)
      thinoxiclayer=(D2m(BoxNumberXY).lt.p_xthick_D2m.and. &
                        (sK4K3(BoxNumberXY)>p_slK4K3))
      if (thinoxiclayer ) then
        !use first order rate description of ammonium model
        lambda=GetInfoFromSet( KNH4(BoxNumberXY), GET, LABDA_1, 11)
        a1 = sK4K3_l* GetInfoFromSet( KNH4(BoxNumberXY), GET, COEFFICIENT,11)
        a2 = sK4K3_l* GetInfoFromSet( KNH4(BoxNumberXY), GET, COEFFICIENT,12)
        a5 = sK4K3_l* GetInfoFromSet( KNH4(BoxNumberXY), GET, COEFFICIENT,15) &
            -ruBPn3(BoxNumberXY)/D1m(BoxNumberXY)/p_poro(BoxNumberXY)
      else
        !use zero order rate description of ammonium model
        rx_any=max(p_slK4K3*K4n(BoxNumberXY),jK4K3n(BoxNumberXY))
        a5= (rx_any -ruBPn3(BoxNumberXY))/ D1m(BoxNumberXY)/p_poro(BoxNumberXY)
      endif
      if (isnan(a5)) then
        if (thinoxiclayer) &
        r = GetInfoFromSet( KNH4(BoxNumberXY), GET, COEFFICIENT, 15)
        write(LOGUNIT,*) 'a5 is Nan sK4K3,r,ruBPn3,',sK4K3_l,r,ruBPn3(BoxNumberXY)
      endif

      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Initialize the set of differential equations giving:
      ! - n. of layers
      ! - n. of coefficients
      ! - layer depths
      ! - environmental conditions (diffusion, porosity and adsorption coeff.)
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


      Dxm=(D1m(BoxNumberXY)+D2m(BoxNumberXY) ) * 0.5
      if ( isanoxic )  then
        KNO3(BoxNumberXY) = InitializeSet( KNO3(BoxNumberXY), 3, 5)
      elseif (thinoxiclayer) then
        KNO3(BoxNumberXY) = InitializeSet( KNO3(BoxNumberXY), 3, 8)
      else
        KNO3(BoxNumberXY) = InitializeSet( KNO3(BoxNumberXY), 3, 6)
      endif

      call DefineSet(KNO3(BoxNumberXY), LAYERS, LAYER1, LAYER2, &
        D1m(BoxNumberXY), Dxm)

      call DefineSet(KNO3(BoxNumberXY), DIFFUSION, FOR_ALL_LAYERS, 0, diff, &
        dummy)

      call DefineSet(KNO3(BoxNumberXY), POROSITY, FOR_ALL_LAYERS, 0, &
        p_poro(BoxNumberXY), dummy)

      call DefineSet( KNO3(BoxNumberXY), ADSORPTION, FOR_ALL_LAYERS, 0, p_p, &
        dummy)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Define coefficients for the steady-state solutions in each layer
      ! General solution of the equilibrium profile:
      ! 1st layer:
      ! N(z) = n11*exp(lambda*z) + n12*exp(-lambda*z) + n13*z^2 + n14*z + n15
      ! 2nd layer:
      ! N(z) = n21*exp(gamma*z) + n22*exp(-gamma*z)
      !    n22 = 0 (boundary condition)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        if (thinoxiclayer) then
          call DefineSet( KNO3(BoxNumberXY), DOUBLE_DEFINE, 11, &
            ZERO_EXPONENTIAL_TERM, lambda, sK4K3_l)

          call DefineSet( KNO3(BoxNumberXY), DOUBLE_DEFINE, 12, &
            ZERO_EXPONENTIAL_TERM, - lambda, sK4K3_l)
        endif

        call DefineSet( KNO3(BoxNumberXY), DEFINE, 13, QUADRATIC_TERM, dummy, &
          dummy)

        call DefineSet(KNO3(BoxNumberXY),DEFINE,14, LINEAR_TERM, dummy,dummy)
        call DefineSet(KNO3(BoxNumberXY),DEFINE,15, CONSTANT_TERM, dummy,dummy)

       call DefineSet(KNO3(BoxNumberXY),DOUBLE_DEFINE,21, EXPONENTIAL_TERM, - &
        gamma_13, sK13G4)

      call DefineSet(KNO3(BoxNumberXY), DOUBLE_DEFINE,22, EXPONENTIAL_TERM, &
        gamma_13, sK13G4)

      call DefineSet(KNO3(BoxNumberXY), DOUBLE_DEFINE,31, EXPONENTIAL_TERM, - &
        gamma_23, sK23G4)

      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Insert other boundary conditions and continuity between layers:
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      call CompleteSet( KNO3(BoxNumberXY), SET_CONTINUITY, FLAG, MASS, dummy)

      call CompleteSet( KNO3(BoxNumberXY), SET_BOUNDARY, LAYER1, &
          EQUATION, ZERO, value=N3n_Ben(BoxNumberXY))

      if (thinoxiclayer) then
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! a11 / (lambda * lambda * diff) = a12 / (lambda * lambda * diff)
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        call FixProportionCoeff(KNO3(BoxNumberXY),12,11,a2,a1)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! a11 / (lambda * lambda * diff) = a15 / (2 * diff)
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        call FixProportionCoeff(KNO3(BoxNumberXY),11,13,a1,a5)
      endif

      if ( InitializeModel== 0) then

        call CompleteSet(KNO3(BoxNumberXY),INPUT_TERM,13, PARAMETER, dummy, &
           value=a5)

        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Calculate for the above defined set of boundary conditions
        ! the steady-state profiles and return the vertically integrated
        ! concentration
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        cK3n = CalculateSet( KNO3(BoxNumberXY), SET_LAYER_INTEGRAL_UNTIL, &
          LAYER1, LAYER3, D2m(BoxNumberXY), ZERO)

        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Calculate the adaptation time to the steady-state profile
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        sx_any=(sK13G4*(D2m(BoxNumberXY)-D1m(BoxNumberXY)) &
            +sK4K3_l*D1m(BoxNumberXY))/D2m(BoxNumberXY)
        Tau  =   CalculateTau(  sx_any,  diff,  p_p,  D2m(BoxNumberXY))

        !=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Estimate the average value of K3n over the actual time step
        ! (transient value).
        ! This value depends on the adaptation time, the actual time step,
        ! the ''old'' value and the ''equilibrium value''
        !=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        cK3n = cK3n+( K3n(BoxNumberXY)- cK3n)* IntegralExp( - LocalDelta/ &
          Tau, DONE)

        !=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Derive the equations for the transient profiles, assuming the same
        ! solution as for the steady-state case and using cK3n as new &
        ! constraint. (r=dummy output)
        !=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        cx_any = CalculateSet( KNO3(BoxNumberXY), ADD, 0, 0, dummy, cK3n)

        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Start calculation of fluxes:
        !
        ! Denitrification:
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        r=D1m(BoxNumberXY); if ( isanoxic) r=ZERO
        jK13G14n = max(ZERO,CalculateFromSet( KNO3(BoxNumberXY), INTEGRAL, &
          RFLUX, r, D2m(BoxNumberXY)))* sK13G4

        if (jK13G14n > ZERO ) then
          cx_any= min(K3n(BoxNumberXY)+jK4K3n(BoxNumberXY)*LocalDelta, &
                    (K16r(BoxNumberXY)/p_qro+jATK16o*LocalDelta)* p_qon_dentri)
          call LimitChange(ANY,jK13G14n,cx_any,max_change_per_step, &
                                                        limit_rate_r)
          limit_rate_r=DONE
        elseif (jK13G14n==ZERO) then
          if (sK13G4 ==ZERO) then
            write(LOGUNIT,*), &
               "Warning denitrification=0 and water-sediment nitrate flux<0"
            call set_warning_for_getm
          endif
        else
          write(LOGUNIT,*) "Negative (or Nan) flux for nitrate"
  !       call PrintSet(KNO3(BoxNumberXY),"Negative (or Nan) flux for nitrate")
          write(LOGUNIT,'(''D1m='',F10.3)') D1m(BoxNumberXY)
          write(LOGUNIT,'(''D2m='',F10.3)') D2m(BoxNumberXY)
          write(LOGUNIT,'(''N3n_Ben/m3 ='',F10.3)') N3n_Ben(BoxNumberXY)
          write(LOGUNIT,'(''nitrate/m2 (K3n)='',F10.3)') K3n(BoxNumberXY)
          write(LOGUNIT,*) 'K16r=',K16r(BoxNumberXY)
          write(LOGUNIT,*) 'rrAto=',rrAto(BoxNumberXY)
          write(LOGUNIT,*) 'jK13G14n=',jK13G14n
          write(LOGUNIT,*) 'jK4K3n=',jK4K3n(BoxNumberXY)
          write(LOGUNIT,*) 'sK13G4,sK23G4=',sK13G4,sK23G4
          write(LOGUNIT,*) 'a5=',a5
          write(LOGUNIT,*) '-a5=', &
            ruBPn3(BoxNumberXY)/D1m(BoxNumberXY)/p_poro(BoxNumberXY)
  !       write(LOGUNIT,*) 'Mp3n=',Mp3n(BoxNumberXY)
  !       write(LOGUNIT,*) 'Mp4n=',Mp4n(BoxNumberXY)
           jK13G14n=ZERO
        endif

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Flux at the water/sediment interface
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        if ( dry_z(BoxNumberXY) <=p_clDxm .and. p_dry_ben) then
          jK3N3n = ZERO
        else
          jK3N3n = CalculateFromSet( KNO3(BoxNumberXY),  &
              DERIVATIVE, RFLUX, ZERO, dummy)
        endif

        if (sK4K3(BoxNumberXY) > p_slK4K3) then
          select case (p_flux_at_deep_end)  
             case(1); r=ZERO;
             case(2); r=min(ZERO,CalculateFromSet(KNO3(BoxNumberXY), &
                               DERIVATIVE,RFLUX,D2m(BoxNumberXY), ZERO))
             case(3); r=CalculateFromSet(KNO3(BoxNumberXY), &
                               DERIVATIVE,RFLUX,D2m(BoxNumberXY), ZERO)
          end select
          jK23K13n(BoxNumberXY)  = r
        else
          jK23K13n(BoxNumberXY)  = ZERO
        endif

        !Check on negative flux: flux from water into sediment
        call LimitChange(NEGATIVE,jK3N3n, &
            N3n_Ben(BoxNumberXY)*Depth_Ben(BoxNumberXY) ,&
            max_change_per_step,limit_rate_0)
        if (jK3N3n.gt.NZERO) then
          !In case of a positive flux to the water colomn,this flux is limited 
          !in case of too fast decrease of the concentration K3n in the sediment
          rx_any=jK3N3n+jK23K13n(BoxNumberXY)
          cx_any=K3n(BoxNumberXY)+ min(ZERO,(jK4K3n(BoxNumberXY) &
           -ruBPn3(BoxNumberXY) -jK13G14n+jK23K13n(BoxNumberXY)) *LocalDelta)
          call LimitChange(ANY,rx_any,cx_any ,max_change_per_step,limit_rate_1)
          jK3N3n=jK3N3n * min(limit_rate_0,limit_rate_1)
        else
          !In case of a negative flux to the water colomn,the denitrification 
          !flux is limited in case of too fast decrease of K3n in the sediment
          jK3N3n=jK3N3n * limit_rate_0
          rx_any=jK13G14n+jK23K13n(BoxNumberXY)
          cx_any=K3n(BoxNumberXY)+ min(ZERO,(jK4K3n(BoxNumberXY) &
              -ruBPn3(BoxNumberXY) -jK3N3n) *LocalDelta)
          call LimitChange(ANY,rx_any,cx_any ,max_change_per_step,limit_rate_1)
          jK13G14n=jK13G14n* limit_rate_1
          jK23K13n(BoxNumberXY)=jK23K13n(BoxNumberXY)* limit_rate_1
        endif

        call addbotflux(ANY,BoxNumberXY,iiBen,ppK3n,iiPel,ppN3n,jK3N3n)
        call flux(          BoxNumberXY,iiBen,ppK3n,      ppG4n,jK13G14n)

        !for output:
        jK3G4n(BoxNumberXY) = jK13G14n

 ! fixed rate ( calculated elsewhere)
        jK23G4n(BoxNumberXY)=CalculateFromSet( KNO3(BoxNumberXY), INTEGRAL, &
                      RFLUX, D2m(BoxNumberXY),p_d_tot_2 )* sK23G4 !*limit_flux
!       if (abs(jK23G4n(BoxNumberXY)).gt.1.0D+10) then
!         write(LOGUNIT,*) "----------------------------------------"
!         write(LOGUNIT,*) 'reBTn',reBTn(BoxNumberXY)
!         write(LOGUNIT,*) 'jK4K3n',jK4K3n(BoxNumberXY)
!         write(LOGUNIT,*) 'sK4K3',sK4K3(BoxNumberXY)
!         write(LOGUNIT,*) 'jK23G4n',jK23G4n(BoxNumberXY)
!         write(LOGUNIT,*) 'jK3G4n',jK3G4n(BoxNumberXY)
!         write(LOGUNIT,*) 'K3n',K3n(BoxNumberXY)
!         write(LOGUNIT,*) 'D1m',D1m(BoxNumberXY)
!         write(LOGUNIT,*) 'D2m',D2m(BoxNumberXY)
!         write(LOGUNIT,*) 'sK13G4',sK13G4
!         call PrintSet(KNO3(BoxNumberXY),'high jK23G4n')
!       endif
        cx_any=K26r(BoxNumberXY)/p_qro *p_qon_dentri
        call LimitChange(ANY,jK23G4n(BoxNumberXY),cx_any, max_change_per_step)

      else
        call CompleteSet( KNO3(BoxNumberXY),INPUT_TERM,13,PARAMETER, dummy, &
           value=a5)
        cx_any = CalculateSet( KNO3(BoxNumberXY), 0, 0,  0, dummy, dummy)
      end if

  end do

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
