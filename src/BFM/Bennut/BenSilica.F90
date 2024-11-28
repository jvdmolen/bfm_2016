#include "DEBUG.h"

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenSilica
!
! DESCRIPTION
!   Description of the diagenitic processes in the sediment
!       Details on the equations and the method used to calculate
!       the equilibrium and transient profiles can be found in
!       Ruardij et al., 1995. Neth. J. Sea Res. 33(3/4):453-483
!
!
!     mode =1 : on basis of environment,biogene silicate and silicate states
!         actual state is calculated and transient profile in case
!         of no changes in the layers
!     mode =2 : all fluxes are calculated  and transitprofiles with
!              desortpion/adsorption in the layers which are shifted.
!     mode=3: mode==1 + mode== 2
!
!
!
! !INTERFACE
  subroutine BenSilicaDynamics(mode)
!
! !USES:

  ! For the following Benthic-states fluxes are defined: K5s, Q6s, D9m
  ! The following Benthic-states are used (NOT in fluxes): D1m, D2m
  ! The following global scalar vars are used: &
  !    NO_BOXES_XY,   &
  !  BoxNumberXY, InitializeModel, LocalDelta
  ! The following Benthic 1-d global boxvars are modified:M5s,KSiO3,KSiO3eq, &
  ! jK25K15s
  ! The following Benthic 1-d global boxvars are used: irrenh, ETW_Ben, &
  ! N5s_Ben, shiftD2m
  ! The following Benthic 1-d global boxpars  are used: p_poro
  ! The following 0-d global parameters are used: p_clD1D2m, p_q10diff, &
  ! p_d_tot

  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO,NZERO,DONE,dummy
  use mem,  ONLY: K5s, K15s,Q6c,Q6s, D6m, D9m, D1m, D2m, D2STATE
  use mem, ONLY:  iiPel,iiBen,ppN5s,ppK5s,ppK15s, ppQ6s, ppD9m, &
    flux,BoxNumberXY,NO_BOXES_XY, LocalDelta,max_change_per_step, &
    InitializeModel, M5s,M15s, KSiO3, Depth_Ben,ruBTs,reBTs,reATs, &
    KSiO3eq,jK25K15s,irrenh,ETW_Ben,N5s_Ben, shiftD2m,dry_z
  use constants, ONLY: LAYER1,LAYER2,LAYER3,LAYER4, SHIFT, STANDARD, &
    EQUATION, SET_LAYER_INTEGRAL, ANY,POSITIVE,NEGATIVE, &
    LAYER2, ADD, INTEGRAL, DERIVATIVE, RFLUX, MASS,EXPONENTIAL_INTEGRAL
  use mem_Param,  ONLY: p_poro, p_clD1D2m, p_q10diff,p_dry_ben, &
                        p_d_tot,p_d_tot_2,p_pK5_ae,p_clDxm
  use mem_Diffusion,ONLY: p_diff=>p_diff_N5
  use botflux,onlY:addbotflux
  use mem_BenSilica
  use mem_BenthicNutrient3, ONLY:sw_special_shift,p_pAn2Ni
  use LimitRates, ONLY:LimitShift,LimitChange
  USE BFM_ERROR_MSG, ONLY: set_warning_for_getm,BFM_ERROR
! use mem,  ONLY: Output2d_1,Output2d_2,Output2d_3, Output2d_4


  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following bennut functions are used:InitializeSet, &
  ! DefineSet, CompleteSet, CalculateSet, CalculateTau, CopySet, &
  ! CalculateFromSet
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use bennut_interface,ONLY:CalculateSet,CalculateTau,CopySet,CalculateFromSet


  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following global functions are used:eTq
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_interface,   ONLY: eTq


  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following sesame functions are used:IntegralExp, insw
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: IntegralExpDist,ExpDist,SetExpDist,IntegralExp, insw,PartQ
!
!
! !AUTHORS
!   Original version by  P. Ruardij
!
! !REVISION_HISTORY
!       September 1999 by M. Vichi     Commented version
!
! COPYING
!
!   Copyright (C) 2006 P. Ruardij & M. Vichi
!   (rua@nioz.nl, vichi@bo.ingv.it)
!
!   This program is free softwaee; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation;
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!EOP
!------------------------------------------------------------------------!
!BOC
!
!
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE
  integer,intent(IN)                    :: mode

  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer     :: i,WithSpecialLayer
  logical     :: dry
  real(RLEN)  :: r,s,h
  real(RLEN)  :: lim_1,lim_2,lim_3,lim_4
  real(RLEN)  :: cD1m,cD2m
  real(RLEN)  :: Dnewm
  real(RLEN)  :: chM5s
  real(RLEN)  :: cM5s
  real(RLEN)  :: Tau
  real(RLEN)  :: alpha,lambda
  real(RLEN)  :: diff
  real(RLEN)  :: M5b0
  real(RLEN)  :: M5b_0_d1
  real(RLEN)  :: M5bD1
  real(RLEN)  :: suD1,suD2,suDn
  real(RLEN)  :: shiftmass
  real(RLEN)  :: jQ6K5s    !minus total flux of Q6 to K5
  real(RLEN)  :: jQ6K5s_1  !minus       flux of Q6 to K5 in first layer
  real(RLEN)  :: jQ6K15s   !minus total flux of Q6 to K5
  real(RLEN)  :: jK5N5s, jK15K5s
  real(RLEN)  :: smQ6
  real(RLEN)  :: sM5s
  real(RLEN)  :: M0
  real(RLEN)  :: R5s,RN5s_Ben,K0s,R0s
  real(RLEN)  :: R15s
  real(RLEN)  :: tom2_K5s,tom2_ae,tom2_an
  real(RLEN)  :: rus
  real(RLEN)  :: mM5s       ! =chM5s- M5s
  real(RLEN)  :: pShift,mShift,lShift,dnm

  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  do BoxNumberXY=1,NO_BOXES_XY

      p_p_ae=p_pK5_ae(BoxNumberXY)
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! nutrient uptake due to primary production-nutrient release by BenPhyto
      ! (only in winter)
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      rus=ruBTs(BoxNumberXY)-reBTs(BoxNumberXY)
      if (abs(rus).gt.12.0) then
        write(LOGUNIT,*) 'BenSilica: rus',rus,ruBTs,reBTs
      endif
      if ( mode.eq.1) rus=1.0D-05

      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Correction due to environmental regulating factors,
      ! saturation value: temperature
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      chM5s = p_chM5s+ p_cvM5s*( eTq( ETW_Ben(BoxNumberXY), p_q10)- DONE)

      cD2m  =   min(  max(  D2m(BoxNumberXY),   p_clD2m),  p_chD2m) 
      cD1m  =   min(  max(  D1m(BoxNumberXY),   p_clD1m),  cD2m- p_clD1D2m)
      select case (sw_special_shift)
        case (0); dnm=     cD2m
        case (1); dnm=    (DONE-p_pAn2Ni) * cD2m +p_pAn2Ni*p_d_tot
      end select

      tom2_ae=( DONE+p_p_ae) * p_poro(BoxNumberXY)
      tom2_an=( DONE+p_p_an) * p_poro(BoxNumberXY)
      tom2_K5s= (tom2_ae*cD2m +tom2_an*(dnm-cD2m))/dnm

      if ( InitializeModel== 1) then
         K5s(BoxNumberXY)= 0.1* (chM5s)*cD2m*tom2_ae
         K15s(BoxNumberXY)=0.9*chM5s*(p_d_tot_2-cD2m)*tom2_an
      endif

      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Correction due to environmental regulating factors,
      ! diffusion coefficient: temperature and bioirrigation
      ! dissolution rate: temperature
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      diff = p_diff* irrenh(BoxNumberXY)* eTq( ETW_Ben(BoxNumberXY), p_q10diff)
      smQ6  =   p_smQ6* eTq(  ETW_Ben(BoxNumberXY),  p_q10diff)

      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate coefficient for the e-folding distribution of the anoxic
      ! mineralization. D9.m is the average penetration depth for biogenic Si
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      alpha =SetExpDist(D9m(BoxNumberXY),Dlm=p_clD9m,Dmm=p_chD2m)

      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate total biogenic silica from m2 --> m3 porewater
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      M5b0=Q6s(BoxNumberXY)/p_poro(BoxNumberXY)/IntegralExpDist(-alpha, p_d_tot_2)

      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Average content of Biogenic silica in the oxic layer
      ! and calculation of the zero-order dissolution term
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      M5b_0_d1  =   M5b0* IntegralExpDist( - alpha,  cD1m)/ cD1m

      if ( InitializeModel ==0 ) then
          mM5s=max(NZERO,chM5s-K5s(BoxNumberXY)/(dnm *tom2_K5s))
      else
          mM5s=0.5 * chM5s
      endif

      R5s =chM5s*dnm*tom2_K5s -K5s(BoxNumberXY)
      R15s=chM5s*tom2_an*(p_d_tot_2-dnm)-K15s(BoxNumberXY)

      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Check if the amount of biogenic silicate is so large that within one
      ! step the dissolution will reach the equilibrium concentration and
      ! no feed back is possible that the profile adapt its self to a new
      ! situation. This is done by limiting the dissolution rates....
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      r=DONE
      M5bD1  =  M5b0* ExpDist( -alpha, cD1m)
      suD1  =   smQ6* M5bD1/ chM5s 
      suD2  =   smQ6* M5bD1/ chM5s*ExpDist(-alpha ,cD2m-cD1m)
      suDn  =   smQ6* M5bD1/ chM5s*ExpDist(-alpha ,dnm-cD1m)

      sM5s= smQ6*M5b_0_d1/chM5s
      M0=smQ6*mM5s/chM5s* M5b_0_d1

      lambda= sqrt(max(0.001D+00,sM5s)/diff)

      if (isnan(lambda)) then
         write(LOGUNIT,*)'Lambda is NAN'
         write(LOGUNIT,*)'mmM5s', mM5s,' smQ6', smQ6
         write(LOGUNIT,*)'chM5s', chM5s
      endif

      if ( R5s.lt.p_lxRs) then
        if ( p_lxRs==ZERO) then
           R5s=0.5*chM5s*tom2_K5s*dnm
           Write(LOGUNIT,*) 'Value of K5s larger than max.equilibirum value'
           Write(LOGUNIT,*) 'K5s reset to:', R5s
           K5s(BoxNumberXY)=R5s
         else
           R5s=p_lxRs*tom2_K5s*dnm
         endif
      endif
      if ( R15s.lt.p_lxRs) then
        if ( p_lxRs==ZERO) then
           R15s=0.5*chM5s*tom2_an*(p_d_tot_2-dnm)
           Write(LOGUNIT,*) 'Value of K15s larger than max.equilibirum value'
           Write(LOGUNIT,*) 'K15s reset to:', R15s
           K15s(BoxNumberXY)=R15s
         else
           R15s=p_lxRs*tom2_an*(p_d_tot_2-dnm)
         endif
      endif
      RN5s_Ben=chM5s-N5s_Ben(BoxNumberXY)

      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate for the above defined set of boundary conditions
      ! the steady-state profiles and return the average concentration
      ! in the oxic and denitrification layers
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      dry=( dry_z(BoxNumberXY) <=p_clDxm .and. p_dry_ben)
      lShift=ZERO;mShift=ZERO;lShift=ZERO
      if ( mode.ne.2) then
        call BenSilicaEquation(KSiO3(BoxNumberXY),InitializeModel,sw_method,&
          .true.,RN5s_Ben, R5s,R15s,cD1m,cD2m,dnm,p_d_tot_2,&
          p_s_ads,pShift,mShift,lShift,p_poro(BoxNumberXY),p_p_ae,p_p_an, &
          diff,lambda,p_s_ads,chM5s,rus,alpha,p_sul,suD1,suD2,suDn)

        if ( InitializeModel== 1) then
          cM5s = CalculateSet( KSiO3(BoxNumberXY), 0, 0, 0, dummy, dummy)
          return
        endif

        cM5s =CalculateSet( KSiO3(BoxNumberXY), SET_LAYER_INTEGRAL, &
                                 LAYER1, LAYER3, dummy, ZERO)/ (dnm*tom2_K5s)
        KSiO3eq(BoxNumberXY)=CopySet(KSiO3(BoxNumberXY), KSiO3eq(BoxNumberXY))

        jQ6K5s=max(ZERO,suD1* CalculateFromSet(KSiO3(BoxNumberXY), &
                   EXPONENTIAL_INTEGRAL, RFLUX, cD1m, cD2m))  &
         +max(ZERO,suD2* CalculateFromSet( KSiO3(BoxNumberXY), &
                   EXPONENTIAL_INTEGRAL, RFLUX, cD2m, dnm))
        r=Q6s(BoxNumberXY)*PartQ(D9m(BoxNumberXY),cD1m,dnm, p_d_tot_2)
        call LimitChange(POSITIVE,jQ6K5s,r,max_change_per_step)

         !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
         ! Start calculation of fluxes controlled  by eq profile
         !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

         !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
         ! Calculate the adaptation time to the steady-state profile
         !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

         r=jQ6K5s/R5s
         Tau  =   CalculateTau(r,  diff,  p_p_ae,  dnm)

         !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
         ! Estimate the average value of M5s over the actual time step
         ! (transient value).
         ! This value depends on the adaptation time, the actual time step,
         ! the ''old'' value and the ''equilibrium value''
         !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

         cM5s = cM5s+( mM5s- cM5s)* IntegralExp(-LocalDelta/Tau,DONE)

         !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
         ! 1.Store equilibrium profile
         ! 2.Derive the equations for the transient profiles, assuming the same
         ! solution as for the steady-state case and using cM5s as new &
         ! constraint.
         !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

         r = CalculateSet( KSiO3(BoxNumberXY),ADD,0,0,dummy,dnm* cM5s*tom2_K5s)
      endif

      if ( mode.le.1 ) then
        ! This option  is meant to calculate at the start of the calculations
        ! for the benthic system dynamics the concentration of silica in
        ! the porewater
        M5s(BoxNumberXY) = max(ZERO, &
         chM5s- max(ZERO,CalculateFromSet( KSiO3eq(BoxNumberXY), &
         INTEGRAL, STANDARD, ZERO, D1m(BoxNumberXY)))/ D1m(BoxNumberXY))

        M15s(BoxNumberXY) = max(ZERO, &
        chM5s- max(ZERO,CalculateFromSet( KSiO3eq(BoxNumberXY), &
         INTEGRAL, STANDARD,  D1m(BoxNumberXY),dnm))/ (dnm-D1m(BoxNumberXY)))

        return
      endif

      Dnewm  = min(max( D2m(BoxNumberXY)+shiftD2m(BoxNumberXY)*LocalDelta, &
                                                         p_clD2m),  p_chD2m)
      lShift=Dnewm-cD2m
      if ( abs(lShift) > 5.0D-07.and. p_p_ae> p_p_an) then
       WithSpecialLayer=sw_special_shift
       if ( lShift .gt.ZERO) then
         mShift=-suD2* CalculateFromSet( KSiO3eq(BoxNumberXY), &
                              EXPONENTIAL_INTEGRAL, RFLUX, cD2m,Dnewm)
         pShift=ZERO
       else
         mshift=-suD1* CalculateFromSet( KSiO3eq(BoxNumberXY), &
                              EXPONENTIAL_INTEGRAL, RFLUX, Dnewm, cD2m)
         !pShift is >0.0
         pShift=- CalculateFromSet( KSiO3eq(BoxNumberXY), &
               SHIFT, LAYER2,cD2m, Dnewm)*max(ZERO,p_p_ae-p_p_an)/(DONE+p_p_ae)
       endif
       call BenSilicaEquation(KSiO3(BoxNumberXY),WithSpecialLayer,sw_method, &
         .true.,RN5s_Ben,R5s,R15s,cD1m,cD2m,dnm,p_d_tot_2,&
         p_s_ads,pShift,mShift,lShift,p_poro(BoxNumberXY),p_p_ae,p_p_an, &
         diff,lambda,p_s_ads,chM5s,rus,alpha,p_sul,suD1,suD2,suDn)

       cM5s = CalculateSet( KSiO3(BoxNumberXY), SET_LAYER_INTEGRAL, LAYER1, &
                                          LAYER4, dummy, ZERO)/ (dnm*tom2_K5s)

       jQ6K5s=max(ZERO,suD1* CalculateFromSet(KSiO3(BoxNumberXY), &
                 EXPONENTIAL_INTEGRAL, RFLUX, cD1m, cD2m+min(ZERO,lShift))) &
         +max(ZERO,suD2* CalculateFromSet( KSiO3(BoxNumberXY), &
                   EXPONENTIAL_INTEGRAL, RFLUX, cD2m+max(ZERO,lShift), dnm))
        r=Q6s(BoxNumberXY)*PartQ(D9m(BoxNumberXY),cD1m,dnm, p_d_tot_2)
        call LimitChange(POSITIVE,jQ6K5s,r,max_change_per_step)

       !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       ! Calculate the adaptation time to the steady-state profile
       !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

       r=jQ6K5s/R5s
       Tau  =   CalculateTau(r,  diff,  p_p_ae,  dnm)

       !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       ! Estimate the average value of M5s over the actual time step
       ! (transient value).
       ! This value depends on the adaptation time, the actual time step,
       ! the ''old'' value and the ''equilibrium value''
       !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

       cM5s = cM5s+( mM5s- cM5s)* IntegralExp(-LocalDelta/Tau,DONE)
       r = CalculateSet( KSiO3(BoxNumberXY),ADD,0,0,dummy,dnm*cM5s*tom2_K5s)

       if (isnan(cM5s)) then
         write(LOGUNIT,*) 'cM5s=',cM5s
         write(LOGUNIT,*) 'mM5s,Tau=',mM5s,Tau
       endif

     else
        WithSpecialLayer=0
     endif

     !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     ! Calculate flux at the sediment/water interface:
     !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      r= CalculateFromSet( KSiO3(BoxNumberXY),INTEGRAL, MASS, ZERO, cD1m)
      K0s= (chM5s-N5s_Ben(BoxNumberXY))* cD1m*tom2_ae-r
      R0s= r+N5s_Ben(BoxNumberXY)* cD1m*tom2_ae

       if (isnan(K0s)) then
         write(LOGUNIT,*) 'K0s=',K0s
         write(LOGUNIT,*) 'chM5s,tom2_ae,r=',chM5s,tom2_ae,r
       endif

     if ( dry) then
       jK5N5s = ZERO
     else
       jK5N5s=-CalculateFromSet(KSiO3(BoxNumberXY),DERIVATIVE,RFLUX,ZERO,dummy)
       if (isnan(jK5N5s)) then
         write(LOGUNIT,*) '1: jK5N5s=',jK5N5s
       endif
       call LimitShift(jK5N5s,N5s_Ben(BoxNumberXY)*Depth_Ben(BoxNumberXY), &
                     K5s(BoxNumberXY)-rus*LocalDelta,max_change_per_step,lim_1)
       if (isnan(jK5N5s)) then
         write(LOGUNIT,*) '2: jK5N5s=',jK5N5s
       endif
       call LimitChange(POSITIVE,jK5N5s,K0s,max_change_per_step,lim_2)
       jK5N5s=min(lim_1,lim_2)*jK5N5s
       if (isnan(jK5N5s)) then
         write(LOGUNIT,*) '2: jK5N5s=',jK5N5s
       endif
     endif

 
     r= CalculateFromSet( KSiO3(BoxNumberXY),INTEGRAL, MASS, ZERO, dnm)
     !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     ! the dissolution fluxes for K5s+K15s
     ! calculate dfissolution in first layer seperately
     ! limit fluxes for to large rates.
     !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

     jQ6K5s_1= max(ZERO,sM5s *CalculateFromSet &
                    (KSiO3(BoxNumberXY),INTEGRAL, RFLUX, ZERO, cD1m))
     r=cD2m+lShift
     jQ6K5s   =  jQ6K5s_1  &
         +max(ZERO,suD1* CalculateFromSet( KSiO3(BoxNumberXY), &
                   EXPONENTIAL_INTEGRAL, RFLUX, cD1m, min(r,cD2m))) &
         -mshift &
         +max(ZERO,suD2* CalculateFromSet( KSiO3(BoxNumberXY), &
                   EXPONENTIAL_INTEGRAL, RFLUX, max(r,cD2m), dnm))
     if (isnan(jQ6K5s).or.isnan(jK5N5s)) then
      write(LOGUNIT,*) 'reBTs=',reBTs(boxNumberXY),' ruBTs=',ruBTs(boxNumberXY)
      write(LOGUNIT,*) 'jK5N5s=',jK5N5s
      write(LOGUNIT,*) 'lim_1, lim_2 ',lim_1, lim_2
      write(LOGUNIT,*) 'N2s_Ben, Depth_Ben ',N5s_Ben(BoxNumberXY),Depth_Ben(BoxNumberXY)
      write(LOGUNIT,*) 'K5s, rus, K0s ',K5s(BoxNumberXY),rus,K0s
      write(LOGUNIT,*) 'jQ6K5s=',jQ6K5s ,' M5b0=',M5b0 ,' Q6s=',Q6s
      write(LOGUNIT,*) 'K5s=',K5s(BoxNumberXY) ,' K15s=',K15s(BoxNumberXY)
      write(LOGUNIT,*) 'suD1=',suD1 ,' suD2=',suD2
      write(LOGUNIT,*) 'sM5s=',sM5s ,' cD1m=',cD1m ,' cD2m=',cD2m
      write(LOGUNIT,*) 'KSiO3=',KSiO3(BoxNumberXY)
      write(LOGUNIT,*) 'jQ6K5s_1=',jQ6K5s_1 ,' alpha=',alpha
      write(LOGUNIT,*) 'dnm=',dnm,'N5s_Ben=',N5s_Ben(BoxNumberXY)
      call PrintSet(KSiO3(BoxNumberXY),"NaN in jK5N5s")
      call BFM_ERROR("BenSilica","")
     endif
     call LimitChange(NEGATIVE,jQ6K5s_1,R0s,max_change_per_step,lim_3)
     call LimitChange(POSITIVE,jQ6K5s,R5s,max_change_per_step,lim_1)
 
     jQ6K15s = max(ZERO,suDn* CalculateFromSet( &
          KSiO3(BoxNumberXY), EXPONENTIAL_INTEGRAL, RFLUX, dnm, p_d_tot_2))
     call LimitChange(POSITIVE,jQ6K15s,R15s ,max_change_per_step,lim_4)
     if (isnan(jQ6K15s)) then
        write(LOGUNIT,*) 'suD1=',suD1 , ' lambda=',lambda
        write(LOGUNIT,*) 'suD2=',suD2 , ' suDn=',suDn
        write(LOGUNIT,*) 'alpha=',alpha , ' dnm=',dnm
     endif

     h= jQ6K5s-jQ6K5s_1
     r=Q6s(BoxNumberXY)*PartQ(D9m(BoxNumberXY),cD1m,dnm, p_d_tot_2)
     call LimitChange(POSITIVE,h,r,max_change_per_step,lim_2)
     h  =h *min(lim_1,lim_2)

     r=Q6s(BoxNumberXY)*PartQ(D9m(BoxNumberXY),ZERO,cD1m, p_d_tot_2)
     call LimitChange(POSITIVE,jQ6K5s_1,r,max_change_per_step,lim_2)
     jQ6K5s_1=jQ6K5s_1* min(lim_1,lim_2,lim_3)
     jQ6K5s  =(h +jQ6K5s_1)

     r=Q6s(BoxNumberXY)*PartQ(D9m(BoxNumberXY),dnm,p_d_tot_2, p_d_tot_2)
     call LimitChange(POSITIVE,jQ6K15s,r,max_change_per_step,lim_2)

     jQ6K15s =jQ6K15s * min(lim_2,lim_4)

     reBTs(BoxNumberXY)=reBTs(BoxNumberXY)+jQ6K5s_1
     reATs(BoxNumberXY)=jQ6K5s-jQ6K5s_1+jQ6K15s
     call flux(BoxNumberXY, iiBen, ppQ6s, ppK5s, jQ6K5s )
     call flux(BoxNumberXY, iiBen, ppQ6s, ppK15s,jQ6K15s )
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! -Calculate flux at cD2m
      ! + Calculate new depth of of cD2m and the flux of silicate related
      !   to this shifting
      ! -limit fluxes for too large changes
      ! -Calculate flux at underside
      ! -limit fluxes for too large changes
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      shiftmass=ZERO; Dnewm=ZERO
      jK15K5s=-CalculateFromSet(KSiO3(BoxNumberXY),DERIVATIVE,RFLUX,dnm,dummy)

      if (isnan(jK15K5s)) then
        write(LOGUNIT,*)  &
               'Calculation flux jK15K5s is NaN, jK15K5s is reset on zero.'
        jK15K5s=ZERO
        write(LOGUNIT,*) 'alpha=',alpha,' jQ6K5s_1=',jQ6K5s_1
        write(LOGUNIT,*) 'M0=',M0,' M5bd1=',M5bD1
        write(LOGUNIT,*) 'chM5s=',chM5s,' mM5s=',mM5s
        write(LOGUNIT,*) 'smQ6=',smQ6,' cD1m=',cD1m
        write(LOGUNIT,*) 'cD2m=',cD2m,' jQ6K5s=',jQ6K5s
        write(LOGUNIT,*) 'jQ6K15s=',jQ6K15s,' R15s=',R15s
        write(LOGUNIT,*) 'K15s=',K15s(BoxNumberXY),' Q6s=',Q6s(BoxNumberXY)
        write(LOGUNIT,*) 'M5s=',K5s(BoxNumberXY)/tom2_ae/cD2m ,&
                         'M15s=',K15s(BoxNumberXY)/tom2_an/(p_d_tot_2-cD2m)
        call PrintSet(KSiO3(BoxNumberXY),'Problem with Silica')
      endif

      Dnewm  =   dnm+ lShift
      shiftmass= CalculateFromSet( KSiO3(BoxNumberXY), SHIFT,  &
                LAYER3+WithSpecialLayer, dnm, Dnewm)/ LocalDelta

      jK15K5s=jK15K5s+ &
               sign(chM5s*tom2_an*abs(lShift),lShift)/LocalDelta-shiftmass
      h=-jK15K5s+shiftmass

      call LimitShift(jK15K5s,K5s(BoxNumberXY),K15s(BoxNumberXY),  &
                                             max_change_per_step,lim_1)

      call LimitShift(h,R5s,R15s ,max_change_per_step,lim_2)
      jK15K5s=min(lim_1,lim_2)*jK15K5s

      call addbotflux(ANY,BoxNumberXY,iiBen,ppK5s,iiPel,ppN5s,jK5N5s)

      call flux(BoxNumberXY, iiBen, ppK15s,ppK5s,  jK15K5s *insw( jK15K5s)  )
      call flux(BoxNumberXY, iiBen, ppK5s, ppK15s,-jK15K5s *insw(-jK15K5s)  )

      r= chM5s- CalculateFromSet( KSiO3(BoxNumberXY), &
                                       EQUATION, STANDARD, p_d_tot_2, dummy)
      i=p_flux_at_deep_end; if ( r< ZERO) i=1
      select case (i)
        case (1); s=   ZERO
        case (2); s = -max(ZERO,CalculateFromSet(KSiO3(BoxNumberXY), &
                                DERIVATIVE, RFLUX, p_d_tot_2, dummy))
        case (3); s = -CalculateFromSet(KSiO3(BoxNumberXY), &
                                DERIVATIVE, RFLUX, p_d_tot_2, dummy)
      end select
      jK25K15s(BoxNumberXY)=s
      call LimitChange(POSITIVE,s,R15s ,max_change_per_step,lim_1)
      call LimitChange(ANY,jK25K15s(BoxNumberXY),K15s(BoxNumberXY) ,&
                                                max_change_per_step,lim_2)

      jK25K15s(BoxNumberXY)=jK25K15s(BoxNumberXY)* min(lim_1,lim_2)
      call flux(BoxNumberXY, iiBen, ppK15s,ppK15s, jK25K15s(BoxNumberXY) )

      r=log(0.5*(DONE+ExpDist(-alpha,cD1m)))/(-alpha)
      s=-jQ6K5s_1*(r- D9m(BoxNumberXY))/(NZERO+Q6s(BoxNumberXY))
      call LimitChange(NEGATIVE,s,abs(D9m(BoxNumberXY)),max_change_per_step)
      call flux(BoxNumberXY,iiBen, ppD9m,ppD9m,s)

      s= (jQ6K5s-jQ6K5s_1)*((dnm+cD1m)*0.5- D9m(BoxNumberXY)) &
                                         /(NZERO+Q6s(BoxNumberXY))
      call flux(BoxNumberXY,iiBen, ppD9m,ppD9m,s)

      s= ((p_d_tot_2+dnm)*0.5D+00-D9m(BoxNumberXY))* (-jQ6K15s) &
                                            /( NZERO+Q6s(BoxNumberXY))
       call flux(BoxNumberXY, iiBen, ppD9m, ppD9m,s)

      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate the pore-water concentrations from the state variables
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      M5s(BoxNumberXY) = max(ZERO, &
               chM5s- max(ZERO,CalculateFromSet( KSiO3(BoxNumberXY), &
               INTEGRAL, STANDARD, ZERO, D1m(BoxNumberXY)))/ D1m(BoxNumberXY))
      M15s(BoxNumberXY) = max(ZERO, &
               chM5s- max(ZERO,CalculateFromSet( KSiO3eq(BoxNumberXY), &
         INTEGRAL, STANDARD,  D1m(BoxNumberXY),dnm))/ (dnm-D1m(BoxNumberXY)))
    enddo
  end subroutine BenSilicaDynamics

  subroutine BenSilicaEquation(KSiO3,mode,option, &
                 dry,N5s,R5s,R15s,D1m,D2m,dnm,d_tot, &
                 p_s_ads,pShift,mShift,lShift,poro,p_ae,p_an, &
                 diff,lambda,sM5s,chM5s,rus,alpha,slu,suD1,suD2,suDn)
  use global_mem, ONLY:RLEN,ZERO,NZERO,DONE,dummy,idummy
  use constants, ONLY: LAYERS, &
    LAYER1,LAYER2,LAYER3,LAYER4, LAYER5,DIFFUSION, FOR_ALL_LAYERS, &
    POROSITY, ADSORPTION, DEFINE, CONSTANT_TERM, SET_CONTINUITY,FLAG, MASS,&
    SET_BOUNDARY,ADD, EQUATION, BESSELI_EXP_TERM, BESSELK_EXP_TERM, &
    EXPONENTIAL_TERM,PARAMETER_DEFINE,QUADRATIC_TERM,LINEAR_TERM, &
    SET_LAYER_INTEGRAL_UNTIL, SET_LAYER_INTEGRAL, DERIVATIVE, &
    INPUT_TERM,STANDARD, MIN_VAL_EXPFUN

  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following bennut functions are used:InitializeSet, &
  ! DefineSet, CompleteSet, CalculateSet, CalculateTau, CalculateFromSet
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use bennut_interface, ONLY: InitializeSet, DefineSet, CompleteSet

  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following sesame functions are used:IntegralExp, insw
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  implicit none


    integer,intent(INOUT)     ::KSiO3
    integer ,intent(IN)       ::mode
    integer ,intent(IN)       ::option
    logical ,intent(IN)       ::dry
    real(RLEN),intent(IN)     ::N5s
    real(RLEN),intent(IN)     ::R5s
    real(RLEN),intent(IN)     ::R15s
    real(RLEN),intent(IN)     ::D1m
    real(RLEN),intent(IN)     ::D2m
    real(RLEN),intent(IN)     ::dnm
    real(RLEN),intent(IN)     ::d_tot
    real(RLEN),intent(IN)     ::p_s_ads
    real(RLEN),intent(IN)     ::pShift
    real(RLEN),intent(IN)     ::mShift
    real(RLEN),intent(IN)     ::lShift
    real(RLEN),intent(IN)     ::poro
    real(RLEN),intent(IN)     ::p_ae
    real(RLEN),intent(IN)     ::p_an
    real(RLEN),intent(IN)     ::diff
    real(RLEN),intent(IN)     ::lambda
    real(RLEN),intent(IN)     ::sM5s
    real(RLEN),intent(IN)     ::chM5s
    real(RLEN),intent(IN)     ::rus
    real(RLEN),intent(IN)     ::alpha
    real(RLEN),intent(IN)     ::slu
    real(RLEN),intent(IN)     ::suD1
    real(RLEN),intent(IN)     ::suD2
    real(RLEN),intent(IN)     ::suDn

      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Initialize the set of differential equations giving:
      ! - n. of layers
      ! - n. of coefficients
      ! - layer depths
      ! - environmental conditions (diffusion, porosity and adsorption coeff.)
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    integer         ::i
    real(RLEN)      :: ds1,ds2,p_ads_3
    real(RLEN)      ::gamma
    real(RLEN)      ::r,s,t          ! help variables

    if ( mode <= 0) then
         p_ads_3=p_an; ds1=D2m;       ds2=D2m
    else if (lShift.lt.ZERO) then
         p_ads_3=p_an; ds1=D2m+lShift;ds2=D2m
    else
         p_ads_3=p_ae; ds1=D2m;       ds2=ds1 +lShift
    endif


    select case (mode)
      case (-1)
        KSiO3 = InitializeSet( KSiO3, 2, 4)
        call  DefineSet( KSiO3, LAYERS, LAYER1, 0,      D1m, dummy)
        call DefineSet(KSiO3,DIFFUSION, LAYER1,LAYER2,diff,diff)
      case (0)
        KSiO3 = InitializeSet( KSiO3, 4, 9)
        call  DefineSet( KSiO3, LAYERS, LAYER1, LAYER2, D1m, D2m)
        call  DefineSet( KSiO3, LAYERS, LAYER3, idummy, dnm, dummy)
        call DefineSet(KSiO3,DIFFUSION, LAYER1,LAYER2,diff,diff)
        call DefineSet(KSiO3,DIFFUSION, LAYER3,LAYER4,diff,diff)
      case (1)
        KSiO3 = InitializeSet( KSiO3, 5, 12)
        call  DefineSet( KSiO3, LAYERS, LAYER1, LAYER2, D1m, ds1)
        call  DefineSet( KSiO3, LAYERS, LAYER3, LAYER4, ds2, dnm)
        call DefineSet(KSiO3,DIFFUSION, LAYER1,LAYER2,diff,diff)
        call DefineSet(KSiO3,DIFFUSION, LAYER3,LAYER4,diff,diff)
        call DefineSet(KSiO3,DIFFUSION, LAYER5,0,diff,dummy)
    end select

    call DefineSet(KSiO3,POROSITY,  FOR_ALL_LAYERS, idummy, poro, dummy)
    call DefineSet(KSiO3,ADSORPTION,LAYER1, LAYER2, p_ae, p_ae)

    select case (mode)
      case (0)
        call DefineSet(KSiO3,ADSORPTION,LAYER3, LAYER4, p_an, p_an)
      case (1)
        call DefineSet(KSiO3,ADSORPTION,LAYER3, LAYER4, p_ads_3, p_an)
        call DefineSet(KSiO3,ADSORPTION,LAYER5, idummy,p_an,dummy)
    end select

    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Define coefficients for the steady-state solutions in each layer
    ! General solution of the equilibrium profile:
    ! C = Ssat - S
    !
    ! 1st layer:C(z) = c13*z^2 + c14*z + c15
    ! 2nd layer:C(z) = c21*I0*exp[-alpha*(z-cD1m)]+ c22*K0*exp[-alpha*(z-cD1m)]
    ! 3nd layer:C(z) = c21*I0*exp[-alpha*(z-cD1m)]+ c22*K0*exp[-alpha*(z-cD1m)]
    ! ( I0 and K0  = modified Bessel functions of 0-order)
    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    gamma=max(MIN_VAL_EXPFUN,sqrt(p_s_ads/diff))
    select case (option)
      case(1,3)
        call DefineSet( KSiO3, DEFINE, 13, QUADRATIC_TERM, dummy, dummy)
        call DefineSet( KSiO3, DEFINE, 14, LINEAR_TERM, dummy, dummy)
        call DefineSet( KSiO3, DEFINE, 15, CONSTANT_TERM, dummy, dummy)
      case(2,4)
        call DefineSet( KSiO3, DEFINE, 11, EXPONENTIAL_TERM, -gamma,  dummy)
        call DefineSet( KSiO3, DEFINE, 12, EXPONENTIAL_TERM,  gamma, dummy)
        call DefineSet( KSiO3, DEFINE, 15, CONSTANT_TERM, dummy, dummy)
    end select

    select case (mode)
      case (-1)
        call DefineSet( KSiO3,PARAMETER_DEFINE, 21, &
                                    BESSELI_EXP_TERM, -alpha, max(slu,suD1))
      case (0)
        call DefineSet( KSiO3, PARAMETER_DEFINE, 21, &
                                    BESSELI_EXP_TERM, -alpha, max(slu,suD1))
        call DefineSet( KSiO3, PARAMETER_DEFINE, 22, &
                                    BESSELK_EXP_TERM, -alpha, max(slu,suD1))
        i=30
      case (1)
        call DefineSet( KSiO3,    PARAMETER_DEFINE, 21, &
                                    BESSELI_EXP_TERM, -alpha, max(slu,suD1))
        call DefineSet( KSiO3, PARAMETER_DEFINE, 22, &
                                    BESSELK_EXP_TERM, -alpha, max(slu,suD1))
        call DefineSet( KSiO3, DEFINE, 31, EXPONENTIAL_TERM, -gamma,  dummy)
        call DefineSet( KSiO3, DEFINE, 32, EXPONENTIAL_TERM,  gamma, dummy)
        call DefineSet( KSiO3, DEFINE, 35, CONSTANT_TERM, dummy, dummy)
        i=40
    end select

    if ( mode.ge.0) then
      call DefineSet( KSiO3, PARAMETER_DEFINE, i+1, &
                                      BESSELI_EXP_TERM, -alpha, max(slu,suD2))
      call DefineSet( KSiO3, PARAMETER_DEFINE, i+2, &
                                      BESSELK_EXP_TERM, -alpha, max(slu,suD2))

      call DefineSet( KSiO3, PARAMETER_DEFINE, i+11, &
                                      BESSELI_EXP_TERM, -alpha, max(slu,suDn))
      call DefineSet( KSiO3, PARAMETER_DEFINE, i+12, &
                                      BESSELK_EXP_TERM, -alpha, max(slu,suDn))
    endif

    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Insert other boundary conditions and continuity between layers:
    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    !2/6/8
    call CompleteSet( KSiO3, SET_CONTINUITY, FLAG, MASS, dummy)
    !3/7/9
      call CompleteSet(KSiO3,SET_BOUNDARY,LAYER1, EQUATION, ZERO, value=N5s)

    select case (mode)
      case(-1)
      case(0)
        !8,9
        call CompleteSet(KSiO3,SET_LAYER_INTEGRAL_UNTIL,LAYER4,LAYER4,d_tot, &
                                                         value=max(NZERO,R15s))
      case (1)
        !10,12
        r=+p_s_ads* pShift/abs(lShift)/(DONE+p_ads_3)/poro+mShift/abs(lShift)
        call CompleteSet(KSiO3,INPUT_TERM,35,STANDARD, dummy, value=r/p_s_ads)

        call CompleteSet(KSiO3,SET_LAYER_INTEGRAL_UNTIL,LAYER5,LAYER5,d_tot, &
                                                        value=max(NZERO,R15s))
    end select
    !5,11,12
    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Explanation last boundary condition:
    ! S:Silica in porewater   C=chM5s-S
    ! dissolution rate: d= sM5s/poro/(p_ae+1)  ,
    !                       desorption rate (mode=2,4) a= p_s_ds*(p_ae+1)
    ! Procedure calculates Mass per/m2 in oxic layer:
    !    (mode 1,2): M_C= (chM5s-S)*poro*D1m*(p_ae+1)
    !    (mode 3,4) :M_S= chM5s*poro_D1m*(p_ae+1) -M_C
    ! Important to know:  dC/dt=-dS/dt
    ! mode ==1 : poro*diff *( dC/dt(0)-dC/dt(D1m)) -d    *M_C +rus =0
    ! mode ==2 : poro*diff *( dC/dt(0)-dC/dt(D1m)) -(d+a)*M_C +rus =0
    ! mode ==3 : poro*diff *(-dS/dt(0)+dS/dt(D1m)) +d    *M_S -rus =0
    ! mode ==4 : poro*diff *(-dS/dt(0)+dS/dt(D1m)) +(d-a)*M_S -rus =0
    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Change sign for options 1 and 2
    i=sign(1,2*option-5)
    t=i
    r=diff*poro*t
    ! No addition for options 1 and 2
    s=0.5*(1-i)
    call CompleteSet( KSiO3,SET_BOUNDARY,LAYER1, DERIVATIVE, D1m, mfac=+r)
    if ( .not.dry) &
    call CompleteSet( KSiO3,ADD,         LAYER1, DERIVATIVE, ZERO,mfac=-r)
    select case (option)
      case (1)
        r=+t*sM5s/poro/(p_ae +DONE)
        s= s* chM5s*sM5s*D1m+t*rus
        call CompleteSet( KSiO3,ADD+SET_LAYER_INTEGRAL,LAYER1, LAYER1, &
                                              dummy, mfac=-r,value=-s)
      case (2)
        r=p_s_ads*(p_ae+DONE)/poro+t*sM5s/poro/(p_ae +DONE)
        s= s * chM5s*(p_s_ads*(p_ae+DONE)+sM5s)*D1m+t*rus
        call CompleteSet( KSiO3,ADD+SET_LAYER_INTEGRAL,LAYER1, LAYER1, &
                                                dummy, mfac=-r,value=-s)
    end select
    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Calculate for the above defined set of boundary conditions
    ! the steady-state profiles and return the average concentration
    ! in the oxic and denitrification layers
    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
return
end subroutine BenSilicaEquation

!BOP
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
