#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenPhosphate
!
! DESCRIPTION
!   Description of the phosphate diagenitic processes in the sediment
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
  subroutine BenPhosphateDynamics(what)
!
! !USES:

  ! For the following Benthic-states fluxes are defined: K1p, K11p, K21p
  ! The following Benthic-states are used (NOT in fluxes): D1m, D2m, D8m
  ! The following global scalar vars are used: &
  ! NO_BOXES_XY, BoxNumberXY, InitializeModel, LocalDelta
  ! The following Benthic 1-d global boxvars are modified : M1p, M11p, M21p, &
  ! KPO4
  ! The following Benthic 1-d global boxvars are used: reBTp, &
  ! reATp, irrenh, ETW_Ben, N1p_Ben, Depth_Ben, shiftD1m, shiftD2m
  ! The following Benthic 1-d global boxpars  are used: p_poro, p_p_ae
  ! The following 0-d global parameters are used: p_d_tot_2, p_clDxm, p_q10diff
  ! The following global constants are used: RLEN,ZERO

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  USE BFM_ERROR_MSG, ONLY: set_warning_for_getm
  use global_mem, ONLY:RLEN,NZERO,ZERO,DONE,LOGUNIT,dummy
  use mem,  ONLY: K1p, K11p, K21p, D1m, D2m, D8m,DH2m,DH3m, D2STATE
  use mem, ONLY: ppN1p, ppK1p, ppK11p, ppK21p, NO_BOXES_XY,   &
     BoxNumberXY, InitializeModel, M1p, M11p, M21p, KPO4,KPO4sh, &
    reBTp, ruBTp, reATp, irrenh, &
    ETW_Ben, N1p_Ben, Depth_Ben, shiftD1m, shiftD2m, & 
    jK31K21p, iiPel, iiBen, flux, LocalDelta,max_change_per_step,dry_z
  use mem,only: Mp1p
  use constants, ONLY: LAYER1, LAYER3 , SET_LAYER_INTEGRAL, &
    STANDARD, ADD, DERIVATIVE, RFLUX, SHIFT, INTEGRAL,ANY,NEGATIVE
  use mem_Param,ONLY:p_poro,p_pK1_ae,p_d_tot,p_d_tot_2,p_clDxm,p_q10diff, &
                     p_dry_ben,combine_anabac
  use mem_Diffusion,ONLY: p_diff=>p_diff_N1
  use botflux,onlY:addbotflux
  use mem_BenPhosphate
  use mem_BenthicNutrient3, ONLY:sw_special_shift,p_pAn2Ni
  use LimitRates, ONLY:LimitShift,LimitChange

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following bennut functions are used:InitializeSet, &
  ! DefineSet, CompleteSet, CalculateSet, CalculateTau, CalculateFromSet
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use bennut_interface, ONLY: CalculateSet, CalculateTau, CalculateFromSet

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following global functions are used:eTq
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_interface,   ONLY: eTq
! use mem,  ONLY: Output2d_1,Output2d_2,Output2d_3, Output2d_4



  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following sesame functions are used:IntegralExp, insw
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: IntegralExpDist,ExpDist,SetExpDist,IntegralExp,insw
!  
!
! !AUTHORS
!   Original version by  P. Ruardij
!
!
! !REVISION_HISTORY
!       September 1999 by M. Vichi !               Commented version
!
! COPYING
!   
!   Copyright (C) 2006 P. Ruardij & M.Vichi
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
  integer,intent(IN)              :: what

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer  :: KSHIFT=0
  real(RLEN)  :: diff
  real(RLEN)  :: jK1N1p
  real(RLEN)  :: jK11K1p
  real(RLEN)  :: jK21K11p
  real(RLEN)  :: seBT
  real(RLEN)  :: zuBT
  real(RLEN)  :: zuD1
  real(RLEN)  :: zuD2
  real(RLEN)  :: alpha,beta
  real(RLEN)  :: cK1p
  real(RLEN)  :: Tau
  real(RLEN)  :: dnm,jflux
  real(RLEN)  :: Dnew
  real(RLEN)  :: r
  real(RLEN)  :: s_ads
  integer  :: mode
  real(RLEN)  :: lShift,pShift,mShift
  logical dry

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  do BoxNumberXY=1,NO_BOXES_XY


    if ( what.eq.0) then
      if ( K1p(BoxNumberXY)<ZERO) then
        r=K1p(BoxNumberXY)
        K1p(BoxNumberXY)=N1p_Ben(BoxNumberXY)*p_poro(BoxNumberXY) * &
                              D1m(BoxNumberXY)*p_pK1_ae(BoxNumberXY)
        K11p(BoxNumberXY)=0.5 * K11p(BoxNumberXY)+ &
                          max(ZERO,0.5*K11p(BoxNumberXY)-K1p(BoxNumberXY)+r)
        write(LOGUNIT,*) "D1m=",D1m(BoxNumberXY)
        write(LOGUNIT,*) "K1p<0.0 , K1p reset to", K1p(BoxNumberXY)
        write(LOGUNIT,*) "K11p reset to", K11p(BoxNumberXY)
        call set_warning_for_getm
      endif
      if ( K21p(BoxNumberXY)<ZERO) then
        r=K11p(BoxNumberXY)/(p_poro(BoxNumberXY) * &
                   (D2m(BoxNumberXY)-D1m(BoxNumberXY))*p_pK1_ae(BoxNumberXY))
        K21p(BoxNumberXY)=0.5*r* (p_poro(BoxNumberXY) * &
                   (p_d_tot-D2m(BoxNumberXY))*p_p_an)
        write(LOGUNIT,*) "K21p reset to", K21p(BoxNumberXY)
        call set_warning_for_getm
      endif
      return
    endif
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate coefficient for the e-folding distribution of the anoxic
      ! mineralization. D8.m is the average penetration depth for P-detritus
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      if ( combine_anabac) then
        alpha  =   SetExpDist(D8m(BoxNumberXY),Dlm=p_clDxm,Dmm=p_d_tot)
        beta=alpha
      else
        alpha  =   SetExpDist(DH2m(BoxNumberXY),Dlm=p_clDxm,Dmm=p_d_tot)
        beta  =    SetExpDist(DH3m(BoxNumberXY),Dlm=p_clDxm,Dmm=p_d_tot)
      endif

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Recalculate Mineralization m2 --> m3 porewater
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Average in the oxic layer:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      r = NZERO+reBTp(BoxNumberXY)/p_poro(BoxNumberXY)/ D1m(BoxNumberXY)
      zuBT=sign(abs(r)+1.079D-6,r)
      seBT = min(p_shxcons,max(p_slxcons,ruBTp(BoxNumberXY) /( &
       (NZERO+K1p(BoxNumberXY))/p_poro(BoxNumberXY)/ &
                                        (DONE+p_pK1_ae(BoxNumberXY)))))


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Anoxic Mineralization at D1.m, using the exponential distribution
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      if ( InitializeModel== 0) then
        !mineralization at sediment water interface
        zuD1 =  reATp(BoxNumberXY)/ p_poro(BoxNumberXY)/ IntegralExpDist( &
          - alpha, p_d_tot_2- D1m(BoxNumberXY))
        zuD1=sign(abs(zuD1)+1.079D-6,zuD1)
      else
        !mineralization at oxygen penetration depth
        zuD1 = max( 1.D-20, reBTp(BoxNumberXY))/ p_poro(BoxNumberXY)/ &
            IntegralExpDist( - alpha, D1m(BoxNumberXY))*ExpDist(-alpha , D1m(BoxNumberXY))
      endif

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Anoxic Mineralization at D2.m
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      dnm=    (DONE-p_pAn2Ni)* D2m(BoxNumberXY) +p_pAn2Ni *p_d_tot 

      !mineralization at sulphide horizon
      zuD2  =   zuD1* ExpDist( - alpha, D2m(BoxNumberXY)- D1m(BoxNumberXY))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Correction due to environmental regulating factors,
      ! diffusion coefficient: temperature and bioirrigation
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      diff = p_diff* irrenh(BoxNumberXY)* eTq( ETW_Ben(BoxNumberXY), p_q10diff)



      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Initialize the set of differential equations giving:
      ! - n. of layers;
      ! - n. of coefficients
      ! - layers depths
      ! - environmental conditions (diffusion, porosity and adsorption coeff.)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      M1p(BoxNumberXY) = K1p(BoxNumberXY)/ p_poro(BoxNumberXY)/( &
        p_pK1_ae(BoxNumberXY)+ DONE)/( D1m(BoxNumberXY))
      s_ads=p_s_ads*M1p(BoxNumberXY)/(p_mK1p+ M1p(BoxNumberXY))
      dry=( dry_z(BoxNumberXY) <=p_clDxm .and. p_dry_ben)
      mode=0;
      call  BenPhosphateEquation(KPO4(BoxNumberXY),-InitializeModel, &
        .FALSE.,N1p_Ben(BoxNumberXY), K1p(BoxNumberXY),K11p(BoxNumberXY), &
        K21p(BoxNumberXY),D1m(BoxNumberXY),D2m(BoxNumberXY),p_d_tot, &
        p_d_tot_2,p_poro(BoxNumberXY),p_pK1_ae(BoxNumberXY),p_p_an,s_ads, &
        pShift,mShift,lShift,dnm, diff, alpha,beta,seBT,zuBT,zuD1,zuD2)


      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate for the above defined set of boundary conditions
      ! the steady-state profiles and return the vertically integrated
      ! concentration in the first layer.
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      cK1p = CalculateSet( KPO4(BoxNumberXY), SET_LAYER_INTEGRAL, LAYER1, &
        LAYER1, ZERO, ZERO)

      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate the adaptation time to the steady-state profile
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      Tau = CalculateTau( s_ads+seBT, diff, p_pK1_ae(BoxNumberXY), &
          D1m(BoxNumberXY))

      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Estimate the average value of K1p over the actual time step
      ! (transient value).
      ! This value depends on the adaptation time, the actual time step,
      ! the ''old'' value and the ''equilibrium value''
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      cK1p = cK1p+( K1p(BoxNumberXY)- cK1p)* IntegralExp( - LocalDelta/ &
             Tau, DONE)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Derive the equations for the transient profiles, assuming the same
      ! solution as for the steady-state case and using cK1p as new &
      ! constraint.
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      r = CalculateSet( KPO4(BoxNumberXY), ADD, 0, 0, ZERO, cK1p)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate the pore-water average concentrations from the state variables
      ! (Diagnostic variables, not used in calculations)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      M11p(BoxNumberXY) =  CalculateFromSet(KPO4(BoxNumberXY),INTEGRAL, &
              STANDARD,D1m(BoxNumberXY),dnm)/(dnm-D1m(BoxNumberXY));
!-------------------------------------------------------------------------
!     if ( isnan(M11p(BoxNumberXY))) then
!        call PrintSet(KPO4(BoxNumberXY),'Nan in M11p')
!        write(LOGUNIT,*) 'N1p_Ben',N1p_Ben(BoxnumberXY)
!        write(LOGUNIT,*) 'pShift',pShift
!        write(LOGUNIT,*) 'mShift',mShift
!        write(LOGUNIT,*) 'lShift',lShift
!        write(LOGUNIT,*) 'dnm',dnm
!        write(LOGUNIT,*)  'diff',diff
!        write(LOGUNIT,*)  'alpha',alpha
!        write(LOGUNIT,*) 'reBTp',reBTp(BoxnumberXY)
!        write(LOGUNIT,*) 'reATp',reATp(BoxnumberXY)
!        write(LOGUNIT,*) 'zuBT',zuBT
!        write(LOGUNIT,*) 'zuD1',zuD1
!        write(LOGUNIT,*) 'zuD2',zuD2
!        call set_warning_for_getm
!     endif
!-------------------------------------------------------------------------
      M21p(BoxNumberXY) = K21p(BoxNumberXY)/ p_poro(BoxNumberXY)/( p_p_an+ &
        DONE)/( p_d_tot_2- dnm )

      if ( InitializeModel== 0) then

         lShift=shiftD2m(BoxNumberXY)*LocalDelta;
         Dnew =   D2m(BoxNumberXY)+ lShift

         pShift=ZERO;
         if ( abs(lShift) > 0.0000005 .and.sw_special_shift==1) then
           if ( lShift<ZERO) pShift =-CalculateFromSet( KPO4(BoxNumberXY),  &
             SHIFT, LAYER3, D2m(BoxNumberXY), Dnew)*  &
            (p_pK1_ae(BoxNumberXY)-p_p_an)/p_pK1_ae(BoxNumberXY)
           mode=1;
           !average 0-order mineralization from D2m to D(2m)new 
           mShift=zuD2*IntegralExpDist(-alpha,lShift) 
!          write(LOGUNIT,*) 'BenPhoshateEquation 1'
           call  BenPhosphateEquation(KPO4sh(BoxNumberXY),1, &
            .FALSE.,N1p_Ben(BoxNumberXY),K1p(BoxNumberXY),K11p(BoxNumberXY), &
            K21p(BoxNumberXY),D1m(BoxNumberXY),D2m(BoxNumberXY),p_d_tot, &
            p_d_tot_2,p_poro(BoxNumberXY),p_pK1_ae(BoxNumberXY),p_p_an,s_ads,&
            pShift,mShift,lShift,dnm, diff, alpha,beta,seBT,zuBT,zuD1,zuD2)
           KSHIFT=KPO4sh(BoxNumberXY)
           cK1p = CalculateSet( KSHIFT, SET_LAYER_INTEGRAL, LAYER1, &
                                                  LAYER1, ZERO, ZERO)

           !=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
           ! Calculate the adaptation time to the steady-state profile
           !=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

           Tau = CalculateTau( s_ads, diff, p_pK1_ae(BoxNumberXY), &
             D1m(BoxNumberXY))

           !=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
           ! Estimate the average value of K1p over the actual time step
           ! (transient value).
           ! This value depends on the adaptation time, the actual time step,
           ! the ''old'' value and the ''equilibrium value''
           !=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

           cK1p = cK1p+( K1p(BoxNumberXY)- cK1p)* IntegralExp( - LocalDelta/ &
             Tau, DONE)

           !=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
           ! Derive the equations for the transient profiles, assuming the same
           ! solution as for the steady-state case and using cK1p as new &
           ! constraint.
           !=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

           r = CalculateSet( KSHIFT, ADD, 0, 0, dummy, cK1p)

         else
           KSHIFT=KPO4(BoxNumberXY)
         endif
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Start calculation of fluxes:
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Vertical fluxes :
        ! There are 2 problems with this model in this version both connected
        ! with the shifting of the layers:
        ! 1. shifting from the oxic+denitrification layer with high
        ! adsorped fraction phosphate to the lower anoxic layer with a very
        !  low percentage of adsorped phosphate.
        ! 2. Too large changes in spring due to large change of D1.m:
        !  This lead sometimes to a calculated phosphate gradient which
        !  has at some depth negative values.
        !
        !  Solution:
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Calculate flux at the sediment/water interface:
        ! Check on; to high fluxes from pelagic and on concisteny of gradient
        ! ( only flux of M1p > N1p!)
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        if ( dry) then
          jK1N1p =  ZERO
        else
          jK1N1p=CalculateFromSet(KPO4(BoxNumberXY),DERIVATIVE,RFLUX,ZERO,ZERO)
        endif
!-------------------------------------------------------------------------
        if (isnan(jK1N1p)) then
          write(LOGUNIT,*) 'Nan in jK1N1p'
          write(LOGUNIT,*) 'Mp1p,ruBTp,reBTp', &
                                              Mp1p,ruBTp,reBTp
          write(LOGUNIT,*) 'seBT,zuBT,zuD1,zuD2:',seBT,zuBT,zuD1,zuD2
         write(LOGUNIT,*) 'reBTp,D1m,p_poro=',reBTp,D1m,p_poro
         write(LOGUNIT,*) 'reATp,D2m,D6m=',reATp,D2m,D8m
          call PrintSet(KPO4(BoxNumberXY),'Nan in M11p')
        endif
!-------------------------------------------------------------------------

        call LimitShift(jK1N1p,N1p_Ben(BoxNumberXY)* Depth_Ben(BoxNumberXY), &
         K1p(BoxNumberXY)+(reBTp(BoxNumberXY)-ruBTp(BoxNumberXY))*LocalDelta, &
         max_change_per_step,r);

        call addbotflux(ANY,BoxNumberXY,iiBen,ppK1p,iiPel,ppN1p,jK1N1p*r)

        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Calculate diffusive flux at the oxic/denitrification interface:
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        !Be sure that concentration K1p >0 ,when flux from water column
        !is limited by small depths of water column.
        jK11K1p = CalculateFromSet( KPO4(BoxNumberXY), DERIVATIVE, RFLUX, &
          D1m(BoxNumberXY), DONE)-min(ZERO,jK1N1p*(DONE-r))
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Calculate new depth of the oxygen horizon
        ! and the flux of phosphate related to this shifting
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        Dnew  =   D1m(BoxNumberXY)+ shiftD1m(BoxNumberXY)* LocalDelta
        jK11K1p =jK11K1p+  CalculateFromSet( KSHIFT, SHIFT, LAYER1, &
          D1m(BoxNumberXY), Dnew)/ LocalDelta

        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! limit flux according to the actual phosphate content in the layer
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        call LimitShift(jK11K1p,max(ZERO,K1p(BoxNumberXY)+(reBTp(BoxNumberXY) &
            -ruBTp(BoxNumberXY)-jK1N1p*r)*LocalDelta),  &
            K11p(BoxNumberXY),max_change_per_step)

        call flux(BoxNumberXY, iiBen, ppK11p, ppK1p,  jK11K1p* insw( jK11K1p) )
        call flux(BoxNumberXY, iiBen, ppK1p, ppK11p,- jK11K1p* insw(-jK11K1p) )
        
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Calculate diffusive flux at the denitrification/anoxic interface:
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        ! flux according gradient
        jK21K11p=CalculateFromSet( KSHIFT, DERIVATIVE, RFLUX,  dnm, DONE)
        Dnew  =   dnm+ lShift
        ! flux according according shift of layers
        jK21K11p =jK21K11p+  CalculateFromSet( KSHIFT, SHIFT, LAYER3+mode, &
          dnm, Dnew)/ LocalDelta

        call LimitShift(jK21K11p,K11p(BoxNumberXY)-jK11K1p*LocalDelta, &
                                        K21p(BoxNumberXY),max_change_per_step)
        call flux(BoxNumberXY,iiBen,ppK21p,ppK11p, jK21K11p* insw( jK21K11p) )
        call flux(BoxNumberXY,iiBen,ppK11p,ppK21p,-jK21K11p* insw(-jK21K11p) )

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Calculate flux at the lower boundary
        ! 1=no flux, 2= only fluxes_downwards (sink), 3,=full_flux
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        select case (p_flux_at_deep_end)  
           case(1); jflux=ZERO
           case(2); jflux=min(ZERO,CalculateFromSet(KPO4(BoxNumberXY), &
                                          DERIVATIVE, RFLUX, p_d_tot_2, ZERO))
           case(3); jflux=CalculateFromSet(KPO4(BoxNumberXY),  &
                                          DERIVATIVE, RFLUX, p_d_tot_2, ZERO)
        end select
 
        call LimitChange(NEGATIVE,jflux,K21p(BoxNumberXY)-jK21K11p*LocalDelta, &
                                                        max_change_per_step)
        jK31K21p(BoxNumberXY) = jflux

        call flux(BoxNumberXY, iiBen, ppK21p, ppK21p, jK31K21p(BoxNumberXY) )

      end if

  end do

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenPhosphate
!
! DESCRIPTION
  subroutine BenPhosphateEquation(KPO4,mode,dry,N1p,K1p,K11p,K21p,D1m,D2m, &
              d_tot_1,d_tot_2,p_poro,p_p_ae,p_p_an,s_ads, &
              pShift,mShift,lShift,dnm, diff, alpha,beta,seBT,zuBT,zuD1,zuD2)
!
! !USES:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO,DONE,LOGUNIT,dummy
  use constants, ONLY: DIST_EXPONENTIAL_TERM, LAYERS,DERIVATIVE, &
    LAYER1,LAYER2,LAYER3,LAYER5,LAYER6 ,DIFFUSION, FOR_ALL_LAYERS, POROSITY, &
    ADSORPTION, LAYER4, DEFINE, LINEAR_TERM, CONSTANT_TERM, SET_CONTINUITY, &
    FLAG, MASS, SET_BOUNDARY, EQUATION, SET_LAYER_INTEGRAL, EXPONENTIAL_TERM, &
    SET_LAYER_INTEGRAL_UNTIL,ADD,INPUT_TERM,PARAMETER, STANDARD, MIN_VAL_EXPFUN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following bennut functions are used:InitializeSet, &
  ! DefineSet, CompleteSet, CalculateSet, CalculateTau, CalculateFromSet
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use bennut_interface, ONLY: InitializeSet, DefineSet, CompleteSet
  use mem_globalfun,   ONLY: ExpDist

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following sesame functions are used:IntegralExp, insw
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!  
!
! !AUTHORS
!   Original version by  P. Ruardij
!
!
! !REVISION_HISTORY
!       September 1999 by M. Vichi !               Commented version
!
! COPYING
!   
!   Copyright (C) 2006 P. Ruardij & M.Vichi
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
  IMPLICIT NONE
  integer,intent(INOUT)     :: KPO4
  integer,intent(IN)     :: mode
  logical,intent(IN)     :: dry
  real(RLEN),intent(IN)  :: pShift
  real(RLEN),intent(IN)  :: mShift
  real(RLEN),intent(IN)  :: lShift
  real(RLEN),intent(IN)  :: dnm
  real(RLEN),intent(IN)  :: alpha,beta
  real(RLEN),intent(IN)  :: diff
  real(RLEN),intent(IN)  :: seBT ! loss due to primary production
  real(RLEN),intent(IN)  :: zuBT
  real(RLEN),intent(IN)  :: zuD1
  real(RLEN),intent(IN)  :: zuD2
  real(RLEN),intent(IN)  :: N1p
  real(RLEN),intent(IN)  :: K1p
  real(RLEN),intent(IN)  :: K11p
  real(RLEN),intent(IN)  :: K21p
  real(RLEN),intent(IN)  :: D1m
  real(RLEN),intent(IN)  :: D2m
  real(RLEN),intent(IN)  :: d_tot_1
  real(RLEN),intent(IN)  :: d_tot_2
  real(RLEN),intent(IN)  :: p_poro
  real(RLEN),intent(IN)  :: p_p_ae
  real(RLEN),intent(IN)  :: p_p_an 
  real(RLEN),intent(IN)  :: s_ads 
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN)  :: r
  real(RLEN)  :: lambda
  integer     :: i
  real(RLEN)  :: dz 
  real(RLEN)  :: ds1,ds2,p_ads_3

      if ( mode <= 0) then
         p_ads_3=p_p_an; ds1=D2m;                         ds2=D2m
      else if (lShift.lt.ZERO) then
         p_ads_3=p_p_an; ds1=D2m+lShift;ds2=D2m
      else
         p_ads_3=p_p_ae; ds1=D2m;                         ds2=ds1 +lShift
      endif

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Subdivide anoxic layer in two sublayers only for the calculation.
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      dz  =  ( d_tot_1+ dnm )* 0.5D+00

      if ( mode .le.0) then
         KPO4 = InitializeSet( KPO4, 5, 14)
         call DefineSet( KPO4, LAYERS, LAYER1,  LAYER2,   D1m, D2m)
         call DefineSet( KPO4, LAYERS, LAYER3,  LAYER4,   dnm, dz)
      else
         KPO4 = InitializeSet( KPO4, 6, 18)
         call DefineSet( KPO4, LAYERS,  LAYER1,  LAYER2, D1m, ds1 )
         call DefineSet( KPO4, LAYERS,  LAYER3,  LAYER4, ds2, dnm)
         call DefineSet( KPO4, LAYERS,  LAYER5,  0     ,  dz,  dummy)
      endif


      call DefineSet( KPO4, DIFFUSION, FOR_ALL_LAYERS, 0, diff, dummy)

      call DefineSet( KPO4, POROSITY, FOR_ALL_LAYERS, 0, p_poro, dummy)

!-------------------------------------------------------------------------
      if ( p_p_ae< ZERO) then
        write(LOGUNIT,*)' Error:BenPhos mode=',mode
        write(LOGUNIT,*)'Error: BenPhos p_p`s=',p_p_ae,p_p_an,p_ads_3
      endif
!-------------------------------------------------------------------------
      call DefineSet( KPO4, ADSORPTION, LAYER1, LAYER2, p_p_ae, p_p_ae)
      call DefineSet( KPO4, ADSORPTION, LAYER3, 0, p_ads_3,dummy)
      if ( mode.le.0) then
         call DefineSet( KPO4, ADSORPTION, LAYER4, LAYER5, p_p_an, p_p_an)
      else
         call DefineSet( KPO4, ADSORPTION, LAYER4, LAYER5, p_p_an, p_p_an)
         call DefineSet( KPO4, ADSORPTION, LAYER6, 0, p_p_an, dummy)
      endif

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Define coefficients for the steady-state solutions in each layer
      ! General solution of the equilibrium profile:
      ! 1st layer:
      ! P(z) = p13*z^2 + p14*z + p15
      ! 2nd layer:
      ! P(z) = p21*exp[-alpha*(z-D1.m)] + p24*z + p25
      ! 3rd layer:
      ! P(z) = p31*exp[-alpha*(z-D2.m)] + p34*z + p35
      ! 4th layer:
      ! P(z) = p41*exp[-alpha*(z-dnm)] + p44*z + p45
      !    p44 = 0 (boundary condition)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      lambda=max(MIN_VAL_EXPFUN,sqrt((s_ads+seBT)/diff)); 
      call DefineSet( KPO4, DEFINE, 11, EXPONENTIAL_TERM, -lambda, dummy)
      call DefineSet( KPO4, DEFINE, 12, EXPONENTIAL_TERM,  lambda, dummy)
      call DefineSet( KPO4, DEFINE, 15, CONSTANT_TERM, dummy, dummy)

      call DefineSet( KPO4, DEFINE, 21, DIST_EXPONENTIAL_TERM, -alpha, dummy)
      call DefineSet( KPO4, DEFINE, 24, LINEAR_TERM, dummy, dummy)
      call DefineSet( KPO4, DEFINE, 25, CONSTANT_TERM, dummy, dummy)

      if ( mode.le.0) then
        i=30
      else
        call DefineSet( KPO4, DEFINE, 31, EXPONENTIAL_TERM, -lambda,  dummy)
        call DefineSet( KPO4, DEFINE, 32, EXPONENTIAL_TERM,  lambda, dummy)
        call DefineSet( KPO4, DEFINE, 34, LINEAR_TERM, dummy, dummy)
        call DefineSet( KPO4, DEFINE, 35, CONSTANT_TERM, dummy, dummy)
        i=40
      endif

      call DefineSet( KPO4, DEFINE, i+1, DIST_EXPONENTIAL_TERM, -beta, dummy)
      call DefineSet( KPO4, DEFINE, i+4, LINEAR_TERM, dummy, dummy)
      call DefineSet( KPO4, DEFINE, i+5, CONSTANT_TERM, dummy, dummy)

      call DefineSet( KPO4, DEFINE, i+11, DIST_EXPONENTIAL_TERM, -beta, dummy)
      call DefineSet( KPO4, DEFINE, i+14, LINEAR_TERM, dummy, dummy)
      call DefineSet( KPO4, DEFINE, i+15, CONSTANT_TERM, dummy, dummy)

      call DefineSet( KPO4, DEFINE, i+21, DIST_EXPONENTIAL_TERM, - beta, dummy)
      call DefineSet( KPO4, DEFINE, i+25, CONSTANT_TERM, dummy, dummy) 

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Insert other boundary conditions and continuity between layers:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      !1-8/1-10
      call CompleteSet( KPO4, SET_CONTINUITY, FLAG, MASS, dummy)

      !9/11
      if (dry) then
        call CompleteSet(KPO4,SET_BOUNDARY, LAYER1, DERIVATIVE,ZERO,value=ZERO)
      else
        call CompleteSet(KPO4,SET_BOUNDARY, LAYER1, EQUATION, ZERO,value=N1p)
      endif


      select case ( mode )
        case ( 0 )
          !10-11
          call CompleteSet(KPO4,SET_LAYER_INTEGRAL, LAYER2,LAYER3,  &
                                                           dummy, value=K11p)
          call CompleteSet( KPO4, SET_LAYER_INTEGRAL_UNTIL, LAYER4, LAYER5,  &
                                                           d_tot_2, value=K21p)

         !12: condition for 30/40 layer....
          call CompleteSet(KPO4,INPUT_TERM,21,PARAMETER,dummy,value=zuD2)
        case ( 1 )
          !12
          r=+s_ads * pShift/abs(lShift)/p_ads_3/p_poro + mShift/abs(lShift);
          call CompleteSet( KPO4, INPUT_TERM, 35, STANDARD, dummy,  &
                                                           value=r/s_ads)

          !13:  next two lines one boundary condition!
          call CompleteSet( KPO4, SET_LAYER_INTEGRAL,  LAYER2, LAYER2, &
                                                             dummy, value=ZERO)
          call CompleteSet( KPO4, ADD+SET_LAYER_INTEGRAL, LAYER4, LAYER4, &
                                                     dummy, value=K11p-pShift)

          !14
          call CompleteSet(KPO4,INPUT_TERM,21,PARAMETER,dummy,value=zuD2)
          !15
          call CompleteSet( KPO4, SET_LAYER_INTEGRAL_UNTIL, LAYER5, LAYER6, &
                                                          d_tot_2, value=K21p)

          r  =   ExpDist( - beta, dnm- ds2)
          call FixProportionCoeff(KPO4,41,51,DONE,r)
        case ( -1 )
          !10-11
          !The mineralization at D1m equal to the  mineralization at D1m under
          !assumuption of that the mineralization distribution in oxic layer is 
          !distributed according detritus distribution alpha 
          ! 12-13
          call CompleteSet( KPO4, INPUT_TERM, 21, PARAMETER, dummy, value=zuD1)
          call CompleteSet( KPO4, INPUT_TERM, 41, PARAMETER, dummy,&
                                             value=zuD2* ExpDist(-alpha,dnm-D2m))
          r  =   ExpDist( - alpha, D2m- D1m)
          call FixProportionCoeff(KPO4,21,31,DONE,r)
      end select

      !13/!17 : condition for fifth/sixth layer....
      r  =   ExpDist( - alpha, dz- dnm)
      call FixProportionCoeff(KPO4,i+11,i+21,DONE,r)

      r=+s_ads * K1p/p_p_ae/p_poro/D1m+zuBT
      call CompleteSet( KPO4, INPUT_TERM, 15, STANDARD, dummy, &
                                                      value=r/(s_ads+seBT))
  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
