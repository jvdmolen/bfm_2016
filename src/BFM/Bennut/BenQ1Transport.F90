#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! !ROUTINE: BenQ1Transport
!
! DESCRIPTION
!   Description of the DOM dynamics in the sediments
!     icycle=0: old way: only C is calculated
!               Assumed is that the distriubtion og N and P are fully coupled
!               to the one of N ( implcitly assuming that internal quoato i
!               between C and nutrients are constant                    
!    icycle=1: separate calculation for the 3 nutreint cucles in R1c
!             For every cycle ( C,N,P) routine is called.
!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine BenQ1TransportCalculate(KQ,mode,R1_Benx,p_qnUc) 
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO,NZERO,DONE,LOGUNIT,dummy
  use mem,only: 
  use mem,  ONLY: D1m, D6m, D2m,DH2m,DH3m
  use mem, ONLY: NO_BOXES_XY, BoxNumberXY,BoxNumber,LocalDelta, &
    max_change_per_step,irrenh, ETW_Ben,&
    iiQu,iiQ1u,iiQ2u,iiQ1,iiQ11,iiQ21,iiC,iiN,dry_z, &
    Source_D2_vector,ppBenUrea,&
    ppBenLabileDetritus, BenUrea,BenLabileDetritus,&
    iiConsumption,iiProduction,PelBoxAbove
  use constants, ONLY: LAYERS, LAYER1, LAYER2,LAYER3,POSITIVE,&
    DIFFUSION, FOR_ALL_LAYERS, POROSITY, ADSORPTION, DEFINE,EXPONENTIAL_TERM, &
    CONSTANT_TERM, ZERO_EXPONENTIAL_TERM, LINEAR_TERM, SET_CONTINUITY, FLAG, &
    MASS, SET_BOUNDARY, DERIVATIVE, SET_LAYER_INTEGRAL, INPUT_TERM, &
    STANDARD, SET_LAYER_INTEGRAL_UNTIL, ADD, EQUATION
  use mem_globalfun, ONLY: IntegralExpDist,ExpDist,SetExpDist,IntegralExp
  use mem_Param,ONLY:p_poro,p_q10diff, p_clDxm, p_d_tot_2,p_d_tot,combine_anabac
  use mem_BenQ1Transport,ONLY: p_slQ1,p_shQ1,p_p, &
                       jQIBTx,jBTQIx,Q1x,Q11x,Q21x
  use mem_BenQ1Transport,ONLY: constituent,iiLOC,iiUrea
!JM is in mem:  use mem_BenthicNutrient3, ONLY:max_change_per_step
  use LimitRates, ONLY:LimitChange
  use mem_Param,only:p_dry_ben !JM is in mem: ,iiN
  use mem_Diffusion,ONLY: p_diff_urea,p_diff_Q1

! use mem,  ONLY: Output2d_1,Output2d_2,Output2d_3, Output2d_4

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following bennut functions are used:InitializeSet, &
  ! DefineSet, CompleteSet, CalculateSet, CalculateTau
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use bennut_interface, ONLY: InitializeSet, DefineSet, CompleteSet, &
    CalculateSet, CalculateTau

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following global functions are used:eTq
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_interface,   ONLY: eTq

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following sesame functions are used:IntegralExp, insw
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun, ONLY: insw
!
! !AUTHORS
!   Original version by  P. Ruardij
!
! !REVISION_HISTORY
!   by Piet Ruardij  *:0 at Sun Dec 04 23:09:55 CET 2005
!	
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
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE
  integer,dimension(NO_BOXES_XY),intent(INOUT)         ::KQ
  integer,intent(IN)                                   ::mode
  real(RLEN),dimension(NO_BOXES_XY),intent(IN),optional::R1_Benx
  real(RLEN),intent(IN),optional                       ::p_qnUc
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  logical     :: dry
  integer     :: ppQ1x,ppQ11x,ppQ21x
  integer     :: ppQ1u,ppQ11u
  real(RLEN)  :: M
  real(RLEN)  :: a15
  real(RLEN)  :: s,r
  real(RLEN)  :: alpha,beta
  real(RLEN)  :: gamma
  real(RLEN)  :: diff,p_diff
  real(RLEN)  :: Tau
  real(RLEN)  :: cQ1c
  real(RLEN)  :: sQ1
  real(RLEN)  :: zuD1
  real(RLEN)  :: rQ11

  external    :: FixProportionCoeff
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! icycle==0 : old model: no exchange with pelagic 
  !           : (Q1n and Q1p are derived form Q1c gradient
  ! reQix     : sink from Q1x ( loss of Q1x == uptake of Q1x by biota)
  ! ruQix     : source for Q1x ( input of Q1x ==excretionof Q1x by biota)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  if (mode.eq.iiN)  then
    constituent=iiUrea
    Q1x=BenUrea(iiQu,iiN);
    Q11x=BenUrea(iiQ1u,iiN);
    ppQ1x=ppBenUrea(iiQu,iiN);
    ppQ11x=ppBenUrea(iiQ1u,iiN);
    p_diff=p_diff_urea
  elseif (mode.eq.0.or.mode.eq.iiC) then
    constituent=iiLOC
    Q1x=BenLabileDetritus(iiQ1,iiC)-BenUrea(iiQu,iiN)/p_qnUc;
    Q11x=BenLabileDetritus(iiQ11,iiC)-BenUrea(iiQ1u,iiN)/p_qnUc;
    ppQ1x=ppBenLabileDetritus(iiQ1,iiC);
    ppQ11x=ppBenLabileDetritus(iiQ11,iiC);
    ppQ1u=ppBenUrea(iiQu,iiN);
    ppQ11u=ppBenUrea(iiQ1u,iiN);
    p_diff=p_diff_Q1
  endif
  jQIBTx(iiQ1,:)=Source_D2_vector(ppQ1x,iiConsumption)
  if (mode.ne.iiN)  &
         jQIBTx(iiQ1,:)=jQIBTx(iiQ1,:)-Source_D2_vector(ppQ1u,iiConsumption)
  jBTQIx(iiQ1,:)=Source_D2_vector(ppQ1x,iiProduction)
  if (mode.ne.iiN)  &
          jBTQIx(iiQ1,:)=jBTQIx(iiQ1,:)-Source_D2_vector(ppQ1u,iiProduction)
  jQIBTx(iiQ11,:)=Source_D2_vector(ppQ11x,iiConsumption)
  if (mode.ne.iiN)  &
          jQIBTx(iiQ11,:)=jQIBTx(iiQ11,:)-Source_D2_vector(ppQ11x,iiConsumption)
  jBTQIx(iiQ11,:)=Source_D2_vector(ppQ11x,iiProduction)
  if (mode.ne.iiN)  &
          jBTQIx(iiQ11,:)=jBTQIx(iiQ11,:)-Source_D2_vector(ppQ11x,iiProduction)
  if ( .not.combine_anabac) then
    if (mode.eq.iiN) then
      Q21x=BenUrea(iiQ2u,iiN);
      ppQ21x=ppBenUrea(iiQ2u,iiN);
    elseif (mode.eq.0.or.mode.eq.iiC) then
      Q21x=BenLabileDetritus(iiQ21,iiC);
      ppQ21x=ppBenLabileDetritus(iiQ21,iiC);
    endif
    jQIBTx(iiQ21,:)=Source_D2_vector(ppQ21x,iiConsumption)
    jBTQIx(iiQ21,:)=Source_D2_vector(ppQ21x,iiProduction)
  else
    jQIBTx(iiQ21,:)=ZERO; jBTQIx(iiQ21,:)=ZERO
  endif

  do BoxNumberXY=1,NO_BOXES_XY
      BoxNumber=PelBoxAbove(BoxNumberXY)
      dry=.false.;
      if (p_dry_ben) dry=( dry_z(BoxNumberXY) <=p_clDxm )

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Get total Net Benthic DOC (Q1.c)
      ! production/consumption in the oxic layer (m2 --> m3 porewater)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      call LimitChange(POSITIVE,jQIBTx(iiQ1, BoxNumberXY),Q1x(BoxNumberXY), &
                                                       max_change_per_step,r)
      sQ1 = min(p_shQ1, &
        max( p_slQ1, r*jQIBTx(iiQ1, BoxNumberXY)/(NZERO+Q1x(BoxNumberXY)) ))
      M  =   jBTQIx(iiQ1,BoxNumberXY)/D1m(BoxNumberXY)/ p_poro(BoxNumberXY)
      a15  =   M/ sQ1

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Correction due to environmental regulating factors,
      ! diffusion coefficient: temperature and bioirrigation
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      diff = p_diff* irrenh(BoxNumberXY)* eTq( ETW_Ben(BoxNumberXY), p_q10diff)
      gamma  =   sqrt(  sQ1/ diff)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Estimate specific bacterial consumption rate in anoxic layers (limited)
      ! Calculate coefficient for the exponential terms of the solution
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      if ( combine_anabac) then
        alpha  =   SetExpDist(D6m(BoxNumberXY),Dlm=p_clDxm,Dmm=p_d_tot)
        beta=alpha
      else
        alpha  =   SetExpDist(DH2m(BoxNumberXY),Dlm=p_clDxm,Dmm=p_d_tot)
        beta  =    SetExpDist(DH3m(BoxNumberXY),Dlm=p_clDxm,Dmm=p_d_tot)
      endif

      rQ11 = ( jBTQIx(iiQ11,BoxNumberXY)- jQIBTx(iiQ11,BoxNumberXY))/ &
                  p_poro(BoxNumberXY)/( p_d_tot_2- D1m(BoxNumberXY))
      zuD1 = sign(rQ11,max( 1.D-20, abs(rQ11)))/ &
                      IntegralExpDist(-alpha,p_d_tot_2-D1m(BoxNumberXY))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Initialize the set of differential equations giving:
      ! - n. of layers;
      ! - n. of coefficients
      ! - layers depths
      ! - environmental conditions (diffusion, porosity and adsorption coeff.)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! in the case of combine_anabaxc we assume that the production of urea
      ! is done only by deposit feeders. Therefor we assume that in the anoxic 
      ! nearly no production will be of urea and therefor it is assumed that the
      ! the concentration is constant.
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      !
      if (combine_anabac ) then
        KQ(BoxNumberXY) = InitializeSet( KQ(BoxNumberXY), 3, 7)
      else
        KQ(BoxNumberXY) = InitializeSet( KQ(BoxNumberXY), 3, 8)
      endif

      call DefineSet(KQ(BoxNumberXY),LAYERS, LAYER1, &
                                 LAYER2, D1m(BoxNumberXY), D2m(BoxNumberXY))
      call DefineSet(KQ(BoxNumberXY),DIFFUSION, FOR_ALL_LAYERS,0, diff, dummy)
      call DefineSet(KQ(BoxNumberXY),POROSITY, FOR_ALL_LAYERS, 0, &
                                                    p_poro(BoxNumberXY), dummy)
      call DefineSet(KQ(BoxNumberXY),ADSORPTION, FOR_ALL_LAYERS,0, p_p, dummy)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Define coefficients for the steady-state solutions in each layer
      ! General solution of the equilibrium profile:
      ! 1st layer:
      ! Q(z) = q13*z^2 + q14*z + q15
      ! 2nd layer:
      ! Q(z) = q21*exp(gamma*z) + q22*exp(-gamma*z)
      ! 3rd layer:
      ! Q(z) = q31*exp(gamma*z) + q32*exp(-gamma*z) `
      !    q32 = 0 (boundary condition)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      call DefineSet(KQ(BoxNumberXY),DEFINE, 11,EXPONENTIAL_TERM, gamma,dummy)
      call DefineSet(KQ(BoxNumberXY),DEFINE, 12,EXPONENTIAL_TERM,-gamma,dummy)
      call DefineSet(KQ(BoxNumberXY),DEFINE, 15, CONSTANT_TERM,   dummy,dummy)

      call DefineSet(KQ(BoxNumberXY),DEFINE, 21,ZERO_EXPONENTIAL_TERM, &
                                                             -alpha,dummy)
      call DefineSet(KQ(BoxNumberXY),DEFINE, 24,LINEAR_TERM,    dummy,dummy)
      call DefineSet(KQ(BoxNumberXY),DEFINE, 25,CONSTANT_TERM,  dummy, dummy)

      if (.not.combine_anabac ) &
         call DefineSet(KQ(BoxNumberXY),DEFINE, 31,ZERO_EXPONENTIAL_TERM, &
                                                              -beta,dummy)
      call DefineSet(KQ(BoxNumberXY),DEFINE, 35,CONSTANT_TERM, dummy,dummy)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Insert other boundary conditions and continuity between layers:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      !4
      call CompleteSet( KQ(BoxNumberXY),SET_CONTINUITY, FLAG, MASS, dummy)
      !5
      if ( constituent==iiLOC .or.dry ) then
        call CompleteSet( KQ(BoxNumberXY),SET_BOUNDARY, LAYER1, DERIVATIVE, &
                                                       ZERO, value=ZERO)
      else 
        call CompleteSet( KQ(BoxNumberXY),SET_BOUNDARY, LAYER1, EQUATION, &
                                      ZERO, value=R1_Benx(BoxNumberXY))
      endif
   
      !6:
      if (combine_anabac ) then
         r  =   ExpDist( - alpha, D2m(BoxNumberXY)- D1m(BoxNumberXY))
!        call FixProportionCoeff(KQ(BoxNumberXY),21,31,DONE,r)
         call CompleteSet( KQ(BoxNumberXY),SET_LAYER_INTEGRAL_UNTIL, LAYER2, &
            LAYER3, p_d_tot_2, value=Q11x(BoxNumberXY))
      else
         call CompleteSet( KQ(BoxNumberXY),SET_LAYER_INTEGRAL, LAYER2, &
            LAYER2, dummy, value=Q11x(BoxNumberXY))
        call CompleteSet( KQ(BoxNumberXY),SET_LAYER_INTEGRAL_UNTIL, LAYER3, &
        LAYER3, p_d_tot_2, value=Q21x(BoxNumberXY))
      endif
      !8:
      if (mode.eq.0 ) then
         call CompleteSet( KQ(BoxNumberXY),SET_LAYER_INTEGRAL, LAYER1, &
            LAYER1, dummy, value=Q1x(BoxNumberXY))
           cQ1c = CalculateSet( KQ(BoxNumberXY), 0, 0, 0, dummy, dummy)
      else
        call CompleteSet( KQ(BoxNumberXY),INPUT_TERM, 15, STANDARD, dummy, &
        value=a15)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Calculate for the above defined set of boundary conditions
        ! the steady-state profiles and return the vertically integrated
        ! concentration
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        cQ1c = CalculateSet( KQ(BoxNumberXY), SET_LAYER_INTEGRAL, &
          LAYER1, LAYER1, dummy, ZERO)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Calculate the adaptation time to the steady-state profile
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        Tau  =   CalculateTau(  sQ1,  diff,  p_p,  D1m(BoxNumberXY))

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Estimate the average value of Q11 over the actual time step
        ! (transient value).
        ! This value depends on the adaptation time, the actual time step,
        ! the ''old'' value and the ''equilibrium value''
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        cQ1c=cQ1c+( Q1x(BoxNumberXY)- cQ1c)*IntegralExp(-LocalDelta/Tau, DONE)

        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! solution as for the steady-state case and using cQ11c as new 
        ! constraint.
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        r  =   CalculateSet(  KQ(BoxNumberXY),  ADD,  0,  0,  dummy,  cQ1c)
      endif
  enddo
  end subroutine BenQ1TransportCalculate

  subroutine BenQ1TransportDynamics(KQ,mode,ppR1x,R1_Benx,p_qnuc) 
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO,NZERO,DONE,LOGUNIT,dummy
  use mem,  ONLY: D3STATE,Q21n, Q21p, Q11n, Q1n, Q11p, Q1p, &
           D1m, D2m
  use mem, ONLY:ppQ21n,ppQ21p, ppQ11n,ppQ11p, ppQ1n,ppQ1p, &
    NO_BOXES_XY, BoxNumberXY,BoxNumber,LocalDelta,max_change_per_step, &
    shiftD1m,shiftD2m,ppBenUrea,ppBenLabileDetritus, Depth_Ben,   &
    iiQu,iiQ1u,iiQ2u,iiQ1,iiQ11,iiQ21,iiC,iiN,iiPel,iiBen,flux,dry_z, &
    PelBoxAbove
  use constants, ONLY: LAYER1, LAYER2,ANY,POSITIVE,NEGATIVE, &
    DERIVATIVE, SHIFT, RFLUX
  use mem_BenQ1Transport,ONLY: jQIBTx,jBTQIx,Q1x,Q11x,Q21x
  use mem_BenQ1Transport,ONLY: constituent,iiLOC,iiUrea
  use mem_Param,  ONLY: p_clDxm, p_d_tot_2,combine_anabac
  use botflux,only:addbotflux
!JM is in mem:  use mem_BenthicNutrient3, ONLY:max_change_per_step
  use LimitRates, ONLY:LimitShift,LimitChange
  use mem_Param,only:p_dry_ben
! use mem,  ONLY: Output2d_1,Output2d_2,Output2d_3, Output2d_4
  use mem_globalfun, ONLY: insw
  use bennut_interface, ONLY: CalculateSet, CalculateFromSet

  IMPLICIT NONE
  integer,dimension(NO_BOXES_XY),intent(IN)            ::KQ
  integer,intent(IN)                                   ::mode
  integer,intent(IN),optional                          ::ppR1x
  real(RLEN),dimension(NO_BOXES_XY),intent(IN),optional::R1_Benx
  real(RLEN),intent(IN),optional                       ::p_qnUc
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  logical     :: dry
  integer     :: ppQ1x,ppQ11x,ppQ21x,iout
  real(RLEN)  :: s,r,corr
  real(RLEN)  :: jQ11Q1x,jQ21Q11x
  real(RLEN)  :: jQ1Q11x,jQ11Q21x
  real(RLEN)  :: jQ1R1x
  real(RLEN)  :: flow
  real(RLEN)  :: Dnew

  if (mode.eq.iiN)  then
    ppQ1x=ppBenUrea(iiQu,iiN);           ppQ11x=ppBenUrea(iiQ1u,iiN);
  elseif (mode.eq.0.or.mode.eq.iiC) then
    ppQ1x=ppBenLabileDetritus(iiQ1,iiC); ppQ11x=ppBenLabileDetritus(iiQ11,iiC);
  endif
  if ( .not.combine_anabac) then
    if (mode.eq.iiN) then
      ppQ21x=ppBenUrea(iiQ2u,iiN);
    elseif (mode.eq.0.or.mode.eq.iiC) then
      ppQ21x=ppBenLabileDetritus(iiQ21,iiC);
    endif
  endif

  corr=DONE;if (present(p_qnUc)) corr=DONE/p_qnuc

  do BoxNumberXY=1,NO_BOXES_XY
    BoxNumber=PelBoxAbove(BoxNumberXY)
    dry=.false.; if (p_dry_ben) dry=( dry_z(BoxNumberXY) <=p_clDxm )

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Start calculation of fluxes:
    ! Calculate new depth of the oxygen horizon and the flux related
    ! to the shifting. Add the flux at D1.m
    !-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    !-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Fluxes at D1m
    !-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    Dnew  =   D1m(BoxNumberXY)+ LocalDelta* shiftD1m(BoxNumberXY)
    ! unit mmolN/m2/d
    flow = CalculateFromSet( KQ(BoxNumberXY), SHIFT, LAYER1, &
    D1m(BoxNumberXY), Dnew)/ LocalDelta

    ! unit mmolN/m2/d
    flow = flow+ CalculateFromSet( KQ(BoxNumberXY), DERIVATIVE, &
    RFLUX, D1m(BoxNumberXY), dummy)

    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Limit for too large fluxes
    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    !Q1x: mode=iiN : unit mmolN/m2  
    !flow mode=iiN :unit mmolN/m2/d  
    call LimitShift(flow,Q1x(BoxNumberXY),Q11x(BoxNumberXY),&
                                                 max_change_per_step)
    !flow mode=iiN :unit mmolN/m2/d  mode=iiC : unit mgC/m2/d
    flow=flow*corr
    jQ1Q11x=-flow*insw(-flow)
    jQ11Q1x= flow*insw( flow)
!   Because for Q1c the fluxes are defind twice, one for R1c and one for 
!   the Cpart in Qunc=Urea
!   we need another flux definition. een flux from->to can only defines one time!   in the code.
    call flux(BoxNumberXY, iiBen, ppQ1x, ppQ1x, flow)
    call flux(BoxNumberXY, iiBen, ppQ11x, ppQ11x,-flow)


    if (isnan(flow)) then
        write(LOGUNIT,*) 'flow',flow
        write(LOGUNIT,*) 'jBTQIx',jBTQIx(iiQ1,:)
        write(LOGUNIT,*) 'Q1x',Q1x(:)
    endif

    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! One of the 2 fluxes between the some constituents is 0!
    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    !avoid to large fluxes leading to negative concentrations

    if (constituent==iiUrea) then
      jQ1R1x=ZERO; 
      if (.not.dry) then
        !jQ1R1x mode=iiN :unit mmolN/m2/d  
        jQ1R1x = CalculateFromSet( KQ(BoxNumberXY), &
                                     DERIVATIVE, RFLUX, ZERO, ZERO)
        !s: mode=iiN : unit mmolN/m2  
        s=R1_BenX(BoxNumberXY)
        !r: mode=iiN : unit mmolN/m2  
        !Q1x: mode=iiN : unit mmolN/m2  
        !jBTQIx,jQIBTx  unit mmolN/m2/d
        r= Q1x(BoxNumberXY)+ &
        (jBTQIx(iiQ1,BoxNumberXY)- jQIBTx(iiQ1,BoxNumberXY))*LocalDelta
        call LimitShift(jQ1R1x,s*Depth_Ben(BoxNumberXY),r,max_change_per_step)
        if ( isnan(jQ1R1x)) then
            call PrintSet(KQ(BoxNumberXY),"Nan flux for urea:jQ1R1x");
        endif
      endif
      flow=jQ1R1x*corr
      call addbotflux(ANY,BoxNumberXY,iiBen,ppQ1x,iiPel,ppR1x,flow)
    else
      r=jQ1Q11x/ Q1x(BoxNumberXY)*  Q1n(BoxNumberXY)  &
         -jQ11Q1x/ Q11x(BoxNumberXY)* Q11n(BoxNumberXY) 
      call LimitChange(POSITIVE,r,Q1n(BoxNumberXY),max_change_per_step)
      call LimitChange(NEGATIVE,r,Q11n(BoxNumberXY),max_change_per_step)
      call flux(BoxNumberXY,iiBen, ppQ1n, ppQ11n, r *insw(r))
      call flux(BoxNumberXY,iiBen, ppQ11n, ppQ1n,-r *insw(-r))

      r=jQ1Q11x/ Q1x(BoxNumberXY)*  Q1p(BoxNumberXY)  &
         -jQ11Q1x/ Q11x(BoxNumberXY)* Q11p(BoxNumberXY) 
      call LimitChange(POSITIVE,r,Q1p(BoxNumberXY),max_change_per_step)
      call LimitChange(NEGATIVE,r,Q11p(BoxNumberXY),max_change_per_step)
      call flux(BoxNumberXY,iiBen, ppQ1p, ppQ11p, r *insw(r))
      call flux(BoxNumberXY,iiBen, ppQ11p, ppQ1p,-r *insw(-r))

      if ( .not.  combine_anabac) then
        Dnew  =   D2m(BoxNumberXY)+ LocalDelta* shiftD2m(BoxNumberXY)
        flow = CalculateFromSet( KQ(BoxNumberXY), SHIFT, LAYER2, &
          D2m(BoxNumberXY), Dnew)/ LocalDelta

        flow = flow+ CalculateFromSet( KQ(BoxNumberXY), DERIVATIVE, &
          RFLUX, D2m(BoxNumberXY), dummy)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Limit for too large fluxes
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        call LimitShift(flow,Q11x(BoxNumberXY),Q21x(BoxNumberXY),&
                                                 max_change_per_step)

        jQ11Q21x=-flow*insw(-flow)
        jQ21Q11x= flow*insw( flow)
        call flux(BoxNumberXY, iiBen, ppQ21x, ppQ11x, jQ21Q11x)
        call flux(BoxNumberXY, iiBen, ppQ11x, ppQ21x, jQ11Q21x)  

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-
        ! One of the 2 fluxes between the some constituents is 0!
        !-==-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        r=jQ11Q21x/ Q11x(BoxNumberXY)*  Q11n(BoxNumberXY)  &
           -jQ21Q11x/ Q21x(BoxNumberXY)* Q21n(BoxNumberXY) 
        call LimitChange(POSITIVE,r,Q11n(BoxNumberXY),max_change_per_step)
        call LimitChange(NEGATIVE,r,Q21n(BoxNumberXY),max_change_per_step)
        call flux(BoxNumberXY,iiBen, ppQ11n, ppQ21n, r *insw(r))
        call flux(BoxNumberXY,iiBen, ppQ21n, ppQ11n,-r *insw(-r))

        r=jQ11Q21x/ Q11x(BoxNumberXY)*  Q11p(BoxNumberXY)  &
         -jQ21Q11x/ Q21x(BoxNumberXY)* Q21p(BoxNumberXY) 
        call LimitChange(POSITIVE,r,Q11p(BoxNumberXY),max_change_per_step)
        call LimitChange(NEGATIVE,r,Q21p(BoxNumberXY),max_change_per_step)
        call flux(BoxNumberXY,iiBen, ppQ11p, ppQ21p, r *insw(r))
        call flux(BoxNumberXY,iiBen, ppQ21p, ppQ11p,-r *insw(-r))
       endif
     endif
   end do
  end subroutine BenQ1TransportDynamics

  subroutine GetPenetrationDepthQ1(mode,nr,l,th,c,out,ierr)
  !get depth where above 50%  of Q11 of Q12 is found
  !This routine is called in BenBac.
  use global_mem, ONLY:RLEN,ZERO
  use mem,only:NO_BOXES_XY
  use mem_BenQ1Transport,only: NUTR,layer,rto,r50
  use mem_globalfun,   ONLY: rtsafe
  implicit none
  integer,INTENT(IN)        :: mode            !1= recalulate profile
  integer,INTENT(IN)        :: nr(NO_BOXES_XY) ! profile sequence number
  integer,INTENT(IN)        :: l               ! layer nr (+LAYER2 or LAYER3)
  REAL(RLEN),INTENT(IN)     :: th(NO_BOXES_XY)! thickness of layer
  REAL(RLEN),INTENT(IN)     :: c(NO_BOXES_XY)  ! conc Q11/Q12 in layer(mgC/m3)
  REAL(RLEN),INTENT(INOUT)  :: out(NO_BOXES_XY)! output depth
  integer,INTENT(OUT)       :: ierr

  integer              ::i,j
  integer     :: kq_zeros(NO_BOXES_XY)    !JM added
  external    ::rtsafe_funcd

!JM added
  interface
    subroutine BenQ1TransportDynamics(KQ,mode,ppR1x,R1_Benx,p_qnuc)
      use global_mem, ONLY:RLEN
      use mem, ONLY:NO_BOXES_XY
      integer,dimension(NO_BOXES_XY),intent(IN)            ::KQ
      integer,intent(IN)                                   ::mode
      integer,intent(IN),optional                          ::ppR1x
      real(RLEN),dimension(NO_BOXES_XY),intent(IN),optional::R1_Benx
      real(RLEN),intent(IN),optional                       ::p_qnUc
    end subroutine
  end interface

  ierr=0
  do j=1,NO_BOXES_XY
    if ( nr(j).gt.0) then
!JM      if (mode.eq.1) call BenQ1TransportDynamics(0,0)
      kq_zeros=0    !JM addded
      if (mode.eq.1) call BenQ1TransportDynamics(kq_zeros,0)
      NUTR=nr(j)
      layer=l
      rto=th(j)
      r50=c(j)*0.5
      i=rtsafe(rtsafe_funcd,ZERO,rto,1.0D-6,out(j))
      ierr=ierr+i
    else
      ierr=1
    endif
  enddo
  return
  end

  subroutine rtsafe_funcd(x,f,df)
  use global_mem, ONLY:RLEN,dummy
  use constants,  ONLY: EQUATION, STANDARD,INTEGRAL,MASS
  use bennut_interface,   ONLY: CalculateFromLayer
  use mem_BenQ1Transport,only: NUTR,layer,rto,r50
  implicit none
  REAL(RLEN),INTENT(IN)   ::x   !x
  REAL(RLEN),INTENT(OUT)  ::f   !values equation at x  
  REAL(RLEN),INTENT(OUT)  ::df  !value fist derivative at x

  f=CalculateFromLayer(NUTR,layer,INTEGRAL,MASS,x,rto)-r50
  df=CalculateFromLayer(NUTR,layer,EQUATION,STANDARD,x,dummy)
  return
  end  
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
