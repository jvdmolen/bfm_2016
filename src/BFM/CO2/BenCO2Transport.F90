#include "INCLUDE.h"
#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenCO2Transport
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
  subroutine BenCO2TransportDynamics
!

#ifdef INCLUDE_BENCO2

! !USES:

  ! For the following Benthic-states fluxes are defined: G13c, G3c
  ! The following Benthic-states are used (NOT in fluxes): D1m, D6m, D2m
  ! The following global scalar vars are used: &
  !    NO_BOXES_XY,  &
  !   BoxNumberXY, InitializeModel, LocalDelta
  ! The following Benthic 1-d global boxvars are modified : KCO2, jbotO3c
  ! The following Benthic 1-d global boxvars got a value: DICae, DICan
  ! The following Benthic 1-d global boxvars are used:  &
  ! irrenh, ETW_Ben, rrATo, O3c_Ben, shiftD1m
  ! The following Benthic 1-d global boxpars  are used: p_poro
  ! The following 0-d global parameters are used: p_d_tot, p_q10diff
  ! The following global constants are used: RLEN,ZERO

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  USE BFM_ERROR_MSG, ONLY: set_warning_for_getm
  use global_mem, ONLY:RLEN,ZERO,LOGUNIT,DONE,dummy
  use mem,  ONLY: G23c,G13c, G3c, D1m, D6m, D2m, D2STATE
  use mem, ONLY: ppG23c, ppG13c, ppG3c, ppO3c, &
    NO_BOXES_XY, jG33G23c,LocalDelta, max_change_per_step,   &
    BoxNumberXY, InitializeModel, KCO2, DICae,Wind,pCO2ae,   &
    ruBTc,reBTc,reATc, irrenh, ETW_Ben,ESW_Ben,ERHO_Ben,O3c_Ben,  &
    DIC_Ben,CO3_Ben,HCO3_Ben,CO2_Ben,shiftD1m,shiftD2m,Depth_Ben, &
    iiPel,iiBen, flux,CO2an,HCO3an,CO3an,DICan,CO2ae,HCO3ae,CO3ae,DICae,dry_z

  use constants,ONLY:LAYERS,LAYER1,LAYER2, DIFFUSION, FOR_ALL_LAYERS, &
    POROSITY, ADSORPTION, DIST_EXPONENTIAL_TERM, DEFINE, QUADRATIC_TERM, &
    LINEAR_TERM, CONSTANT_TERM, SET_CONTINUITY, FLAG, MASS, SET_BOUNDARY, &
    EXPONENTIAL_TERM, EQUATION, INPUT_TERM, STANDARD,PARAMETER, &
    SET_LAYER_INTEGRAL_UNTIL, LAYER3, SET_LAYER_INTEGRAL, &
    ADD, DERIVATIVE, RFLUX, SHIFT, LAYER4,LAYER5,MW_C,SEC_PER_DAY,ANY
  use mem_Param,  ONLY: p_poro, p_d_tot,p_d_tot_2, p_q10diff,p_dry_ben, &
                            p_clDxm
  use mem_Diffusion,ONLY: p_diff=>p_diff_O3c
  use mem_BenCO2Transport
  use botflux,ONLY:addbotflux
  use LimitRates, ONLY:LimitShift,LimitChange
  use mem_CO2,only:pCO2_air

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following bennut functions are used:&
  ! InitializeSet, DefineSet, CompleteSet, CalculateSet, CalculateTau, &
  ! CalculateFromSet
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use bennut_interface, ONLY: InitializeSet, DefineSet, &
    CompleteSet, CalculateSet, CalculateTau, CalculateFromSet

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following global functions are used:eTq
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_interface,   ONLY: eTq

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following sesame functions are used:IntegralExp
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,ONLY: IntegralExpDist,ExpDist,SetExpDist,IntegralExp,insw

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
  real(RLEN)  :: r
  real(RLEN)  :: alpha,beta,lambda
  real(RLEN)  :: diff1,diff0
  real(RLEN)  :: jflux
  real(RLEN)  :: Tau
  real(RLEN)  :: cG3c
  real(RLEN)  :: zu
  real(RLEN)  :: zuD1
  real(RLEN)  :: zuD2
  real(RLEN)  :: jG13G3c
  real(RLEN)  :: jG23G13c
  real(RLEN)  :: jG3O3c
  real(RLEN)  :: Dnew
  real(RLEN)  :: Dx
  real(RLEN)  :: Dy
  real(RLEN)  :: sO3c
  real(RLEN)  :: mdiff1,mdiff0,mdiffp
  real(RLEN)  :: DIC13,DIC23
  logical     :: first_order

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  do BoxNumberXY=1,NO_BOXES_XY

    if  (p_diff==ZERO) then
     !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     ! Calc Diffusion at actual temperture according Zeebe,RE(2011)
     ! Zeebe,RE(2011),On the molecular diffusion coefficients of dissolved 
     ! CO2,HCO2-and  HCO3-- and there dependence on istopic mass
     ! Geochimica et Cosmochimica Acta 75,2483-2948
     !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      call  CalcDiffusionCO2(ppG3c,ETW_Ben(BoxNumberXY),CO2_Ben(BoxNumberXY), &
        HCO3_Ben(BoxNumberXY),CO3_Ben(BoxNumberXY),DIC_Ben(BoxNumberXY),mdiffp)
      call  CalcDiffusionCO2(ppG3c,ETW_Ben(BoxNumberXY),CO2ae(BoxNumberXY), &
        HCO3ae(BoxNumberXY), CO3ae(BoxNumberXY),DICae(BoxNumberXY),mdiff0)
      call  CalcDiffusionCO2(ppG3c,ETW_Ben(BoxNumberXY),CO2an(BoxNumberXY), &
          HCO3an(BoxNumberXY),CO3an(BoxNumberXY),DICan(BoxNumberXY),mdiff1)
      mdiff1=(mdiff0+mdiff1)*0.5
      mdiff0=(mdiffp+mdiff0)*0.5
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Correction for porosity and bioirrigation
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      diff0 = mdiff0*SEC_PER_DAY* irrenh(BoxNumberXY)* p_poro(BoxNumberXY)
      diff1 = mdiff1*SEC_PER_DAY* irrenh(BoxNumberXY)* p_poro(BoxNumberXY)
    else
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! (Classic) Temperature Correction (diffusion)
      ! Correction for porosity and bioirrigation
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      diff0=p_diff*SEC_PER_DAY* irrenh(BoxNumberXY)* p_poro(BoxNumberXY) *& 
                                      eTq( ETW_Ben(BoxNumberXY), p_q10diff)
      diff1=diff0
    endif

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Check on using first order or zero-order equation for DIC-uptake 
    ! used for primary production process
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    sO3c=ruBTc(BoxNumberXY)/G3c(BoxNumberXY)
    first_order=sO3c.gt.sqrt(diff0)

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Calculate the belonging concentrations for the state variables
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    DIC23 = G23c(BoxNumberXY)/MW_C/ p_poro(BoxNumberXY)/( &
      p_p+ DONE)/( p_d_tot_2- D2m(BoxNumberXY))
    DIC13 = G13c(BoxNumberXY)/MW_C/ p_poro(BoxNumberXY)/( &
      p_p+ DONE)/( D2m(BoxNumberXY)- D1m(BoxNumberXY))

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Recalculate Mineralization mmo/m2 --> mgC/m3 porewater
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    if ( first_order) then
      zu = reBTc(BoxNumberXY) / D1m(BoxNumberXY)/ p_poro(BoxNumberXY)
      lambda=sqrt(sO3c/diff0)
    else
      zu = (reBTc(BoxNumberXY)-ruBTc(BoxNumberXY)) &
                                   / D1m(BoxNumberXY)/ p_poro(BoxNumberXY)
    endif

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Assume Negative Exponential Distribution of Part.Carb. according D6.m
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    alpha  =   SetExpDist(D6m(BoxNumberXY),Dlm=p_clDxm,Dmm=p_d_tot)
    beta=alpha

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Recalculate Mineralization m2 --> m3 porewater
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ! Anoxic Mineralization at D1.m of Q6; (mgC /m3/d)

    zuD1 = ( reATc(BoxNumberXY))/ p_poro(BoxNumberXY)/ IntegralExpDist( &
      -alpha, p_d_tot_2- D1m(BoxNumberXY))
   
    zuD2  =   zuD1* ExpDist( - alpha, D2m(BoxNumberXY)- D1m(BoxNumberXY))

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Initialize and input physical boundaries and forcing:
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    KCO2(BoxNumberXY)  =   InitializeSet(  KCO2(BoxNumberXY),  5,  14)
    Dx=(D1m(BoxNumberXY) +D2m(BoxNumberXY)) * 0.5
    Dy=(D2m(BoxNumberXY) +p_d_tot) * 0.5

    call DefineSet(KCO2(BoxNumberXY),LAYERS,LAYER1,LAYER2, &
                                         D1m(BoxNumberXY),Dx)
    call DefineSet(KCO2(BoxNumberXY),LAYERS,LAYER3,LAYER4, &
                                         D2m(BoxNumberXY),Dy)

    call DefineSet(KCO2(BoxNumberXY),DIFFUSION, LAYER1, LAYER2, diff0,diff1)
    call DefineSet(KCO2(BoxNumberXY),DIFFUSION, LAYER3, LAYER4, diff1,diff1)
    call DefineSet(KCO2(BoxNumberXY),DIFFUSION, LAYER5, 0, diff1,dummy)

    call DefineSet(KCO2(BoxNumberXY),POROSITY, FOR_ALL_LAYERS, 0, &
                                         p_poro(BoxNumberXY), dummy)

    call DefineSet(KCO2(BoxNumberXY),ADSORPTION,FOR_ALL_LAYERS,0, p_p,dummy)

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Give particular solution for all layers:
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    if ( first_order) then
      call DefineSet(KCO2(BoxNumberXY),DEFINE,11, &
                                EXPONENTIAL_TERM, lambda,dummy)
      call DefineSet(KCO2(BoxNumberXY),DEFINE,12, &
                                EXPONENTIAL_TERM,-lambda,dummy)
      call DefineSet(KCO2(BoxNumberXY),DEFINE,15,CONSTANT_TERM, dummy,dummy)
    else
      call DefineSet(KCO2(BoxNumberXY),DEFINE,13,QUADRATIC_TERM, dummy,dummy)
      call DefineSet(KCO2(BoxNumberXY),DEFINE,14,LINEAR_TERM, dummy,dummy)
      call DefineSet(KCO2(BoxNumberXY),DEFINE,15,CONSTANT_TERM, dummy,dummy)
    endif

    call DefineSet( KCO2(BoxNumberXY), DEFINE,21,DIST_EXPONENTIAL_TERM, &
                                                            -alpha, dummy)
    call DefineSet( KCO2(BoxNumberXY), DEFINE,24,LINEAR_TERM,  dummy,dummy)
    call DefineSet( KCO2(BoxNumberXY), DEFINE,25,CONSTANT_TERM,dummy,dummy)

    call DefineSet( KCO2(BoxNumberXY), DEFINE,31,DIST_EXPONENTIAL_TERM, &
                                                           -alpha, dummy)
    call DefineSet( KCO2(BoxNumberXY), DEFINE,34,LINEAR_TERM, dummy,dummy)
    call DefineSet( KCO2(BoxNumberXY), DEFINE,35,CONSTANT_TERM, dummy,dummy)

    call DefineSet( KCO2(BoxNumberXY), DEFINE,41,DIST_EXPONENTIAL_TERM, &
                                                            -beta, dummy)
    call DefineSet( KCO2(BoxNumberXY), DEFINE,44,LINEAR_TERM, dummy,dummy)
    call DefineSet( KCO2(BoxNumberXY), DEFINE,45,CONSTANT_TERM, dummy,dummy)

    call DefineSet( KCO2(BoxNumberXY), DEFINE,51,DIST_EXPONENTIAL_TERM, &
                                                             -beta, dummy)
    call DefineSet( KCO2(BoxNumberXY), DEFINE,55,CONSTANT_TERM, dummy,dummy)

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Insert boundary conditions:
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    !1-8 boundary conditions:
    call CompleteSet( KCO2(BoxNumberXY), SET_CONTINUITY, FLAG, MASS, &
      dummy, dummy)

    !9th boundary condition:
    if ( dry_z(BoxNumberXY) <=p_clDxm .and. p_dry_ben) then
      call SurfaceCO2Diffusion(0,ETW_Ben(BoxNumberXY),ESW_Ben(BoxNumberXY),& 
         ERHO_Ben(BoxNumberXY),Wind,D1m(BoxNumberXY), dummy,r)
      r= -(pCO2_air-pCO2ae(BoxNumberXY))/D1m(BoxNumberXY)* &
                                        r/p_poro(BoxNumberXY)
      call CompleteSet( KCO2(BoxNumberXY), SET_BOUNDARY, LAYER1, &
                                         DERIVATIVE, ZERO, value=r)
    else
       call CompleteSet( KCO2(BoxNumberXY), SET_BOUNDARY, LAYER1, &
                        EQUATION, ZERO, value=O3c_Ben(BoxNumberXY))
    endif

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! a11 / (labda * labda * diff) = a12 / (labda * labda * diff)
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    !10th boundary condition:
    r  =   ExpDist( - alpha, Dx - D1m(BoxNumberXY))
    call FixProportionCoeff(KCO2(BoxNumberXY),21,31,DONE,r)

    !11-12 boundary condition:
    select case ( InitializeModel)
      case ( 0 )
        call CompleteSet(KCO2(BoxNumberXY), SET_LAYER_INTEGRAL, &
          LAYER2, LAYER3, dummy, G13c(BoxNumberXY))
        call CompleteSet(KCO2(BoxNumberXY), SET_LAYER_INTEGRAL_UNTIL, &
          LAYER4, LAYER5, p_d_tot_2, G23c(BoxNumberXY))
      case ( 1 )
        call CompleteSet(KCO2(BoxNumberXY), INPUT_TERM, 21,  &
                                          PARAMETER, dummy, zuD1)
        r  =   ExpDist( - alpha, D2m(BoxNumberXY)-Dx )
        call FixProportionCoeff(KCO2(BoxNumberXY),31,41,DONE,r)
    end select

    !13th boundary condition:
    r  =   ExpDist( - beta, Dy- D2m(BoxNumberXY))
    call FixProportionCoeff(KCO2(BoxNumberXY),41,51,DONE,r)

    !14th boundary condition:
    if (first_order) then
       call CompleteSet(KCO2(BoxNumberXY), INPUT_TERM, 15,  &
                                  STANDARD, dummy, zu/sO3c)
    else
       call CompleteSet(KCO2(BoxNumberXY), INPUT_TERM, 13,  &
                                          PARAMETER, dummy, zu)
    endif

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! 1. Calculate for above defined set of boundary conditions
    !      the gradient of the nutrient
    ! 2. Replace last condition by an alternative new one
    ! 3. Calculate the value belongin by the alternative condition with
    !  the gradient calculated under 1.
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    cG3c = CalculateSet( KCO2(BoxNumberXY), SET_LAYER_INTEGRAL, LAYER1, &
      LAYER1, dummy, ZERO)

    if ( InitializeModel== 0) then
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate adaptation time absolute and Delta() relative to Tau:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      Tau  =   CalculateTau(  ZERO,  diff0,  p_p,  D1m(BoxNumberXY))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Estimate the average value of K4n during the next time step:
      ! Hence this value depend as well as on adaptation time,
      ! ''old'' value, and on ''equilibrium value''
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      cG3c = cG3c+( G3c(BoxNumberXY)- cG3c)* IntegralExp( - LocalDelta/ &
        Tau, DONE)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Recalulate gradient, now using cG3c
      ! Complete first the alternative condition by inputting the value
      ! for G3c (cG3c)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      r = CalculateSet( KCO2(BoxNumberXY), ADD, 0, 0, dummy, cG3c)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! flux at D1.m
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      jG13G3c = CalculateFromSet( KCO2(BoxNumberXY), DERIVATIVE, &
        RFLUX, D1m(BoxNumberXY), dummy)

      Dnew  =   D1m(BoxNumberXY)+ LocalDelta* shiftD1m(BoxNumberXY)
      jG13G3c = jG13G3c+ CalculateFromSet( KCO2(BoxNumberXY), SHIFT, &
        LAYER1, D1m(BoxNumberXY), Dnew)/ LocalDelta

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Damp for too large fluxes
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      call LimitShift(jG13G3c,G3c(BoxNumberXY),G13c(BoxNumberXY), &
                                                         max_change_per_step)
      call flux(BoxNumberXY,iiBen,ppG13c, ppG3c,  jG13G3c* insw(  jG13G3c) )
      call flux(BoxNumberXY,iiBen,ppG3c, ppG13c,- jG13G3c* insw( -jG13G3c) )

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! flux at D2m
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      jG23G13c = CalculateFromSet( KCO2(BoxNumberXY), DERIVATIVE, RFLUX, &
        D2m(BoxNumberXY), dummy)

      Dnew  =   D2m(BoxNumberXY)+ LocalDelta* shiftD2m(BoxNumberXY)
      jG23G13c = jG23G13c+ CalculateFromSet( KCO2(BoxNumberXY), SHIFT, &
        LAYER3, D2m(BoxNumberXY), Dnew)/ LocalDelta

      jG23G13c = jG23G13c - zuD2* p_poro(BoxNumberXY)* &
                    IntegralExpDist(-alpha, p_d_tot_2- D2m(BoxNumberXY))

      call LimitShift(jG23G13c,G13c(BoxNumberXY),G23c(BoxNumberXY), &
                                                 max_change_per_step)

      call flux(BoxNumberXY,iiBen,ppG23c, ppG13c,  jG23G13c*insw( jG23G13c))
      call flux(BoxNumberXY,iiBen,ppG13c, ppG23c,- jG23G13c*insw(-jG23G13c))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! flux at underside
      ! 0=no flux, 1= only fluxes_downwards (sink), 2,=full_flux
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       select case (p_flux_at_deep_end)  
         case(0); jflux=ZERO
         case(1); jflux=min(ZERO, CalculateFromSet(KCO2(BoxNumberXY),  &
                                 DERIVATIVE, RFLUX, p_d_tot_2, ZERO))
         case(2); jflux=CalculateFromSet(KCO2(BoxNumberXY),  &
                                  DERIVATIVE, RFLUX, p_d_tot_2, ZERO)
      end select

      jG33G23c(BoxNumberXY) = jflux
      call LimitChange(ANY,jG33G23c(BoxNumberXY),G23c(BoxNumberXY), &
                                                max_change_per_step)
      call flux(BoxNumberXY,iiBen, ppG23c, ppG23c, jG33G23c(BoxNumberXY) )

      jG3O3c = CalculateFromSet( KCO2(BoxNumberXY), DERIVATIVE, &
        RFLUX, ZERO, ZERO)

      call LimitShift(jG3O3c,O3c_Ben(BoxNumberXY)*Depth_Ben(BoxNumberXY) ,&
                                      G3c(BoxNumberXY),max_change_per_step)

        call addbotflux(ANY,BoxNumberXY,iiBen,ppG3c,iiPel,ppO3c,jG3O3c)
    end if

end do
end 
subroutine CalcDiffusionCO2(mode,ETW,CO2,HCO3,CO3,DIC,diff)
  use global_mem, ONLY:RLEN,LOGUNIT,DONE,ZERO_KELVIN
  use constants,only:SEC_PER_DAY
  use mem,only:ppG3c,ppG3h
  implicit none
  integer,intent(IN)     :: mode
  real(RLEN),intent(IN)  :: ETW,CO2,HCO3,CO3,DIC
  real(RLEN),intent(OUT) :: diff

  select case (mode)
    case (ppG3c)
     diff=( CO2*14.6836D-9*((ETW-ZERO_KELVIN)/217.2056D+00-DONE)**1.9970D+00 &
       +   HCO3*7.0158D-9*((ETW-ZERO_KELVIN)/204.0282D+00-DONE)**2.3974D+00 &
       +   CO3*5.4468D-9*((ETW-ZERO_KELVIN)/210.2626D+00-DONE)**2.1929D+00)/DIC
    case (ppG3h)
     diff= 7.0158D-9*((ETW-ZERO_KELVIN)/204.0282D+00-DONE)**2.3974D+00 
   end select

end

#endif

 
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
