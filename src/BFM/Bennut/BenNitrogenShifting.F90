#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenNitrogenShifting
!
! DESCRIPTION
!   Description of shifting of dissolved N (amm and nitrate) between
!       layers
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine BenNitrogenShiftingDynamics
!
! !USES:

  ! For the following Benthic-states fluxes are defined: K14n, K4n, K24n, K3n
  ! The following Benthic-states are used (NOT in fluxes): D1m, D2m
  ! The following global scalar vars are used:NO_BOXES_XY,BoxNumberXY,LocalDelta
  ! The following Benthic 1-d global boxvars are used: shiftD1m, KNH4, &
  ! shiftD2m, KNO3
  ! The following 0-d global parameters are used: p_d_tot
  ! The following global constants are used: RLEN,ZERO

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,LOGUNIT,ZERO,DONE,NZERO,dummy
  use mem,  ONLY: K4n, K14n, K24n, K3n,K13n,K23n, K6r,K16r,K26r, D1m, D2m, D2STATE
  use mem, ONLY: ppK4n, ppK14n, ppK24n, ppK3n,ppK13n,ppK23n, &
    ppK6r,ppK16r,ppK26r, &
    BoxNumberXY, NO_BOXES_XY,LocalDelta,max_change_per_step, &
     shiftD1m, KNH4p,KRED, shiftD2m,jK3G4n,M3n,source, &
    KNO3, jK34K24n, jK13K3n, jK23K13n, iiBen, flux,source
  use mem,only:reBTn,ruBTn,ruBPn3,jbotN4n,jK4K3n,jK36K26r
  use constants,ONLY: SHIFT,MASS,INTEGRAL,DERIVATIVE, RFLUX, LAYER1, &
           LAYER2,LAYER3, ANY, NEGATIVE
  use mem_Param,ONLY: p_d_tot_2,combine_anabac
  use mem_BenAmmonium,ONLY: p_flux_K24_at_deep_end=>p_flux_at_deep_end
  use mem_BenAnoxic,ONLY: p_flux_K26_at_deep_end=>p_flux_at_deep_end
  use mem_BenNitrate, ONLY: p_chK3n,sw_method
  use mem_BenDenitriDepth,ONLY:p_clM3n
  use BFM_ERROR_MSG, ONLY: set_warning_for_getm,BFM_ERROR

  use LimitRates, ONLY:LimitShift,LimitChange,insw

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following bennut functions are used:CalculateFromSet
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use bennut_interface,   ONLY: CalculateFromSet


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following sesame functions are used:insw, IntegralExp
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: insw
! use mem,  ONLY: Output2d_1,Output2d_2,Output2d_3, Output2d_4
!  
!
! !AUTHORS
!   Original version by  P. Ruardij
!
!
!
! !REVISION_HISTORY
!   April 15, 1994 by EGM Embsen and P Ruardij:
!               Created a new version of the this process
!               so that it can be used with OpenSESAME.
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
  real(RLEN)  :: Dnew
  real(RLEN)  :: shiftmass
  real(RLEN)  :: jK14K4n,jK24K14n,lK13n,lK23n,jK16K6r,jK26K16r
  real(RLEN)  :: cx_any,rx_any

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  do BoxNumberXY=1,NO_BOXES_XY

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Ammonium Fluxes at the oxic/denitrification boundary
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      Dnew  =   D1m(BoxNumberXY)+ LocalDelta* shiftD1m(BoxNumberXY)

      ! Calculate mass shifted in upwards direction:
      shiftmass = CalculateFromSet( KNH4p(BoxNumberXY), SHIFT, LAYER1, &
        D1m(BoxNumberXY), Dnew)/ LocalDelta

      jK14K4n = CalculateFromSet( KNH4p(BoxNumberXY), DERIVATIVE, RFLUX, &
        D1m(BoxNumberXY), ZERO)+ shiftmass


      call LimitShift(jK14K4n,max(ZERO,K4n(BoxNumberXY)+ &
           source(iiBen,BoxNumberXY,ppK4n)*LocalDelta), &
            K14n(BoxNumberXY),max_change_per_step)
      call flux(BoxNumberXY, iiBen, ppK14n, ppK4n,   jK14K4n* insw(jK14K4n) )
      call flux(BoxNumberXY, iiBen, ppK4n, ppK14n, - jK14K4n* insw(-jK14K4n) )
   
      cx_any= source(iiBen,BoxNumberXY,ppK4n)*LocalDelta
      if ( -cx_any .gt.K4n(BoxNumberXY) ) then
       write(LOGUNIT,*) "K4n,reBTn,-ruBTn,-ruBPn3,-jK4K3n,jbotN4n,jK14K4n:", &
         K4n(BoxNumberXY)/LocalDelta,reBTn(BoxNumberXY), -ruBTn(BoxnumberXY), &
         -ruBPn3(BoxNumberXY),-jK4K3n(BoxNumberXY),jbotN4n(BoxNumberXY),jK14K4n
       call set_warning_for_getm
      endif

      jK16K6r=CalculateFromSet( KRED(BoxNumberXY), &
                                   DERIVATIVE, RFLUX, D1m(BoxNumberXY), dummy)

      ! Calculate mass shifted in upwards direction:
      jK16K6r = jK16K6r+  CalculateFromSet( KRED(BoxNumberXY), SHIFT, LAYER1, & 
                                            D1m(BoxNumberXY), Dnew)/ LocalDelta
      call LimitShift(jK16K6r,max(ZERO,K6r(BoxNumberXY)+ &
           source(iiBen,BoxNumberXY,ppK6r)*LocalDelta), &
            K16r(BoxNumberXY),max_change_per_step)

      call flux(BoxNumberXY, iiBen, ppK16r, ppK6r,   jK16K6r* insw(  jK16K6r) )
      call flux(BoxNumberXY, iiBen, ppK6r,  ppK16r,- jK16K6r* insw( -jK16K6r) )




      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! All the nutrient mineralization source term in the anoxic layer
      ! has been added to K14.n in BenBacDynamics
      ! However in the model this layer is subdivided and hence a partition
      ! flux is here calculated according to the exponential distribution.
      !
      ! Calculate coefficient for the e-folding distribution of the anoxic
      ! mineralization. D7.m is the average penetration depth for N-detritus
      !                 +
      ! Anoxic Mineralization at D1.m, using the exponential distribution
      !                 +
      !          Anoxic Mineralization at D2.m
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Ammonium Fluxes at the denitrification/anoxic boundary
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       
      Dnew  =   D2m(BoxNumberXY)+ LocalDelta* shiftD2m(BoxNumberXY)

      ! Calculate mass shifted in upwards direction:
      jK24K14n = CalculateFromSet( KNH4p(BoxNumberXY), SHIFT, LAYER2, &
                                 D2m(BoxNumberXY), Dnew)/ LocalDelta  &
      + CalculateFromSet( KNH4p(BoxNumberXY), DERIVATIVE, RFLUX, &
                                                D2m(BoxNumberXY), ZERO)

      call LimitShift(jK24K14n,K14n(BoxNumberXY),K24n(BoxNumberXY),&
                                                  max_change_per_step)

      call flux(BoxNumberXY, iiBen, ppK24n, ppK14n, jK24K14n* insw( jK24K14n) )
      call flux(BoxNumberXY, iiBen, ppK14n, ppK24n,-jK24K14n* insw(-jK24K14n) )

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Flux at the lower boundary of ammonium 
      ! 1=no flux, 2= only fluxes_downwards (sink),3=full_flux
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      select case (p_flux_K24_at_deep_end) 
        case(1); rx_any=ZERO;
        case(2); rx_any=min(ZERO,CalculateFromSet(KNH4p(BoxNumberXY), &
                                      DERIVATIVE,RFLUX,p_d_tot_2, ZERO))
        case(3); rx_any=CalculateFromSet(KNH4p(BoxNumberXY), &
                                      DERIVATIVE,RFLUX,p_d_tot_2, ZERO)
      end select
 
      jK34K24n(BoxNumberXY)=rx_any

      call flux(BoxNumberXY, iiBen, ppK24n, ppK24n, jK34K24n(BoxNumberXY) )

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Flux at the lower boundary of nitrate 
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        shiftmass = CalculateFromSet( KNO3(BoxNumberXY), SHIFT, LAYER3, &
            D2m(BoxNumberXY), Dnew)/ LocalDelta
        !if sign of shift mass is different from the shift in the 
        !dentrification depth, the calculate profile in unrealistic 
        ! the may happen when nitrifation input is very low.
        shiftmass=shiftmass * insw(shiftmass *shiftD2m(BoxNumberXY))

        !the nitrate profile will be high near the source and low a near 
        !the denitrifation depth due to denitrfication in the 
        !denitrifcation layer. Therefor it can be assumed that the shift 
        !in biomass will be always smaller that total conentration per m3 
        !this is protection to avoid to use results from unrealistic 
        !calculated profiles 

        lK13n = CalculateFromSet( KNO3(BoxNumberXY), INTEGRAL, MASS, &
            D1m(BoxNumberXY), Dnew)
        lK23n = CalculateFromSet( KNO3(BoxNumberXY), INTEGRAL, MASS, &
             Dnew,p_d_tot_2)
        if (sw_method==2) then
          call LimitChange(NEGATIVE,shiftmass,K13n(BoxNumberXY), &
             max_change_per_step)
        elseif (jK3G4n(BoxNumberXY).gt.NZERO.or.M3n(BoxNumberXY).gt.p_clM3n &
             .and.lK13n>NZERO.and.abs(shiftmass).gt.NZERO ) then
            rx_any=abs(K3n(BoxNumberXY)/LocalDelta)
            shiftmass=max(-rx_any,min(rx_any,shiftmass)) 
        else
            shiftmass=ZERO
        endif

        !In case of low nitrification profile is not mainly deterimined 
        !by nitrifcation input
        !An again a chance on unrealitic profiles: therfore we do not take 
        !in account the flux at the denitrification depth

        ! 1=no flux, 2= only fluxes_downwards (sink), 3,=full_flux
        jK23K13n(BoxNumberXY)=jK23K13n(BoxNumberXY) +shiftmass
        if (K3n(BoxNumberXY).gt.p_chK3n) then
          jK23K13n(BoxNumberXY)=jK23K13n(BoxNumberXY) &
                 -(K3n(BoxNumberXY) -p_chK3n)/LocalDelta
          write(LOGUNIT,*) "K3n is too high: K3n is kept on ",p_chK3n
          call set_warning_for_getm
        endif
        cx_any=K3n(BoxNumberXY)+Source(iiBen,BoxNumberXY,ppK3n)*LocalDelta
        call LimitChange(ANY,jK23K13n(BoxNumberXY),cx_any,max_change_per_step)
        call flux(BoxNumberXY, iiBen, ppK3n, ppK3n,  jK23K13n(BoxNumberXY) ) 
   
        jK13K3n(BoxNumberXY)=CalculateFromSet(KNO3(BoxNumberXY), &
                               DERIVATIVE,RFLUX,D1m(BoxNumberXY), ZERO)


        if (sw_method ==2) then
          call flux(BoxNumberXY,iiBen,ppK13n,ppK13n, &
                                  (lK13n-K13n(BoxNumberXY))/LocalDelta)
          call flux(BoxNumberXY,iiBen,ppK23n,ppK23n, &
                                  (lK23n-K23n(BoxNumberXY))/LocalDelta)
        endif


        jK26K16r = CalculateFromSet( KRED(BoxNumberXY), SHIFT, LAYER3, &
                              D2m(BoxNumberXY), Dnew)/ LocalDelta    &
                +CalculateFromSet( KRED(BoxNumberXY), &
                                    DERIVATIVE, RFLUX, D2m(BoxNumberXY), dummy)
        call LimitShift(jK26K16r,K16r(BoxNumberXY),K26r(BoxNumberXY),&
                                                  max_change_per_step)
      call flux(BoxNumberXY, iiBen, ppK26r, ppK16r, jK26K16r* insw( jK26K16r) )
      call flux(BoxNumberXY, iiBen, ppK16r, ppK26r,-jK26K16r* insw(-jK26K16r) )
     !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! flux at underside of K6r
      ! p_flux_at_deep_end 1=no flux,2= only fluxes_downwards(sink),3=full_flux
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      select case (p_flux_K26_at_deep_end) 
        case (1); rx_any = ZERO
        case (2); rx_any = min(ZERO,CalculateFromSet(KRED(BoxNumberXY), &
                                         DERIVATIVE, RFLUX, p_d_tot_2, dummy))
        case (3); rx_any = CalculateFromSet(KRED(BoxNumberXY), &
                                          DERIVATIVE, RFLUX, p_d_tot_2, dummy)
      end select
      jK36K26r(BoxNumberXY)=rx_any 
      call flux(BoxNumberXY, iiBen, ppK26r, ppK26r,  jK36K26r(BoxNumberXY))



  enddo
  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
