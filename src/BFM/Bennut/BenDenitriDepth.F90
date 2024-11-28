#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenDenitriDepth
!
! DESCRIPTION
!   !

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine BenDenitriDepthDynamics
!
! !USES:

  ! For the following Benthic-states fluxes are defined: D2m
  ! The following Benthic-states are used (NOT in fluxes): K3n,K4n,D1m
  ! The following global scalar vars are used: &
  !    NO_BOXES_XY,   &
  !  BoxNumberXY, InitializeModel
  ! The following Benthic 1-d global boxvars are modified : shiftD1m,shiftD2m
  ! The following Benthic 1-d global boxvars  are used: KNO3E, KNH4
  ! The following 0-d global parameters are used: p_d_tot, p_clD1D2m
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO,LOGUNIT,NZERO,DONE,dummy
  use mem,  ONLY: D2m, D1m,D6m, D2STATE,rrATo
  use mem, ONLY: ppD2m, jK3G4n,NO_BOXES_XY,BoxNumberXY, LocalDelta, &
    InitializeModel, shiftD1m, shiftD2m, KNO3,KNH4,KRED,iiBen, flux,&
    max_change_per_step, jK4K3n,jbotN3n,K4n
  use constants,ONLY:EQUATION,STANDARD, GET,LABDA_1,LABDA_2,ONE_PER_DAY,ANY, &
                                                 LABDA_2,DERIVATIVE,RFLUX
  use mem_Param,  ONLY: p_d_tot, p_clD1D2m,p_qon_dentri,p_qro,p_poro, &
                                                         p_clDxm,p_d_tot_2
  use mem_BenDenitriDepth
  use mem_BenAmmonium,only: p_slK4K3
  use LimitRates, ONLY:LimitChange

!  use mem,  ONLY: Output2d_1,Output2d_2,Output2d_3, Output2d_4


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following bennut functions are used:CalculateFromSet, GetInfoFromSet
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use bennut_interface,   ONLY: CalculateFromSet, GetInfoFromSet
  use mem_globalfun,   ONLY: IntegralExpDist,SetExpDist
!  
!
! !AUTHORS
!   Original version by  P. Ruardij
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
! integer  :: control
  integer     ::i
  real(RLEN)  :: Dxm
  real(RLEN)  :: Dzm
  real(RLEN)  :: D2mNew
  real(RLEN)  :: M3n_D1m,M3n_Dxm,M3n_D2m
  real(RLEN)  :: M6r_D1m,M6r_D2m
  real(RLEN)  :: sK3G4,sK4K3
  real(RLEN)  :: jlK3G4n
  real(RLEN)  :: lambda
  real(RLEN)  :: q_jK3G4_rrATo
  real(RLEN)  :: clD1D2m,r,zuD1,alpha

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! integer, external  :: PrintSet
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  do BoxNumberXY=1,NO_BOXES_XY
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Set minimum distance of D2m  to D1m accroding nexct rule:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       
    clD1D2m=min(0.03D+00,max(D1m(BoxNumberXY),p_clD1D2m))

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Calculate the depth of the denitrification layer on basis of the nitrate 
    ! at oxygen penetration depth.
    ! This is done only when enough input is of nitrate by nitrification
    ! In this case only is the profile determined by the large source of nitrate 
    ! in the oxic layer
    ! If not the dentriffication depth is dtermined by parmeters which set the 
    ! denitrification depth on a fixed distance of the oxygen penetration depth
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    alpha=SetExpDist(D6m(BoxNumberXY),Dlm=p_clDxm,Dmm=p_d_tot)
    zuD1= p_qro*rrATo(BoxNumberXY)/p_poro(BoxNumberXY)/IntegralExpDist(-alpha, &
             p_d_tot_2- D1m(BoxNumberXY))
    r= zuD1*IntegralExpDist(-alpha,D2m(BoxNumberXY)-D1m(BoxNumberXY)) &
                                                       *p_poro(BoxNumberXY) &
     +CalculateFromSet(KRED(BoxNumberXY),DERIVATIVE,RFLUX,D2m(BoxNumberXY),dummy) &
     - CalculateFromSet( KRED(BoxNumberXY), DERIVATIVE, RFLUX, D1m(BoxNumberXY), dummy)
     q_jK3G4_rrATo= jK3G4n(BoxNumberXY)/p_qon_dentri/(NZERO+r/p_qro)

    D2mnew=-999.0;
    !Use first order nitrification coefficeint from ammonium model
    !This is limited for very low and high values.
    sK4K3=GetInfoFromSet( KNH4(BoxNumberXY), GET, LABDA_2, 11)

    r=sK4K3*K4n(BoxNumberXY)-min(ZERO,jbotN3n(BoxnumberXY))
    if (r >p_slK4K3*K4n(BoxNumberXY)) then

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate concentration of nitrate in porewater in M3n:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      M3n_D1m = CalculateFromSet( KNO3(BoxNumberXY), EQUATION, &
        STANDARD, D1m(BoxNumberXY), dummy)
      Dxm=D1m(BoxNumberXY)+D2m(BoxNumberXY) * 0.5
      M3n_Dxm= CalculateFromSet( KNO3(BoxNumberXY), EQUATION, &
        STANDARD, Dxm, dummy)
      M3n_D2m= CalculateFromSet( KNO3(BoxNumberXY), EQUATION, &
        STANDARD, D2m(BoxNumberXY), dummy)
      M6r_D1m = CalculateFromSet( KRED(BoxNumberXY), EQUATION, &
        STANDARD, D1m(BoxNumberXY), dummy)
      M6r_D2m= CalculateFromSet( KRED(BoxNumberXY), EQUATION, &
        STANDARD, D2m(BoxNumberXY), dummy)

      select case ( M3n_D1m<= ZERO)

      case( .TRUE. )

!       ! Do nothing give a warning
!       ! Let if D1m increase let D2m increase too to avoid zero thickness of 
!       ! denitrifiaction layer
!       if ( InitializeModel ==0 ) then
!         write(LOGUNIT,'(''D1m='',F12.3,'' D2m='',F12.3)')  &
!                                      D1m(BoxNumberXY),D2m(BoxNumberXY)
!         write(LOGUNIT,'(''D6m='',F12.3,'' D7m='',F12.3)')  &
!                                      D6m(BoxNumberXY),D7m(BoxNumberXY)
!         write(LOGUNIT,'(''K3n='',F12.3,'' K4n='',F12.3)')  &
!                                      K3n(BoxNumberXY),K4n(BoxNumberXY)
!         write(LOGUNIT,'(''K14n='',F12.3,'' K24n='',F12.3)')  &
!                                      K14n(BoxNumberXY),K24n(BoxNumberXY)
!         write(LOGUNIT,'(''reBTn='',F12.3,'' reATn='',F12.3)')  &
!                                      reBTn(BoxNumberXY),reATn(BoxNumberXY)
!         write(LOGUNIT,'(''fluxN3='',F12.3,'' fluxK4n='',F12.3)')  &
!                                    jbotN3n(BoxNumberXY),jbotN4n(BoxNumberXY)
!         write(LOGUNIT,'(''M3n='',F12.3,'' M4n='',F12.3)')  &
!                                      M3n(BoxNumberXY),M4n(BoxNumberXY)
!         write(LOGUNIT,'(''N3n='',F12.3,'' N4n='',F12.3)') &
!                                    N3n_Ben(BoxNumberXY),N4n_Ben(BoxNumberXY)
!         write(LOGUNIT,'(''K6r='',F12.3,'' sK3G4='',F12.3)') &
!        K6r(BoxNumberXY),GetInfoFromSet( KNO3(BoxNumberXY), GET, LABDA_2, 21)  
!        write(LOGUNIT,'(''jK3G4n='',F12.3,'' jK4K3n='',F12.3)') &
!                                  jK3G4n(BoxNumberXY),jK4K3n(BoxNumberXY)
!         call PrintSet(  KNH4(BoxNumberXY) "concentration nitrate on D1m < 0")
!         call PrintSet(  KNO3(BoxNumberXY) "concentration nitrate on D1m < 0")
!      endif

        D2mNew =D2m(BoxNumberXY)
        shiftD2m(BoxNumberXY) = max( min( p_d_tot- &
          p_clD1D2m, D2mNew), D1m(BoxNumberXY)+ clD1D2m)- D2m(BoxNumberXY)

      case( .FALSE. )

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! According solution nitrate concentration decreases
        ! exponentially. Calculate depth at where below the 
        ! denitrification /m2 is  below a (small) fixed number. 
        ! Use this new depth as the (uncorrected) new denitrification depth
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        sK3G4= GetInfoFromSet(KNO3(BoxNumberXY), GET, LABDA_2, 31,dummy,dummy) 
        ! nagative vale!
        lambda= GetInfoFromSet(KNO3(BoxNumberXY), GET, LABDA_1, 31,dummy,dummy) 
        ! Calculate the denitrification rate which will below the new D2m 
        ! (10% of total)

        jlK3G4n=(jK3G4n(BoxNumberXY)+sK3G4*M3n_D2m/(-lambda)) * p_pK3G4n

        ! Calculate the depth under which the nitrate concentration is 
        ! below 0.1 mmol/m3
        !Solve the equation:
        !(p_clM3n)=M3n_Dxm*e(-d*lamba)
        r=+log(p_clM3n/M3n_Dxm)/lambda

        !Calculate the depth under which 10% of the denitrification take place
        D2mNew= log(jlK3G4n*(-lambda)/(sK3G4* M3n_Dxm))/(lambda) &
                                     *min(DONE,q_jK3G4_rrATo) 
        ! The new sulphide horizon is the minimum of the 2 conditions
        D2mnew=max(ZERO,min(D2mnew,r))+Dxm
        Dzm=p_d_tot-p_clD1D2m

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! 1. Calculate uncorrected shift of D2.m
        ! 2. limit shift incase of D2mNew moves in the direction of D1m
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        shiftD2m(BoxNumberXY) = (D2mNew - D2m(BoxNumberXY))*ONE_PER_DAY
        shiftD2m(BoxNumberXY) = shiftD2m(BoxNumberXY) *  &
                              p_mD2m /(p_mD2m+abs(shiftD2m(BoxNumberXY) ))

        shiftD2m(BoxNumberXY)=shiftD2m(BoxnumberXY)*D2m(BoxNumberXY)/ &
                (D2m(BoxNumberXY)+abs(shiftD2m(BoxNumberXY)))
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Correct by damping the change of D2m in case large changes:
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        if ( InitializeModel== 0) then

          r=max(ZERO,Dzm-D2mnew)
          shiftD2m(BoxNumberXY) = shiftD2m(BoxNumberXY)* r/(Dzm+r)

          call LimitChange(ANY,shiftD2m(BoxNumberXY),D2m(BoxNumberXY), &
                                                   max_change_per_step)
         
          shiftD2m(BoxNumberXY)=min(Dzm, &
            max(shiftD2m(BoxNumberXY)+D2m(BoxNumberXY),&
            clD1D2m+D1m(BoxNumberXY)+shiftD1m(BoxNumberXY))) -D2m(BoxNumberXY)

           call flux(BoxNumberXY, iiBen, ppD2m, ppD2m, shiftD2m(BoxNumberXY) )

        end if
    end select

   else
     Dzm=p_d_tot-p_clD1D2m
     shiftD2m(BoxNumberXY)=min(Dzm,max((D1m(BoxNumberXY)-D2m(BoxNumberXY))* &
          (DONE-q_jK3G4_rrATo)+D1m(BoxNumberXY),&
          clD1D2m+D1m(BoxNumberXY)+shiftD1m(BoxNumberXY))) -D2m(BoxNumberXY)
     shiftD2m(BoxNumberXY) = shiftD2m(BoxNumberXY) *  &
                        p_mD2m /(p_mD2m+abs(shiftD2m(BoxNumberXY) ))
     call flux(BoxNumberXY, iiBen, ppD2m, ppD2m, shiftD2m(BoxNumberXY) )
   endif
  end do

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
