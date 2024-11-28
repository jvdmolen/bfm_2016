#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenAmmonium
!
! DESCRIPTION
!   Description of the diagenetic ammonium processes in the sediment
!   Details on the equations and the method used to calculate
!   the equilibrium and transient profiles can be found in
!   Ruardij et al., 1995. Neth. J. Sea Res. 33(3/4):453-483
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine BenAmmoniumDynamics
!
! !USES:

  ! For the following Benthic-states fluxes are defined: K4n, K3n, G2o
  ! The following Benthic-states are used (NOT in fluxes): D1m, K24n, &
  ! D7m
  ! The following Benthic 1-d global boxvars are modified: M4n, KNH4, jG2K3o
  ! The following Benthic 1-d global boxvars got a value: M14n, M24n
  ! The following Benthic 1-d global boxvars are used: reBTn,  &
  ! irrenh, ETW_Ben, N4n_Ben, Depth_Ben
  ! The following Benthic 1-d global boxpars  are used: p_poro
  ! The following 0-d global box parametes are used: p_d_tot, p_clDxm, &
  ! p_q10diff, p_qon_nitri

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  USE BFM_ERROR_MSG, ONLY: set_warning_for_getm
  use global_mem,ONLY:RLEN,ZERO,NZERO,DONE,dummy
  use mem,       ONLY: K4n, G2o, D1m, K14n, D2m, K24n, D7m,DH2m,DH3m,D2STATE
  use mem,       ONLY: ppN4n, ppK4n, ppK3n, ppG2o, NO_BOXES_XY, &
    BoxNumberXY, InitializeModel, LocalDelta,max_change_per_step,dry_z, &
    M4n,KNH4,KNH4p,jKuK4n,jG2K3o,M14n,M24n, reBTn, ruBPn3, ruBTn, &
    irrenh, ETW_Ben,jK4K3n, N4n_Ben, Depth_Ben, iiBen,iiPel,iiHN, flux,sK4K3
  use constants, ONLY: LAYER1, ADD, ANY, POSITIVE, &
                 INTEGRAL, RFLUX, DERIVATIVE,SET_LAYER_INTEGRAL
  use mem_Param,ONLY: p_poro, p_d_tot_2,p_d_tot,p_clDxm,p_q10diff,p_qon_nitri, &
                 combine_anabac, CalcBenBacteria,p_dry_ben
  use mem_Diffusion,ONLY: p_diff=>p_diff_N4
  use LimitRates,ONLY:LimitShift,LimitChange
  use botflux,   ONLY:addbotflux
  use mem_BenAmmonium
#ifdef INCLUDE_BENCO2
  use mem_CO2,ONLY:p_qhK4K3n
#endif
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following bennut functions are used:InitializeSet, &
  ! DefineSet, CompleteSet, CalculateSet, CalculateTau, CalculateFromSet
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use bennut_interface, ONLY: CalculateSet, CalculateTau, CalculateFromSet, &
                              CopySet

  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! The following global functions are used:eTq
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use global_interface,   ONLY: eTq

  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following sesame functions are used:IntegralExp
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: SetExpDist,IntegralExp
!  
!
! !AUTHORS
!   Original version by  P. Ruardij
!
!
! !REVISION_HISTORY
!   September 1999 by M. Vichi !               Commented version
!
!
! COPYING
!   
!   Copyright (C) 2006 P. Ruardij, the BFM team 
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
  real(RLEN)  :: rdummy,r
  real(RLEN)  :: zuBT
  real(RLEN)  :: alpha,beta
  real(RLEN)  :: diff
  real(RLEN)  :: Tau
  real(RLEN)  :: labda
  real(RLEN)  :: a15
  real(RLEN)  :: cK4n
  real(RLEN)  :: jK4N4n,jK14K4n
  real(RLEN)  :: cO2
  real(RLEN)  :: scK4K3
  real(RLEN)  :: senK4n !relative loss of N4n per unit N4n due to update by Phytoplankton
  real(RLEN)  :: renK4n !loss of N4n due to update by Phytoplankton
  logical     :: dry

  external    :: EquationBenAmmonium

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  do BoxNumberXY=1,NO_BOXES_XY
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate the pore-water average concentrations from the state variables
      ! (Diagnostic variables, not used in calculations)
      ! Calculate pore-water oxygen concentration in the oxic layer
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      M4n(BoxNumberXY) = K4n(BoxNumberXY)/(p_poro(BoxNumberXY) &
                                               *(p_p+DONE)* D1m(BoxNumberXY))
      M14n(BoxNumberXY) = K14n(BoxNumberXY)/ p_poro(BoxNumberXY)/ &
                        ( p_p+ DONE)/( D2m(BoxNumberXY)- D1m(BoxNumberXY))
      M24n(BoxNumberXY) = K24n(BoxNumberXY)/ p_poro(BoxNumberXY)/ &
                        ( p_p+ DONE)/( p_d_tot_2- D2m(BoxNumberXY))

      if (K4n(BoxNumberXY)<0) then
        K4n(BoxNumberXY)=max(NZERO,N4n_Ben(BoxNumberXY)+M14n(BoxNumberXY)) &
          *0.5D+00*D1m(BoxnumberXY) *(p_p+DONE)*p_poro(BoxNumberXY)
        write(LOGUNIT,*) 'K4n is reset to',K4n(BoxNumberXY)
        call set_warning_for_getm()
      endif
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate coefficient for the e-folding distribution of the anoxic
      ! mineralization. D7.m is the average penetration depth for N-detritus
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      if ( combine_anabac) then
        alpha  =   SetExpDist(D7m(BoxNumberXY),Dlm=p_clDxm,Dmm=p_d_tot)
        beta=alpha
      else
        alpha  =   SetExpDist(DH2m(BoxNumberXY),Dlm=p_clDxm,Dmm=p_d_tot)
        beta  =    SetExpDist(DH3m(BoxNumberXY),Dlm=p_clDxm,Dmm=p_d_tot)
      endif

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Recalculate Mineralization m2 --> m3 porewater
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        ! uptake of N by heterotrophic bacteria (shortage), nitrifying bacteria 
        ! & benthic diatoms (prim.prod),corrected for nitrate uptake by diatoms
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Netto Average in the oxic layer:
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

        zuBT  =   (reBTn(BoxNumberXY)+jKuK4n(BoxNumberXY)) &
                                      / p_poro(BoxNumberXY)/ D1m(BoxNumberXY)
        zuBT=max(zuBT,1.079D-6)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Correction due to environmental regulating factors,
      ! diffusion coefficient: temperature and bioirrigation
      ! nitrification rate: temperature and oxygen
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      diff = p_diff* irrenh(BoxNumberXY)* &
                                eTq( ETW_Ben(BoxNumberXY), p_q10diff)

      if ( .not.CalcBenBacteria(iiHN) ) then
        sK4K3(BoxNumberXY) = p_sK4K3* eTq( ETW_Ben(BoxNumberXY), p_q10) 

        cO2 =max(1.D-20,G2o(BoxNumberXY)/ p_poro(BoxNumberXY)/ D1m(BoxNumberXY))
        if (InitializeModel == 0 ) &
          sK4K3(BoxNumberXY)=sK4K3(BoxNumberXY) * &
            cO2/( cO2+ p_clO2o)* M4n(BoxNumberXY)/( M4n(BoxNumberXY)+ p_clM4n)
      endif
      
      if ( isnan(sK4K3(BoxNumberXY))) then
        write(LOGUNIT,*) 'nitrification rate=NAN'
      endif

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Calculate the uptake by diatoms of ammonium ruBTn-ruBPN3
      ! ruBTn =ammonium uptake by Benthic diatoms +uptake by bacteria
      !        (only) to  fill the internal buffer in the cell.
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      renK4n= ruBTn(BoxNumberXY)- ruBPn3(BoxNumberXY)
      senK4n=ZERO;if (renK4n.gt.ZERO) senK4n=renK4n/ p_poro(BoxNumberXY)/ &
                    D1m(BoxNumberXY)/max(NZERO,M4n(BoxNumberXY)) 

      dry=.false.
      if (p_dry_ben) dry=( dry_z(BoxNumberXY) <=p_clDxm )

      select case ( InitializeModel)
        case ( 0 )
          !the flux is determined on basis of the total change  of NH4
          ! = nitrification + uptake by diatom + uptake by bacteria 
          scK4K3  =   min(p_shK4K3,max(p_slK4K3,sK4K3(BoxNumberXY)+senK4n))
          labda  =   sqrt( scK4K3 / diff)
          ! Calculate coefficient of the zero order term
          a15  =   zuBT/scK4K3 
          call EquationBenAmmonium(.FALSE.,KNH4p(BoxNumberXY),InitializeModel, &
          D1m(BoxNumberXY),D2m(BoxNumberXY),diff,p_poro(BoxNumberXY),p_p,&
          labda,scK4K3,alpha,beta,a15,reBTn(BoxNumberXY),p_d_tot_2 ,&
          N4n_Ben(BoxNumberXY),K14n(BoxNumberXY),K24n(BoxNumberXY) )
          cK4n = CalculateSet( KNH4p(BoxNumberXY), SET_LAYER_INTEGRAL, LAYER1, &
            LAYER1, dummy, ZERO)
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          ! Calculate for the above defined set of boundary conditions
          ! the steady-state profiles and return the vertically integrated
          ! concentration in the first layer.
          !
          ! Calculate the adaptation time to the steady-state profile
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

          Tau  =   CalculateTau(  scK4K3,  diff,  p_p,  D1m(BoxNumberXY))

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          ! Estimate the average value of K4n over the actual time step
          ! (transient value).
          ! This value depends on the adaptation time, the actual time step,
          ! the ''old'' value and the ''equilibrium value''
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
          cK4n = cK4n+( K4n(BoxNumberXY)- cK4n)* IntegralExp( - LocalDelta/ &
            Tau, DONE)
  
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          ! Derive the equations for the transient profiles, assuming the same
          ! solution as for the steady-state case and using cK4n as new &
          ! constraint.
          !-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
          rdummy = CalculateSet( KNH4p(BoxNumberXY), ADD, 0, 0, dummy, cK4n)
  
          !-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          ! Start calculation of fluxes:
          !
          ! Next calculation is done in BenNitrogenShifting which
          ! include all fluxes between layers!
          ! All the nutrient mineralization source term in the anoxic layer
          ! has been added to K14.n in BenBacDynamics
          ! However in the model this layer is subdivided and hence a partition
          ! flux is here calculated according to the exponential distribution.
          !
          ! Calculate Nitrification flux in the first layer and the related
          ! oxygen consumption flux:
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

          if ( p_sK4K3 >ZERO .or.  (.not.CalcBenBacteria(iiHN)) ) then 
            jK4K3n(BoxNumberXY)=max(ZERO, sK4K3(BoxNumberXY)* K4n(BoxNumberXY))
            call LimitChange(POSITIVE,jK4K3n(BoxNumberXY),K4n(BoxNumberXY),&
                                                           max_change_per_step)
          endif
          call flux(BoxNumberXY, iiBen, ppK4n, ppK3n, jK4K3n(BoxNumberXY) )

          jG2K3o(BoxNumberXY)  =   jK4K3n(BoxNumberXY)* p_qon_nitri
          call flux(BoxNumberXY, iiBen, ppG2o, ppG2o, -( jG2K3o(BoxNumberXY)) )

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          ! 1.Estimation of the Vertical fluxes from set of transient solutions:
          !
          ! 2.Avoid to large fluxes leading to negative concentrations
          !   Fluxes between the layers are calculated in NitrogenBenShifting
          !   However we caculate here the flux from layer2->layer1 to have
          !   a better correction for the flux jK4N4n
          !-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

          jK14K4n = CalculateFromSet( KNH4p(BoxNumberXY), DERIVATIVE, RFLUX, &
                                                        D1m(BoxNumberXY), ZERO)
          jK4N4n = CalculateFromSet( KNH4p(BoxNumberXY), &
                                         DERIVATIVE, RFLUX, ZERO, ZERO)
          r=K4n(BoxNumberXY)+ LocalDelta*  &
           (reBTn(BoxNumberXY)-renK4n-jK4K3n(BoxNumberXY)+jK14K4n)
          call LimitShift(jK4N4n,N4n_Ben(BoxNumberXY)*Depth_Ben(BoxNumberXY), &
            r, max_change_per_step)

          call addbotflux(ANY,BoxNumberXY,iiBen,ppK4n,iiPel,ppN4n,jK4N4n)

          !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          !  All transport between layers are done in BenNitrogenShifting:
          !
          !  jK14K4n = CalculateFromSet(KNH4, DERIVATIVE, RFLUX, D1.m, 0.0);
          !  jK24K14n= CalculateFromSet(KNH4, DERIVATIVE, RFLUX, D2.m, 0.0);
          !  jK34K24n= CalculateFromSet(KNH4, DERIVATIVE, RFLUX, p_d_tot, 0.0);
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
          if ( senK4n < NZERO ) then
             !neglectable loss other  the nitrification
             KNH4(BoxNumberXY) = CopySet(KNH4p(BoxNumberXY),KNH4(BoxNumberXY))
          else 
            !the nitrate production as calculated by the ammonium model is used
            !as input for the nitrate model and all other loss terms other than
            ! nitrification has to be excluded
            scK4K3  =   min(p_shK4K3,max(p_slK4K3,sK4K3(BoxNumberXY)))
            labda  =   sqrt( abs(scK4K3) / diff)
            a15  =   zuBT/scK4K3 
            call EquationBenAmmonium(.FALSE.,KNH4(BoxNumberXY),InitializeModel,&
             D1m(BoxNumberXY),D2m(BoxNumberXY),diff,p_poro(BoxNumberXY),p_p,&
             labda,scK4K3,alpha,beta,a15,reBTn(BoxNumberXY),p_d_tot_2 ,&
             N4n_Ben(BoxNumberXY),K14n(BoxNumberXY),K24n(BoxNumberXY) )
             cK4n=CalculateSet(KNH4(BoxNumberXY), SET_LAYER_INTEGRAL, LAYER1, &
                                           LAYER1, dummy, ZERO)
             cK4n = cK4n+( K4n(BoxNumberXY)- cK4n)* IntegralExp(-LocalDelta/ &
                         Tau, DONE)
             rdummy = CalculateSet( KNH4(BoxNumberXY), ADD, 0, 0, dummy, cK4n)
          endif
        case ( 1 )    ! Determination of Initial condition:
          scK4K3  =   min(p_shK4K3,max(p_slK4K3,sK4K3(BoxNumberXY)))
          labda  =   sqrt( scK4K3 / diff)
          a15  =   zuBT/scK4K3 
          call EquationBenAmmonium(.false.,KNH4(BoxNumberXY),InitializeModel, &
             D1m(BoxNumberXY),D2m(BoxNumberXY),diff,p_poro(BoxNumberXY),p_p,&
             labda,scK4K3,alpha,beta,a15,reBTn(BoxNumberXY),p_d_tot_2 ,&
              N4n_Ben(BoxNumberXY),ZERO,ZERO )

          rdummy = CalculateSet( KNH4(BoxNumberXY), 0, 0,  0, dummy, dummy)

          if (dry .and. (.not.CalcBenBacteria(iiHN)) ) then
             jK4K3n(BoxNumberXY) =  ZERO
          else
             jK4K3n(BoxNumberXY) = sK4K3(BoxNumberXY)*  &
             CalculateFromSet( KNH4(BoxNumberXY), INTEGRAL, RFLUX, ZERO, &
                                                       D1m(BoxNumberXY))
          endif
          jG2K3o(BoxNumberXY)  =   jK4K3n(BoxNumberXY)* p_qon_nitri

          if (jK4K3n(BoxNumberXY) < ZERO) then
            write(LOGUNIT,'(''Negative value calculated for K4n'')') 
            write(LOGUNIT,'(''Minrealisation (reBTT) ='',F10.3)') &
                                                      reBTn(BoxNumberXY)
            write(LOGUNIT,'(''Temperature (ETW) ='',F10.3)')ETW_Ben(BoxNumberXY)
            write(LOGUNIT,'(''Amommonium at interface (N4n_Ben) ='',F10.3)') &
                                                      N4N_Ben(BoxNumberXY)
            call PrintSet(KNH4(BoxNumberXY), &
                              "Error in Initialization of ammonium gradient")
          endif

      end select

  end do

  end
  subroutine EquationBenAmmonium(dry,KNH4,mode,D1m,D2m,diff,p_poro,p_p, &
                   labda,scK4K3,alpha,beta,a15,reBTn,p_d_tot,N4n,K14n,K24n )
  use global_mem, ONLY:RLEN,ZERO,dummy,DONE
  use constants, ONLY: LAYERS, LAYER1, LAYER2,LAYER3,DEFINE,   &
  DIFFUSION, FOR_ALL_LAYERS, POROSITY, ADSORPTION, DOUBLE_DEFINE, &
  EXPONENTIAL_TERM, CONSTANT_TERM, DIST_EXPONENTIAL_TERM,  &
  LINEAR_TERM, SET_CONTINUITY, FLAG, MASS, SET_BOUNDARY, EQUATION,&
  SET_LAYER_INTEGRAL, SET_LAYER_INTEGRAL_UNTIL, INPUT_TERM,  &
  PARAMETER, STANDARD,DERIVATIVE

  use mem_globalfun,   ONLY: IntegralExpDist,ExpDist
  use bennut_interface, ONLY: InitializeSet, DefineSet, CompleteSet

  implicit none
  logical,intent(IN)      :: dry
  integer,intent(INOUT)   :: KNH4
  integer,intent(IN)      :: mode
  REAL(RLEN),intent(IN)   :: D1m
  REAL(RLEN),intent(IN)   :: D2m
  REAL(RLEN),intent(IN)   :: diff
  REAL(RLEN),intent(IN)   :: p_poro
  REAL(RLEN),intent(IN)   :: p_p
                  
  REAL(RLEN),intent(IN)   :: labda
  REAL(RLEN),intent(IN)   :: scK4K3
  REAL(RLEN),intent(IN)   :: alpha,beta
  REAL(RLEN),intent(IN)   :: a15
  REAL(RLEN),intent(IN)   :: reBTn
  REAL(RLEN),intent(IN)   :: p_d_tot
  REAL(RLEN),intent(IN)   :: N4n
  REAL(RLEN),intent(IN)   :: K14n
  REAL(RLEN),intent(IN)   :: K24n 

  REAL(RLEN)              ::r
  external                ::FixProportionCoeff
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Initialize the set of differential equations giving:
  ! - n. of layers;
  ! - n. of coefficients
  ! - layer depths
  ! - environmental conditions (diffusion, p_porosity and adsorption coeff.)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      KNH4 = InitializeSet( KNH4, 3, 8)

      call DefineSet( KNH4, LAYERS, LAYER1, LAYER2, D1m, D2m)
      call DefineSet( KNH4, DIFFUSION, FOR_ALL_LAYERS, 0, diff, dummy)
      call DefineSet( KNH4, POROSITY, FOR_ALL_LAYERS, 0, p_poro, dummy)
      call DefineSet( KNH4, ADSORPTION, FOR_ALL_LAYERS,0, p_p, dummy)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Define coefficients for the steady-state solutions in each layer
      ! General solution of the equilibrium profile:
      ! 1st layer:
      ! A(z) = a11*exp(labda*z) + a12*exp(-labda*z) + a13*z^2 + a14*z + a15
      ! 2nd layer:
      ! A(z) = a21*exp[-alpha*(z-D1.m)] + a24*z + a25
      ! 3rd layer:
      ! A(z) = a31*exp[-alpha*(z-D2.m)] + a34*z + a35
      !    a34 = 0 (boundary condition)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      call DefineSet(KNH4,DOUBLE_DEFINE,11,EXPONENTIAL_TERM,labda,scK4K3)
      call DefineSet(KNH4,DOUBLE_DEFINE, 12, EXPONENTIAL_TERM, &
                                             -labda, scK4K3)
      call DefineSet(KNH4,DEFINE, 15,CONSTANT_TERM, dummy, dummy)

      call DefineSet(KNH4,DEFINE, 21,DIST_EXPONENTIAL_TERM, &
                                                   -alpha, dummy)
      call DefineSet(KNH4,DEFINE, 24,LINEAR_TERM, dummy,dummy)
      call DefineSet(KNH4,DEFINE, 25,CONSTANT_TERM,dummy,dummy)

      call DefineSet( KNH4, DEFINE,31, DIST_EXPONENTIAL_TERM, &
                                                 -beta, dummy)
      call DefineSet( KNH4, DEFINE,35, CONSTANT_TERM, dummy, dummy)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Insert other boundary conditions and continuity between layers:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      !4/6
      call CompleteSet( KNH4, SET_CONTINUITY, FLAG, MASS, dummy)
      !5/7
      if (dry) then
        call CompleteSet(KNH4,SET_BOUNDARY,LAYER1,DERIVATIVE,ZERO,value=ZERO)
      else
        call CompleteSet(KNH4,SET_BOUNDARY,LAYER1,EQUATION,ZERO,value=N4n)
      endif

      select case (mode)
         case(0)
           !6/8
           call CompleteSet( KNH4, SET_LAYER_INTEGRAL, LAYER2, &
            LAYER2, dummy, value=K14n)
           !7/9
           call CompleteSet( KNH4, SET_LAYER_INTEGRAL_UNTIL, &
            LAYER3, LAYER3, p_d_tot, value=K24n)
           !8
           call CompleteSet( KNH4, INPUT_TERM, 15, STANDARD, &
            dummy,value=a15)
         case(1)           !6 
          r  =   ExpDist( - alpha, D2m- D1m)
          call FixProportionCoeff(KNH4,21,31,DONE,r)
          !7
          call CompleteSet(KNH4,INPUT_TERM,15, STANDARD, dummy, &
            value=a15)
          !8
          r= reBTn / p_poro/   &
            IntegralExpDist( - alpha, D1m) *ExpDist(-alpha , D1m)
          call CompleteSet( KNH4, INPUT_TERM, 21, PARAMETER, &
            dummy, value=r)
      end select
      return
      end


!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
