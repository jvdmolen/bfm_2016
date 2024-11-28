#include "INCLUDE.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Param
!
! DESCRIPTION
!   List of global model parameters 
!      (global variables which can be changed during the model initialization

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  MODULE mem_Param

!
! !USES:

  USE global_mem
  USE constants
  USE mem, ONLY: iiPhytoPlankton, iiMesoZooPlankton, &
    iiMicroZooPlankton, iiBenOrganisms, iiBenBacteria, &
    iiBenthicPhosphate, iiBenthicAmmonium,iiBenPhyto

!  
!
! !AUTHORS
!   mfstep/ERSEM team
!
! !REVISION_HISTORY
!   --------
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
  ! Default all is public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  PUBLIC

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Global Model Parameters
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ! 0d-parameter used in pelagic submodel
  real(RLEN)   :: &
      p_small=1.0D-80  ,  &
      p_q10diff=1.49  ,  &  ! Temperature-dependency porewater diffusion
      p_qro=0.5  ,  &  ! stoichiometry O2-->S2-
      p_qon_dentri=1.25  ,  &  ! stoichiometry O2-->N denitrification
      p_qon_nitri=2.0  ,  &  ! stoichiometry O2-->N nitrification
      p_clDxm=0.001, &  ! minimal value of D?.m for calculation of the alpha
      p_1d_poro=0.40 !porosity value used when running in 1D
  !  alpha is used in expo.func and values > 1/0.001 leadt

  ! 0d-parameter used in pelagic submodel
  logical   :: &
      CalcPelagicFlag=.TRUE.  ! Switch for pelagic system
  ! 0d-parameter used in pelagic submodel
  integer   :: &
      CalcBenthicFlag=3  ! Switch for benthic system

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !  Allocate the logical flags for switch on the LFG
  !  Initialise to TRUE (overwritten by the namelist values)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  logical   :: CalcPhytoPlankton(iiPhytoPlankton) = .TRUE.
  logical   :: CalcMicroZooPlankton(iiMicroZooPlankton) = .TRUE.
  logical   :: CalcMesoZooPlankton(iiMesoZooPlankton) = .TRUE.
  logical   :: CalcBenOrganisms(iiBenOrganisms) = .TRUE.
  logical   :: CalcSuspensionFeeders(iiBenOrganisms) = .TRUE.
  logical   :: CalcBenBacteria(iiBenBacteria) = .TRUE.
  logical   :: CalcBacteria = .TRUE.
  logical   :: CalcBenPhyto(iiBenPhyto) = .FALSE.
  logical   :: CalcMacroPhyto = .FALSE.
  logical   :: combine_anabac = .TRUE.

  logical   :: &
      CalcPelChemistry=.TRUE.  ,  &  !
      AssignPelBenFluxesInBFMFlag=.TRUE.  ,  &  ! Switches to make choice to define boundary
      AssignAirPelFluxesInBFMFlag=.TRUE.        ! fluxes in physical of biological model
  !      %dim2D%

  ! 1d-parameter used in benthic submodel
  real(RLEN),public,dimension(:),allocatable   :: &
      p_pK1_ae  ,  &
      p_pK5_ae  ,  &
      p_poro
  ! 0d-parameter used in pelagic submodel
  integer   :: &
      check_fixed_quota=0  , &
      use_function_for_eps0=1
  ! 0d-parameter used in pelagic submodel
  real(RLEN)   :: &
      XLatitude=54.0  ,  &  ! Latitude
      p_PAR=0.50  ! Photosynthetically available radiation
  ! 0d-parameter used in pelagic submodel
  integer   :: &
      ChlLightFlag=2  ,  &  ! Switch between light prop.(=1) or Chla.(=2) as a state
      LightForcingFlag=1  ! Switch between instantaneous light and day light average
  ! 1d-parameter used in pelagic submodel
  
  ! 0d-parameter used in pelagic submodel
  real(RLEN)   :: &
      p_eps0=0.04  ,  &  ! Background extinction (abiotic)
      p_epsESS=0.04e-3  ,  &  ! Inorg. suspended matter extinction coeff. (abiotic)
      p_InitSink=100.0  ,  &  ! parameter to Initialize BenthicSInk var.
      p_d_tot=0.30  ,  &  ! m # maximal Thickness of  of D2m
      p_d_tot_2=0.35  ,  &  ! m # Thickness of modelled benthic sediment layers
      p_clD1D2m=0.01  ,  &  ! m # minimum distancebetween D1m and D2m
      p_pe_R1c=0.60  ,  &  ! Fraction of excretion going to PLOC of phyto
      p_pe_R1n=0.72  ,  &  ! Fraction of excretion going to PLOC of phyto
      p_pe_R1p=0.832  ,  &  ! Fraction of excretion going to PLOC of phyto
      p_peZ_R1c=0.40  ,  &  ! Fraction of excretion going to PLOC of phyto
      p_peZ_R1n=0.72  ,  &  ! Fraction of excretion going to PLOC of phyto
      p_peZ_R1p=0.832  ,  &  ! Fraction of excretion going to PLOC of phyto
      p_pe_R1s=0.06  ,  &
      p_epsChla=10.0e-3  ,  &  ! Chla-contribution to extinction
      p_epsR6=0.1e-3, &  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      p_xeff_an=0.5, &    ! relative efficieny of anaerobic processes compared to aerobic
      p_mdTdz=0.01 ! if the max(DTdz) between all layers>0.01 the mixlayer depths are
                   ! calculated in CalulateThermoVars; otherwise all three vars which
                   ! describe the these depths (BMLD,TMld,SdTdzTha) are set on ZERO. 

   
     real(RLEN) :: &
        p_max_state_change=0.1, &        !max allowed rate of change (used in 
                              ! predictor/corrector routine LimitNutrientUptake
        p_max_state_rate=0.1, &
        p_qpPhc=0.002, &      ! quota where above phophatase take place
        p_qnUlc=0.025, &      ! quota where above urea uptake take place 
        p_qR6cQ9x=0.4, &      ! C-detritus coupled to inorganic Silt.
        p_xdensilt_k=2.798e+6 ! density of silt ( kg/m3)  
   real(RLEN) :: p_surface_1d=1.0e6 ! for modelling macrophyt a size of the 
                                    ! area is needed. It is only when the 
                                    ! macrophytmodel is appled in a 1d-mode 
                                    !(=column model). In case of application in
                                    ! 3d-model the surface is derived from the 
                                    !bathymetry
   integer  ::p_check_track=0 ! parameter works only whene tracking is active.
            ! p_check_track>0 check if content of tracked state variable is <= state variable
            ! tracking take place at start of Ecology routine (p_check_track=1)
            ! and/or tracking take place at start of Track routine (p_check_track=+10)
            ! and/or tracking take place at end of Track routine (p_check_track=+100)

  logical  ::p_dry_ben=.FALSE.
  logical  ::nan_check=.FALSE.
  logical  ::mass_conservation_check=.FALSE.
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")

  public InitParam
  
  contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitParam()
  use mem
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  namelist /Param_parameters/ p_small, p_q10diff, p_qro, p_qon_dentri,  &
    p_qon_nitri,p_clDxm,p_1d_poro, CalcPelagicFlag, CalcBenthicFlag, &
    CalcPhytoPlankton,CalcMicroZooPlankton, CalcPelChemistry, &
    CalcMesoZooPlankton,CalcBenOrganisms,CalcSuspensionFeeders, &
    CalcBenBacteria, CalcBacteria, CalcBenPhyto, CalcMacroPhyto, &
    combine_anabac, AssignPelBenFluxesInBFMFlag, AssignAirPelFluxesInBFMFlag, &
    p_PAR, ChlLightFlag, LightForcingFlag, use_function_for_eps0,p_eps0, & 
    p_epsESS, p_InitSink, p_d_tot,p_d_tot_2, p_clD1D2m, p_pe_R1c, p_pe_R1n, &
    p_pe_R1p, p_pe_R1s, p_peZ_R1c, p_peZ_R1n,p_peZ_R1p,   &
    p_epsChla, p_epsR6,p_xeff_an,p_max_state_change,p_max_state_rate, &
    p_qpPhc, p_qnUlc,p_qR6cQ9x, p_xdensilt_k, p_mdTdz,p_surface_1d, &
    check_fixed_quota,p_dry_ben,p_check_track,nan_check,mass_conservation_check
    integer :: i
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !BEGIN compute

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
   write(LOGUNIT,*) "#  Reading Param parameters.."
   open(NMLUNIT,file='Param.nml',status='old',action='read',err=100)
   read(NMLUNIT,nml=Param_parameters,err=101)
   close(NMLUNIT)
   write(LOGUNIT,*) "#  Namelist is:"
   write(LOGUNIT,nml=Param_parameters)
   !larvae of bentic of suspension feeders exist only if young filterdeeders are modelled!
   CalcMesoZooPlankton(3)=CalcMesoZooPlankton(3) &
                                       .and. CalcSuspensionFeeders(iiYy3)
   CalcBenBacteria(3)=CalcBenBacteria(3).and.(.not.combine_anabac) 
#ifndef INCLUDE_MACROPHYT
!JM   if ( CalcMacroPhyto==1) then 
   if ( CalcMacroPhyto .eqv. .true.) then 
        write(LOGUNIT,*)'Parameter MacroPhyto is reset to ''zero'''
        write(LOGUNIT,*)'Sub model for Macrophyto is NOT included  in this version!'
   endif
!JM   CalcMacrophyto=0;
   CalcMacrophyto=.false.;
#endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !END compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"InitParam.f90","Param.nml")
101 call error_msg_prn(NML_READ,"InitParam.f90","Param_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end  subroutine InitParam

  end module
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
