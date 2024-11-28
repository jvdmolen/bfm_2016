!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Phaeo
!
! DESCRIPTION
!   Parameter values for the Phaeocystis functinal group
!   +Integer which describe the different actions
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  module mem_Phaeo
!
! !USES:

  use global_mem
  use mem,  ONLY: NO_BOXES

!  
!
! !AUTHORS
!   the ERSEM group, Piet Ruardij
!
!
!
! !REVISION_HISTORY
!   !
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
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Default all is public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  public

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Phaeo PARAMETERS (read from nml)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !
  !  ---------------- Physiological parameters -----------------
  !
  ! Parameters valid only for one functional group/all functional groups:
  real(RLEN)  :: p_sP2P6 ! specific rate in which P2 transfer themselves to Phaeocystis colonies
  real(RLEN)  :: p_steep_init !steepness of the increase of mortality at unfavorable circumstances.
  real(RLEN)  :: p_cmP2P6n ! dissolved N concentration where below a optimal transfer from P2 to P6 is
  real(RLEN)  :: p_cmP2P6p ! dissolved P concentration where below a optimal transfer from P2 to P6 is.
  real(RLEN)  :: p_wP6c ! weight of one Phaeocystis colonial cell
  real(RLEN)  :: p_steep_mort !steepness of the increase of mortality at unfavorable circumstances.
  real(RLEN)  :: p_rsP6m  !sedimentation rate which become stress for pheao colonies.
  real(RLEN)  :: p_qnR3c !NC quotum in colony but outside cell 
  real(RLEN)  :: p_qpR3c !PC quotum in colony but outside cell 
  real(RLEN)  :: p_cuY3a(2) ! size range of Phaeo.as food for FilterFeeders (cell number in colony)
  real(RLEN)  :: p_cuYy3a(2) !size range of Phaeocystis as food for Young FilterF. (cell number in colony)
  real(RLEN)  :: p_cuZ5a(2) ! size range of Phaeocystis as food for Ciliates (number of cells in colony)
  real(RLEN)  :: p_cuZ4a(2) ! size range of Phaeoc. as food for Onniv. MesozooP. (cell number in colony)
  real(RLEN)  :: p_smxInC   ! mortality of celss in colony
  real(RLEN)  :: p_chnxInC  ! nr cells in colony at which rate p_smxInCa result in halvwe a rate 
  real(RLEN)  :: p_mxDepthm ! Depth at which mortality is equal p_sdmo due to mechanical stress
  real (RLEN) ::p_mSal ! Salinity where below mortality of Phaeo is large 
  integer     :: sw_select=1!if  Phaeo as food is smaller or large than according selection in p_cuY3,p_cuZ5a,p_cuZ4a
                            ! biomass of Pcc which is grazing-rate=input/pcu(1 or 2) in stead of rate=input/SizeCol
!JM  integer     :: sw_detloc=1 !TEP is distributed over LOCand det (sw_loc=0)
  integer     :: sw_detloc=0 !TEP is distributed over LOCand det (sw_loc=0)
                             !all TEP is put in R6 ( sw_loc=1)
                             !all TEP is put in R1 ( sw_loc=-1)
                              
  integer,parameter        :: CALC_MORTALITY=1
  integer,parameter        :: NEW_COLONIES=2
  integer,parameter        :: CALC_FOOD_FILTERFEEDER=3
  integer,parameter        :: CALC_FOOD_MESOZOO=4
  integer,parameter        :: CALC_FOOD_MICROZOO=5
  integer,parameter        :: CALC_FOOD_YOUNG_FILTERFEEDER=6
  integer,parameter        :: CALC_GRAZING_FILTERFEEDER=-3
  integer,parameter        :: CALC_GRAZING_MESOZOO=-4
  integer,parameter        :: CALC_GRAZING_MICROZOO=-5
  integer,parameter        :: CALC_GRAZING_YOUNG_FILTERFEEDER=-6
  integer,parameter        :: CALC_LIMIT_FILTERCAP=8
  integer,parameter        :: CALC_SEDIMENTATION=9
  integer,parameter        :: CALC_REL_RADIUS=10
  integer,parameter        :: CALC_REL_PHOSPHATE_UPTAKE=11
  integer,parameter        :: CALC_REL_NITRATE_UPTAKE=13
  integer,parameter        :: CALC_REL_AMMONIUM_UPTAKE=14
  integer,parameter        :: CALC_REL_UREA_UPTAKE=15
  integer,parameter        :: CALC_NET_NITROGEN_UPTAKE=16
  integer,parameter        :: CALC_NET_PHOSPHATE_UPTAKE=17
  integer,parameter        :: LIMIT_BOUND=18
  integer,parameter        :: CALC_MORTALITY_CELLS_IN_COLONY=20
  integer,parameter        :: FORCE_MORTALITY=21
  integer,parameter        :: COLONY_DEGRADATION=22
  integer,parameter        :: CALC_LOC_DET_FLUX=23
  integer,parameter        :: CALC_TEMPERATURE_DEPENDENCE=25

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")

  public InitPhaeo
  contains

 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitPhaeo()

! default value:
!JM  sw_detloc=0;
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  namelist /Phaeo_parameters/ p_sP2P6,p_steep_init,p_cmP2P6n,p_cmP2P6p,&
       p_wP6c,p_steep_mort,p_rsP6m, p_qnR3c,p_qpR3c, &
       p_cuY3a,p_cuYy3a,p_cuZ4a,p_cuZ5a,p_smxInC,p_chnxInC,p_mxDepthm, &
       p_mSal,sw_select, sw_detloc

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !BEGIN compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Default values the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   p_qnR3c =ZERO
   p_mxDepthm=DONE


  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

   write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
   write(LOGUNIT,*) "#  Reading Phaeo parameters.."
   open(NMLUNIT,file='Phaeo.nml',status='old',action='read',err=100)
    read(NMLUNIT,nml=Phaeo_parameters,err=101)
    close(NMLUNIT)
    write(LOGUNIT,*) "#  Namelist is:"
    write(LOGUNIT,nml=Phaeo_parameters)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !END compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"InitPhaeo.f90","Phaeo.nml")
101 call error_msg_prn(NML_READ,"InitPhaeo.f90","Phaeo_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end  subroutine InitPhaeo
  end module mem_Phaeo
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
