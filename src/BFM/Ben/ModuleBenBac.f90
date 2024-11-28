!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenBac
!
! DESCRIPTION
!   The parameter value file for benthic bacteria (H1,H2)   
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  module mem_BenBac
!
! !USES:

  use global_mem
  use mem,  ONLY: iiH3

!  
!
! !AUTHORS
!   ERSEM group, HBB   
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
  ! Default all is public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! BenBac PARAMETERS (read from nml)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN)  :: p_qnc  ! Optimal Internal N quota
  real(RLEN)  :: p_qpc  ! Optimal Internal P quota
  real(RLEN)  :: p_qlnc  ! Optimal Internal N quota
  real(RLEN)  :: p_qlpc  ! Optimal Internal P quota
  real(RLEN)  :: p_q10  ! Q10
  real(RLEN)  :: p_sum  ! Potential uptake rate       (1/d)
  real(RLEN)  :: p_suhR6  ! Specific (high) uptake rate (1/d)
  real(RLEN)  :: p_sulR6  ! Specific (slow) uptake rate (1/d)
  real(RLEN)  :: p_sulR1 ! Specific uptake rate of Q1 (1/d)
  real(RLEN)  :: p_suhR1 ! Specific uptake rate of Q1 with optimal nutrient content(1/d)
  real(RLEN)  :: p_suR2 ! Specific uptake rate of Q2 
  real(RLEN)  :: p_cuR6np ! Prefferential uptake of nutrients of Q6
  real(RLEN)  :: p_pur  ! Fraction of uptake respired
  real(RLEN)  :: p_srr  ! Specific respiration        (1/d)
  real(RLEN)  :: p_qun  ! max. uptake of KIn (m3pw/mmoln/d)
  real(RLEN)  :: p_qup  ! max. uptake of KIp (m3pw/mmoln/d)
  real(RLEN)  :: p_lN3N4n ![NH4+] at which [NO3-]uptake are equal (BFM)
  real(RLEN)  :: p_lureaN4n ![NH4+] at which urea uptake is equal to [NH4+] uptake

  real(RLEN)  :: p_sd  ! Specific nutrient stress mortality (1/d)
  real(RLEN)  :: p_sd2  ! Density Dependent mortality  (1/d)
  real(RLEN)  :: p_thdo !internal rel.quotum where nutrient stress occur

  real(RLEN)  :: p_clO2o(1:iiH3)  ! limitation at low oxygen (p_clO3o>0 : aerobic, otherwise anoxic)
  integer     ::sw_an=0
  integer  :: p_iK4(1:iiH3)  ! BenthicAmmonium(p_IK4) =1,2 --> K4n K14n
  integer  :: p_iK1(1:iiH3)  ! BenthicPhosphate(p_IK1) =1,2 --> K1p K11p
  integer  :: p_iQ1(1:iiH3)  ! BenDetritus(p_IQ1) =1,2 --> Q1 Q11
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")

  public InitBenBac
  contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitBenBac()

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  namelist /BenBac_parameters/ p_qpc, p_qlpc, p_qnc, p_qlnc, p_q10, &
    p_sum, p_suhR6, p_sulR6, p_suhR1,p_sulR1,p_suR2,p_cuR6np, p_pur,&
    p_srr, p_qun, p_qup, p_lN3N4n,p_lureaN4n,p_sd,p_sd2,p_thdo,  &
    p_clO2o,sw_an, p_iK4, p_iK1, p_iQ1
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !BEGIN compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

   p_lureaN4n=1000.
   write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
   write(LOGUNIT,*) "#  Reading BenBac parameters.."
open(NMLUNIT,file='BenBac.nml',status='old',action='read',err=100)
    read(NMLUNIT,nml=BenBac_parameters,err=101)
    close(NMLUNIT)
    write(LOGUNIT,*) "#  Namelist is:"
    write(LOGUNIT,nml=BenBac_parameters)
    if ( sw_an<0.or.sw_an.gt.1) stop 'sw_an outside range 0-1'  
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !END compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"InitBenBac.f90","BenBac.nml")
101 call error_msg_prn(NML_READ,"InitBenBac.f90","BenBac_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end  subroutine InitBenBac
  end module mem_BenBac
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
