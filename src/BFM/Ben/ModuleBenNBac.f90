!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenNbac
!
! DESCRIPTION
!   The parameter value file for benthicnitrifying  bacteria (H3)
!
!

!   This file is generated directly from OpenSesame model code, using a code
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  module mem_BenNBac
!
! !USES:

  use global_mem

!
!
! !AUTHORS
!   ERSEM group, HBB
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
!   the Free Software Foundation
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
  ! BenNBac PARAMETERS (read from nml)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  real(RLEN)  :: p_qunB  ! affinity for  N (m3/mgrC/d)
  real(RLEN)  :: p_qupB  ! affinity for  N (m3/mgrC/d)
  real(RLEN)  :: p_suB  ! max uptake rate (-)
  real(RLEN)  :: p_qunA  ! affinity for  N (m3/mgrC/d)
  real(RLEN)  :: p_qupA   ! max. uptake of KIp (m3/mmolP/d)
  real(RLEN)  :: p_suA  ! max uptake rate (-)
  real(RLEN)  :: p_qnNc  ! reverse of growth yield (mol N/mgC oxidixed (1/(1.32 grC/molN))
  real(RLEN)  :: p_pra  !cost for being active as a fraction of max.growthrate
  real(RLEN)  :: p_srA  ! maintenance
  real(RLEN)  :: p_srB  ! maintenance
  real(RLEN)  :: p_clO2o  ! half saturation constant (15 mmol O2/m3))
  real(RLEN)  :: p_clHCO3  ! half saturation constant (500 mmol HCO3/m3))
  real(RLEN)  :: p_qnc   ! Optimal Internal N quota
  real(RLEN)  :: p_qpc   ! Optimal Internal P quota
  real(RLEN)  :: p_qlnc   ! Optimal Internal N quota
  real(RLEN)  :: p_qlpc   ! Optimal Internal P quota
  real(RLEN)  :: p_q10   ! Q10
  real(RLEN)  :: p_sd   ! nutrient stress mortality          (1/d)
  real(RLEN)  :: p_sd2   ! density dependent mortality          (1/d)
  real(RLEN)  :: p_thdo
  real(RLEN)  :: p_lureaN4an ![NH4+] at which urea uptake is equal to [NH4+] uptake
  real(RLEN)  :: p_lureaN4bn ![NH4+] at which urea uptake is equal to [NH4+] uptake
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")

  public InitBenNBac
  contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitBenNBac()

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  namelist /BenNBac_parameters/p_qunB,p_qupB,p_suB, p_qunA,p_qupA,p_suA, &
    p_qnNc, p_pra,p_srA,p_srB, p_clO2o, p_clHCO3, p_qpc,p_qlpc, p_qnc,p_qlnc, p_q10, &
    p_sd,p_sd2,p_thdo,p_lureaN4an,p_lureaN4bn

  !BEGIN compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
   write(LOGUNIT,*) "#  Reading BenNBac parameters.."
open(NMLUNIT,file='BenNBac.nml',status='old',action='read',err=100)
    read(NMLUNIT,nml=BenNBac_parameters,err=101)
    close(NMLUNIT)
    write(LOGUNIT,*) "#  Namelist is:"
    write(LOGUNIT,nml=BenNBac_parameters)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !END compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"InitBenNBac.f90","BenNBac.nml")
101 call error_msg_prn(NML_READ,"InitBenNBac.f90","BenNBac_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end  subroutine InitBenNBac
  end module mem_BenNBac
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
