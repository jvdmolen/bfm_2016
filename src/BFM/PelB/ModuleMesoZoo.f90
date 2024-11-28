!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: MesoZoo
!
! DESCRIPTION
!   This file contains the parameter values for the mesozooplankton submodel.
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  module mem_MesoZoo
!
! !USES:

  use global_mem
  use mem,  ONLY: iiMesoZooPlankton, iiPhytoPlankton,iiMicroZooPlankton

!  
!
! !AUTHORS
!   N. Broekhuizen and A.D. Bryant, ERSEM group
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
  ! Default all is public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! MesoZoo PARAMETERS (read from nml)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !
  !  The variable Z3 represents carnivorous mesozooplankton and Z4 represents
  !  omnivorous mesozooplankton.
  !
  real(RLEN)  :: p_q10(iiMesoZooPlankton)  ! Doubling temperature
  real(RLEN)  :: p_srm(iiMesoZooPlankton)  ! RestRespiration rate at 10 degrees C
  real(RLEN)  :: p_srs(iiMesoZooPlankton)  ! Respiration rate at 10 degrees C
  real(RLEN)  :: p_sum(iiMesoZooPlankton)  ! Maximal productivity at 10 degrees C
  real(RLEN)  :: p_sd(iiMesoZooPlankton)  ! Background natural mortality
  real(RLEN)  :: p_vum(iiMesoZooPlankton)  ! Specific search volume
  real(RLEN)  :: p_pur(iiMesoZooPlankton)  ! activity respitation
  real(RLEN)  :: p_peuPI(iiMesoZooPlankton)  ! Faeces production
  real(RLEN)  :: p_peuR3(iiMesoZooPlankton)  ! Faeces production
  real(RLEN)  :: p_peuMIZ(iiMesoZooPlankton)  ! Faeces production
  real(RLEN)  :: p_peuMEZ(iiMesoZooPlankton)  ! Faeces production
  real(RLEN)  :: p_smd(iiMesoZooPlankton)  ! Fractional density-dependent mortality
  real(RLEN)  :: p_qpc(iiMesoZooPlankton)  ! Quotum phosphate
  real(RLEN)  :: p_qnc(iiMesoZooPlankton)  ! Quotum nitrate
  real(RLEN)  :: p_puPI(iiMesoZooPlankton,iiPhytoPlankton)
  real(RLEN)  :: p_puMIZ(iiMesoZooPlankton,iiMicroZooPlankton)
  real(RLEN)  :: p_puMEZ(iiMesoZooPlankton,iiMesoZooPlankton)
  real(RLEN)  :: p_clO2o(iiMesoZooPlankton)  ! Low oxygen
  integer     :: p_sw_nut_part(iiMesoZooPlankton) !sw to partioning of nutrients  between uptake/detritus
  integer     :: p_sw_faecelpell_sed(iiMesoZooPlankton) !sw to seton/off fast sedimetation of faecelpellets
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")

  public InitMesoZoo
  contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitMesoZoo()

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  namelist /MesoZoo_parameters/ p_q10, p_srm,p_srs, p_puPI, p_puMIZ, p_puMEZ, p_sd, &
    p_sum, p_vum, p_pur, p_peuPI,p_peuR3,p_peuMIZ,p_peuMEZ, &
    p_smd,p_qpc, p_qnc, p_clO2o, p_sw_nut_part,p_sw_faecelpell_sed

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !BEGIN compute

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
   write(LOGUNIT,*) "#  Reading MesoZoo parameters.."
open(NMLUNIT,file='MesoZoo.nml',status='old',action='read',err=100)
    read(NMLUNIT,nml=MesoZoo_parameters,err=101)
    close(NMLUNIT)
    write(LOGUNIT,*) "#  Namelist is:"
    write(LOGUNIT,nml=MesoZoo_parameters)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !END compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"InitMesoZoo.f90","MesoZoo.nml")
101 call error_msg_prn(NML_READ,"InitMesoZoo.f90","MesoZoo_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end  subroutine InitMesoZoo
  end module mem_MesoZoo
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
