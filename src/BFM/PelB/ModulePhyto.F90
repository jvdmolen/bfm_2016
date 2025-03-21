#include "DEBUG.h"
#include "INCLUDE.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Phyto
!
! DESCRIPTION
!   Parameter values for the phytoplankton groups
!
!

!   This file is generated directly from OpenSesame model code, using a code
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  module mem_Phyto
!
! !USES:

  use global_mem
  use mem,  ONLY: iiPhytoPlankton,NO_BOXES
  use mem_Param,ONLY:CalcPhytoPlankton

!
!
! !AUTHORS
!   the ERSEM group, Marcello Vichi, JWB, HBB
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
  ! Phyto PARAMETERS (read from nml)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !
  !  ---------------- Physiological parameters -----------------
  !
  real(RLEN)  :: p_q10(iiPhytoPlankton)   ! Doubling temperature
  real(RLEN)  :: p_sum(iiPhytoPlankton)   ! Maximal productivity at 10 degrees C
  real(RLEN)  :: p_srm(iiPhytoPlankton)   ! Respiration rate at 10 degrees C
  real(RLEN)  :: p_srs(iiPhytoPlankton)   ! Respiration rate at 10 degrees C
  real(RLEN)  :: p_ses(iiPhytoPlankton)   ! Excretion rate at 10 degrees C
  real(RLEN)  :: p_sdmo(iiPhytoPlankton)  ! Max.specific nutrient-stress lysis rate
  real(RLEN)  :: p_thdo(iiPhytoPlankton)  ! Half value for nutrient stress lysis
  real(RLEN)  :: p_sd2(iiPhytoPlankton)   ! Desnity dependent mortality
  real(RLEN)  :: p_pu_ea(iiPhytoPlankton) ! Fraction of pp excreted as PLOC/PDET
  real(RLEN)  :: p_pu_ra(iiPhytoPlankton) ! Activity respiration rate
  !
  !  ---------------- Nutrient parameters in phytoplankton -----------------
  !
  integer     :: p_iRI(iiPhytoPlankton)  ! excretion as sugars or as TEP
  real(RLEN)  :: p_qnlc(iiPhytoPlankton) !miniumum quoata
  real(RLEN)  :: p_lqnlc(iiPhytoPlankton) !quoata where below limiting of chla systhere
  real(RLEN)  :: p_qnRc(iiPhytoPlankton) !optimal quoata
  real(RLEN)  :: p_xqn(iiPhytoPlankton)
  real(RLEN)  :: p_qplc(iiPhytoPlankton)
  real(RLEN)  :: p_qpRc(iiPhytoPlankton)
  real(RLEN)  :: p_xqp(iiPhytoPlankton)
  real(RLEN)  :: p_qslc(iiPhytoPlankton)  ! Minimum quotum Si in PI
  real(RLEN)  :: p_qsRc(iiPhytoPlankton)  ! Reference quotum Si in PI
  real(RLEN)  :: p_xqs(iiPhytoPlankton)
  real(RLEN)  :: p_qun(iiPhytoPlankton)
  real(RLEN)  :: p_qup(iiPhytoPlankton)
  real(RLEN)  :: p_qus(iiPhytoPlankton)  ! affinity of PI for Si
  real(RLEN)  :: p_lN3N4n(iiPhytoPlankton)   ! NH4-concentration where limitation of NO3 uptake to 0.5
  real(RLEN)  :: p_lureaN4n(iiPhytoPlankton)! NH4-concentration where limitation of urea uptake to 0.5
  real(RLEN)  :: p_lN1(iiPhytoPlankton)
  real(RLEN)  :: p_chPs(iiPhytoPlankton)
  real(RLEN)  :: p_Ke(iiPhytoPlankton)  ! Initial slope P-I curve
  real(RLEN)  :: p_clO2o(iiPhytoPlankton)

  !
  !  ------------- Chlorophyll parameters -----------
  !  skel: Skeletonema costatum pav: Pavlova lutheri
  !  syn: Synechoccus sp. (significant alpha decrease with irradiance)
  !  gyr: Gyrodinium sp. iso: Isochrysis galbana
  !              skel     iso      syn      gyr
  real(RLEN)  :: p_sdchl_l(iiPhytoPlankton)  ! Specific turnover rate for Chla [d-1] in light
  real(RLEN)  :: p_sdchl_d(iiPhytoPlankton)  ! Specific turnover rate for Chla [d-1] in dark
  ! p_qchlc =    0.03,    0.025,   0.1,     0.02    # Maximum quotum Chla:C
  !             +-0.024  +-0.001  +-0.003  +-0.004
  !  Thalassiosira sp. [0.05+-0.01]
  !              skel     pav      syn      gyr
  !   p_qchlc  = 0.05,      0.03,      0.07,      0.02 # Maximum quotum Chla:C
  real(RLEN)   :: p_qchlc(iiPhytoPlankton)  ! Fixed/Maximum quotum Chla:C dependent on ChlLightFlag [mg Chla (mg C)-1]
  !                                          %p_qchlc%   0.05, 0.03,  0.07,  0.02
  real(RLEN)   :: p_qlPlc(iiPhytoPlankton)  ! Fixed/Maximum quotum Chla:C dependent on ChlLightFlag [mg Chla (mg C)-1]
  real(RLEN)   :: p_xdiv(iiPhytoPlankton)  ! max number of division per day
  real(RLEN)   :: p_xsize_m(iiPhytoPlankton) ! radius cell size in m
  real(RLEN)   :: p_xsize_c(iiPhytoPlankton) ! weight cell n mgC
  !                                         %p_qchlc%   0.05, 0.03,  0.07,  0.02
  integer      :: p_xgeneric(iiPhytoPlankton)  ! generic/specific
  real(RLEN)   :: p_EIR(iiPhytoPlankton)  ! ligt where above bleaching appears
  integer      :: TestStatus(iiPhytoPlankton) ! test of phytoplankton group can be exluded in 1d-model during run.
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Parameters valid only for one functional group/all functional groups:
  real        :: p_qnR2c=ZERO
  integer     :: p_limnut  ! switch for nut. limitation (Liebig is default)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  logical       :: CalcPhytoCopy(iiPhytoPlankton)
  integer       :: sw_old_CalcPhyto(iiPhytoPLankton)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")

  public InitPhyto,test_Phyto_status
  contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine InitPhyto()

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  implicit none
  integer    ::i,j
  real(RLEN),dimension(iiPhytoPlankton)   ::r

  namelist /Phyto_parameters/ p_q10, p_sum, p_srs,p_srm, p_ses,p_sdmo, p_sd2, &
    p_pu_ea, p_pu_ra, p_qnlc,p_lqnlc, p_qplc, p_qslc, p_qnRc, p_qpRc, p_qsRc, &
    p_qun, p_qup, p_qus, p_xqn, p_xqp, p_xqs, p_thdo, &
    p_lN3N4n,p_lureaN4n,p_lN1,p_chPs, p_iRI,p_kE, &
    p_sdchl_l, p_sdchl_d,p_clO2o,p_qchlc,p_qlPlc,p_xdiv, p_xgeneric, &
    p_xsize_m,p_xsize_c,p_EIR, &
    TestStatus, p_limnut,p_qnR2c
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !BEGIN compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Open the namelist file(s)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    p_EIR=300.0
    p_xgeneric=1
    TestStatus=0
write(LOGUNIT,*) "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
   write(LOGUNIT,*) "#  Reading Phyto parameters.."
open(NMLUNIT,file='Phyto.nml',status='old',action='read',err=100)
    read(NMLUNIT,nml=Phyto_parameters,err=101)
    close(NMLUNIT)
    do i=1,iiPhytoPlankton
       CalcPhytoCopy(i)=CalcPhytoPlankton(i)
       if ( p_qnRc(i) < p_lqnlc(i) .or. &
            p_lqnlc(i) < p_qnlc(i) .or. &
            p_qpRc(i) < p_qplc(i) .or. &
            p_qsRc(i) < p_qslc(i)  ) then
         write(LOGUNIT,*)  &
          'Error in the quoata parameters: a minimum quota is larger than the optimum'
         goto 101
       endif
       j= abs(p_xgeneric(i))
       if (j >3 .or.j<0) then
         write(LOGUNIT,*)  'Error in p_xgeneric: value outside range -1...1'
         write(LOGUNIT,*) ' 0= stress to temperature'
         write(LOGUNIT,*) ' 1= adaptive to temperature'
         write(LOGUNIT,*) ' 2= stress'
         write(LOGUNIT,*) ' 3=more stressa'
         goto 101
       endif
    enddo
    write(LOGUNIT,*) "#  Namelist is:"
    write(LOGUNIT,nml=Phyto_parameters)
    r=p_sum/p_qChlc/p_Ke/86400
    write(LOGUNIT,*)'alpha=(p_sum/p_qChlc/p_Ke/86400) (gC/(gChla d) umol/(m2 d)'
    write(LOGUNIT,'(100G24.14)') r
    r=r/24.0*86400
    write(LOGUNIT,*)'alpha=(p_sum/p_qChlc/p_Ke/24 )(gC/(gChla h) umol/(m2 s)'
    write(LOGUNIT,'(100G24.14)') r
    write(LOGUNIT,*)'Max.Spec.Rate of Photosynthesis P^(chl)_m=(p_sum/p_qchl/24) (gC/gChla /h)'
    r=p_sum/p_qchlc/24
    write(LOGUNIT,'(100G24.14)') r
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !END compute
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  return
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Error Messages
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
100 call error_msg_prn(NML_OPEN,"InitPhyto.f90","Phyto.nml")
101 call error_msg_prn(NML_READ,"InitPhyto.f90","Phyto_parameters")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end  subroutine InitPhyto
!-------------------------------------------------------------------------!
!BOP
!
! !ROUTINE: test_Phyto_status
!
! DESCRIPTION
!   With routine processes of or more Phytoplankton groups  can be set off
!   in case of too low numbers.
!   The reason is that a very low numbers the diffusion routine works in a way
!   that no mass conservation is kept, leading to a laot of warnings when running
!   the model.
!
!
! !INTERFACE
  subroutine test_Phyto_status
!
!
! !AUTHORS
!   the ERSEM group, Marcello Vichi, JWB, HBB
!
!
!
! !REVISION_HISTORY
!   !
!
!EOP
!-------------------------------------------------------------------------!
  use mem,ONLY: PhytoPlankton,iiPhytoPlankton,ppPhytoPlankton,iiC,iiM,iiS, &
      iiN,iiP,Depth,D3SOURCE,ppR3c,ppPcc,ppN1p,ppN3n,ppN4n,ppN5s,ppO3c, &
      ppR1c,ppR2c,ppR6c,ppR1n,ppR1p,ppR6s,iiP6,iiPel, &
      flux_vector,ResetSource_D3_vector,CoupledtoBDc,sw_CalcPhyto,D3SINK
  use global_mem,only:LOGUNIT,ZERO
#ifdef BFM_GOTM
  use bfm_output,only:var_names
#else
  use api_bfm,only: var_names
#endif

  use BFM_ERROR_MSG, ONLY: set_warning_for_getm

     implicit none
     real(RLEN)        :: scalar
     real(RLEN),dimension(NO_BOXES)        :: rzero
     integer           :: i,j,l,m,k

     rzero=ZERO
     do i=1,iiPhytoPlankton
       m=0
       sw_old_CalcPhyto(i)=sw_CalcPhyto(1,i)
       k=sw_CalcPhyto(1,i)
       j=ppPhytoPlankton(i,iiC)
       if ( CalcPhytoCopy(i)) then
         scalar=sum(PhytoPlankton(i,iiC)*Depth)
         if (k ==1) then
          !To overcome problems with the vertical diffusion calculation in gotm
          ! we stop calculating processes for this phytoplankton type for very
          ! low concentrations, with exception of the phytoplankton type is
          ! coupled to a benthic phytoplankton type which is "active"
          ! See further at ModuleBenPhyto.F90
           l=CoupledtoBDc(i)
           if (scalar.lt.1.0D-6.and.(l==0.or.l>iiPhytoPlankton) &
                                        .and.TestStatus(i)==1)then
!JM              write(LOGUNIT,*) 'Warning ',trim(var_names(j)), &
!JM              'processes are excluded in the calculation:concentration<1.d-06'
            !set all modelled constituents for this functional group on
!JM             call set_warning_for_getm
             sw_CalcPhyto(1,i)=0
             m=1
             k=0
             CalcPhytoPlankton(i)=.false.
           endif
         elseif( scalar.gt.1.0D-5) then
            k=1
            m=1
            sw_CalcPhyto(1,i)=1
            CalcPhytoPlankton(i)=.true.
            write(LOGUNIT,*) 'Warning ',trim(var_names(j)), &
           'processes are (re)included in the calculation:concentration>1.d-05'
            call set_warning_for_getm
            CalcPhytoPlankton(i)=.true.
         endif
       endif
       if ( m.eq.1.or.k.ne.sw_old_CalcPhyto(i)) then
         l=1
         do while (j.gt.0.and.l<=iiM)
            call ResetSource_D3_vector(j)
            if (l.eq.iiC) then
              if (ppO3c.gt.0) then
                call flux_vector(iiPel,ppO3c,j,rzero)
                call flux_vector(iiPel,j,ppO3c,rzero)
              endif
              call flux_vector(iiPel,j,ppR2c,rzero)
              call flux_vector(iiPel,j,ppR1c,rzero)
              call flux_vector(iiPel,j,ppR6c,rzero)
            endif
            if (l.eq.iiP)then
              call flux_vector(iiPel,ppN1p,j,rzero)
              call flux_vector(iiPel,j,ppN1p,rzero)
              call flux_vector(iiPel,ppR1p,j,rzero)
              call flux_vector(iiPel,j,ppR1p,rzero)
            endif
            if (l.eq.iiN)then
              call flux_vector(iiPel,ppN3n,j,rzero)
              call flux_vector(iiPel,j,ppN3n,rzero)
              call flux_vector(iiPel,ppR1n,j,rzero)
              call flux_vector(iiPel,j,ppR1n,rzero)
              call flux_vector(iiPel,ppN4n,j,rzero)
              call flux_vector(iiPel,j,ppN4n,rzero)
            endif
            if (l.eq.iiS)then
              call flux_vector(iiPel,ppN5s,j,rzero)
              call flux_vector(iiPel,j,ppN5s,rzero)
              call flux_vector(iiPel,j,ppR6s,rzero)
            endif
            l=l+1
            j=ppPhytoPlankton(i,l)
         enddo
         if (i==iiP6 ) then
             call ResetSource_D3_vector(ppR3c)
             call ResetSource_D3_vector(ppPcc)
         endif
       endif
     enddo
!    write(LOGUNIT,*) 'CalcPhyto;',CalcPhytoPlankton
!    endif
  end subroutine test_Phyto_status

  end module mem_Phyto
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
