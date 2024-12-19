#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: MicroZoo
!
! DESCRIPTION
!   !

!   This file is generated directly from OpenSesame model code, using a code
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine MicroZooDynamics(zoo,  ppzooc, ppzoon, ppzoop)
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO,NZERO,DONE
  use constants,ONLY:p_qnUc,MW_C
#ifdef NOPOINTERS
  use mem,  ONLY: D3STATE
#else
  use mem, ONLY: D3STATE, B1c,Bac,R3c,PhytoPlankton, MicroZooPlankton
#endif
  use mem, ONLY: ppBac,ppB1c, ppB1n, ppB1p, ppO2o, ppO3c, ppR6c, &
    ppR1n, ppR6n, ppR1p, ppR6p, ppR1c, ppR1n, ppN1p, ppPhytoPlankton, &
    ppMicroZooPlankton,Depth, NO_BOXES,flR3R2c,flBaZTc,jPTZTn,jPTZTp,jPLO3c, &
    flPIR6s, ETW, eO2mO2, qnB1c, qpB1c, qn_mz,qp_mz, jnetMiZc,jZIR6n,jZIR6p, &
    jZIDIn,jZIDIp,jrrMic, iiPhytoPlankton, iiMicroZooPlankton, &
    iiC,iiN,iiP,iiL,iiS, iiPel, iiP6, flux_vector,fixed_quota_flux_vector
  use mem_Param,ONLY: p_pe_R1c,p_pe_R1n,p_pe_R1p, &
     p_peZ_R1c,p_peZ_R1n,p_peZ_R1p, check_fixed_quota,CalcPhytoPlankton
  use mem_MicroZoo
  use mem_PelBac,only: p_version_PelBac=>p_version
  use mem_Phaeo, ONLY:CALC_FOOD_MICROZOO,CALC_GRAZING_MICROZOO

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following vector functions are used:eTq_vector, MM_vector
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: eTq_vector, MM_vector,insw_vector
  use global_interface,ONLY:PhaeocystisCalc

! use mem,  ONLY: Output2d_1,Output2d_2,Output2d_3, Output2d_4,iiZ5,Source_D3_vector

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

! !INPUT:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer,intent(IN)  :: zoo
  integer,intent(IN) :: ppzooc
  integer,intent(IN) :: ppzoon
  integer,intent(IN) :: ppzoop

!
!
! !AUTHORS
!   ERSEM group, Hanneke Baretta-Bekker
!
!
!
! !REVISION_HISTORY
!   by Piet Ruardij at Thu Mar 16 08:34:04 CET 2006
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
  ! Set up Local Variable for copy of state var. object
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),dimension(:),pointer :: zooc
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer  :: i,j,iout
  real(RLEN),dimension(NO_BOXES)  :: put_u
  real(RLEN),dimension(NO_BOXES)  :: et
  real(RLEN),dimension(NO_BOXES)  :: eO2
  real(RLEN),dimension(NO_BOXES)  :: ru_xfs_n,ru_xfs_p
  real(RLEN),dimension(NO_BOXES)  :: rumc,rumn,rump
  real(RLEN),dimension(NO_BOXES)  :: rugc,rugn,rugp
  real(RLEN),dimension(NO_BOXES)  :: runc,runn,runp
  real(RLEN),dimension(NO_BOXES)  :: efood
  real(RLEN),dimension(NO_BOXES)  :: rrsc
  real(RLEN),dimension(NO_BOXES)  :: rrac
  real(RLEN),dimension(NO_BOXES)  :: rea1c, rea6c
  real(RLEN),dimension(NO_BOXES)  :: rea1n, rea6n
  real(RLEN),dimension(NO_BOXES)  :: rea1p, rea6p
  real(RLEN),dimension(NO_BOXES)  :: rdc
  real(RLEN),dimension(NO_BOXES)  :: ruB1c
  real(RLEN),dimension(NO_BOXES)  :: ruPIc
  real(RLEN),dimension(NO_BOXES)  :: ruZIc
  real(RLEN),dimension(NO_BOXES)  :: rumR3c
  real(RLEN),dimension(NO_BOXES)  :: rumB1c
  real(RLEN),dimension(NO_BOXES)  :: rumBac
  real(RLEN),dimension(NO_BOXES,iiPhytoPlankton)  :: rumPIc
  real(RLEN),dimension(NO_BOXES,iiPhytoPlankton)  :: food
  real(RLEN),dimension(NO_BOXES,iiPhytoPlankton)  :: suPI
  real(RLEN),dimension(NO_BOXES,iiMicroZooPlankton)  :: rumZIc
  real(RLEN),dimension(NO_BOXES)  :: rr1c,rr6c
  real(RLEN),dimension(NO_BOXES)  :: rr1n,rr6n
  real(RLEN),dimension(NO_BOXES)  :: rr1p,rr6p
  real(RLEN),dimension(NO_BOXES)  :: renc,renn,renp
  real(RLEN),dimension(NO_BOXES)  :: rx_any,cx_any,sx_any
  real(RLEN),dimension(NO_BOXES)  :: tfluxc
  real(RLEN),dimension(NO_BOXES)  :: tfluxn
  real(RLEN),dimension(NO_BOXES)  :: tfluxp
  real(RLEN),dimension(:),pointer :: phytoc,phytonps,micro
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Copy  state var. object in local var
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  zooc => D3STATE(:,ppzooc)

  tfluxc=ZERO
  tfluxn=ZERO
  tfluxp=ZERO

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Temperature effect
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  et  =   eTq_vector( ETW(:),  p_q10(zoo))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Oxygen limitation
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! If there is saturation (eO2mO2i==1) eO2 is 1!
  eO2  = min(DONE,(DONE+ p_chro(zoo))  * MM_vector(  eO2mO2(:),   p_chro(zoo)))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Available food, etc...
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rumB1c=ZERO
  if ( p_version_PelBac==5) then
    cx_any=max(ZERO,B1c(:)-Bac(:))
    rumBac= p_suBa(zoo)* Bac * Bac  /( Bac + p_minfood(zoo))
    rumB1c=rumB1c+rumBac
  else
    cx_any=B1c(:)
  endif
  rumB1c= rumB1c + p_suB1(zoo)* cx_any * cx_any /( cx_any + p_minfood(zoo))
  rumc    =   rumB1c
  rumn    =   rumB1c* qnB1c(:)
  rump    =   rumB1c* qpB1c(:)

  food=ZERO
  do i = 1 , iiPhytoPlankton
    j=i;if (p_type(i)>0 ) j=p_type(i)
    sx_any=p_suPI(zoo,i)
    !if phyto is NOT phaeocystis output (r) is equal to input (r)
    !input and output=r input:p_suPI
    call PhaeocystisCalc(CALC_FOOD_MICROZOO,i,sx_any,sx_any,p_suPI(zoo,i))
    suPI(:,i)=sx_any
    phytoc => PhytoPlankton(i,iiC)
    food(:,j)=food(:,j)+sx_any*phytoc
  enddo
  do i = 1 , iiPhytoPlankton
    if (p_type(i)>0 ) food(:,i)=food(:,p_type(i))
  enddo
  do i = 1 , iiPhytoPlankton
    if ( p_suPI(zoo,i) .gt.ZERO.and.CalcPhytoPlankton(i)) then
      phytoc => PhytoPlankton(i,iiC)
      rumPIc(:,i) = suPI(:,i)* phytoc
      rumPIc(:,i) =rumPIc(:,i) * food(:,i)/(food(:,i) + p_minfood(zoo))
      rumc  =   rumc+ rumPIc(:,i)
      phytonps=>PhytoPlankton(i,iiN)
      rumn  =   rumn+ rumPIc(:,i)* phytonps/(NZERO+phytoc)
      phytonps=>PhytoPlankton(i,iiP)
      rump  =   rump+ rumPIc(:,i)* phytonps/(NZERO+phytoc)
      if ( i==iiP6) then
      !TEP material in colony`
        where (insw_vector(phytoc) >NZERO )
         rumR3c =  rumPIc(:,i) *R3c/(NZERO+ phytoc)
         rumc= rumc+ rumR3c
        endwhere
      endif
    else
      rumPIc(:,i) =  ZERO
    endif
 enddo

  do i = 1 , ( iiMicroZooPlankton)
    micro =>MicroZooPlankton(i,iiC)
    rumZIc(:,i) = p_suZI(zoo,i)* micro* micro/( micro+  p_minfood(zoo))
    rumc  =   rumc+ rumZIc(:,i)
    rumn  =   rumn+ rumZIc(:,i)* qn_mz(:,i)
    rump  =   rump+ rumZIc(:,i)* qp_mz(:,i)
  end do

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Uptake
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rumc=rumc
  efood  =   MM_vector(  rumc,  p_chuc(zoo))
  rugc  =   p_sum(zoo)* et* zooc* efood* eO2

  rrsc  =   p_srs(zoo)* et* zooc
  rugc =rugc *insw_vector(rugc*(DONE-p_pu_ea(zoo)-p_pu_ra(zoo))-rrsc)
  put_u  =   rugc/(NZERO+rumc)
  call findnan(put_u,NO_BOXES,iout)
  if ( iout>0) write(LOGUNIT,*) 'MicroZoo put_u=NAN', &
          iout,rumc(iout),rugc(iout),R3c(iout),B1c(iout)

  rugc=ZERO;rea1c=ZERO;rea6c=ZERO;
  rugn=ZERO;rea1n=ZERO;rea6n=ZERO;
  rugp=ZERO;rea1p=ZERO;rea6p=ZERO;
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Fluxes into microzooplankton
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ruB1c  =   put_u* rumB1c
  if ( p_version_PelBac==5) then
    call flux_vector( iiPel, ppBac,ppBac, -put_u* rumBac )
    flBaZTc=flBaZTc+put_u*rumBac
  endif
  rugc  =rugc+   ruB1c
  ru_xfs_n  =   ruB1c* qnB1c(:)
  rugn  =rugn+   ru_xfs_n
  ru_xfs_p  =   ruB1c* qpB1c(:)
  rugp  = rugp+  ru_xfs_p
  iout= fixed_quota_flux_vector( check_fixed_quota,iiPel, ppzooc,ppB1c,ppzooc, &
                                                                 ruB1c ,tfluxc)
  iout= fixed_quota_flux_vector( check_fixed_quota,iiPel, ppzoon,ppB1n,ppzoon, &
                                                       ru_xfs_n,tfluxn)
  iout= fixed_quota_flux_vector( check_fixed_quota,iiPel, ppzoop,ppB1p,ppzoop, &
                                                       ru_xfs_p,tfluxp)
  rx_any=   ruB1c*min(p_pu_ea(zoo),DONE-p_peZ_R1c)
  rea6c =   rea6c+rx_any
  rea1c =   rea1c+ruB1c*p_pu_ea(zoo)-rx_any
  rx_any=   ru_xfs_n*min(p_pu_ea(zoo),DONE-p_peZ_R1n)
  rea6n =   rea6n+rx_any
  rea1n =   rea1n+ru_xfs_n*p_pu_ea(zoo)-rx_any
  rx_any=   ru_xfs_p*min(p_pu_ea(zoo),DONE-p_peZ_R1p)
  rea6p =   rea6p+rx_any
  rea1p =   rea1p+ ru_xfs_p*p_pu_ea(zoo)-rx_any

  do i = 1 , iiPhytoPlankton
    if (  p_suPI(zoo,i) > ZERO.and.CalcPhytoPlankton(i) ) then
      ruPIc  =   put_u* rumPIc(:,i)
      !if phyto is NOT phaeocystis output (ruPIc) is equal to input (ruPIc)
      !rx_any=dummy output ruPIc=input
      call PhaeocystisCalc(CALC_GRAZING_MICROZOO,i,rx_any,ruPIc,p_suPI(zoo,i))
      rugc  =   rugc+ ruPIc
      phytoc => PhytoPlankton(i,iiC);phytonps=>PhytoPlankton(i,iiN)
      ru_xfs_n=ruPIc* phytonps/(NZERO+phytoc)
      rugn  =   rugn+ ru_xfs_n
      phytonps=>PhytoPlankton(i,iiP)
      ru_xfs_p=ruPIc* phytonps/(NZERO+phytoc)
      rugp  =   rugp+ ru_xfs_p
      iout= fixed_quota_flux_vector( check_fixed_quota,iiPel, ppzooc, &
                           ppPhytoPlankton(i,iiC),ppzooc, ruPIc ,tfluxc)
      iout= fixed_quota_flux_vector( check_fixed_quota,iiPel, ppzoon, &
                           ppPhytoPlankton(i,iiN), ppzoon, ru_xfs_n,tfluxn)
      iout= fixed_quota_flux_vector( check_fixed_quota,iiPel, ppzoop, &
                    ppPhytoPlankton(i,iiP),ppzoop, ru_xfs_p,tfluxp)
      jPTZTn(1)=jPTZTn(1)+ sum(ru_xfs_n*Depth)
      jPTZTp(1)=jPTZTp(1)+ sum(ru_xfs_p*Depth)
      ! Chl is transferred to the sink
      phytonps=>PhytoPlankton(i,iiL)
      rx_any=ruPIc* phytonps/(NZERO+phytoc)
      call flux_vector(iiPel,ppPhytoPlankton(i,iiL),ppPhytoPlankton(i,iiL), -rx_any)
      j=ppPhytoPlankton(i,iiS)
      if ( j.gt.0 ) then
        ! P1s is directly transferred to R6s
        ! PhytoPlankton[i].s -> R6.s = ruPIc * qsPc[i]
        phytonps=>PhytoPlankton(i,iiS)
        flPIR6s(:,i)  =   flPIR6s(:,i)+ ruPIc* phytonps/(NZERO+phytoc)
      end if
      rx_any = ruPIc*min(p_pu_ea(zoo),DONE-p_pe_R1c)
      rea6c  = rea6c +  rx_any
      rea1c  = rea1c +  ruPIc*p_pu_ea(zoo)-rx_any
      rx_any = ru_xfs_n*min(p_pu_ea(zoo),DONE-p_pe_R1n)
      rea6n  = rea6n +  rx_any
      rea1n  = rea1n +  ru_xfs_n*p_pu_ea(zoo)-rx_any
      rx_any = ru_xfs_p*min(p_pu_ea(zoo),DONE-p_pe_R1p)
      rea6p  = rea6p +  rx_any
      rea1p  = rea1p +  ru_xfs_p*p_pu_ea(zoo)-rx_any
      if ( i==iiP6) then
        !assumed is that all R3 not selected as food.
        flR3R2c=flR3R2c+ruPIc/(NZERO+rumPIc(:,i)) * rumR3c
        call findnan(flR3R2c,NO_BOXES,iout)
        if ( iout>0) write(LOGUNIT,*) 'MicroZoo flR3R2=NAN', &
                 iout,rumR3c(iout),ruPIc(iout),rumPIc(i,iout),put_u(iout)
      endif
   else
     rx_any=ZERO
     call flux_vector(iiPel,ppPhytoPlankton(i,iiC),ppzooc, rx_any)
     if (ppzoon>0) call flux_vector(iiPel,ppPhytoPlankton(i,iiN),ppzoon, rx_any)
     if (ppzoop>0) call flux_vector(iiPel,ppPhytoPlankton(i,iiP),ppzoop, rx_any)
    endif

  end do


  do i = 1 , ( iiMicroZooPlankton)
    ruZIc  =   put_u* rumZIc(:,i)
    rugc  =   rugc+ ruZIc
    ru_xfs_n  =ruZIc* qn_mz(:,i)
    rugn  =   rugn+ ru_xfs_n
    ru_xfs_p  =ruZIc* qp_mz(:,i)
    rugp  =   rugp+ ru_xfs_p
    ! intra-group predation is not computed
    if ( i/= zoo) then
      iout= fixed_quota_flux_vector( check_fixed_quota, &
                iiPel,ppzooc,ppMicroZooPlankton(i,iiC),ppzooc, ruZIc,tfluxc )
      iout= fixed_quota_flux_vector( check_fixed_quota, &
       iiPel,ppzoon,ppMicroZooPlankton(i,iiN), ppzoon,ru_xfs_n,tfluxn)
      iout= fixed_quota_flux_vector( check_fixed_quota, &
       iiPel,ppzoop,ppMicroZooPlankton(i,iiP), ppzoop,ru_xfs_p,tfluxp)
    end if

    rx_any = ruZIc*min(p_pu_ea(zoo),DONE-p_peZ_R1c)
    rea6c  = rea6c +  rx_any
    rea1c  = rea1c +  ruZIc*p_pu_ea(zoo)-rx_any
    rx_any = ru_xfs_n*min(p_pu_ea(zoo),DONE-p_peZ_R1n)
    rea6n  = rea6n +  rx_any
    rea1n  = rea1n +  ru_xfs_n*p_pu_ea(zoo)-rx_any
    rx_any = ru_xfs_p*min(p_pu_ea(zoo),DONE-p_peZ_R1p)
    rea6p  = rea6p +  rx_any
    rea1p  = rea1p +  ru_xfs_p*p_pu_ea(zoo)-rx_any
    jnetMiZc(1)=jnetMiZc(1)-sum(Depth(:)*ruZIc)
  end do

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Rest respiration
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  jrrMic(1)= jrrMic(1)+sum(rrsc*Depth)


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !       Fluxes from microzooplankton
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    !  activity, total respiration fluxes
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    rrac  =   rugc* p_pu_ra(zoo)
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Dissolved nutrient dynamics
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    runc  =   max(  ZERO,  rugc-rea1c-rea6c-rrac)
    runn  =   max(  ZERO,  rugn-rea1n-rea6n+ rrsc* qn_mz(zoo, :))
    runp  =   max(  ZERO,  rugp-rea1p-rea6p+ rrsc* qp_mz(zoo, :))

    ! Carbon which will be excreted/respired when internal quota (N/C) are 
    ! below optimum. renc has te be renc >=rrsc
    renc  =   max(ZERO,(-runn/( NZERO+ runc)/ p_qnc(zoo)+DONE),   &
                    (-runp/( NZERO+ runc)/ p_qpc(zoo)+DONE))* runc

    ! take in consideration that a (small) part of runC is used to make urea
    runc=runc-renc
    renn=max(ZERO,runn-p_qnc(zoo)*runc)/(DONE-p_qnc(zoo)/p_qnUc)
    runc=runc-renn/p_qnUc
    renp=max(ZERO,runp-p_qpc(zoo)*runc)

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Mortality (rdc) + Excretion (reac)
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    rdc  =  (( DONE- eO2)* p_sdo(zoo)+ p_sd(zoo))* zooc

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Fluxes due to mortality and excretion
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    rr1c  =   rea1c + rdc*p_peZ_R1c
    rr6c  =   rea6c + rdc*(DONE-p_peZ_R1c)

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Organic Nitrogen dynamics
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    rr1n  =   rea1n + rdc* qn_mz(:,zoo)* p_peZ_R1n
    rr6n  =   rea6n + rdc* qn_mz(:,zoo)* (DONE-p_peZ_R1n)

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Organic Phosphorus dynamics
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    rr1p  =   rea1p + rdc* qp_mz(:,zoo)* p_peZ_R1p
    rr6p  =   rea6p + rdc* qp_mz(:,zoo)* (DONE-p_peZ_R1p)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Fluxes
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  iout= fixed_quota_flux_vector( check_fixed_quota,iiPel, ppzooc,ppzooc,ppO3c, &
                                                  rrac+rrsc,tfluxc )
  call flux_vector( iiPel, ppO2o,ppO2o,- (rrac+rrsc)/ MW_C )
  jPLO3c(1)=jPLO3c(1)+sum((rrac+rrsc)*Depth)

  iout= fixed_quota_flux_vector(check_fixed_quota,iiPel,ppzooc,ppzooc,ppR6c, &
                                              rr6c+renc,tfluxc)
  iout= fixed_quota_flux_vector(check_fixed_quota,iiPel,ppzoon,ppzoon,ppR6n, &
                                             rr6n ,      tfluxn)
  iout= fixed_quota_flux_vector(check_fixed_quota,iiPel,ppzoop,ppzoop,ppR6p, &
                                              rr6p ,      tfluxp)

  iout= fixed_quota_flux_vector(check_fixed_quota,iiPel,ppzooc,ppzooc,ppR1c, &
                                              rr1c+renn/p_qnUc ,tfluxc)
  iout= fixed_quota_flux_vector(check_fixed_quota,iiPel,ppzoon,ppzoon,ppR1n, &
                                              rr1n+ renn       ,tfluxn)
  iout= fixed_quota_flux_vector(check_fixed_quota,iiPel,ppzoop,ppzoop,ppR1p, &
                                              rr1p            ,tfluxp)
  iout= fixed_quota_flux_vector(check_fixed_quota,iiPel,ppzoop,ppzoop,ppN1p, &
                                              renp ,tfluxp)

  rx_any=rugc-rrac-rea1c-rea6c
  jnetMiZc(1)=jnetMiZc(1)+sum(Depth(:)*rx_any)
  jZIDIn(1)= jZIDIn(1)+sum(Depth(:)*(rr1n+renn))
  jZIDIp(1)= jZIDIp(1)+sum(Depth(:)*(rr1p+renp))
  jZIR6n(1)= jZIR6n(1)+sum(Depth(:)*rr6n)
  jZIR6p(1)= jZIR6p(1)+sum(Depth(:)*rr6p)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! test if assumed fixed-quota fits
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rx_any=tfluxc*p_qnc(zoo)
  iout= fixed_quota_flux_vector( check_fixed_quota, &
                             -iiN,0,0,0,rx_any,tfluxn,"MicroZoo")
  rx_any=tfluxc*p_qpc(zoo)
  iout= fixed_quota_flux_vector( check_fixed_quota, &
                              -iiP,0,0,0,rx_any,tfluxp,"MicroZoo")

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
