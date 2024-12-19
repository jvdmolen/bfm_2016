#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: MesoZoo
!
! DESCRIPTION
!   This submodel describes the carbon dynamics and associated
!    nutrient dynamics in carnivorous mesozooplankton (represented
!    by the state variable Z3) and in omnivorous zooplankton (in
!    the model known as Z4).
!
!    sparate handling of LOC and Detritus
!

!   This file is generated directly from OpenSesame model code, using a code
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine MesoZooDynamics(zoo,  ppzooc, ppzoon, ppzoop)
!
! !USES:

  ! For the following Pelagic-states fluxes are defined: O2o, N1p, N4n, R6c, &
  ! R6p, R6n
  ! For the following Pelagic-group-states fluxes are &
  ! defined: PhytoPlankton, MicroZooPlankton, MesoZooPlankton
  ! The following Pelagic 1-d global boxvars are modified : flPIR6s
  ! The following Pelagic 1-d global boxvars  are used: ETW
  ! The following Pelagic 2-d global boxvars are used: &
  ! qn_mz, qp_mz, qnZc, qpZc
  ! The following groupmember vars are used: iiPhytoPlankton, &
  ! iiMicroZooPlankton, iiMesoZooPlankton
  ! The following constituent constants  are used: iiC, iiN, iiP, iiL
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO,DONE,NZERO
#ifdef NOPOINTERS
  use mem,  ONLY: D3STATE
#else
  use mem, ONLY: D3STATE, O2o,R3c, &
    PhytoPlankton, MicroZooPlankton, MesoZooPlankton
#endif
  use mem,ONLY: ppO2o, ppO3c, ppN1p, ppR1c, ppR1n,ppR1p, ppRZc, ppR6c, ppR6p, &
    ppR6n, iiC, iiN, iiP, iiL,iiS,iiZ4,iiZ3,iiZ2,iiP6,ppR3c,iiYy3, &
    iiPhytoPlankton, iiMicroZooPlankton, iiMesoZooPlankton,  &
    ppPhytoPlankton,ppMicroZooPlankton,ppMesoZooPlankton, NO_BOXES,&
    qn_mz, qp_mz, qnZc, qpZc, jnetYy3c, jnetMeZc,jmY3c,flPIR6s, ETW,  &
    jPelFishInput, iiPel, flux_vector,fixed_quota_flux_vector,Depth, &
    pMIupZ4,sediMeZ,jrrMec,jZIR6n,jZIR6p,jZIDIn,jZIDIp,jPTZTn,jPTZTp
  use mem_Param,  ONLY: check_fixed_quota,p_pe_R1c, p_pe_R1n, p_pe_R1p, &
                   p_peZ_R1c, p_peZ_R1n, p_peZ_R1p, &
                   CalcMesoZooPlankton,CalcMicroZooPlankton,CalcPhytoPlankton
  use constants,ONLY: p_qnUc,MW_C
  use mem_Phaeo,ONLY:CALC_FOOD_MESOZOO,CALC_GRAZING_MESOZOO
  use mem_MesoZoo

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following vector functions are used:eTq_vector, MM_vector
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: eTq_vector, MM_vector
  use global_interface,ONLY:PhaeocystisCalc,TopPredLosses

!use mem,  ONLY: Output2d_1,Output2d_2,Output2d_3, Output2d_4


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
!   N. Broekhuizen and A.D. Bryant, ERSEM group
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
  ! Set up Local Variable for copy of state var. object
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),dimension(:),pointer :: zooc
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer  :: i,iout
  integer  :: j
  real(RLEN),dimension(NO_BOXES)  :: put_u,put_zu
  real(RLEN),dimension(NO_BOXES)  :: cmuc,cmuac,cmuzc
  real(RLEN),dimension(NO_BOXES)  :: ceuc
  real(RLEN),dimension(NO_BOXES)  :: eo
  real(RLEN),dimension(NO_BOXES)  :: et
  real(RLEN),dimension(NO_BOXES)  :: rrsc,rrsn,rrsp
  real(RLEN),dimension(NO_BOXES)  :: rrac
  real(RLEN),dimension(NO_BOXES)  :: rugc,rugn,rugp
  real(RLEN),dimension(NO_BOXES)  :: psloppy
  real(RLEN),dimension(NO_BOXES)  :: rea1c, rea1n, rea1p
  real(RLEN),dimension(NO_BOXES)  :: rea6c, rea6n, rea6p
  real(RLEN),dimension(NO_BOXES)  :: reslc, resln, reslp
  real(RLEN),dimension(NO_BOXES)  :: ru_xfs_c, ru_xfs_n, ru_xfs_p
  real(RLEN),dimension(NO_BOXES)  :: rx_any,px_any
  real(RLEN),dimension(NO_BOXES)  :: flZIR1c, flZIR1n, flZIR1p
  real(RLEN),dimension(NO_BOXES)  :: rmc,rmn,rmp
  real(RLEN),dimension(NO_BOXES)  :: sm
  real(RLEN),dimension(NO_BOXES)  :: rmdc
  real(RLEN),dimension(NO_BOXES)  :: runc,runn,runp
  real(RLEN),dimension(NO_BOXES,iiPhytoPlankton)     :: cumPIc
  real(RLEN),dimension(NO_BOXES)  :: cumR3c
  real(RLEN),dimension(NO_BOXES,iiMicroZooPlankton)  :: cumMIZc
  real(RLEN),dimension(NO_BOXES,iiMesoZooPlankton)   :: cumMEZc
  real(RLEN),dimension(NO_BOXES)  :: ruPIc
  real(RLEN),dimension(NO_BOXES)  :: ruMIZc
  real(RLEN),dimension(NO_BOXES)  :: ruMEZc
  real(RLEN),dimension(NO_BOXES)  :: rq6c,rq6n,rq6p
  real(RLEN),dimension(NO_BOXES)  :: renc,renn,renp
  real(RLEN),dimension(NO_BOXES)  :: tfluxc,tfluxn,tfluxp
  real(RLEN),dimension(NO_BOXES)  :: net
  real(RLEN),dimension(NO_BOXES)  :: pvum
  real(RLEN),dimension(NO_BOXES)  :: limit_grazing
  real(RLEN),dimension(NO_BOXES)  :: rZ2Z3c
  real(RLEN),dimension(:),pointer :: phytoc,phytonps
  
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Copy  state var. object in local var
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  zooc => D3STATE(:,ppzooc)

  tfluxc=ZERO
  tfluxn=ZERO
  tfluxp=ZERO

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !Physiological temperature response
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  eo  =   MM_vector(  max(NZERO,O2o(:)),  p_clO2o(zoo))
  et  =   eTq_vector(  ETW(:),  p_q10(zoo))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate total food
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! rua: food uptake generalized. Addition of new FG become much more simpler!

  cmuc  = ZERO
  cmuac = ZERO
  cmuzc = ZERO
  ceuc  = ZERO
  cumR3c= ZERO
  do i = 1 , ( iiPhytoPlankton)
    if (CalcPhytoPlankton(i)) then
     phytoc=> PhytoPlankton(i,iiC)
     cumPIc(:,i)  =   p_puPI(zoo,i)* phytoc

     !output and input cumPIc input:p_puPI
     call PhaeocystisCalc(CALC_FOOD_MESOZOO,i,cumPIc(:,i),cumPIc(:,i), &
                                                           p_puPI(zoo,i))
     cmuc =  cmuc+ cumPIc(:,i)
     cmuac=  cmuac+ cumPIc(:,i)
     ceuc =  ceuc+ cumPIc(:,i) *p_peuPI(zoo)
     if ( i==iiP6) then
      !TEP material in colony`
       cumR3c=  cumPIc(:,i) *R3c/(NZERO+phytoc)
       cmuc  = cmuc+ cumR3c
       cmuac =   cmuac+ cumR3c
       ceuc  =   ceuc+ cumR3c *p_peuR3(zoo)
     endif
   else
     cumPIc(:,i) = ZERO
   endif
  end do
  do i = 1 , ( iiMicroZooPlankton)
    if (CalcMicroZooPlankton(i).and.p_puMiZ(zoo,i)>ZERO) then
      cumMIZc(:,i)  =   p_puMIZ(zoo,i)* MicroZooPlankton(i,iiC)
      cmuc =   cmuc +cumMIZc(:,i)
      cmuzc=   cmuzc+cumMIZc(:,i)
      ceuc =   ceuc +cumMIZc(:,i) *p_peuMIZ(zoo)
    else
      cumMIZc(:,i) = ZERO
    endif
  end do

  do i = 1 , ( iiMesoZooPlankton)
    if (CalcMesoZooPlankton(i).and.p_puMEZ(zoo,i)>ZERO) then
      cumMEZc(:,i) =   p_puMEZ(zoo,i)* MesoZooPlankton(i,iiC)
      cmuc =   cmuc +cumMEZc(:,i)
      cmuzc=   cmuzc+cumMEZc(:,i)
      ceuc =   ceuc +cumMEZc(:,i) *p_peuMEZ(zoo)
   else
      cumMEZc(:,i) = ZERO
   endif
  end do

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate total food uptake
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   if ( p_srm(zoo) > ZERO) then
     !efficiency 1-excreted part of food - avticity respiration
     px_any=DONE- (cmuc-ceuc)/(NZERO+cmuc) - p_pur(zoo)
     px_any=p_sum(zoo) *( px_any /(p_srs(zoo)) - DONE/(NZERO+ p_vum(zoo) *cmuc))
     pvum=min(DONE,max(1.0D-6,px_any))
   else
     pvum=DONE
   endif

  rugc  =eo*et* p_sum(zoo)* MM_vector(  pvum* p_vum(zoo)* cmuc,p_sum(zoo))* zooc

! foodsaturation=ZERO;
! where (rugc.gt.ZERO) &
! foodsaturation=max(ZERO,eo*et*zooc*(p_vum(zoo)*cmuc -p_sum(zoo)))/(rugc)
  psloppy=ZERO;
! where (cmuac.gt.ZERO) &
! !Calculate the food uptake part of phyto plankton whicnb will  neglectged
! ! for uptake and will be detritus.
! psloppy=foodsaturation/(NZERO+rugc)*cmuc/cmuac

  put_u = rugc/ ( NZERO + cmuc)
  put_zu=put_u


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Total Gross Uptake,excretion,sloppy feeding
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rugc  =   ZERO
  rugn  =   ZERO
  rugp  =   ZERO
  rea1c  =   ZERO;rea6c  =   ZERO; reslc  =   ZERO
  rea1n  =   ZERO;rea6n  =   ZERO; resln  =   ZERO
  rea1p  =   ZERO;rea6p  =   ZERO; reslp  =   ZERO

  do i = 1 , iiPhytoPlankton
    phytoc=> PhytoPlankton(i,iiC)
    if (CalcPhytoPlankton(i).and.p_puPI(zoo,i) > ZERO ) then
      ruPIc  =   (put_u +psloppy)* cumPIc(:,i)
      !output:  ( =rx_any)  input ruPIc input:p_puPI
      call PhaeocystisCalc(CALC_GRAZING_MESOZOO,i,rx_any,ruPIc,p_puPI(zoo,i))
      rugc  =   rugc+ ruPIc
      reslc =  reslc+ ruPIc           *psloppy
      rx_any = ruPIc*min(p_peuPI(zoo),DONE-p_pe_R1c)
      rea6c  = rea6c +  rx_any
      rea1c  = rea1c +  ruPIc*p_peuPI(zoo)-rx_any


      iout= fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzooc, &
                          ppPhytoPlankton(i,iiC),ppzooc, ruPIc ,tfluxc)
      phytonps=> PhytoPlankton(i,iiN)
      ru_xfs_n= ruPIc* phytonps/(NZERO+phytoc)
      rugn  =   rugn+ ru_xfs_n
      resln =  resln+ ru_xfs_n               *psloppy
      rx_any = ru_xfs_n*min(p_peuPI(zoo),DONE-p_pe_R1n)
      rea6n  = rea6n +  rx_any
      rea1n  = rea1n +  ru_xfs_n*p_peuPI(zoo)-rx_any
      iout= fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzoon, &
                ppPhytoPlankton(i,iiN),ppzoon, ru_xfs_n,tfluxn )
      jPTZTn(1)=jPTZTn(1)+ sum(rugn*Depth)
      phytonps=> PhytoPlankton(i,iiP)
      ru_xfs_p= ruPIc* phytonps/(NZERO+phytoc)
      rugp  =   rugp+ ru_xfs_p
      reslp =  reslp+ ru_xfs_p            *psloppy
      rx_any = ru_xfs_p*min(p_peuPI(zoo),DONE-p_pe_R1p)
      rea6p  = rea6p +  rx_any
      rea1p  = rea1p +  ru_xfs_p*p_peuPI(zoo)-rx_any
      jPTZTp(1)=jPTZTp(1)+ sum(rugp*Depth)
      iout= fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzoop, &
              ppPhytoPlankton(i,iiP),ppzoop, ru_xfs_p,tfluxp )
      ! Chl is transferred to the sink
      phytonps=> PhytoPlankton(i,iiL)
      rx_any= ruPIc* phytonps/(NZERO+phytoc)
      call flux_vector(iiPel, &
              ppPhytoPlankton(i,iiL), ppPhytoPlankton(i,iiL),-rx_any )
      ! PIs is directly transferred to R6s
      j=ppPhytoPlankton(i,iiS)
      if (j>0) then
         phytonps=> PhytoPlankton(i,iiS)
         flPIR6s(:,i)  =   flPIR6s(:,i)+ ruPIc* phytonps/(NZERO+phytoc)
      endif
      if ( i==iiP6) then
        ru_xfs_c=ruPIc/(NZERO+cumPIc(:,i)) * cumR3c
        iout= fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzooc, &
                       ppR3c,ppzooc, ru_xfs_c ,tfluxc)
        rugc  =   rugc+ ru_xfs_c
        rea6c  =  rea6c+ ru_xfs_c *p_peuR3(zoo)
        reslc =  reslc+ ru_xfs_c  *psloppy
      endif
    else
       rx_any=ZERO
       call flux_vector(iiPel,ppPhytoPlankton(i,iiC),ppzooc, rx_any)
       if (ppzoon> 0) &
           call flux_vector(iiPel,ppPhytoPlankton(i,iiN),ppzoon, rx_any)
       if (ppzoop> 0) &
           call flux_vector(iiPel,ppPhytoPlankton(i,iiP),ppzoop, rx_any)
       if (i==iiP6) call flux_vector(iiPel,ppR3c,ppzooc, rx_any)
    endif
  end do
  if ( zoo== iiZ4) pMIupZ4=rugc

  do i = 1 , ( iiMicroZooPlankton)
    if (CalcMicroZooPlankton(i).and.p_puMiZ(zoo,i)>ZERO) then
      ruMIZc  =   put_zu* cumMIZc(:,i)
      rugc  =   rugc+ ruMIZc
      ru_xfs_n  =   + ruMIZc* qn_mz(:,i)
      rugn  =   rugn+ ru_xfs_n
      ru_xfs_p  =   + ruMIZc* qp_mz(:,i)
      rugp  =   rugp+ ru_xfs_p
      rx_any = ruMIZc*min(p_peuMIZ(zoo),DONE-p_peZ_R1c)
      rea6c  = rea6c +  rx_any
      rea1c  = rea1c +  ruMIZc*p_peuMIZ(zoo)-rx_any
      rx_any = ru_xfs_n*min(p_peuMIZ(zoo),DONE-p_peZ_R1n)
      rea6n  = rea6n +  rx_any
      rea1n  = rea1n +  ru_xfs_n*p_peuMIZ(zoo)-rx_any
      rx_any = ru_xfs_p*min(p_peuMIZ(zoo),DONE-p_peZ_R1p)
      rea6p  = rea6p +  rx_any
      rea1p  = rea1p +  ru_xfs_p*p_peuMIZ(zoo)-rx_any

    else
      ruMIZc=ZERO;ru_xfs_n=ZERO;ru_xfs_p=ZERO;
    endif
    iout= fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzooc, &
        ppMicroZooPlankton(i,iiC),ppzooc, ruMIZc,tfluxc )
    iout= fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzoon, &
        ppMicroZooPlankton(i,iiN),ppzoon, ru_xfs_n,tfluxn)
    iout= fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzoop, &
        ppMicroZooPlankton(i,iiP),ppzoop, ru_xfs_p,tfluxp)
  end do

  if ( zoo== iiZ4) pMIupZ4=(rugc-pMIupZ4)/(NZERO+put_u *cmuc)

  rZ2Z3c=ZERO
  do i = 1 , ( iiMesoZooPlankton)
    if (CalcMesoZooPlankton(i).and.p_puMEZ(zoo,i)>ZERO) then
      ruMEZc  =   put_zu* cumMEZc(:,i)
      ! intra-group predation is not computed
      if ( zoo.eq.i.and.zoo.eq.iiZ3)  &
            jPelFishInput(1)= jPelFishInput(1) + sum(ruMEZc*Depth)

      rugc  =   rugc+ ruMEZc
      ru_xfs_n  =   + ruMEZc* qnZc(:,i)
      rugn  =   rugn+ ru_xfs_n
      ru_xfs_p  =   + ruMEZc* qpZc(:,i)
      rugp  =   rugp+ ru_xfs_p
      rx_any = ruMEZc*min(p_peuMEZ(zoo),DONE-p_peZ_R1c)
      rea6c  = rea6c +  rx_any
      rea1c  = rea1c +  ruMEZc*p_peuMEZ(zoo)-rx_any
      rx_any = ru_xfs_n*min(p_peuMEZ(zoo),DONE-p_peZ_R1n)
      rea6n  = rea6n +  rx_any
      rea1n  = rea1n +  ru_xfs_n*p_peuMEZ(zoo)-rx_any
      rx_any = ru_xfs_p*min(p_peuMEZ(zoo),DONE-p_peZ_R1p)
      rea6p  = rea6p +  rx_any
      rea1p  = rea1p +  ru_xfs_p*p_peuMEZ(zoo)-rx_any
      if ( zoo.eq.iiZ3.and.i.eq.iiZ2 ) rZ2Z3c=ruMEZc-rea1c-rea6c
      if (i.eq.iiZ2) jmY3c(1,iiYy3)=jmY3c(1,iiYy3)+sum(ruMEZc*Depth)
    else
      ruMEZc=ZERO;ru_xfs_n=ZERO;ru_xfs_p=ZERO;
    endif
    if ( i/= zoo) then
       iout= fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzooc, &
                    ppMesoZooPlankton(i,iiC),ppzooc, ruMEZc, tfluxc )
       iout= fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzoon, &
            ppMesoZooPlankton(i,iiN),ppzoon, ru_xfs_n , tfluxn)
       iout= fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzoop, &
            ppMesoZooPlankton(i,iiP),ppzoop, ru_xfs_p , tfluxp)
    end if
    if (i.ne.iiZ3) jnetMeZc(1)=jnetMeZc(1)-sum(Depth(:)*ruMEZc)
  end do


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Respiration and basal metabolism
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rrac  =   p_pur(zoo)* rugc

  rrsc  =   max(pvum*p_srs(zoo),p_srm(zoo))* et*eo* zooc
  rrsn  =   rrsc * qnZc(:,zoo)
  rrsp  =   rrsc * qpZc(:,zoo)

  jrrMec(1)= jrrMec(1)+sum(rrsc*Depth)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Assimilated material
  ! Respectively Carbon, Nitrogen and Phosphorus
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  runc  =    max(ZERO, rugc -rea1c -rea6c -reslc -rrac)
  runn  =    max(ZERO, rugn -rea1n -rea6n -resln +rrsn )
  runp  =    max(ZERO, rugp -rea1p -rea6p -reslp +rrsp )

! ! calculate how C will be excreted
  renc  =   max(ZERO,(-runn/( NZERO+ runc)/ p_qnc(zoo)+DONE),   &
                    (-runp/( NZERO+ runc)/ p_qpc(zoo)+DONE))* runc

  !Correct excretion of renn for the fact that for ecretion some C is needed
  ! for excretion as urea.
  runc=runc-renc
  renn=max(ZERO,runn-p_qnc(zoo)*runc) /(DONE-p_qnc(zoo)/p_qnUc)
  runc=runc-renn/p_qnUc
  renp=max(ZERO,runp-p_qpc(zoo)*runc)
  call findnan(renp,NO_BOXES,iout)
  rx_any =   rea6c+  renc
  if (p_sw_faecelpell_sed(zoo)==1)  &
                    call flux_vector( iiPel, ppRZc,ppRZc, rx_any)
  if ( iout>0) then
         write(logunit,*) 'MesoZoo :NAN in renp layer',iout
         write(LOGUNIT,*) 'zooc=',zooc(iout)
         write(LOGUNIT,*) 'renn=',renn(iout)
         write(LOGUNIT,*) 'runp=',runp(iout)
         write(LOGUNIT,*) 'renc=',renc(iout)
         write(LOGUNIT,*) 'runc=',runc(iout)
         write(LOGUNIT,*) 'rugc=',rugc(iout)
         write(LOGUNIT,*) 'rea1c,rea6c=',rea1c(iout),rea6c(iout)
         write(LOGUNIT,*) 'reslc=',reslc(iout)
         write(LOGUNIT,*) 'rrac=',rrac(iout)
         write(LOGUNIT,*) 'put_u=',put_u(iout)
         write(LOGUNIT,*) 'cmuc=',cmuc(iout)
         write(LOGUNIT,*) 'et*eo=',et(iout)*eo(iout)
  endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Natural mortality + low oxygen mortality
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rmc  =   (p_sd(zoo)+ p_srs(zoo)*(DONE-eo))* zooc *et
  rmn  =   (p_sd(zoo)+ p_srs(zoo)*(DONE-eo))* zooc *et * qnZc(:,zoo)
  rmp  =   (p_sd(zoo)+ p_srs(zoo)*(DONE-eo))* zooc *et*  qpZc(:,zoo)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Density dependent mortality
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  sm   =   p_smd(zoo)* zooc
  rmdc  =   sm * zooc
! rmdn  =   rmdc * qnZc(:,zoo)
! rmdp  =   rmdc * qpZc(:,zoo)

  if (zoo.eq.iiZ2) jmY3c(1,iiYy3)=jmY3c(1,iiYy3)+sum((rmc+rmdc)*Depth)

  jPelFishInput(1)= jPelFishInput(1) + sum((rmdc+rmc)*Depth)

  flZIR1c=ZERO
  flZIR1n=ZERO
  flZIR1p=ZERO
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! sloppy feeding
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
  flZIR1c=flZIR1c+rmc*p_peZ_R1c+reslc*p_pe_R1c+rea1c
  flZIR1n=flZIR1n+rmn*p_peZ_R1n+resln*p_pe_R1n+rea1n
  flZIR1p=flZIR1p+rmp*p_peZ_R1p+reslp*p_pe_R1p+rea1p

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Fluxes for eliminated excess nutrients
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rq6c  =   rmc*(DONE-p_peZ_R1c)+ rea6c+ reslc*(DONE-p_pe_R1c)
  rq6n  =   rmn*(DONE-p_peZ_R1n)+ rea6n+ resln*(DONE-p_pe_R1n)
  rq6p  =   rmp*(DONE-p_peZ_R1p)+ rea6p+ reslp*(DONE-p_pe_R1p)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! prepare flow statements
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  !nitrogen is excreted as urea
  !add organic exretion to loc production due to sloppy feeding
  flZIR1n=flZIR1n+renn
  flZIR1c=flZIR1c+renn/p_qnUc

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! flow statements
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   rx_any=rq6c
   call TopPredLosses(iiPel,NO_BOXES,sm,0.6D+00, &
        zooc,zooc*qnZc(:,zoo),zooc*qpZc(:,zoo), &
        rrac,flZIR1n,flZIR1c,renp,rq6c,rq6n,rq6p)
   call flux_vector( iiPel, ppRZc,ppRZc, rq6c-rx_any)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! flow statements
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! jPLO3c(1)=jPLO3c(1)+sum((rrac+rrsc)*Depth)
! Output2d_3(1)=Output2d_3(1)+rrac(NO_BOXES)+rrsc(NO_BOXES)

  call flux_vector( iiPel, ppO2o,ppO2o,- (rrac+rrsc)/ MW_C )
  iout= fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzooc, &
                             ppzooc,ppO3c,rrac+rrsc,tfluxc )
  iout= fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzooc, &
                             ppzooc,ppR1c, flZIR1c,tfluxc)
  iout= fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzoon, &
                             ppzoon,ppR1n, flZIR1n ,tfluxn)
  iout= fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzoop, &
                             ppzoop,ppR1p, flZIR1p ,tfluxp)
  !phosphate is excreted as inorganic dissolved
  iout= fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzoop, &
                             ppzoop,ppN1p, renp ,tfluxp)

  iout= fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzooc, &
                             ppzooc,ppR6c, rq6c+renc ,tfluxc )
  iout= fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzoon, &
                             ppzoon,ppR6n, rq6n ,tfluxn)
  iout= fixed_quota_flux_vector( check_fixed_quota, iiPel, ppzoop, &
                             ppzoop,ppR6p, rq6p ,tfluxp)

  jZIDIn= jZIDIn+sum(Depth(:)*renn)  ! output as urea!
  jZIDIp= jZIDIp+sum(Depth(:)*renp)

  jZIR6n= jZIR6n+sum(Depth(:)*rq6n)
  jZIR6p= jZIR6p+sum(Depth(:)*rq6p)

  net=rugc-rrac-rea1c-rea6c-renc
  if ( zoo.ne.iiZ2) then
    jnetMeZc(1)=jnetMeZc(1)+sum(Depth(:)*net)
    !Correction for filterfeeder larvae eaten by carnivorous microzooplankton
    rx_any= net*rZ2Z3c/(NZERO+rugc-rea1c-rea6c)
    jnetYy3c(1)=jnetYy3c(1)-sum(Depth(:)*rx_any)
  else
    !feeding of benthic larvae is added to the filterfeeder production.
    jnetYy3c(1)=jnetYy3c(1)+sum(Depth(:)*net)
    !decreases the filterfeeder larvae will hesitate
    !to go back to benthos in case of presence of much (pelagic) food...
    px_any=DONE/(DONE+p_sum(zoo)/(NZERO+p_vum(zoo) *cmuc))
    sediMeZ(:,zoo)=sediMeZ(:,zoo)*min(DONE,px_any)
  endif

  renn=tfluxc*p_qnc(zoo)
  iout= fixed_quota_flux_vector( check_fixed_quota,-iiN,0,0,0,renn,tfluxn,"Mesozoo")
  if (iout>0) write(LOGUNIT,*) 'zooc=',zooc
  renp=tfluxc*p_qpc(zoo)
  iout= fixed_quota_flux_vector( check_fixed_quota,-iiP,0,0,0,renp,tfluxp,"Mesozoo")
  if (iout>0) write(LOGUNIT,*) 'zooc=',zooc

  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
