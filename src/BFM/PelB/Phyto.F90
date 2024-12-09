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
!   This process describes the dynamics of all phytoplankton
!    groups in the ERSEM model. The differences in behaviour
!    are expressed by differences in parameter-values only.
!
!
!

!   This file is generated directly from OpenSesame model code, using a code
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine PhytoDynamics(phyto, ppphytoc, ppphyton, ppphytop, ppphytos, &
    ppphytol)
!
! !USES:

  ! For the following Pelagic-states fluxes are defined: R6c, O2o, R2c, &
  ! N3n, N4n, N1p, R1n, R6n, R1p, R6p, N5s
  ! The following global scalar vars are used: SUNQ, ThereIsLight
  ! The following Pelagic 1-d global boxvars are modified : flPIR6s
  ! The following Pelagic 1-d global boxvars  are used: ETW, EIR, xEPS, Depth
  ! The following Pelagic 2-d global boxvars are modified : eiPII
  ! The following Pelagic 2-d global boxvars got a value: sunPI
  ! The following Pelagic 2-d global boxvars  are used: qpPc, qnPc, qlPc
  ! The following 0-d global parameters are used: &
  ! ChlLightFlag, LightForcingFlag
  ! The following global constants are used: RLEN
  ! The following constants are used: HOURS_PER_DAY

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO,DONE,NZERO
#ifdef NOPOINTERS
  use mem,  ONLY: D3STATE
#else
  use mem, ONLY: D3STATE, O2o,N5s,N1p,N3n,N4n,N5s
#endif
  use mem, ONLY:ppO2o,ppR1c,ppR1n,ppR1p,ppR2c,ppR2n,ppR3c,ppR6c, &
    ppN1p,ppN3n,ppN4n,ppN5s,ppN6r, ppO3c,iiR1,iiR2,iiR3,iiP6,iiP1, iiPel,   &
    SUNQ,ThereIsLight, ETW,EIR, xEPS,Nun,Nup, &
    jnetPIc,jrrPTc, rnetPTc, fr_lim_PI_n,fr_lim_PI_p,&
    Depth, eiPI, sunPI, sugPI, sdoPI, qpPc,qnPc,qsPc,qlPc,NO_BOXES, &
    rml,flux_vector,sourcesink_flux_vector,Source_D3_vector,&
    flPIR6n,flPIR1n,iNPI,max_change_per_step,flN3N4n, &
    flPIN4n,flPIR6p,flPIR1p,flPIR6s,flR1O3c,flR1N4n,fl_xgrazing_PIc, &
    LocalDelta,dry_z,jPLO3c,iiConsumption,iiC,max_rate_per_step, &
    ppMesoZooPlankton,iiMesoZooPlankton,ppMicroZooPlankton,iiMicroZooPlankton

#ifdef INCLUDE_PELCO2
  use mem,  ONLY: HCO3
#endif
  use constants,  ONLY: MW_C,HOURS_PER_DAY,p_qnUc,POSITIVE
  use mem_Param,  ONLY: p_pe_R1c, p_pe_R1n, p_pe_R1p,ChlLightFlag, &
                        LightForcingFlag, p_qon_nitri,p_dry_ben,p_qro,p_xeff_an

  use mem_Phyto
  use LimitRates, ONLY:LimitChange_vector,DoubleLimitChange_vector
  use mem_Phaeo,ONLY:CALC_MORTALITY,COLONY_DEGRADATION,CALC_NET_NITROGEN_UPTAKE, &
    CALC_MORTALITY_CELLS_IN_COLONY,CALC_LOC_DET_FLUX, &
    CALC_NET_PHOSPHATE_UPTAKE,CALC_REL_AMMONIUM_UPTAKE, &
    CALC_REL_NITRATE_UPTAKE,CALC_REL_UREA_UPTAKE,CALC_REL_PHOSPHATE_UPTAKE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following vector functions are used:eTq_vector, MM_vector, insw_vector
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: eTq_vector,  insw_vector,exp_limit
  use global_interface,ONLY:PhaeocystisCalc
  use SourceFunctions,only: Source_D3_withgroup
  use botflux,only:getbotflux_3d

! use mem,  ONLY: Output2d_1,Output2d_2,Output2d_3, Output2d_4

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

! !INPUT:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer,intent(IN)  :: phyto
  integer,intent(IN) :: ppphytoc
  integer,intent(IN) :: ppphyton
  integer,intent(IN) :: ppphytop
  integer,intent(IN) :: ppphytos
  integer,intent(IN) :: ppphytol

!
!
! !AUTHORS
!   ERSEM group + J.G. Baretta-Bekker + W.Ebenhoeh
!     P. Ruardij (NIOZ)
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
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Set up Local Variable for copy of state var. object
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),dimension(NO_BOXES):: phytoc
  real(RLEN),dimension(NO_BOXES):: phyton
  real(RLEN),dimension(NO_BOXES):: phytop
  real(RLEN),dimension(NO_BOXES):: phytos
  real(RLEN),dimension(NO_BOXES):: phytol
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer                         :: silica_control   ,iout
  real(RLEN)                      :: avPhyto,dry_pel,rscalar
  real(RLEN)                      :: sudc !max growth rate related to max.number of divions perday.
  real(RLEN),parameter            :: low_phytoc=1.0D-10
  real(RLEN),dimension(NO_BOXES)  :: s,ex_sw
  real(RLEN),dimension(NO_BOXES)  :: qx_any,px_any   !any quotum,part
  real(RLEN),dimension(NO_BOXES)  :: cx_any   !any concentration
  real(RLEN),dimension(NO_BOXES)  :: rx_any,sx_any !any rate,specific rate
  real(RLEN),dimension(NO_BOXES)  :: dl   !day length
  real(RLEN),dimension(NO_BOXES)  :: eo,et
  real(RLEN),dimension(NO_BOXES)  :: sumc
  real(RLEN),dimension(NO_BOXES)  :: sea
  real(RLEN),dimension(NO_BOXES)  :: set
  real(RLEN),dimension(NO_BOXES)  :: ses
  real(RLEN),dimension(NO_BOXES)  :: sdo
  real(RLEN),dimension(NO_BOXES)  :: sx_main,sx_act
  real(RLEN),dimension(NO_BOXES)  :: rugc
  real(RLEN),dimension(NO_BOXES)  :: sra,srs,srt
  real(RLEN),dimension(NO_BOXES)  :: slc
  real(RLEN),dimension(NO_BOXES)  :: rufc
  real(RLEN),dimension(NO_BOXES)  :: rupp,rupn,rups
  real(RLEN),dimension(NO_BOXES)  :: rbufp,rbufn,rbufs
  real(RLEN),dimension(NO_BOXES)  :: rumn,rump,rums,rum5s
  real(RLEN),dimension(NO_BOXES)  :: ru_xCol_n,ru_xCol_p
  real(RLEN),dimension(NO_BOXES)  :: misn,misp,miss,sx_grazing
  real(RLEN),dimension(NO_BOXES)  :: luxn,luxp,luxs
  real(RLEN),dimension(NO_BOXES)  :: tN
  real(RLEN),dimension(NO_BOXES)  :: iNIn,iN1p,iNIs
  real(RLEN),dimension(NO_BOXES)  :: qnoc,qpoc,qsoc
  real(RLEN),dimension(NO_BOXES)  :: eN5s
  real(RLEN),dimension(NO_BOXES)  :: rrtc,rrmc,rrx_an_c
  real(RLEN),dimension(NO_BOXES)  :: rr2c
  real(RLEN),dimension(NO_BOXES)  :: rr1c,rr1n,rr1p,reUc,reUn
  real(RLEN),dimension(NO_BOXES)  :: rr6c,rr6n,rr6p,rr6s
  real(RLEN),dimension(NO_BOXES)  :: runc,runn,runp,runs,run5s
  real(RLEN),dimension(NO_BOXES)  :: runun,run3n,run4n,runup,run1p
  real(RLEN),dimension(NO_BOXES)  :: rum4n,rum3n,rumun,rum1p,rumup
  real(RLEN),dimension(NO_BOXES)  :: Irr
  real(RLEN),dimension(NO_BOXES)  :: rho_Chl
  real(RLEN),dimension(NO_BOXES)  :: rate_Chl
  real(RLEN),dimension(NO_BOXES)  :: new_Chl
  real(RLEN),dimension(NO_BOXES)  :: flChydrate,rr3c,leftover
  real(RLEN),dimension(NO_BOXES)  :: x_light_limitation
  real(RLEN),dimension(NO_BOXES)  :: ex_limit
  real(RLEN),dimension(NO_BOXES)  :: sx_stress_degradation
  real(RLEN),dimension(NO_BOXES)  :: cqun3,cquR1n,cqup,lim_qu

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  silica_control =0 : no silica component present in cell
  !  silica_control =1 : external regulation of silica limitation & limitation
  !                       of carbon fixation under silica depletion
  !  silica_control =2 : internal regulation of silica limitation & excretion
  !                      of fixed carbon under nutrient stress
  !                      Process description based on:
  !                Growth physiology and fate of diatoms in the ocean: a review
  !                G.Sarthou, K.R. Timmermans, S. Blain, & P. Treguer
  !                JSR 53 (2005) 25-42
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Copy  state var. object in local var
  !  A small concentration is to avoid that one of the state vars. --->0.0
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  phytoc = D3STATE(:,ppphytoc)
  dry_pel=DONE;if (p_dry_ben) dry_pel=dry_z(1)

  avPhyto=sum(phytoc*Depth)/(NZERO+sum(Depth))
  if (p_sum(phyto)==ZERO .or. avPhyto < low_phytoc) then
    rx_any=ZERO
    call sourcesink_flux_vector( iiPel, ppO3c,ppphytoc, rx_any )  !reset on zero
    call sourcesink_flux_vector( iiPel, ppphytoc,ppO3c, rx_any )  !reset on zero
    call flux_vector( iiPel, ppphytoc,ppR1c, rx_any )  !reset on zero
    call flux_vector( iiPel, ppphytoc,ppR2c, rx_any )  !reset on zero
    call flux_vector( iiPel, ppphyton,ppR2n, rx_any )  !reset on zero
    call flux_vector( iiPel, ppphytoc,ppR3c, rx_any )  !reset on zero
    call flux_vector( iiPel, ppN3n,ppphyton, rx_any )  !reset on zero
    call flux_vector( iiPel, ppN4n,ppphyton, rx_any )  !reset on zero
    call flux_vector( iiPel, ppR1n,ppphyton, rx_any )  !reset on zero
    call flux_vector(iiPel, ppphyton,ppN4n, rx_any)  !reset on zero
    call flux_vector( iiPel, ppN1p,ppphytop, rx_any )  !reset on zero
    call flux_vector( iiPel, ppR1p,ppphytop, rx_any )  !reset on zero
    call flux_vector(iiPel, ppphytop,ppN1p, rx_any )  !reset on zero
    if ( ppphytos>0) then
      call flux_vector( iiPel, ppN5s,ppphytos, rx_any)  !reset on zero
      call flux_vector( iiPel, ppphytos,ppN5s, rx_any)  !reset on zero
      call flux_vector( iiPel, ppphytoc,ppR6c, rx_any )  !reset on zero
    endif
    return
  endif

  phyton = D3STATE(:,ppphyton)
  phytop = D3STATE(:,ppphytop)
  phytol = D3STATE(:,ppphytol)
  if ( ppphytos > 0 )  phytos = D3STATE(:,ppphytos)

  silica_control=0
  if ( p_qus(phyto) > ZERO )  then
    silica_control=2
    if ( p_xqs(phyto)>ZERO ) silica_control=3
  elseif ( p_chPs(phyto) > ZERO ) then
    silica_control=1
  endif
  fl_xgrazing_PIc(phyto,:)=Source_D3_withgroup(ppphytoc, &
                ppMesoZooPlankton,iiMesoZooPlankton,iiC,iiConsumption) + &
               Source_D3_withgroup(ppphytoc, &
                ppMicroZooPlankton,iiMicroZooPlankton,iiC,iiConsumption) &
                -getbotflux_3d(ppphytoc,Depth)
   sx_grazing=fl_xgrazing_PIc(phyto,:)/(NZERO+phytoc)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient limitation (intracellular) N, P
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  iN1p = max( NZERO, ( qpPc(phyto,:)&
                             - p_qplc(phyto))/( p_qpRc(phyto)- p_qplc(phyto)))
  iNIn = max( NZERO, ( qnPc(phyto,:)&
                            - p_qnlc(phyto))/( p_qnRc(phyto)- p_qnlc(phyto)))
  iNIs=DONE;
  if ( silica_control >= 2 ) iNIs = max( NZERO, ( qsPc(phyto,:)&
                            - p_qslc(phyto))/( p_qsRc(phyto)- p_qslc(phyto)))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Phytoplankton growth is limited by nitrogen and phosphorus
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    select case ( p_limnut)
      case ( 0 ) ; iNPI(phyto,:)  =   (iN1p* iNIn)**(0.5D+00)
      case ( 1 ) ; iNPI(phyto,:)  =   min(  iN1p,  iNIn)  ! Liebig rule
      case ( 2 ) ; iNPI(phyto,:)  =   2.0D+00/( DONE/ iN1p+ DONE/ iNIn)
    end select

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient limitation due to intra- extracellular silicate.
  ! eN5s limit externally nutrient limitation.
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !calculate internal quota
  eN5s  =   DONE
  if ( silica_control > 0 ) then
    select case (silica_control)
      case(1)
         eN5s = min( eN5s,N5s/(N5s + p_chPs(phyto)))
         tN=min(DONE,iNPI(phyto,:),eN5s)
      case(2,3)
        tN=min(DONE,iNIs,iNPI(phyto,:))
      end select
  else
      tN=iNPI(phyto,:)
  endif

  ex_limit=2.0* p_thdo(phyto)* (DONE/( tN+ p_thdo(phyto)))
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Temperature response of Phytoplankton
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   et  =   eTq_vector(  ETW,  p_q10(phyto))


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Photosynthesis (Irradiance EIR is in uE m-2 s-1, Irr is mid-layer EIR)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  eiPI(phyto,:)=DONE
  if ( ChlLightFlag== 2) then
    Irr = ZERO !in case of a tidal flat production when dry_pel=0
    px_any= Depth*dry_pel*xEPS
    where (px_any.gt.ZERO) Irr = max( NZERO, EIR*( DONE- exp(-px_any))/ px_any)
    eiPI(phyto,:) = max(ZERO,  &
        DONE- exp_limit(- qlPc(phyto, :)/ p_qchlc(phyto)/ p_Ke(phyto)* Irr))
  end if

  select case ( LightForcingFlag)
    case ( 1 ) ;sumc=p_sum(phyto)* et* eiPI(phyto,:)
    case ( 2 ) ;sumc=p_sum(phyto)* et* eiPI(phyto,:)*( SUNQ/HOURS_PER_DAY)
    case ( 3 ) ;sumc=p_sum(phyto)* et* eiPI(phyto,:)* ThereIsLight
  end select


#ifdef INCLUDE_PELCO2
  rx_any = sumc*phytoc/MW_C
  !2ways to limit Phyto
  call LimitChange_vector(POSITIVE,rx_any,HCO3,max_change_per_step,px_any)
  where (HCO3.gt.ZERO)! (HCO3==0 PH calculation is not succeeded).
    sumc=sumc*min(px_any,HCO3/(5.0D+00+HCO3))
  end where
#endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Activity processes
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  sea  =   sumc* p_pu_ea(phyto)
  sra  =   p_pu_ra(phyto)*max(ZERO, sumc- sea)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Rest processes
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ses  =   p_ses(phyto)*et            ! activity +rest excretion
  srs  =   et* p_srs(phyto)

  jrrPTc(1)=jrrPTc(1)+sum(srs*phytoc*Depth)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Potential-Net prim prod. (mgC /m3/d):
  !  a. limit prim prod at high temperature with the (maximum) number of
  !     divisions per day by make use of the length of the light period.
  !  c. assume that activity respiration increases with temperature and
  !     is NOT limited by the number of divisions per day.
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  dl=(SUNQ/ HOURS_PER_DAY)
  sx_any  =max( ZERO, dl*( sumc- sea -sra)-ses -srs)  ! net production
  sudc=p_xdiv(phyto)*log(2.0D+00)
  px_any=DONE
  where (sx_any.gt.ZERO) px_any=min(DONE,sudc/(NZERO+sx_any))
  sumc=px_any*sumc; sea=px_any*sea
  sra=sra*px_any**p_xgeneric(phyto)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Lysis and excretion
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  !density dependent mortality:only used for phytoplankton which is not beeaten
  sdo  =   p_sd2(phyto)* phytoc *phytoc
  !mortality due too low chlorophyl content
  sdo=sdo +p_sdmo(phyto)* exp(-(qlPc(phyto,:) +p_qlPlc(phyto))/p_qlPlc(phyto))
  !low oxygen  dependent mortality

   rr3c=ZERO
   cx_any=DONE
   if ( phyto .eq.iiP6) then
     !Calculate mortality add this to sdo
     ! !output in sx_any. no input: px_any is a dummy var.
     px_any=ZERO
     call PhaeocystisCalc(CALC_MORTALITY,phyto,sx_any,px_any,p_sdmo(phyto))
     call LimitChange_vector(POSITIVE,sx_any,cx_any,max_change_per_step)
     ! remove colonies by decreasing the biomass Pcc
     ! !input in sx_any. no output: px_any is a dummy var.
     call PhaeocystisCalc(COLONY_DEGRADATION,phyto,px_any,sx_any,ZERO)
     ! from R3->R2, P6->R1n P6c->R1p, P56->R1c are calculated in this call
     call PhaeocystisCalc(CALC_LOC_DET_FLUX,phyto,px_any,sx_any,ZERO)
     sdoPI(phyto,:)=sx_any

     call PhaeocystisCalc(CALC_MORTALITY_CELLS_IN_COLONY,phyto, &
                                            sx_any,et,1.0D+00)
     call LimitChange_vector(POSITIVE,sx_any,cx_any,max_change_per_step)
     sdo=sdo+sdoPI(phyto,:)+sx_any
   else
! if p_iRI(phyto)==1 we deal with "normal phytoplankton"
! which under go at nutrient stress autolyis which may be caused by viral lysis.
     sdo  = sdo + ex_limit * p_sdmo(phyto)  &
        * phytoc(:)/(phytoc(:)+0.01D+00)  ! nutr. -stress lysis
     call LimitChange_vector(POSITIVE,sdo,cx_any,max_change_per_step)
   endif
   rml=rml+sdo *phytoc

  !-=-=0-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Apportioning over R1 and R6:
  ! Cell lysis generates both DOM and POM.
  ! The nutr.-depleted cell has a nutrient-carbon ratio equal to p_q?lc.
  ! Assuming that this structural part is not easily degradable,
  ! at least a fraction equal to the minimum quota is released as POM.
  ! Therefore, nutrients (and C) in the structural part go to R6.
  ! Mortality in P6 is defined as colony breakdown
  ! in case of P6 sdo stands for mortality rate  of cells in colony:
  ! this leads only to flux of P6 to R3c
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if (phyto.ne.iiP6) sdoPI(phyto,:)=sdo
  rx_any=sdoPI(phyto,:)*phytoc
  rr1c  =  p_pe_R1c* rx_any
  rr6c  =  (DONE-p_pe_R1c)* rx_any

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! total processes
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  srt  =   sra+ srs     ! total specific respiration
  set  =   sea+ ses     ! total specific excretion
  slc  =   sea+ sra+ sdoPI(phyto,:)  ! specific loss terms

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Correct for mortality : in order to model an effective
  ! mortality at the end of a bloom it is necessary to make all phyto (Phaeo)
  ! which will die during next time-step inactive, otherwise the mortality
  ! is counteracted by the new production
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  px_any=DONE-sdo*LocalDelta
  sumc=sumc*px_any

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Flows of realized gross Production and respiration
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rugc=sumc * phytoc  ! gross production
  call sourcesink_flux_vector( iiPel, ppO3c,ppphytoc, rugc )  ! source/sink.c
  call FindNaNInRates(iiPel,ppO3c,'Phyto:after rugc')
  call findnan(rugc,NO_BOXES,iout)
  if ( iout.gt.0) then
    write(LOGUNIT,*)'phyto_nr,phytoc,phyton,phytop,sdo,iNIn,iN1p,qlPc,tN,xEps',&
         phyto,phytoc(iout), phyton(iout),phytop(iout), &
         sdo(iout),iNIn(iout),iN1p(iout),qlPc(phyto,iout),tN(iout),xEps(iout) 
  endif
  call flux_vector( iiPel, ppO2o,ppO2o, rugc/ MW_C ) 
  rrtc  =   srt* phytoc  ! total actual respiration
  rrmc  =   srs* phytoc  ! total rest respiration
  jPLO3c(1)=jPLO3c(1)+sum(srs*phytoc*Depth)


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of Potential-Net prim prod. (mgC /m3/d)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  !Only anoxic mortality of respiration is larger than net production
  ! eo=1 : no mortality  eo<1 : mortality
  px_any=max(ZERO,-(sumc-sea -ses-srs)/(ses+srs))
  eo  =(DONE-px_any)+px_any*( DONE-exp(-(O2o+p_qon_nitri*N3n)/p_clO2o(phyto)))

  rufc  =   max(  NZERO, ( sumc- slc)* phytoc)  ! net production
  sugPI(phyto,:)  =   rufc/( NZERO+ phytoc)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient Uptake: calculate maximal uptake of N,P
  ! Check if C-fixation is larger to make of all C new biomass
  ! If not increase excretion of C
  ! NZERO is added in case that parameter p_... and the nutrient are both 0.
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rscalar=NZERO+p_lN3N4n(phyto)
  cqun3 = rscalar/( rscalar+N4n)
  rscalar=NZERO+p_lureaN4n(phyto)
  cquR1n = rscalar/( rscalar+N4n)
  rscalar=NZERO+p_lN1(phyto)
  cqup  = rscalar/( rscalar+ N1p)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! The nutrient uptake is controlled by the fr_lim_PI. This variable stands
  ! for the fraction of nutrients which is avilable for  uptake by this
  ! functional group. The factor is deterimined by comparing the potential
  ! net growth rates (in mgC/m3/d) of all functional groups which take up
  ! nutrients  (See PelB/LimitNutrientUptake.F90).In this way nutrient uptake
  ! is coupled to the net-growth, however without excluding the effect of
  ! nutrient affinity defined by the p_qu[nps] parameters.
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  px_any=fr_lim_PI_n(phyto,:)
  !----------------potential N-uptake-------------------------------
  lim_qu= DONE
  ! max pot. uptake of N4
  call PhaeocystisCalc(CALC_REL_AMMONIUM_UPTAKE,phyto,lim_qu,N4n, &
                                                             p_qun(phyto))
  rum4n=   max(ZERO,p_qun(phyto)* N4n* phytoc*et* lim_qu)
  call DoubleLimitChange_vector(POSITIVE,rum4n,N4n,px_any,max_change_per_step)
  !max pot.uptake of N3
  call PhaeocystisCalc(CALC_REL_NITRATE_UPTAKE,phyto,lim_qu,N3n, &
                                                             p_qun(phyto))
  rum3n=max(ZERO,p_qun(phyto)* N3n* phytoc*et *lim_qu*cqun3 )
  call DoubleLimitChange_vector(POSITIVE,rum3n,N3n,px_any,max_change_per_step)
  !max pot. uptake of R1n
  call PhaeocystisCalc(CALC_REL_UREA_UPTAKE,phyto,lim_qu,Nun, &
                                                           p_qun(phyto))
  rumun=max(ZERO,p_qun(phyto)*Nun* phytoc*et *lim_qu*cquR1n )
  call DoubleLimitChange_vector(POSITIVE,rumun,Nun,px_any,max_change_per_step)
  !----------------potential P-uptake-------------------------------
  px_any=fr_lim_PI_p(phyto,:)
  lim_qu=DONE
  call PhaeocystisCalc(CALC_REL_PHOSPHATE_UPTAKE,phyto,lim_qu,N1p, &
                                                               p_qup(phyto))
  rum1p=  max(ZERO,p_qup(phyto)* N1p* phytoc*et *lim_qu )
  call DoubleLimitChange_vector(POSITIVE,rum1p,N1p,px_any,max_change_per_step)
!---------------------------------------------------------------------------Y
      call findnan(rum1p,NO_BOXES,iout)
      if ( iout.gt.0) then
        write(LOGUNIT,*),'After LimitChange Nan in rum1p',rum1p(iout), &
        lim_qu(iout),phytoc(iout),px_any(iout),fr_lim_PI_n(phyto,iout)
      endif
!---------------------------------------------------------------------------Y

  ! max pot. uptake of org. P
  ! output=lim_qu ( affinity depends on size of colonies)
  lim_qu=DONE
  rumup=  max(ZERO,p_qup(phyto)*Nup* phytoc*et*lim_qu*cqup)
  call DoubleLimitChange_vector(POSITIVE,rumup,Nup,px_any,max_change_per_step)
  !----------------potential Si-uptake-------------------------------
  if ( silica_control > 0 )  then
    rum5s  =  max(ZERO,p_qus(phyto)*N5s* phytoc*et)
    call LimitChange_vector(POSITIVE,rum5s,N5s,max_change_per_step)
  endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !  Determintation of optimum quotum..
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
     qnoc=p_xqn(phyto)*p_qnRc(phyto)
     qpoc=p_xqp(phyto)*p_qpRc(phyto)
     if ( silica_control>= 2) qsoc=p_xqs(phyto)*p_qsRc(phyto)
!    px_any=min(qnPc(phyto,:)/p_qnRc(phyto),qpPc(phyto,:)/p_qpRc(phyto))
!    if ( silica_control>= 2) px_any=min(px_any,qsPc(phyto,:)/p_qsRc(phyto))
!    qnoc=p_qnRc(phyto)*min(max(p_qnlc(phyto)/p_qnRc(phyto), &
!                                       (px_any+p_xqn(phyto))*0.5),p_xqn(phyto))
!    qpoc=p_qpRc(phyto)*min(max(p_qplc(phyto)/p_qpRc(phyto), &
!                                       (px_any+p_xqp(phyto))*0.5),p_xqp(phyto))
!    if ( silica_control>= 2) qsoc=p_qsRc(phyto)*min( &
!       max(p_qslc(phyto)/p_qsRc(phyto),(px_any+p_xqs(phyto))*0.5),p_xqs(phyto))


  !in cx_any N present in cell buffer
  !For Phaeocystis it means that N which is outside the cell but in the colony
  ! is not seen as part of the buffer.
  cx_any= phyton- phytoc*p_lqnlc(phyto)
! rbufn=max(ZERO, rufc)* p_xqn(phyto)* p_qnRc(phyto)
  rbufn=max(ZERO, rufc)* qnoc
  call LimitChange_vector(POSITIVE,rbufn,cx_any,max_rate_per_step)

  ! N uptake brom 4 sources: nitrate, ammonium, organic N, from buffer incell.
  ex_sw=insw_vector(cx_any)
  rumn  =   rum3n+ rum4n  +rumun +rbufn*ex_sw

!---------------------------------------------------------------------------Y
  call findnan(rumn,NO_BOXES,iout)
  if (iout.gt.0) then
    write(LOGUNIT,*)  &
     'Phyto:rumn NaN',phytoc(iout),rum3n(iout),rum4n(iout),rumun(iout)
    write(LOGUNIT,*)'Phyto:nr,N3n,N4n,Nun:',N3n(iout),N4n(iout),Nun(iout)
    write(LOGUNIT,*)'Phyto:nr,rbufn,qnoc,ex_sw:',rbufn(iout),qnoc(iout),ex_sw(iout)
  endif
!---------------------------------------------------------------------------Y

  ! P uptake brom 3 sources: phosphate, organic P, from buffer in cell.
  !For Phaeocystis it means that P which is outside the cell but in the colony
  ! is not seen as part of the buffer.
  cx_any = phytoc*qpPc(phyto,:)- phytoc*p_qplc(phyto)
  rbufp = max(ZERO,rufc)* qpoc
  call LimitChange_vector(POSITIVE,rbufp,cx_any, max_rate_per_step)
  ex_sw=insw_vector(cx_any)
  rump  = rum1p +rumup+rbufp*ex_sw  ! max pot. uptake of NI

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Check if which fraction C-fixation can be used for new biomass
  ! by checking the potential nutrient avilability
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  runc = min( rufc, rumn/p_qnlc(phyto),  &
                                 rump/p_qplc(phyto))

! ! Si uptake brom 2 sources: dissolved silicate, from diatom-skeleton
  if ( silica_control  >=  2) then
    cx_any = phytos- phytoc*p_qslc(phyto)
!   rbufs =max(ZERO, rufc)*  p_xqs(phyto)* p_qsRc(phyto)
    rbufs =max(ZERO, rufc)* qsoc
    call LimitChange_vector(POSITIVE,rbufs,cx_any,max_rate_per_step)
    ex_sw=insw_vector(cx_any)
    rums = rum5s+rbufs*ex_sw
    runc = min( runc, rums/p_qslc(phyto))
  endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! C flows ; phyto --> detritus
  ! a. activity excretion of sugars:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  flChydrate  =   set* phytoc

 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 ! b. stress exretion:all C which can not be used for growth as carbo-hydrates:
 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  runc  =   min(rufc,max( runc,ZERO))


  leftover = rufc-runc
  !distribute too much fixed carbon over stress respiration and excretion
  flChydrate  =   flChydrate+ leftover


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! C flows: calculate fluxes
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rr2c=ZERO
  select case (p_iRI(phyto))
     case (iiR1) ; rr1c =flChydrate + rr1c
     case (iiR2) ; rr2c= flChydrate
     case (iiR3)
          rr3c=rr3c+flChydrate
          if ( phyto.ne.iiP6) stop 'error in definition of rr3c'
  end select

!---------------------------------------------------------------------
  call findnan(rr2c,NO_BOXES,iout)
  if ( iout.gt.0) then
         write(LOGUNIT,*),'Phyto:Nan in phytoc,rumup' &
         ,EIR(iout),xEPS(iout),runc(iout),rufc(iout)
  endif
!---------------------------------------------------------------------

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Apparent Net prim prod. (mgC /m3/d)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  sunPI(phyto,:)  =   runc/( NZERO+ phytoc)
  sx_main=srs+ses
  sx_act=max(sx_main,sumc)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient dynamics: NITROGEN
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! Intracellular missing amount of N, ruXCol is only>0 for P6(Phaeocystis)
  call PhaeocystisCalc(CALC_NET_NITROGEN_UPTAKE,phyto,ru_xCol_n,flChydrate,ZERO)
  !Nutrient uptake always >=0, corrrection of content via misn
  qx_any=qnoc
  !In rupn an estimation is made about the size of the net growth on next
  !timestep. To avoid large changes during the night we assume only net growth
  !during the day and therefor the loss terms are correctd for daylength,
  !implictily assuming that grazing during the night is the same during the day.
  rupn= &
      max(ZERO,ru_xCol_n+qnoc*(sunPI(phyto,:)-(sx_grazing+sx_main)/dl)*phytoc)
  qx_any=qnoc- qnPc(phyto,:)
  ex_sw=insw_vector(qx_any)
  misn=(sx_act*qx_any-qnPc(phyto,:)*sx_main*(DONE-ex_sw))*phytoc
  luxn=sx_act*(p_xqn(phyto)*p_qnRc(phyto)-qnoc)*phytoc
  runn  =   min( rumn,rupn+misn )  ! actual uptake of NI

  !Calculate the uptake of the different sources.
  ! (uptake from the cell itself changes automatically internal buffer)
  ex_sw  =   insw_vector(runn)
  px_any=min(DONE,runn*ex_sw/( NZERO+ rumn))  
  run3n  =   px_any* rum3n
  run4n  =   px_any* rum4n
  runun  =   px_any* rumun

  call FindNaNInRates(iiPel,ppN4n,'Before phyto->N4n')
  call flux_vector( iiPel, ppN3n,ppphyton, run3n )  ! source/sink.n
  call flux_vector( iiPel, ppN4n,ppphyton, run4n )  ! source/sink.n
  call flux_vector( iiPel, ppR1n,ppphyton, runun )  ! source/sink.n
  call flux_vector( iiPel, ppR1c,ppphytoc, runun/p_qnUc )  ! source/sink.n
  ! oxygen which freed by phytoplankton when they take up nitrate
  call flux_vector( iiPel, ppO2o,ppO2o,run3n*p_qon_nitri )  ! source/sink.o
  !correction in winter situation of negative growth to keep nutrients below
  ! max.ratio
  runn=min( rumn,rupn+misn+luxn )  ! actual uptake of NI
  rx_any=-runn*(DONE-insw_vector(runn))
! reun=ZERO; reUc=ZERO;
  reUn=rx_any;reUc=reUn/p_qnUc
! call flux_vector(iiPel, ppphyton,ppN4n,rx_any)
  call FindNaNInRates(iiPel,ppN4n,'After phyto->N4n')
  !urea dynamics assuming that breakdown happen at outside phytoplankton cel
  !Collect rates: (See PelGlobal and PelChem)
  flPIN4n=flPIn4n- runn*( DONE-ex_sw)  !only for output  !only for output

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nuttrient dynamics: PHOSPHORUS
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Intracellular missing amount of P, ruXCol is only>0 for P6(Phaeocystis)
  call PhaeocystisCalc(CALC_NET_PHOSPHATE_UPTAKE,phyto,ru_xCol_p, &
                                                           flChydrate,ZERO)

  rupp= &
     max(ZERO,ru_xCol_p+qpoc *(sunPI(phyto,:)-(sx_grazing+sx_main)/dl)*phytoc)
  !Nutrient uptake always >=0, corrrection of content via rn_xlux_p
  qx_any=qpoc - qpPc(phyto,:)
  ex_sw=insw_vector(qx_any)
  misp=(sx_any*qx_any-qpPc(phyto,:) *sx_main*(DONE-ex_sw))*phytoc
  luxp=sx_any*(p_xqp(phyto)*p_qpRc(phyto)-qpoc)*phytoc

  runp  =   min( rump,rupp+ misp)  ! actual uptake

  !Calculate the uptake of the different sources.
  ! (uptake form the cell itself changes automatically internal buffer)
  ex_sw  =   insw_vector(runp)
  px_any=min(DONE,runp*ex_sw/( NZERO+ rump))  ! actual uptake of Nn
  run1p  =   px_any* rum1p
  runup  =   px_any* rumup  ! actual uptake of Nn

  call flux_vector( iiPel, ppN1p,ppphytop, run1p )  ! source/sink.n
  call flux_vector( iiPel, ppR1p,ppphytop, runup )  ! source/sink.n

  !correction in winter for negative growth to keep nutrients below max.ratio
  runp  =   min( rump,rupp+ misp+luxp)  ! actual uptake
  ex_sw=insw_vector(runp)
  call flux_vector(iiPel,ppphytop,ppN1p,-runp*(DONE-ex_sw)) !source/sink.

  call FindNaNInRates(iiPel,ppN1p,'After ppN1p-fluxes')

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Excretion of N and P to PON and POP
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rx_any=sdoPI(phyto,:)*qnPc(phyto,:)*phytoc
  s=phyton-p_qnlc(phyto)*low_phytoc
  ! extreme cases no flux at lower below :
  call LimitChange_vector(POSITIVE,rx_any,s,max_change_per_step)
  rr1n  =   p_pe_R1n* rx_any
  rr6n  =   rx_any- rr1n

  rx_any=sdoPI(phyto,:)*qpPc(phyto,:)*phytoc
  s=phytop-p_qplc(phyto)*low_phytoc
  ! extreme cases no flux at lower below :
  call LimitChange_vector(POSITIVE,rx_any,s,max_change_per_step)
  rr1p  =   p_pe_R1p* rx_any
  rr6p  =   rx_any- rr1p


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient dynamics: SILICATE
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  if ( silica_control > 0 )  then
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Losses of Si
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    rr6s  =   sdoPI(phyto,:)* phytos  ! Lysis, particulate

    ! Collect first all fluxes of P-->silica
    flPIR6s(phyto,:)  =   flPIR6s(phyto,:)+ rr6s
    select case (silica_control)

    case (1)
     runs = max(ZERO, p_qsRc(phyto) * runc );          ! net uptake
     call flux_vector( iiPel, ppN5s,ppphytos, runs)  ! source/sink.c

    case (2,3)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      !  Nutrient uptake
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Si uptake based on C uptake, in qx_any maximum quotum
      rups= &
          max(ZERO,qsoc*(sunPI(phyto,:)-(sx_grazing+sx_main)/dl))*phytoc
      ! miss is used to correct for too much silica in the cells,
      ! in q_ra diff between max and acutual
      qx_any=qsoc- qsPc(phyto,:)
      ex_sw=insw_vector(qx_any)
      miss= (sx_act*qx_any-qsPc(phyto,:)*sx_main*(DONE-ex_sw))*phytoc
      luxs=sx_act*(p_xqs(phyto)*p_qsRc(phyto)-qsoc)*phytoc
      runs  =  min( rums,rups+miss ) ! actual uptake
      ex_sw=insw_vector(runs)
      !Calculate the uptake of the different sources.
      ! (uptake form the cell itself changes automatically internal buffer)
      ! source/sink.s:
      run5s  =   ex_sw* runs* rum5s/( NZERO+ rums)  ! actual uptake of Nn
      call flux_vector( iiPel, ppN5s,ppphytos, run5s)
      runs  =  min( rums,rups+miss+luxs ) ! actual uptake
      ex_sw=insw_vector(runs)
      call flux_vector( iiPel, ppphytos,ppN5s,  &
                                      - runs*(DONE-ex_sw))  ! source/sink.s
    end select
  endif

  if ( ChlLightFlag== 2) then
    x_light_limitation=min(DONE,p_EIR(phyto)/(NZERO+EIR))
    qx_any=p_lqnlc(phyto)*x_light_limitation
    ! due to (too much) light and N-limitation
    ! incase of extreme limitation:
    ! DONE/NZERO:keep x_stress_degradation below the maximum real number...
    px_any=max(ZERO,DONE- &
          (qx_any+qnPc(phyto,:)-p_qnlc(phyto))/ (p_lqnlc(phyto)-p_qnlc(phyto)))
    sx_stress_degradation=p_sdchl_l(phyto)*px_any
    cx_any=DONE
    call LimitChange_vector(POSITIVE,sx_stress_degradation, &
                                              cx_any,max_change_per_step)
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Chl-a synthesis and photoacclimation
    ! Carbon fixation is limited at high temperatures leading to a lower
    ! growth rate.
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    where (Irr> ZERO)
      rho_Chl = p_qchlc(phyto)* p_qchlc(phyto)* eiPI(phyto,:) &
          /qlPc(phyto,:)* p_Ke(phyto)/ Irr* x_light_limitation
    elsewhere
       rho_Chl=ZERO
    endwhere

! total synthesis, only when there is net production (runc > 0)
    new_Chl = rho_Chl*(( sumc- sra)* phytoc-flChydrate)
    !determine net growth without loss
    sx_any=sumc-slc+sdo-ses-srs
    !rate=new -toomuch (especially in winter)  -degradation

    rate_Chl=new_Chl- sdo* phytol-et*p_sdchl_d(phyto)* phytol
!     + min(ZERO,sx_any )*  max( ZERO, phytol- p_qchlc(phyto)* phytoc)

    ! calculate total rate of change of Chla in phyto
    rate_Chl=rate_Chl-phytol * sx_stress_degradation
    call flux_vector( iiPel, ppphytol,ppphytol, rate_Chl )

    ! calulate loss of N due to stress degradation due oxygen radicals mgC/m2/d)
    ! it only happen in the light and when N become limiting.
    rx_any= sx_stress_degradation*max(ZERO,phyton-phytoc*p_qnlc(phyto))
!--------------------------------------------------------------------------
    call findnan(rx_any,NO_BOXES,iout)
    if (iout.gt.0) then
      write(LOGUNIT,*) 'Chl degradation rate is NaN phyto,qnPc,phytoc,sx_any', &
                        phyto,qnPc(phyto,iout),phytoc(iout),sx_any(iout)
    endif
!--------------------------------------------------------------------------
!   where( rx_any>ZERO)
!   ! distribute degradation produces over R1n and R6n
!   rr6n=rr6n+(DONE-p_pe_R1n) * rx_any
!   rr1n=rr1n+p_pe_R1n * rx_any
!   ! calulate loss of C due to stress degradation due oxygen radicals mgC/m2/d)
!   rx_any =rx_any/min(p_xqn(phyto)*p_qnRc(phyto), &
!                              max(p_qnlc(phyto),qnPc(phyto,:)))
!   ! distribute degradation produces over R1c and R6c
!   rr6c=rr6c+(DONE-p_pe_R1c) * rx_any
!   rr1c=rr1c+p_pe_R1c * rx_any
!   ! calulate loss of C due to stress degradation due oxygen radicals mgC/m2/d)
!   rx_any=rx_any*qpPc(phyto,:)
!   ! distribute degradation produces over R1c and R6c
!   rr6p=rr6p+(DONE-p_pe_R1p) * rx_any
!   rr1p=rr1p+p_pe_R1p * rx_any
!   endwhere
  endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Respiration flows + (only for CO2) release of CO2 due to Urea->NH4
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !reverse nitrification
  ! maximal part of respiration where probably is not enough O2 present
  rrx_an_c=rrmc*(DONE-eo)
  ! part of the respiration where reverse nitrification is used...
  call flux_vector( iiPel, ppO2o,ppO2o,-( rrtc-rrx_an_c)/ MW_C ) ! source/sink.o
  px_any=p_qon_nitri*N3n/(NZERO+O2o+p_qon_nitri*N3n)
  rx_any=rrx_an_c*px_any/MW_C
  !subtract
  !add to flux reverse nitrification
  call LimitChange_vector(POSITIVE,rx_any,N3n*p_qon_nitri, &
                                                  max_change_per_step,s)
  flN3N4n=flN3N4n+rx_any*px_any*s
  !unit mol02/m3
  rx_any=rrx_an_c*(DONE-px_any*s/p_xeff_an)/MW_C
  !unit-> mol S2-
  call flux_vector( iiPel, ppN6r,ppN6r,+rx_any/p_xeff_an*p_qro )
  rx_any=rrx_an_c*(px_any*s+(DONE-px_any*s))/p_xeff_an
  call sourcesink_flux_vector( iiPel, ppphytoc,ppO3c, rrtc-rrx_an_c+rx_any)
  call findnega(rrtc,NO_BOXES,iout)
  if ( iout.gt.0) then
     write(LOGUNIT,*) 'rrtc,srs,phytc=',rrtc(iout),srs(iout),phytoc(iout)
  endif
  call FindNaNInRates(iiPel,ppphytoc,'After phyto->ppO3c')
  call FindNaNInRates(iiPel,ppO3c,'Phyto:after rrtc')

  select case (phyto.eq.iiP6)
    case (.true.)
   ! for phaeo the colony mortality is included in sdoPI(phyto,:). The total
   ! mortality in sdo,the difference is the mortality of cells in the colony.
      px_any=sdoPI(phyto,:)/(NZERO+sdo)
      !rr1 and rr6c are produced by mortality of the cell.
      rx_any=(rr1c+rr6c)*(DONE-px_any)+rr3c
      call flux_vector( iiPel, ppphytoc,ppR3c, rx_any )  ! source/sink.c
    case (.false.)
      px_any=DONE
  end select

  call flux_vector( iiPel, ppphytoc,ppR1c, rr1c*px_any +reUc )  ! source/sink.c
  call flux_vector( iiPel, ppphytoc,ppR6c, rr6c*px_any )  ! source/sink.c
  call flux_vector( iiPel, ppphytoc,ppR2c, rr2c)  ! source/sink.c
  rx_any=rr2c* &
    max(ZERO,min(DONE,qnPc(phyto,:)-p_qnlc(phyto))/(p_qnRc(phyto)-p_qnlc(phyto))*p_qnR2c)
  call flux_vector( iiPel, ppphyton,ppR2n, rx_any)  ! source/sink.c
  call FindNaNInRates(iiPel,ppphytoc,'After phyto->R6c')

  flPIR6n(phyto,:)=flPIR6n(phyto,:) + rr6n*px_any !fluxes define in PelChem
  flPIR6p(phyto,:)=flPIR6p(phyto,:) + rr6p*px_any
  flPIR1n(phyto,:)=flPIR1n(phyto,:) + rr1n*px_any+reUn
  flPIR1p(phyto,:)=flPIR1p(phyto,:) + rr1p*px_any

  rx_any=(sumc-sra)*phytoc-(flChydrate-ses*phytoc)
  rnetPTc=rnetPTc+rx_any
  jnetPIc(phyto,1)=jnetPIc(phyto,1)+sum(rx_any*Depth)
! End of computation section for process PhytoDynamics


  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
