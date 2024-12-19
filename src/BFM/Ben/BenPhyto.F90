#include "DEBUG.h"
#include "INCLUDE.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenPhyto
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
  subroutine BenPhytoDynamics(bp, ppbpc, ppbpn, ppbpp, ppbps, &
    ppbpl)
!
! !USES:

  ! For the following Pelagic-states fluxes are defined: Q6c, G2o, Q2c, &
  ! K3n, R1n, Q6n, R1p, Q6p
  ! The following global scalar vars are used: SUNQ, ThereIsLight
  ! The following Pelagic 1-d global boxvars are modified : jBIQ6s
  ! The following Pelagic 1-d global boxvars  are used: ETW_Ben, EIR_Ben
  ! The following Pelagic 2-d global boxvars got a value: sunBI
  ! LightForcingFlag
  ! The following global constants are used: RLEN
  ! The following constants are used: HOURS_PER_DAY

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO,NZERO,DONE
#ifdef NOPOINTERS
  use mem,  ONLY: D2STATE
#else
  use mem, ONLY: D2STATE,Dfm,Dcm,Dlm,D1m,Q6c,Q6n,Q6p,Q6s,&
                         D6m,D7m,D8m,D9m,Q6c,Q1p,Q1c,G2o,G3c
#endif
  use mem, ONLY:NO_BOXES_XY,ppG2o, ppQ1c,ppQ2n,ppQ2c,ppQ6c,ppQ6n,ppQ6p,&
    ppK3n,ppK4n,ppG3c, ppK1p, ppQ1n,ppQun, ppQ1p, ppK5s,ppQ6s,ppK16r, &
    ppDfm,ppDcm,ppDlm,ppD6m,ppD7m,ppD8m,ppD9m, iiR1,iiR2,iiBen,iiC,iiL, &
    jnetBPTc, Source_D2_vector,iiConsumption,iiPhytoPlankton, &
    sunBI, sugBI, flux_vector,sourcesink_flux_vector,shiftDlm, &
    ruBPc,ruBPn,ruBPp,ruBPs,fr_lim_BPI_n,fr_lim_BPI_p,xEPS_Sedi,&
    jPIQ6s,LocalDelta,max_change_per_step,max_rate_per_step, &
    ETW_Ben,EIR_Ben,SUNQ,ThereIsLight,CoupledtoBDc,turenh,pxturinD1, &
    Kp1p,K1p,Kp5s,K5s,Kp3n,K3n,Kp4n,K4n,Qpun,Qun,ResetSource_D2_Vector, &
    ppBenOrganisms,iiBenOrganisms,iiConsumption
#ifdef INCLUDE_BENCO2
  use mem,  ONLY: HCO3ae
#endif
  use constants,  ONLY: HOURS_PER_DAY,INTEGRAL,AVERAGE, &
         MW_C,p_qnUc,POSITIVE,NEGATIVE,ANY
  use mem_Param,  ONLY: p_pe_R1c, p_pe_R1n, p_pe_R1p, &
                 p_pK1_ae,p_pK5_ae,LightForcingFlag, p_poro,p_epsChla,&
                 p_d_tot,p_epsR6,p_epsESS,p_qpPhc,p_qon_nitri,p_qro,p_xeff_an

  use mem_BenPhyto
  use mem_Phyto
  use mem_Bioturbation,ONLY:p_Etur,p_cturm
  use mem_BenAmmonium, ONLY: p_pK4=>p_p
  use mem_BenNitrate,  ONLY: p_pK3=>p_p
  use SourceFunctions,only: Source_D2_withgroup

  use BFM_ERROR_MSG, ONLY: set_warning_for_getm,BFM_ERROR

  use LimitRates, ONLY:LimitChange_vector,DoubleLimitChange_vector

! use mem,  ONLY: Output2d_1,Output2d_2,Output2d_3, Output2d_4

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following vector functions are used:eTq_vector, insw_vector
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: eTq_vector, insw_vector,exp_limit,PartQ_vector

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

! !INPUT:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer,intent(IN)  :: bp
  integer,intent(IN) :: ppbpc
  integer,intent(IN) :: ppbpn
  integer,intent(IN) :: ppbpp
  integer,intent(IN) :: ppbps
  integer,intent(IN) :: ppbpl

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
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Set up Local Variable for copy of state var. object
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),dimension(NO_BOXES_XY):: bphytoc
  real(RLEN),dimension(NO_BOXES_XY):: bphyton
  real(RLEN),dimension(NO_BOXES_XY):: bphytop
  real(RLEN),dimension(NO_BOXES_XY):: bphytos
  real(RLEN),dimension(NO_BOXES_XY):: bphytol
  real(RLEN),dimension(NO_BOXES_XY):: phytoc
  real(RLEN),dimension(NO_BOXES_XY):: phyton
  real(RLEN),dimension(NO_BOXES_XY):: phytop
  real(RLEN),dimension(NO_BOXES_XY):: phytos
  real(RLEN),dimension(NO_BOXES_XY):: phytol
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer                         :: silica_control,nrphyto,i,iout
  real(RLEN)                      :: scalar2
  real(RLEN)                      :: sudc !max growth rate related to max.number of divions perday.
  real(RLEN),parameter            :: low_phytoc=1.0D-10
  real(RLEN),parameter            :: p_optimum=0.90
  real(RLEN),dimension(NO_BOXES_XY) :: ex_sw  !answer of insw_vector ( 0 or 1)
  real(RLEN),dimension(NO_BOXES_XY) :: cx_any,sx_any,rx_any,qx_any 
  real(RLEN),dimension(NO_BOXES_XY) :: mx_any,rx_any_m,px_any !any depth,part
  real(RLEN),dimension(NO_BOXES_XY) :: dl   !day length
  real(RLEN),dimension(NO_BOXES_XY) :: silt
  real(RLEN),dimension(NO_BOXES_XY) :: et,eo,ea
  real(RLEN),dimension(NO_BOXES_XY) :: sumc,sumc0
  real(RLEN),dimension(NO_BOXES_XY) :: sea ,set, ses
  real(RLEN),dimension(NO_BOXES_XY) :: sdo
  real(RLEN),dimension(NO_BOXES_XY) :: sx_main,sx_act,sx_grazing,px_nett,px_affl
  real(RLEN),dimension(NO_BOXES_XY) :: rugc
  real(RLEN),dimension(NO_BOXES_XY) :: sra, srs , srt
  real(RLEN),dimension(NO_BOXES_XY) :: slc
  real(RLEN),dimension(NO_BOXES_XY) :: run_xbc_c,run_xbn_c,runc,runl
  real(RLEN),dimension(NO_BOXES_XY) :: rupn,rupp, rups
  real(RLEN),dimension(NO_BOXES_XY) :: rumn,rump,rums,rum5s
  real(RLEN),dimension(NO_BOXES_XY) :: misn, misp, miss
  real(RLEN),dimension(NO_BOXES_XY) :: luxn, luxp, luxs
  real(RLEN),dimension(NO_BOXES_XY) :: tN
  real(RLEN),dimension(NO_BOXES_XY) :: iN
  real(RLEN),dimension(NO_BOXES_XY) :: iNIp,iNIn,iNIs
  real(RLEN),dimension(NO_BOXES_XY) :: qnoc,qpoc,qsoc
  real(RLEN),dimension(NO_BOXES_XY) :: eN5s
  real(RLEN),dimension(NO_BOXES_XY) :: rrc,rr_dark_m2_c,rrt_m2_c
  real(RLEN),dimension(NO_BOXES_XY) :: rr2c
  real(RLEN),dimension(NO_BOXES_XY) :: rr1c,rr1n,rr1p
  real(RLEN),dimension(NO_BOXES_XY) :: rr6c,rr6n,rr6p,rr6s
  real(RLEN),dimension(NO_BOXES_XY) :: runn,run3n,run4n,runun,ruQuc
  real(RLEN),dimension(NO_BOXES_XY) :: runp,run1p,runup
  real(RLEN),dimension(NO_BOXES_XY) :: runs,run5s
  real(RLEN),dimension(NO_BOXES_XY) :: rbufn,rbufp,rbufs
  real(RLEN),dimension(NO_BOXES_XY) :: xIrr_level,xIrr_av
  real(RLEN),dimension(NO_BOXES_XY) :: rho_Chl
  real(RLEN),dimension(NO_BOXES_XY) :: rate_Chl
  real(RLEN),dimension(NO_BOXES_XY) :: new_Chl
  real(RLEN),dimension(NO_BOXES_XY) :: x_light_limitation
  real(RLEN),dimension(NO_BOXES_XY) :: sx_stress_degradation
  real(RLEN),dimension(NO_BOXES_XY) :: flChydrate,leftover
  real(RLEN),dimension(NO_BOXES_XY) :: rum3n,rumun,rum4n,rum1p,rumup
  real(RLEN),dimension(NO_BOXES_XY) :: M1p,M5s,M3n,M4n,Mun,Mup
  real(RLEN),dimension(NO_BOXES_XY) :: cqun3,cqup,cquR1n
  real(RLEN),dimension(NO_BOXES_XY) :: m3tm2,px_l_aerob_light,px_c_aerob
  real(RLEN),dimension(NO_BOXES_XY) :: jenn,jenp,jens
  real(RLEN),dimension(NO_BOXES_XY) :: eiBPI
  real(RLEN),dimension(NO_BOXES_XY) :: lim_D1,lim_turm
  real(RLEN),dimension(NO_BOXES_XY) :: d0,d1,st
  real(RLEN),dimension(NO_BOXES_XY) :: step,p,xEPS_l,eiBPI_l,sumc_l
  real(RLEN),dimension(NO_BOXES_XY) :: px_c_light,px_l_light,eh_xphyto
  real(RLEN),dimension(NO_BOXES_XY) :: ZERO_vector
  real(RLEN),dimension(NO_BOXES_XY) :: qpBDc,qnBDc,qsBDc,qlBDc,qpBLc,qnBLc,qsBLc, qlBLc
  real(RLEN),dimension(NO_BOXES_XY) :: dx_c_mid_light_m,dx_l_mid_light_m
  real(RLEN),external  ::GetDelta

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! silica_control =0 : no silica component present in cell
  ! silica_control =1 : external regulation of silica limitation & limitation of
  !                      carbon fixation under silica depletion
  ! silica_control =2 : internal regulation of silica limitation & excretion
  !                      of fixed carbon under nutrient stress
  !                      Process description based on:
  !                      Growth physiology and fate of diatoms in the ocean: a review
  !                      G.Sarthou, K.R. Timmermans, S. Blain, & P. Treguer
  !                      JSR 53 (2005) 25-42
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


  nrphyto=p_useparams(bp)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Copy  state var. object in local var
  !  A small concentration is to avoid that one of the state vars. --->0.0
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !biomass per square meter
  bphytoc = D2STATE(:,ppbpc)
  bphyton = D2STATE(:,ppbpn)
  bphytop = D2STATE(:,ppbpp)
  bphytol = D2STATE(:,ppbpl)

  if ( ppbps > 0 )  then
     bphytos=D2STATE(:,ppbps)
  !-----------------------------------------------------------------------
     if ( bphytos(1) < ZERO) then
        cx_any=bphytos
        bphytos=p_qsRc(nrphyto)*bphytoc
        D2STATE(:,ppbps)=bphytos
       write(LOGUNIT,*)'Benphytoc:BP1s is negative, reset to appropiate value'
       write(LOGUNIT,*) 'bphytos=',bphytos,'bphytoc==',bphytoc
        call set_warning_for_getm
     endif
  !-----------------------------------------------------------------------
  !JM-----------------------------------------------------------------------
     if ( bphytoc(1) < ZERO) then
!        cx_any=bphytoc
       write(LOGUNIT,*)'BenPhyto: BP1c is negative, small number added'
       write(LOGUNIT,*)'bphytoc==',bphytoc
       bphytoc=bphytoc+NZERO
       D2STATE(:,ppbpc)=bphytoc
!       call set_warning_for_getm
     endif
  !-----------------------------------------------------------------------
  endif

  ! if CoupledtoBDc(nrphyto)>iiPhytoPlankton ther is already tested tat
  ! the concentration  is lower then 1.0e-5
  if( CoupledtoBDc(nrphyto)>iiPhytoPlankton) then
!   write(LOGUNIT,*) "BenPhyto Reset"
    rx_any=ZERO
    call ResetSource_D2_Vector(ppbpc)
    call ResetSource_D2_Vector(ppbpn)
    call ResetSource_D2_Vector(ppbpp)
    call ResetSource_D2_Vector(ppbps)
    call ResetSource_D2_Vector(ppbpl)
    call sourcesink_flux_vector(iiBen, ppG3c,ppbpc,rx_any) ! source/sink.c
    call sourcesink_flux_vector(iiBen, ppbpc,ppG3c, rx_any )  ! source/sink.c
    call flux_vector(iiBen, ppbpc,ppQ1c, rx_any )  ! source/sink.c
    call flux_vector(iiBen, ppbpc,ppQ2n, rx_any )  ! source/sink.c
    call flux_vector(iiBen, ppbpc,ppQ2c, rx_any )  ! source/sink.c
    call flux_vector(iiBen, ppK4n,ppbpn, rx_any )  ! source/sink.n
    call flux_vector(iiBen, ppK3n,ppbpn, rx_any )  ! source/sink.n
    call flux_vector(iiBen, ppQun,ppbpn, rx_any )  ! source/sink.n
    call flux_vector(iiBen, ppK1p,ppbpp, rx_any )  ! source/sink.n
    call flux_vector(iiBen, ppQ1p,ppbpp, rx_any )  ! source/sink.n
    call flux_vector(iiBen, ppbpc,ppQ6c, rx_any )  ! source/sink.c
    call flux_vector(iiBen, ppbpn,ppQ6n, rx_any )  ! source/sink.c
    call flux_vector(iiBen, ppbpp,ppQ6p, rx_any )  ! source/sink.c
    call flux_vector(iiBen, ppbpn,ppQ1n, rx_any )  ! source/sink.c
    call flux_vector(iiBen, ppbpp,ppQ1p, rx_any )  ! source/sink.c
    call flux_vector(iiBen, ppbpn,ppK4n,rx_any)  ! source/sink.n
    call flux_vector(iiBen, ppbpp,ppK1p,rx_any)  ! source/sink.p
    if ( ppbps > 0 )  then
       call flux_vector(iiBen, ppbps,ppQ6s, rx_any)  ! source/sink.c
       call flux_vector(iiBen, ppK5s,ppbps, rx_any)  ! source/sink.c
       call flux_vector(iiBen, ppbps,ppK5s,rx_any)  ! source/sink.p
    endif
    if(maxval(Dfm) .gt.1.0D-5.or.maxval(Dcm).gt.1.0D-5) then
      rx_any_m=0.5*Dfm;cx_any=Dfm-1.0D-5
      call LimitChange_vector(POSITIVE,rx_any_m,cx_any,max_change_per_step)
      call flux_vector(iiBen, ppDfm,ppDfm,-rx_any_m)  ! source/sink.p
      rx_any_m=0.5*Dcm;cx_any=Dcm-1.0D-5
      call LimitChange_vector(POSITIVE,rx_any_m,cx_any,max_change_per_step)
      call flux_vector(iiBen, ppDcm,ppDcm,-rx_any_m)  ! source/sink.p
    endif
    if (maxval(Dlm) .gt.1.0D-5) then
    rx_any_m=0.5* Dlm;cx_any=Dlm-1.0D-5
    call LimitChange_vector(POSITIVE,rx_any_m,cx_any,max_change_per_step)
    call flux_vector(iiBen, ppDlm,ppDlm,-rx_any_m)  ! source/sink.p
    endif
    shiftDlm=ZERO
!   Output2d_1=ZERO; Output2d_2=ZERO;Output2d_3=ZERO; Output2d_4=ZERO;
    return
  endif

!   write(LOGUNIT,*) "BenPhyto Calculated"
  silica_control=0
  if ( p_qus(nrphyto) > ZERO )  then
    silica_control=2
    if ( p_xqs(nrphyto)>ZERO ) silica_control=3
  elseif ( p_chPs(nrphyto) > ZERO ) then
    silica_control=1
  endif

   ZERO_vector=ZERO
   sx_grazing=Source_D2_withgroup(ppbpc, &
       ppBenOrganisms,iiBenOrganisms,iiC,iiConsumption)/(NZERO+bphytoc)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Temperature response of Phytoplankton
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   et  =   eTq_vector(  ETW_Ben(:),  p_q10(nrphyto))

  !Based on number of Colijn the kd in the sediment in tidal area is 5500
  !a  t a porosity  0.46
  ! %silt which is present in the solution.
  call CouplingBioSilt(1,silt);silt=silt *p_px_silt_in_pw(bp)
  m3tm2=Dlm*p_poro
  !Calculate the fraction of the benthic diatoms above the compensation depth
  ! pxinlight(c/l)= fraction benphyto in the light layer
  px_c_light=CalculateBenPhyto_vector(iiC,INTEGRAL,ZERO_vector,Dlm)
  px_l_light=CalculateBenPhyto_vector(iiL,INTEGRAL,ZERO_vector,Dlm)
  px_c_aerob=CalculateBenPhyto_vector(iiC,INTEGRAL,ZERO_vector,D1m)
  dx_c_mid_light_m= &
              CalculateBenPhyto_vector(iiC,AVERAGE,ZERO_vector,Dlm)
  dx_l_mid_light_m= &
              CalculateBenPhyto_vector(iiL,AVERAGE,ZERO_vector,Dlm)
  !mass of Bentic diatoms per m2 in the light
  phytoc=px_c_light*bphytoc/m3tm2
  phyton=px_l_light*bphyton/m3tm2
  phytop=px_l_light*bphytop/m3tm2
  if ( ppbps > 0 ) phytos=px_l_light*bphytos/m3tm2
  phytol=px_l_light*bphytol/m3tm2

  px_any=DONE-px_l_light
  qnBDc(:)=max(ZERO,px_any*bphyton/(NZERO+(DONE-px_c_light)*bphytoc))
  qpBDc(:)=max(ZERO,px_any*bphytop/(NZERO+(DONE-px_c_light)*bphytoc))
  qlBDc(:)=max(ZERO,px_any*bphytol/(NZERO+(DONE-px_c_light)*bphytoc))
  qsBDc(:)=max(ZERO,px_any*bphytos/(NZERO+(DONE-px_c_light)*bphytoc))
  qnBLc(:)=max(ZERO,phyton/(NZERO+phytoc))
  qpBLc(:)=max(ZERO,phytop/(NZERO+phytoc))
  qlBLc(:)=max(ZERO,phytol/(NZERO+phytoc))
  qsBLc(:)=max(ZERO,phytos/(NZERO+phytoc))
  !maximal density diatoms in top layer if eh_xphyto>1 new input benthic
  ! diatoms will be redistributed over larger depth
  eh_xphyto=max(DONE,px_c_light*bphytoc/p_poro/p_hBPc(bp) )

  M1p=max(ZERO,Kp1p/p_poro/Dlm/(DONE+p_pK1_ae))
  M5s=max(ZERO,Kp5s/p_poro/Dlm/(DONE+p_pK5_ae))
  M3n=max(ZERO,Kp3n/p_poro/Dlm/(DONE+p_pK3))
  M4n=max(ZERO,Kp4n/p_poro/Dlm/(DONE+p_pK4))
  Mun=Qpun/p_poro/Dlm
  if (isnan(Kp4n(1))) then
     write(LOGUNIT,*) 'Kp4n is NaN',KP1p,Kp5s,Kp3n,Dlm
  endif
  cx_any=min(p_qpPhc*Q1c,Q1p(:))
  Mup(:) = max(ZERO,Q1p(:)-cx_any)/p_poro/D1m


  scalar2=NZERO+p_lN3N4n(nrphyto)
  cqun3 = scalar2/( NZERO+scalar2+M4n)
  scalar2=NZERO+p_lureaN4n(nrphyto)
  cquR1n = scalar2/( NZERO+scalar2+M4n)
  scalar2=NZERO+p_lN1(nrphyto)
  cqup  = scalar2/( NZERO+scalar2+ M1p)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient limitation (intracellular) N, P
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  iNIp = min( DONE, max( NZERO, ( qpBLc(:)&
                       - p_qplc(nrphyto))/( p_qpRc(nrphyto)- p_qplc(nrphyto))))
  iNIn = min( DONE, max( NZERO, ( qnBLc(:)&
                      - p_qnlc(nrphyto))/( p_qnRc(nrphyto)- p_qnlc(nrphyto))))
  iNIs = min( DONE, max( NZERO, ( qsBLc(:) &
                     - p_qslc(nrphyto))/( p_qsRc(nrphyto)- p_qslc(nrphyto))))


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Phytoplankton growth is limited by nitrogen and phosphorus
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  select case ( p_limnut)
    case ( 0 ) ; iN  =   (iNIp* iNIn)**(0.5D+00)  ! geometric mean
    case ( 1 ) ; iN  =   min(  iNIp,  iNIn)  ! Liebig rule
    case ( 2 ) ; iN  =   2.0D+00/( DONE/ iNIp+ DONE/ iNIn)  ! combined
  end select

  ! tN controls sedimentation of phytoplankton
  tN= iN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient limitation due to intra- extracellular silicate.
  ! eN5s limit externally nutrient limitation.
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !calculate internal quota
  eN5s  =   DONE
  if ( silica_control > 0 ) then
    select case (silica_control)
      case(1)
         eN5s = min( eN5s,M5s/(M5s + p_chPs(nrphyto)))
         tN=min(iN,eN5s)
      case(2,3)
        tN=min(iNIs,iN)
      end select
  endif

 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 !Photosynthesis(Irradiance EIR_Ben in uE m-2 s-1, xIrr_av is mid-layer EIR_Ben)
 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  xIrr_level=max(NZERO,EIR_Ben(:)) !light at sediment water interface
  if ( sw_precision.eq.0) then
    xEPS_Sedi  = silt * p_epsESS + phytol*p_epsChla  &
          +Q6c*PartQ_vector(D6m(:),ZERO_vector,Dlm,p_d_tot)/Dlm*p_epsR6/p_poro
    xIrr_av =xIrr_level / xEPS_Sedi/Dlm*( DONE- exp_limit(-xEPS_Sedi* Dlm ))
    eiBPI = max(ZERO, DONE- &
      exp_limit(- qlBLc(:)/ p_qchlc(nrphyto)/ p_Ke(nrphyto)* xIrr_av))

    select case ( LightForcingFlag)
      case (1); sumc = p_sum(nrphyto)* et* eiBPI   *  eN5s
      case (2); sumc = p_sum(nrphyto)* et* eiBPI*( SUNQ/ HOURS_PER_DAY) * eN5s
      case (3); sumc = p_sum(nrphyto)* et* eiBPI* ThereIsLight * eN5s
    end select
    sumc=max(ZERO,sumc)
  else
    xEPS_Sedi=ZERO;eiBPI=ZERO;sumc=ZERO;p=ZERO
    d0=ZERO
    st=CalculateBenPhyto_vector(iiC,INTEGRAL,ZERO_vector,Dlm)
    step=st/sw_precision
    do i=1,sw_precision
      p=p+step
      d1=CalculateBenPhyto_vector(iiL,INTEGRAL,ZERO_vector,Dlm,inverse=p)
      px_any=CalculateBenPhyto_vector(iiL,INTEGRAL,d0,d1)
      xEPS_l  = silt * p_epsESS +px_any* bphytol/(d1-d0)/p_poro*p_epsChla  &
          +Q6c*PartQ_vector(D6m(:),d0,d1,p_d_tot)/(d1-d0)*p_epsR6/p_poro
      xIrr_av=xIrr_level / xEPS_l/(d1-d0)*( DONE- exp_limit(-xEPS_l* (d1-d0) ))
      eiBPI_l = max(ZERO, DONE- &
               exp_limit(- qlBLc(:)/ p_qchlc(nrphyto)/ p_Ke(nrphyto)* xIrr_av))
      select case ( LightForcingFlag)
        case (1); sumc_l= p_sum(nrphyto)*et*eiBPI_l*eN5s
        case (2); sumc_l= p_sum(nrphyto)*et*eiBPI_l*eN5s*(SUNQ/HOURS_PER_DAY)
        case (3); sumc_l= p_sum(nrphyto)*et*eiBPI_l*eN5s* ThereIsLight
      end select
      xEPS_Sedi=xEPS_Sedi+xEPS_l
      eiBPI=eiBPI+eiBPI_l
      sumc=sumc+sumc_l
      if ( i< sw_precision) then
        xIrr_level =xIrr_level * exp_limit(-xEPS_l* (d1-d0) )
        d0=d1
      endif
    enddo
    xEPS_Sedi=xEPS_Sedi/sw_precision
    eiBPI=eiBPI/sw_precision
    sumc=(max(ZERO,sumc/sw_precision))
  endif
  if (xEPS_Sedi(1).le.NZERO) then
    write(LOGUNIT,*) silt,p_epsESS , bphytol,Q6c
  endif

  !-----------------------------------------------------------------------
  if ( isnan(sumc(1))) then
   write(LOGUNIT,*) 'BP1: sumc,et,eiBPI_l,eN5s, phytoc,', &
                                        sumc,et,eiBPI_l,eN5s,phytoc
   write(LOGUNIT,*) 'BP1: d0,d1,xEPS_l,xIrr_level,EIR_Ben', &
                                        d0,d1,xEPS_l,xIrr_level,EIR_ben
  endif
  !-----------------------------------------------------------------------

  !density dependent mort.: only used for phytoplankton which is not grazed
  sdo  =   ZERO

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Lysis and excretion
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! select case (p_iRI(nrphyto))
!  case(1)
     sdo  = sdo +  2.0* ( p_thdo(nrphyto)/ &
     ( tN+ p_thdo(nrphyto)))* p_sdmo(nrphyto)  
!  case(2)
!    sdo  = sdo +  2.0* ( p_thdo(nrphyto)/ &
!    ( min(iNIs,iNIn)+ p_thdo(nrphyto)))* p_sdmo(nrphyto)  ! nutr. -stress lysis
! end select

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Rest processes
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ses  =   p_ses(nrphyto)*et            ! activity +rest excretion
  srs  =   p_srs(nrphyto)*et

  rx_any=(sumc*(DONE-p_pu_ea(nrphyto))-srs)*phytoc*m3tm2
  call LimitChange_vector(POSITIVE,rx_any,G3c, max_rate_per_step,p)
#ifdef INCLUDE_BENCO2
   !2 ways to limit Phyto
   rx_any  = rx_any /MW_C
   cx_any=max(ZERO,HCO3ae)*D1m/p_poro
   call LimitChange_vector(POSITIVE,rx_any,cx_any,max_change_per_step,px_any)
   sumc=sumc*min(p,px_any,cx_any/(DONE+cx_any))       ! gross production
#endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Activity processes
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  sea  =   sumc* p_pu_ea(nrphyto)
  sra  =   p_pu_ra(nrphyto)*max(ZERO, sumc- sea)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Potential-Net prim prod. (mgC /m3/d):
  !  a. limit prim prod at high temperature with the (maximum) number of
  !     divisions per day by make use of the length of the light period.
  !  c. assume that activity respiration increases with temperature and
  !     is NOT limited by the number of divisions per day.
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  dl=(SUNQ/ HOURS_PER_DAY)
  sx_any  =max(ZERO, dl*( sumc- sea-sra)-ses-srs)  ! net production
  sudc=p_xdiv(nrphyto)*log(2.0D+00)
  px_any=DONE
  where (sx_any.gt.NZERO) px_any=min(DONE,sudc/(NZERO+sx_any))
  sumc=px_any*sumc; sea=px_any*sea
  sra=sra*px_any**p_xgeneric(nrphyto)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Correct for mortality : in order to model an effective
  ! mortality at the end of the bloom it is necessary to make all Phaeo
  ! which will die during next time-step inactive, otherwise the mortality
  ! is counteracted by the new production
  px_any=DONE-min(max_change_per_step,sdo)*LocalDelta
  sumc0=sumc
  sumc=px_any*sumc

  !-=-=0-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Apportioning over R1 and Q6:
  ! Cell lysis generates both DOM and POM.
  ! The nutr.-depleted cell has a nutrient-carbon ratio equal to p_q?lc.
  ! Assuming that this structural part is not easily degradable,
  ! at least a fraction equal to the minimum quota is released as POM.
  ! Therefore, nutrients (and C) in the structural part go to Q6.
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rx_any=sdo*phytoc
  rr1c  =  p_pe_R1c* rx_any
  rr6c  =  (DONE-p_pe_R1c)* rx_any


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !test on (temporary) anarobic circumstances in the oxic layer during the night
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ex_sw=insw_vector(sumc-ses-sea-sra-srs)
  call CalcZeroOrderOxInLayer(G2o,D1m,ZERO_vector,Dlm,cx_any)
  eo=ex_sw+ (DONE-ex_sw)* &
       (DONE-exp(-(cx_any+Kp3n*p_qon_nitri)/p_poro/D1m/p_clO2o(nrphyto)))
  ea=p_xeff_an+(DONE-p_xeff_an)*eo
  srs=srs/ea

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! total processes
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  set  =   sea+ ses     ! total specific excretion
  srt  =   sra+ srs     ! total specific respiration
  slc  =   sea+ sra+ sdo  ! specific loss terms

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Production, productivity and respiration flows
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  px_l_aerob_light=DONE; where (Dlm>D1m) px_l_aerob_light=px_c_aerob
  rugc=sumc * phytoc  ! gross production
  call sourcesink_flux_vector(iiBen, ppG3c,ppbpc, rugc*m3tm2 )
  call flux_vector(iiBen, ppG2o,ppG2o, px_l_aerob_light*rugc*m3tm2/ MW_C )
  call flux_vector(iiBen, ppK16r,ppK16r,  &
                        -p_qro*(DONE-px_l_aerob_light)*rugc*m3tm2/ MW_C )
  ruBPc=ruBPc+rugc*m3tm2

  rrc  =   srt* phytoc  ! total actual respiration

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Potential-Net prim prod. (mgC /m3/dd0
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  run_xbc_c  =   max(  NZERO, ( sumc- slc)* phytoc)  ! net production
  sugBI(:,bp)  =   run_xbc_c/( NZERO+ phytoc)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !  Determintation of optimum quotum..
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  qnoc=p_xqn(nrphyto)*p_qnRc(nrphyto)
  qpoc=p_xqp(nrphyto)*p_qpRc(nrphyto)
  if ( silica_control>= 2) qsoc=p_xqs(nrphyto)*p_qsRc(nrphyto)

! px_any=min(qnBLc(:)/p_qnRc(nrphyto),qpBLc(:)/p_qpRc(nrphyto))
! if ( silica_control>= 2) px_any=min(px_any,qsBLc(:)/p_qsRc(nrphyto))
! qnoc=p_qnRc(nrphyto)*min(max(p_qnlc(nrphyto)/p_qnRc(nrphyto), &
!                          (px_any+p_xqn(nrphyto))*0.5),p_xqn(nrphyto))
! qpoc=p_qpRc(nrphyto)*min(max(p_qplc(nrphyto)/p_qpRc(nrphyto), &
!                          (px_any+p_xqp(nrphyto))*0.5),p_xqp(nrphyto))
! if ( silica_control>= 2) qsoc=p_qsRc(nrphyto)*min( &
!       max(p_qslc(nrphyto)/p_qsRc(nrphyto),(px_any+p_xqs(nrphyto))*0.5),p_xqs(nrphyto))


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient Uptake: calculate maximal uptake of N,P
  ! Check if C-fixation is larger to make of all C new biomass
  ! If not increase excretion of C
  ! In BenLimitnutrientuptake at forehand the availability of nutrients  for
  ! each group is calculated on bases of affinity (p_qu's) and the biomass
  ! In case of nutrient limitations with this predictor-corrector
  ! nutrients are distributed such that the chance on negative values
  ! of the nutrient concentrations is 0.
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  px_any=fr_lim_BPI_n(:,bp)*max(DONE,K4n/(NZERO+Kp4n))
  rum4n =         p_qun(nrphyto)* phytoc *M4n
  call DoubleLimitChange_vector(POSITIVE,rum4n,min(Kp4n,K4n)/m3tm2,px_any, &
                                                       max_change_per_step)
  px_any=fr_lim_BPI_n(:,bp)*max(DONE,K3n/(NZERO+Kp3n))
  rum3n =  cqun3 *p_qun(nrphyto)* phytoc *M3n
  call DoubleLimitChange_vector(POSITIVE,rum3n,min(Kp3n,K3n)/m3tm2,px_any, &
                                                       max_change_per_step)
  px_any=fr_lim_BPI_n(:,bp)*max(DONE,Qpun/(NZERO+Qun))
  rumun =  cquR1n*p_qun(nrphyto)* phytoc *Mun
  call DoubleLimitChange_vector(POSITIVE,rumun,min(Qpun,Qun)/m3tm2,px_any, &
                                                       max_change_per_step)
  px_any=fr_lim_BPI_p(:,bp)*max(DONE,K1p/(NZERO+Kp1p))
  rum1p =         p_qup(nrphyto)* phytoc *M1p
  call DoubleLimitChange_vector(POSITIVE,rum1p,min(Kp1p,K1p)/m3tm2,px_any, &
                                                       max_change_per_step)
  rumup =  cqup  *p_qup(nrphyto)* phytoc *Mup
  call DoubleLimitChange_vector(POSITIVE,rumup,Q1p,px_any,max_change_per_step)
  if ( silica_control  >=  2) then
    rum5s =         p_qus(nrphyto)* phytoc *M5s
    ! no-predictor-corrector value needed for silicate: there is only 1 consumer
    cx_any=min(Kp5s,K5s)/m3tm2
    call LimitChange_vector(POSITIVE,rum5s,cx_any,max_change_per_step)
  endif

  rbufn=max(ZERO,run_xbc_c) *p_xqn(nrphyto) *p_qnRc(nrphyto)
  call LimitChange_vector(POSITIVE,rbufn,phyton, max_rate_per_step)
  ! N uptake brom 4 sources: nitrate, ammonium, organic N, from buffer incell.
  ex_sw=insw_vector(phyton- phytoc*p_qnlc(nrphyto))
  rumn(:)  =   rum4n(:) +rum3n(:) + rumun(:) +rbufn(:)*ex_sw

  ! P uptake brom 3 sources: phosphate, organic P, from buffer in cell.
  rbufp=run_xbc_c*p_xqp(nrphyto)*p_qpRc(nrphyto)
  call LimitChange_vector(POSITIVE,rbufp,phytop, max_rate_per_step)
  ex_sw=insw_vector(phytop- phytoc*p_qplc(nrphyto))
  rump(:)= rum1p(:)+rumup(:)+ex_sw*rbufp  ! max pot. uptake of NI

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Check if which fraction C-fixation can be used for new biomass
  ! by checking the potential nutrient avilability
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  run_xbn_c(:) = min(run_xbc_c,rumn/p_qnlc(nrphyto), rump/p_qplc(nrphyto))

  ! Si uptake brom 2 sources: dissolved silicate, from diatom-skeleton
  if ( silica_control  >=  2) then
    rbufs =min(ZERO,run_xbc_c )*p_xqs(nrphyto) *p_qsRc(nrphyto)
    call LimitChange_vector(POSITIVE,rbufs,phytos, max_rate_per_step)
    ex_sw=insw_vector(phytos- phytoc*p_qslc(nrphyto))
    rums= rum5s+ex_sw*rbufs
    run_xbn_c(:) = min( run_xbn_c(:), rums/ p_qslc(nrphyto))
  endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! C flows ; phyto --> detritus
  ! a. activity excretion of sugars:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  flChydrate  =   set* phytoc

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! b.stress:All C which cannot be used for growth is excreted :
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  runc  =   max(ZERO,min(run_xbc_c,run_xbn_c))
  leftover = run_xbc_c-runc
  ! distribute too much fixed carbon over stress respiration and excretion
  flChydrate  =   flChydrate+ leftover

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! C flows: calculate fluxes
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rr2c=ZERO
  select case (p_iRI(nrphyto))
     case (iiR1) ; rr1c =flChydrate + rr1c
     case (iiR2) ; rr2c= flChydrate
  end select

  call flux_vector( iiBen, ppbpc,ppQ1c, rr1c*m3tm2 )  ! source/sink.c
  call flux_vector( iiBen, ppbpc,ppQ2c, rr2c*m3tm2 )  ! source/sink.c
  rx_any=max(ZERO,rr2c* min(DONE,qnBLc(:)-p_qnlc(nrphyto)) &
                       /(p_qnRc(nrphyto)-p_qnlc(nrphyto))*p_qnR2c)
  call flux_vector( iiBen, ppbpn,ppQ2n, rx_any*m3tm2 )  ! source/sink.c

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Apparent Net prim prod. (mgC /m3/d)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  sunBI(:,bp)  =   runc/( NZERO+ phytoc)
  sx_main=ses+srs
  sx_act=max(sx_main,sumc)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate compensation depth
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  xIrr_level = max( NZERO, EIR_Ben) !light at sediment water interface
  !Determine average net production.
  mx_any=ZERO;px_nett=ZERO;px_affl=ZERO
  where (xIrr_level.gt.NZERO)
    !Calc net production
    !use "potential production"
    sx_any=p_sum(nrphyto)* et * eN5s
    !Calc losses of new production : activity respiration + excretion
    px_nett=max(NZERO,sumc-sra-sea)/(NZERO+sumc)
    !affinity for light
    px_affl=- qlBLc(:)/ p_qchlc(nrphyto)/ p_Ke(nrphyto)* xIrr_level
    px_any=max(ZERO,DONE- (sx_main/(sx_any*px_nett)))
    mx_any=max(ZERO,-log((log(px_any)/px_affl))/xEPS_Sedi)
  endwhere

  shiftDlm=(mx_any-Dlm)*insw_vector(mx_any)
  call LimitChange_vector(NEGATIVE,shiftDlm,Dlm,max_change_per_step)
  call flux_vector(iiBen,ppDlm,ppDlm,shiftDlm)
  !-----------------------------------------------------------------------
  if ( isnan(shiftDlm(1)).or.abs(shiftDlm(1))*LocalDelta.gt.p_d_tot) then
    write(LOGUNIT,*)'Nan ShiftDlm---------------------------'
   write(LOGUNIT,*)'mx_any,px_any,px_affl,xEPS_Sedi,sx_main,sx_xany,px_nett',&
                     mx_any,px_any,px_affl,xEPS_Sedi,sx_main,sx_any,px_nett
    write(LOGUNIT,*)'eN5s,sdo,xIrr_level,EIR_Ben',en5s,sdo,xIrr_level,EIR_Ben
    write(LOGUNIT,*)'runc,run_xbc_c,eo,ea',runc,run_xbc_c,eo,ea
    write(LOGUNIT,*),'phyto[cnpsl]',phytoc,phyton,phytop,phytos,phytol
    write(LOGUNIT,*),'bphyto[cnpsl]',bphytoc,bphyton,bphytop,bphytos,bphytol
    write(LOGUNIT,*)'silt,Q6,PartQ,Dlm', &
               silt ,Q6c,PartQ_vector(D6m(:),ZERO_vector,Dlm,p_d_tot),Dlm
  endif
  !-----------------------------------------------------------------------

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient dynamics: NITROGEN
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  qx_any = p_xqn(nrphyto)* p_qnRc(nrphyto)
  rupn = max(ZERO,qnoc* (runc -(sx_grazing+sx_main)*phytoc/dl) )
  qx_any = qnoc- qnBLc(:)
  ex_sw=insw_vector(qx_any)
  misn= (sx_act*qx_any - (DONE-ex_sw)*sx_main*qnBLc(:))*phytoc
  luxn=sx_act*(p_xqn(nrphyto)*p_qnRc(nrphyto)-qnoc)*phytoc

  runn  =   min( rumn, rupn +misn) !actual uptake of NI

  !Calculate the uptake of the different sources.uptake form the cell
  ! itself changes automatically internal buffer

  ex_sw  =   insw_vector(runn)
  run4n  =   ex_sw* runn* rum4n(:)/( NZERO+ rumn)  ! actual uptake of Nn
  run3n  =   ex_sw* runn* rum3n(:)/( NZERO+ rumn)  ! actual uptake of Nn
  runun  =   ex_sw* runn* rumun(:)/( NZERO+ rumn)  ! actual uptake of Nn
  !-----------------------------------------------------------------------
  if (isnan(run3n(1)).or.isnan(run4n(1))) then
     write(LOGUNIT,*) 'run3n or run4n=NaN run3n,run4n,rumun,runn,rumn:', &
                                       run3n,run4n,rumun,runn,rumn
     write(LOGUNIT,*) 'M3n M4n Mun,K4n,kP4n,K3n,Kp3n=',  &
                                         M3n,M4n,Mun,K4n,Kp4n,K3n,Kp3n
     write(LOGUNIT,*) 'BP1c,BP1n,BP1p,BP1s=',bphytoc,bphyton,bphytop,bphytos
     write(LOGUNIT,*) 'Dfm,Dlm=',Dfm,Dlm,phytoc,phyton,phytop
  endif
  !-----------------------------------------------------------------------

  call flux_vector( iiBen, ppK4n, ppbpn,  run4n*m3tm2 )  ! source/sink.n
  call flux_vector( iiBen, ppK3n, ppbpn,  run3n*m3tm2 )  ! source/sink.n
  ! oxygen which freed by phytoplankton when they take up nitrate
  call flux_vector( iiBen, ppG2o,ppG2o,run3n*p_qon_nitri*m3tm2 )  !source/sink.o
  call flux_vector( iiBen, ppQun, ppbpn,  runun*m3tm2 )  ! source/sink.n
  !CO2 production urease: ureae ->NH4 +CO2
  ruQuc=runun/p_qnUc
  call flux_vector( iiBen, ppQ1c,ppbpc,ruQuc*m3tm2)

  ruBPn =ruBPn  +ex_sw*(run3n+run4n+runun)*m3tm2
  !  correction in winter situation of negative growth to keep
  !  nutrients below max.ratio
  runn  =   min( rumn,rupn+ misn+luxn)  ! actual uptake
  ex_sw=insw_vector(runn)

  jenn=- runn*(DONE-ex_sw)*m3tm2  ! source/sink.n

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nuttrient dynamics: PHOSPHORUS
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Intracellular missing amount of P
  qx_any = p_xqp(nrphyto)* p_qpRc(nrphyto)
  rupp = max(ZERO,qpoc* (runc-(sx_grazing+sx_main)*phytoc/dl) )
  qx_any = qpoc- qpBLc(:)
  ex_sw=insw_vector(qx_any)
  misp =(sx_act*qx_any -(DONE-ex_sw)*sx_main*qpBLc(:))* phytoc
  luxp=sx_act*(p_xqp(nrphyto)*p_qpRc(nrphyto)-qpoc)*phytoc
  ex_sw=insw_vector(misp)

  runp  =   min( rump, rupp +misp)

  !Calculate the uptake of the different sources. (uptake form the
  ! cell itself changes automatically internal buffer)
  ex_sw  =   insw_vector(runp)
  run1p  =   ex_sw* runp* rum1p(:)/( NZERO+ rump)  ! actual uptake of Nn
  runup  =   ex_sw* runp* rumup(:)/( NZERO+ rump)  ! actual uptake of Nn

  call flux_vector( iiBen, ppK1p,ppbpp, run1p*m3tm2 )  ! source/sink.n
  call flux_vector( iiBen, ppQ1p,ppbpp, runup*m3tm2 )  ! source/sink.n
  ruBPp=ruBPp+(run1p+runup)*m3tm2
  !-----------------------------------------------------------------------
  if (isnan(ruBPp(1))) then
     write(LOGUNIT,*) 'BP:ruBPp is Nan, run1p,runup',run1p,runup
  endif
  !-----------------------------------------------------------------------

  !correction in winter situation of negative
  ! growth to keep nutrients below max.ratio
  runp  =   min( rump,rupp+ misp+luxp)  ! actual uptake
  ex_sw=insw_vector(runp)

  jenp= - runp* (DONE-ex_sw)*m3tm2  ! source/sink.p

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Excretion of N and P to PON and POP
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rx_any=sdo*qnBLc(:)*phytoc;cx_any=phyton-p_qnlc(nrphyto)*low_phytoc
  ! extreme cases no flux at lower below :
  call LimitChange_vector(POSITIVE,rx_any,cx_any,max_change_per_step)
  rr1n  =   p_pe_R1n* rx_any
  rr6n  =   rx_any- rr1n

  rx_any=sdo*qpBLc(:)*phytoc;cx_any=phytop-p_qplc(nrphyto)*low_phytoc
  ! extreme cases no flux at lower below :
  call LimitChange_vector(POSITIVE,rx_any,cx_any,max_change_per_step)
  rr1p  =   p_pe_R1p* rx_any
  rr6p  =   rx_any- rr1p

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient dynamics: SILICATE
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rr6s=ZERO
  jens=ZERO
  if ( silica_control > 0 )  then
    select case (silica_control)
    case (1)
     runs = max(ZERO, p_qsRc(nrphyto) * runc );          ! net uptake
     call flux_vector( iiBen, ppK5s,ppbps, runs*m3tm2)  ! source/sink.c
     ruBPs=ruBPs+runs*m3tm2
    case (2,3)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      !  Silia uptake according Droop: Si uptake based on C uptake
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      qx_any = p_xqs(nrphyto)* p_qsRc(nrphyto)
      rups = max(ZERO,p_optimum*qsoc* (runc-(sx_grazing+sx_main)*phytoc/dl) )
      !all silica can only be used to make new biomass!
      ! miss is here only used to correct for too much silica in the cells
      qx_any = qsoc- qsBLc(:)
      ex_sw=insw_vector(qx_any)
      miss=(sx_act*qx_any - (DONE-ex_sw)*sx_main*qsBLc(:)) * phytoc
      luxs=sx_act*(p_xqs(nrphyto)*p_qsRc(nrphyto)-qsoc)*phytoc
      runs  =   min( rums, rups +miss)
      !Calculate the uptake of the different sources.
      !(uptake from the cell itself changes automatically internal buffer)
      ! source/sink.c
      ex_sw=insw_vector(runs)
      run5s  =   ex_sw* runs* rum5s/( NZERO+ rums)
      call flux_vector( iiBen, ppK5s,ppbps, ex_sw*run5s*m3tm2 )
      ruBPs=ruBPs+runs*m3tm2*ex_sw
      !
      !
      runs  =   min( rums,rups+ miss+luxs)  ! actual uptake
      ex_sw=insw_vector(runs)
      jens = -runs*(DONE-ex_sw)*m3tm2  ! source/sink.c
    end select

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Losses of Si
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    rx_any=sdo*qsBLc(:)*phytoc;cx_any=phytos-p_qslc(nrphyto)*low_phytoc
    ! extreme cases no flux at lower below :
    call LimitChange_vector(POSITIVE,rx_any,cx_any,max_change_per_step)
    rr6s  = rx_any    ! Lysis, particulate
  endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate the changes in the distribution due to new production-respiration
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !Calculate the depth where the modulus of the  BDIA is found
  !Calculate the fraction between actual mass and maximal mass
  !Correct Dfm for primary production
  rx_any_m=(dx_c_mid_light_m*eh_xphyto-Dfm) *runc*m3tm2/(NZERO+bphytoc)
  call LimitChange_vector(NEGATIVE,rx_any_m,Dfm,max_change_per_step)
  call flux_vector(iiBen, ppDfm,ppDfm,rx_any_m)
  call FindNaNInRates(iiBen,ppDfm,'Correction for production')
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !respiration  in the production layer (oxygen units)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Respiration flows : respiration
  rrt_m2_c=rrc*m3tm2
  ! flux above the minimum of compensation depth and oxygen penetration depth
  ! flux in denitrification layer but above the compensation depth
  call flux_vector(iiBen,ppK16r,ppK16r, &
                              p_qro*(DONE-px_l_aerob_light)*rrt_m2_c )
  call flux_vector(iiBen,ppG2o,ppG2o,-px_l_aerob_light* rrt_m2_c/MW_C )
  !only changes in Dfm due to respiration other activity respitation is already
  !included in the changes due to netto production
  rx_any_m=-(dx_c_mid_light_m-Dfm) *(srs*phytoc*m3tm2/bphytoc)
  call LimitChange_vector(NEGATIVE,rx_any_m,Dfm,max_change_per_step)
  call flux_vector(iiBen, ppDfm,ppDfm,rx_any_m)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !respiration below the production layer (oxygen units)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  cx_any=(bphytoc-phytoc*m3tm2)
  ! in p_any part of the dark respiration which take place in the aerobic part.
  where (D1m>Dlm)
     rx_any=srs*cx_any*(px_c_aerob-px_c_light)
     ! flux in oxic layer below the compensation depth
     mx_any= CalculateBenPhyto_vector(iiC,AVERAGE,Dlm,D1m)
     rx_any_m=-(mx_any-Dfm) *(rx_any/bphytoc)
     rrt_m2_c=rrt_m2_c+rx_any
  elsewhere
     rx_any_m=ZERO;rx_any=ZERO
  endwhere
  call flux_vector(iiBen,ppG2o,ppG2o,-rx_any/MW_C)
  call LimitChange_vector(NEGATIVE,rx_any_m,Dfm,max_change_per_step)
  call flux_vector(iiBen, ppDfm,ppDfm,rx_any_m)
  ! on case that the production layers finish in the anerobic layer
  ! the oxygen produced by prim.prod >= respiration.
  ! flux in the anaerobic layer
  rr_dark_m2_c=srs*ea/p_xeff_an*bphytoc*(DONE-px_l_aerob_light) &
       *min(p_xqn(nrphyto),max(DONE/p_xqn(nrPhyto),p_qchlc(nrphyto)/qlBDc(:)))
  mx_any=p_d_tot
  mx_any= CalculateBenPhyto_vector(iiC,AVERAGE,D1m,mx_any)
  rx_any_m=-(mx_any-Dfm)*(rr_dark_m2_c/bphytoc)
  call flux_vector(iiBen,ppK16r,ppK16r,p_qro*rr_dark_m2_c/MW_C)
  call LimitChange_vector(NEGATIVE,rx_any_m,Dfm,max_change_per_step)
  call flux_vector(iiBen, ppDfm,ppDfm,rx_any_m)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !total respiration plus CO2 freed by urease
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call sourcesink_flux_vector( iiBen, ppbpc,ppG3c,  &
                                         rrt_m2_c+rr_dark_m2_c+ruQuc*m3tm2)
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! Chl-a synthesis and photoacclimation
! C- fixation is limited at high temperatures leading to a lower growth rate.
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   x_light_limitation=min(DONE,p_EIR(nrphyto)/(NZERO+EIR_Ben))
   qx_any=p_lqnlc(nrphyto)*max(ZERO,x_light_limitation)
   px_any=max(ZERO,DONE- &
      (qx_any+qnBLc(:)-p_qnlc(nrphyto))/ (p_lqnlc(nrphyto)-p_qnlc(nrphyto)))
   sx_stress_degradation=p_sdchl_l(nrphyto)*px_any
   cx_any=DONE
   call LimitChange_vector(POSITIVE,sx_stress_degradation, &
                                              cx_any,max_change_per_step)
  where (xIrr_level> ZERO) !light at sediment water interface
    rho_Chl = p_qchlc(nrphyto)* p_qchlc(nrphyto)* eiBPI &
      * phytoc/(phytol+ NZERO) * p_Ke(nrphyto)/ xIrr_level*x_light_limitation
  elsewhere
    rho_Chl=ZERO
  endwhere

  new_Chl = rho_Chl*((sumc-sra)* phytoc-flChydrate)
  !determine net growth without loss
  runl=new_Chl - (sdo+sx_stress_degradation)*phytol
  rate_Chl=runl*m3tm2-et*p_sdchl_d(nrphyto) * bphytol

  rx_any_m=(dx_l_mid_light_m*eh_xphyto-Dcm)*runl*m3tm2/(NZERO+bphytol)
  call LimitChange_vector(NEGATIVE,rx_any_m,Dcm,max_change_per_step)
  !Correct Dcm for primary production
  call flux_vector(iiBen,ppDcm,ppDcm,rx_any_m)
! OUtput2d_1=Source_D2_vector(ppDcm,0)

  call LimitChange_vector(NEGATIVE,rate_Chl,bphytol,max_change_per_step)
  !calculate net rate per m2
  call flux_vector( iiBen, ppbpl,ppbpl,rate_Chl )


!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! Define total fluxes from benthic diatom  to detritus
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call flux_vector( iiBen, ppbpc,ppQ6c, rr6c*m3tm2 )  ! source/sink.c
  call flux_vector( iiBen, ppbpn,ppQ6n, rr6n*m3tm2 )  ! source/sink.c
  call flux_vector( iiBen, ppbpp,ppQ6p, rr6p*m3tm2 )  ! source/sink.c
  jPIQ6s(:,bp)  =   jPIQ6s(:,bp)      + rr6s*m3tm2

  call flux_vector( iiBen, ppbpn,ppQ1n, rr1n*m3tm2 )  ! source/sink.c
  call flux_vector( iiBen, ppbpp,ppQ1p, rr1p*m3tm2 )  ! source/sink.c


!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! respiration:
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!remove nutrients freed by respiration of diatoms which are longer in the dark.
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! processes below D1m----------------------------------------------------
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  jenn=jenn+rr_dark_m2_c*max(ZERO,qnBDc(:)-p_xqn(nrphyto)*p_qnRc(nrphyto))
  jenp=jenp+rr_dark_m2_c*max(ZERO,qpBDc(:)-p_xqp(nrphyto)*p_qpRc(nrphyto))
  if ( silica_control > 0 )&
    jens=jens+rr_dark_m2_c*max(ZERO,qsBDc(:)-p_xqs(nrphyto)*p_qsRc(nrphyto))
  mx_any=p_d_tot
  !calculate Chlc/quotum in anerobic part of sediment
  qx_any=bphytol*CalculateBenPhyto_vector(iiL,INTEGRAL,max(Dlm,D1m),mx_any)/ &
   (NZERO+bphytoc* CalculateBenPhyto_vector(iiC,INTEGRAL,max(Dlm,D1m),mx_any))
  !in rx_any change in Chl concentration degradated due to dark respiration
  rx_any=rr_dark_m2_c*max(ZERO,qx_any-p_xqn(nrphyto)*p_qchlc(nrphyto))

  call flux_vector(iiBen, ppbpn,ppK4n,jenn)  ! source/sink.n
  call flux_vector(iiBen, ppbpp,ppK1p,jenp)  ! source/sink.p
  call flux_vector(iiBen, ppbps,ppK5s,jens)
  call flux_vector(iiBen, ppbpl,ppbpl,-rx_any)
  mx_any=p_d_tot
  mx_any= CalculateBenPhyto_vector(iiL,AVERAGE,max(Dlm,D1m),mx_any)
  rx_any_m=-(mx_any-Dcm) *(rx_any/bphytol)
  call LimitChange_vector(NEGATIVE,rx_any_m,Dcm,max_change_per_step)
  call flux_vector(iiBen, ppDcm,ppDcm,rx_any_m)
! OUtput2d_2=Source_D2_vector(ppDcm,0)
  mx_any=p_d_tot
  mx_any= CalculateBenPhyto_vector(iiC,AVERAGE,max(Dlm,D1m),mx_any)
  rx_any=rr_dark_m2_c*(DONE-DONE/p_xeff_an)
  rx_any_m=-(mx_any-Dfm) *(rx_any/bphytoc)
  call LimitChange_vector(NEGATIVE,rx_any_m,Dfm,max_change_per_step)
  call flux_vector(iiBen, ppDfm,ppDfm,rx_any_m)

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! Calculate the changes in the distribution of part.detritus due to  mortality
! of benthic diatoms.
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rx_any=(dx_c_mid_light_m- D6m(:))*(rr6c*m3tm2)/(NZERO+ Q6c(:))
  call LimitChange_vector(NEGATIVE,rx_any,D6m,max_change_per_step)
  call flux_vector(iiBen, ppD6m,ppD6m,rx_any)
  rx_any=(dx_l_mid_light_m- D7m(:))*(rr6n*m3tm2)/(NZERO+ Q6n(:))
  call LimitChange_vector(NEGATIVE,rx_any,D7m,max_change_per_step)
  call flux_vector(iiBen, ppD7m,ppD7m,rx_any)
  rx_any=(dx_l_mid_light_m- D8m(:))*(rr6p*m3tm2)/(NZERO+ Q6p(:))
  call LimitChange_vector(NEGATIVE,rx_any,D8m,max_change_per_step)
  call flux_vector(iiBen, ppD8m,ppD8m,rx_any)
  rx_any=(dx_l_mid_light_m- D9m(:))*(rr6s*m3tm2)/(NZERO+ Q6s(:))
  call LimitChange_vector(NEGATIVE,rx_any,D9m,max_change_per_step)
  call flux_vector(iiBen, ppD9m,ppD9m,rx_any)

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! Bioturbation of denthic diatoms due activity of macrobenthos
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  mx_any=p_cturm
! lim_D1=D1m/(D1+Dfm);lim_turm=p_cturm/(p_cturm+Dfm)
  lim_D1=1-exp(-D1m/Dfm);lim_turm=1-exp(-p_cturm/Dfm)
  rx_any_m=p_Etur*turenh(:)/Dfm * (pxturinD1 * px_c_aerob *lim_D1 + &
    (DONE-pxturinD1)*lim_turm*CalculateBenPhyto_vector(iiC,INTEGRAL,D1m,mx_any))
  call LimitChange_vector(NEGATIVE,rx_any_m,Dfm,max_change_per_step)
  call flux_vector(iiBen, ppDfm,ppDfm,rx_any_m)
! lim_D1=D1m/(D1+Dcm);lim_turm=p_cturm/(p_cturm+Dcm)
  lim_D1=1-exp(-D1m/Dcm);lim_turm=1-exp(-p_cturm/Dcm)
  rx_any_m=p_Etur* turenh(:) /Dcm * (pxturinD1 *lim_D1 * &
      CalculateBenPhyto_vector(iiL,INTEGRAL,ZERO_vector,D1m) +lim_turm * &
      (DONE-pxturinD1) * CalculateBenPhyto_vector(iiL,INTEGRAL,D1m,mx_any))
  call LimitChange_vector(NEGATIVE,rx_any_m,Dcm,max_change_per_step)
  call flux_vector(iiBen, ppDcm,ppDcm,rx_any_m)
! OUtput2d_3=Source_D2_vector(ppDcm,0)
  call FindNaNInRates(iiBen,ppDfm,'Correction for bioturbation')
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! Calculate the effect on distribution of benthic diatoms by migration of cells
! in the direction of the light
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !running upward
  ! in cases of high concentration near surface running upward is limited.
  rx_any_m=- p_xladm(bp)*insw_vector(EIR_Ben)/max(DONE,eh_xphyto)* &
                                                   min(DONE,Dfm/(NZERO+Dlm))
  !Correct Dfm for transport in upward direction
  call LimitChange_vector(NEGATIVE,rx_any_m,Dfm,max_change_per_step)
  call flux_vector(iiBen, ppDfm,ppDfm,rx_any_m)
  !Correct Dfm for transport in upward direction
  rx_any_m=- p_xladm(bp)*insw_vector(EIR_Ben)/max(DONE,eh_xphyto) &
                                                   *min(DONE,Dcm/(NZERO+Dlm))
  call LimitChange_vector(NEGATIVE,rx_any_m,Dcm,max_change_per_step)
  call flux_vector(iiBen, ppDcm,ppDcm,rx_any_m)
  call FindNaNInRates(iiBen,ppDcm,'Correction for migration')
! OUtput2d_4=Source_D2_vector(ppDcm,0)


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! MODEL  BFM - Biogeochemical Flux Model version 2.50
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  jnetBPTc=jnetBPTc+((sumc-sra-sdo)*phytoc-flChydrate)*m3tm2

  end
  !BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
