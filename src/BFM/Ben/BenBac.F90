#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: BenBac
!
! DESCRIPTION
!   !    This submodel describes the carbon dynamics and associated
!    nutrient dynamics in benthic bacteria (represented
!    by state variables H1-H2)
!
!

!   This file is generated directly from OpenSesame model code, using a code
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine BenBacDynamics(hx,  pphxc, pphxn, pphxp)
!
! !USES:

  ! For the following Benthic-states fluxes are defined: Q6c, Q6n, Q6p, G2o, &
  ! K16r, D6m, D7m, D8m
  ! The following Benthic-states are used (NOT in fluxes): D1m, D2m
  ! For the following Benthic-group-states fluxes are defined: &
  ! BenLabileDetritus, BenthicAmmonium, BenthicPhosphate
  ! The following Benthic 1-d global boxvars are modified :  &
  ! The following Benthic 1-d global boxvars  are used: ETW_Ben
  ! The following Benthic 2-d global boxvars got a value: ruQic,reQic
  ! The following groupmember vars  are used: iiH1, iiH2
  ! The following constituent constants  are used: iiC, iiN, iiP
  ! The following 0-d global parameters are used: p_d_tot, p_small, &
  ! p_peZ_R1c, p_peZ_R1n, p_peZ_R1p, p_qro
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,DONE,NZERO,ZERO
#ifdef NOPOINTERS
  use mem,  ONLY: D2STATE
#else
  use mem, ONLY: D2STATE, Q6c, Q6n, Q6p, D6m,K1p,K11p,K21p,K3n,K4n,K14n,K24n, &
    D7m, D8m, D1m, D2m, G2o,Q11c,Q21c,Q2c,Q2n,K13n,K23n,K6r
#endif
  use mem,ONLY:ppQ6c,ppQ6n,ppQ6p, ppG3c,ppG13c,ppG23c, ppG2o, ppK6r,ppK16r,ppK26r, &
    ppD6m,ppD7m,ppD8m,ppQ2n,ppQ2c,ppK3n,ppK21p,ppK24n,NO_BOXES_XY, &
    ppBenUrea,ppBenLabileDetritus,ppBenthicAmmonium, &
    jcrrBTo, flux_vector,max_change_per_step, max_rate_per_step, &
    iiC,iiN, iiP, iiBen,iiH1,iiH2,iiH3,BenUrea,BenLabileDetritus,&
    ppBenthicPhosphate, ETW_Ben,G2_xavail_o, sourcesink_flux_vector, &
    KQ1c,DH2m,DH3m,ppDH2m,ppDH3m,jnetHIc,fr_lim_HI_n,fr_lim_HI_p,fr_lim_HI_o, &
    iiConsumption,LocalDelta,ppBenOrganisms,iiBenOrganisms

#ifdef INCLUDE_BENCO2
    use mem,ONLY:pHAn,phDn
#endif
  use mem_Param, ONLY: p_d_tot, p_peZ_R1c, p_peZ_R1n, p_peZ_R1p, p_qro,p_xeff_an
  use mem_Param, ONLY: p_poro,p_p_ads_K1_ae=>p_pK1_ae,p_qpPhc,combine_anabac
  use mem_BenAmmonium, ONLY:p_p_ads_K4=>p_p
  use mem_BenPhosphate, ONLY:p_p_ads_K1_an=>p_p_an
  use mem_BenBac
  use constants,  ONLY:p_qnUc,MW_C,LAYER2,LAYER3,POSITIVE,NEGATIVE
  use LimitRates, ONLY:LimitChange_vector,DoubleLimitChange_vector
  use SourceFunctions,only: Source_D2_withgroup

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following vector functions are used:eTq_vector, MM_vector, &
  ! PartQ_vector, insw_vector
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun, ONLY: eTq_vector, MM_vector, PartQ_vector, insw_vector
  use global_interface, only: FindMidPointExp_vector

! use mem,  ONLY: Output2d_1,Output2d_2,Output2d_3, Output2d_4

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

! !INPUT:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer,intent(IN)  :: hx
  integer,intent(IN) :: pphxc
  integer,intent(IN) :: pphxn
  integer,intent(IN) :: pphxp

!
!
! !AUTHORS
!   W. Ebenhoh and C. Kohlmeier
!
!
!
! !REVISION_HISTORY
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
  real(RLEN),dimension(NO_BOXES_XY) :: hxc
  real(RLEN),dimension(NO_BOXES_XY) :: hxn
  real(RLEN),dimension(NO_BOXES_XY) :: hxp
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer                            :: iout
  real(RLEN),dimension(NO_BOXES_XY)  :: clm
  real(RLEN),dimension(NO_BOXES_XY)  :: cm
  real(RLEN),dimension(NO_BOXES_XY)  :: chm
  real(RLEN),dimension(NO_BOXES_XY)  :: cmm
  real(RLEN),dimension(NO_BOXES_XY)  :: et,ex_sw,cx_any,rx_any,qx_any,px_any
  real(RLEN),dimension(NO_BOXES_XY)  :: r_xgrazing_c,ex_any
  real(RLEN),dimension(NO_BOXES_XY)  :: puhQ6,iN
  real(RLEN),dimension(NO_BOXES_XY)  :: eo,ea,epH
  real(RLEN),dimension(NO_BOXES_XY)  :: availQ6_c,availQ6_p,availQ6_n
  real(RLEN),dimension(NO_BOXES_XY)  :: suQ1
  real(RLEN),dimension(NO_BOXES_XY)  :: ruQ1c,ruQ2c,ruQ2n
  real(RLEN),dimension(NO_BOXES_XY)  :: puhQ1,puhQ2
  real(RLEN),dimension(NO_BOXES_XY)  :: ruQ1n,ruQun,ruQ1p,ruK3n
  real(RLEN),dimension(NO_BOXES_XY)  :: ruQ6c,ruQ6n,ruQ6p
  real(RLEN),dimension(NO_BOXES_XY)  :: rumn,rum4n,rum3n,rumun,rump,rumup,rum1p
  real(RLEN),dimension(NO_BOXES_XY)  :: rrmc,rrac,rrtc
  real(RLEN),dimension(NO_BOXES_XY)  :: sm
  real(RLEN),dimension(NO_BOXES_XY)  :: misn,misp
  real(RLEN),dimension(NO_BOXES_XY)  :: renn,renp
  real(RLEN),dimension(NO_BOXES_XY)  :: reR7c
  real(RLEN),dimension(NO_BOXES_XY)  :: jHIQ6c,jHIQ6n,jHIQ6p,jHIKIp
  real(RLEN),dimension(NO_BOXES_XY)  :: rmnuc,rmnun,rmnip,rmHc,rmHn,rmHp
  real(RLEN),dimension(NO_BOXES_XY)  :: qnQ1c,qpQ1c
  real(RLEN),dimension(NO_BOXES_XY)  :: qnHc,qpHc
  real(RLEN),dimension(NO_BOXES_XY)  :: qnQ6c,qpQ6c
  real(RLEN),dimension(NO_BOXES_XY)  :: rutc
  real(RLEN),dimension(NO_BOXES_XY)  :: rumc
  real(RLEN),dimension(NO_BOXES_XY)  :: rugc,rufc,runc
  real(RLEN),dimension(NO_BOXES_XY)  :: t,s,r,xeff
  real(RLEN),dimension(NO_BOXES_XY)  :: rupp,rupn
  real(RLEN),dimension(NO_BOXES_XY)  :: cQ1c,cQ2c,cQ2n,uQ1c
  real(RLEN),dimension(NO_BOXES_XY)  :: cK3n,cK4n,cK1p,cQun,cQ1p
  real(RLEN),dimension(NO_BOXES_XY)  :: to_pm3pw
  real(RLEN),dimension(NO_BOXES_XY)  :: M4n,M3n,Mun,M1p,Mup
  real(RLEN),dimension(NO_BOXES_XY)  :: jNIBIn,jNIBIp
  real(RLEN),dimension(NO_BOXES_XY)  :: limit_oxygen
  real(RLEN),dimension(NO_BOXES_XY)  :: cquN3,cquR1n,pxlayer2,pxlayer3

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  Copy  state var. object in local var
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  hxc = D2STATE(:,pphxc)
  hxn = D2STATE(:,pphxn)
  hxp = D2STATE(:,pphxp)

  qnHc=max(ZERO,hxn/(NZERO+hxc))
  qpHc=max(ZERO,hxp/(NZERO+hxc))
  call findlarge(qnHc,NO_BOXES_XY,4.D+00*p_qnc,iout)
!  if (iout.gt.0) then
!    write(LOGUNIT,*) 'BenBac quotum N/C too high' 
!    write(LOGUNIT,*) 'BenBac hx,qnHchxc,hxn,hxp',hx,qnHc,hxc,hxn,hxp 
!  endif

  jHIQ6c  = ZERO
  jHIQ6n  = ZERO
  jHIQ6p  = ZERO
  jHIKIp  = ZERO

  r_xgrazing_c=Source_D2_withgroup(pphxc, &
       ppBenOrganisms,iiBenOrganisms,iiC,iiConsumption)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Assign functional group-dependent parameters:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  select case ( hx)
    case ( iiH1 )
      clm  =   ZERO
      cm  =   D1m(:)
      chm  =   D1m(:)
      cmm  =   D1m(:)/ 2.0D+00
      to_pm3pw =   DONE/(chm-clm)/p_poro
      cK1p=K1p
      cK4n=K4n
      cK3n=max(ZERO,K3n-K13n-K23n)
      M1p =   cK1p/(p_p_ads_K1_ae+DONE)*to_pm3pw
      cQ2c=Q2c
      cQ2n=Q2n
    case ( iiH2 )
      if ( combine_anabac ) then
        clm  =   D1m(:)
        cm  =   D2m(:)
        chm  =   p_d_tot
        to_pm3pw=DONE/(chm-clm)/p_poro
        pxlayer2=PartQ_vector(  D6m(:),  clm,  cm,  p_d_tot)
        pxlayer3=PartQ_vector(  D6m(:),  cm,  chm,  p_d_tot)
        cx_any=pxlayer2+pxlayer3
        pxlayer2=pxlayer2/(NZERO+cx_any);pxlayer3=pxlayer3/(ZERO+cx_any)
        cK1p= (pxlayer2*max(ZERO,K11p)+pxlayer3*max(ZERO,K21p))
        cK4n=max(NZERO,pxlayer2*K14n) +max(NZERO,pxlayer3*K24n)
        cK3n=min(K3n,K13n+K23n)
        M1p=max(ZERO,pxlayer2*K11p/(DONE+p_p_ads_K1_ae) &
            +pxlayer3*K21p/(DONE+p_p_ads_K1_an)) *to_pm3pw
        cQ2c=ZERO;cQ2n=ZERO
      else
        clm  =   D1m(:)
        cm  =   (D2m(:)+D1m(:))*0.5
        chm  =   D2m(:)
        cmm  =   (D2m(:)+D1m(:))*0.5
        to_pm3pw=DONE/(chm-clm)/p_poro
        cK1p= K11p
        cK4n=K14n
        cK3n=K13n
        M1p=max(ZERO,cK1p/(DONE+p_p_ads_K1_ae)*to_pm3pw)
        cQ2c=ZERO;cQ2n=ZERO
      endif
    case ( iiH3 )
        if ( combine_anabac) STOP 'wrong use of combine_anabac'
        clm  =   D2m(:)
        cm  =   (D2m(:)+p_d_tot)*0.5
        chm  =   p_d_tot
        cmm  =   (D2m(:)+p_d_tot)*0.5
        to_pm3pw=DONE/(p_d_tot-clm)/p_poro
        cK1p= K21p
        cK4n=K24n
        cK3n=K23n
        M1p=max(ZERO,cK1p/(DONE+p_p_ads_K1_an)*to_pm3pw)
        cQ2c=ZERO;cQ2n=ZERO
  end select
  cQun=max(ZERO,BenUrea(p_iQ1(hx),iiN))
  Mun =   cQun*to_pm3pw
  uQ1c=BenLabileDetritus(p_iQ1(hx),iiC)
  cQ1c=max(ZERO,uQ1c-cQun/p_qnUc)
  M3n = max(ZERO,cK3n)*to_pm3pw
  M4n = cK4n/(p_p_ads_K4+DONE)*to_pm3pw
  cx_any=BenLabileDetritus(p_iQ1(hx),iiN)
  qnQ1c=max(ZERO,cx_any/(cQ1c+NZERO))
  cQ1p=BenLabileDetritus(p_iQ1(hx),iiP)
  qpQ1c=max(ZERO,cQ1p/(cQ1c+NZERO))
  cx_any=min(p_qpPhc*cQ1c,cQ1p(:))
  Mup(:) =  max(ZERO,min(cQ1p(:),cQ1c(:)*p_qpPhc)-cx_any)*to_pm3pw


  ! Determine where the median is of Q6 in the range from clm to D1m

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Physiological temperature response
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  et  =   eTq_vector(  ETW_Ben(:),  p_q10)
  select case (hx)
     case (iiH1)
          eo = DONE-exp(- G2o(:)/D1m/p_poro / p_clO2o)
          ea=p_xeff_an+(DONE-p_xeff_an)*eo
     case (iiH2)
          eo = DONE
          ea=p_xeff_an
  end select

#ifdef INCLUDE_BENCO2
  epH=DONE
  if (hx==iiH2) then
     epH=min(DONE,0.1*exp(-4.0D+00*(pHAn-9.0D+00)))
  endif
#endif


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate total food
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Detritus (if eaten): calculate available amount
  ! and add it to the total amount of food
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  availQ6_c  =   Q6c(:)* PartQ_vector(  D6m(:),  clm,  chm,  p_d_tot)
  availQ6_n  =   Q6n(:)* PartQ_vector(  D7m(:),  clm,  chm,  p_d_tot)
  availQ6_p  =   Q6p(:)* PartQ_vector(  D8m(:),  clm,  chm,  p_d_tot)

  qnQ6c  =   availQ6_n* p_cuR6np/( availQ6_c+ NZERO)
  qpQ6c  =   availQ6_p* p_cuR6np/( availQ6_c+ NZERO)
 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Growth is controlled by quality of detritus (N and P content):
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  xeff=DONE-p_pur
  iN=min((qnHc-p_qlnc)/(p_qnc-p_qlnc), &
                                 (qpHc-p_qlpc)/(p_qpc-p_qlpc))
  r=max(ZERO, iN)
  puhQ1=min(DONE,max(ZERO,DONE-min(abs(qnQ1c-xeff*p_qnc)/(xeff*p_qnc), &
               abs(qpQ1c-xeff*p_qpc)/(xeff*p_qpc)), r))
  puhQ2=max(ZERO,DONE-min(abs(qnHc-xeff*p_qnc)/(xeff*p_qnc), &
                            abs(qpHc-xeff*p_qpc)/(xeff*p_qpc)),iN)

  puhQ6=max(ZERO,DONE-max(abs(qnQ6c-xeff*p_qnc)/(xeff*p_qnc), &
               abs(qpQ6c-xeff*p_qpc)/(xeff*p_qpc)))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Total substrate availability:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ruQ6c  =max(ZERO, ea*et*availQ6_c* ( p_suhR6* puhQ6+ p_sulR6*(DONE-puhQ6)))

  suQ1= puhQ1*p_suhR1+(DONE-puhQ1)*p_sulR1
  ruQ1c=max(ZERO,ea*et*suQ1* cQ1c)
  ruQ2c=max(ZERO,ea*et*puhQ2*p_suR2* cQ2c)

  rutc  =   ruQ6c+ ruQ2c+ ruQ1c

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Potential uptake by bacteria
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  rumc  =   p_sum* et* epH* hxc

  ! Actual uptake by bacteria
  rugc  = min(rumc, rutc)
  ! Maintenance respiration
  rrmc  =   p_srr* hxc* et

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Limitation of uptake by shortage of oxygen
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  limit_oxygen=DONE
  if ( hx==iiH1) then
    rx_any= eo*(rugc* p_pur/MW_C +rrmc)
    call DoubleLimitChange_vector(POSITIVE,rx_any,G2_xavail_o, &
          fr_lim_HI_o(hx,:),max_change_per_step, limit_oxygen)
    rugc=rugc*((DONE-eo)+eo*limit_oxygen)
    eo=eo*limit_oxygen/((DONE-eo)+eo*limit_oxygen)
  endif

  rufc  =   rugc*(DONE-p_pur)-rrmc

  !if net-gain is less than the costs for the rest-respration the bacteria 
  !stop all activity.
  ex_sw=insw_vector(rufc)
  rufc=rufc*ex_sw
  rugc=rugc*ex_sw

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of respiration:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  rrac= rugc*p_pur
  rrtc=rrmc+rrac

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Carbon fluxes into bacteria
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ruQ6c = rugc* ruQ6c/ (NZERO + rutc )
  ruQ2c = rugc* ruQ2c/ (NZERO + rutc )
  ruQ1c = rugc* ruQ1c/ (NZERO + rutc )

  call LimitChange_vector(POSITIVE,ruQ1c,cQ1c,max_change_per_step)
  call LimitChange_vector(POSITIVE,ruQ2c,cQ2c,max_change_per_step)
  call flux_vector( iiBen,ppQ6c,pphxc, ruQ6c )
  call flux_vector( iiBen,ppQ2c,pphxc, ruQ2c )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Nutrient fluxes into bacteria from carbon fluxes
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ruQ6n  =   max(ZERO,ruQ6c* qnQ6c)
  ruQ6p  =   max(ZERO,ruQ6c* qpQ6c)
  ruQ1n  =   ruQ1c* qnQ1c
  ruQ1p  =   ruQ1c* qpQ1c
  ruQ2n  =   max(ZERO,ruQ2c*cQ2n/(NZERO+cQ2c) )

  call LimitChange_vector(POSITIVE,ruQ6n,Q6n,max_change_per_step)
  call LimitChange_vector(POSITIVE,ruQ6p,Q6p,max_change_per_step)
  call flux_vector( iiBen,ppQ6n,pphxn, ruQ6n )
  call flux_vector( iiBen,ppQ6p,pphxp, ruQ6p )
  call flux_vector( iiBen,ppQ2n,pphxn, ruQ2n )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of maximal nutrient uptake  on basis of  nutrient affinity
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  cquN3 =max(ZERO, p_lN3N4n/(NZERO+ p_lN3N4n+ M4n))
  cquR1n= p_lureaN4n/( NZERO + p_lureaN4n+ M4n)
  cK1p=max(ZERO,cK1p)
  cK4n=max(ZERO,cK4n)
  cK3n=max(ZERO,cK3n)

  t=fr_lim_HI_n(hx,:)
  rum4n = max(ZERO,p_qun       *hxc *M4n)
  call DoubleLimitChange_vector(POSITIVE,rum4n,cK4n,t,max_change_per_step)
  rum3n = max(ZERO,p_qun*cquN3*hxc *M3n)
  call DoubleLimitChange_vector(POSITIVE,rum3n,cK3n,t,max_change_per_step)
  rumun = max(ZERO,p_qun*cquR1n*hxc *Mun)
  call DoubleLimitChange_vector(POSITIVE,rumun,cQun,t,max_change_per_step,s)
  rx_any=rumun/p_qnUc
  call DoubleLimitChange_vector(POSITIVE,rx_any,uQ1c,t,max_change_per_step,t)
  rumun=rumun*min(s,t)
  t=fr_lim_HI_p(hx,:)
  rum1p =  max(ZERO,p_qup       *hxc *M1p)
  call DoubleLimitChange_vector(POSITIVE,rum1p,cK1p,t,max_change_per_step)
  rumup = Mup *et* p_qup *hxc
  call DoubleLimitChange_vector(POSITIVE,rumup,cQ1p,t,max_change_per_step,s)
  rx_any=rumup/p_qpPhc
  call DoubleLimitChange_vector(POSITIVE,rx_any,cQ1c,t,max_change_per_step,t)
  rumup=rumup*min(s,t)


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of potential nutrient uptake on basis of internal quotum and
  ! uptake
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  qx_any=qnHc-p_qlnc
  ex_sw=insw_vector(qx_any)
  !a too high nutrient content of food is related to the food upake     (r< 0)
  !a too low nutrient content of food is related to net growth          (r> 0)
  rx_any=rufc*p_qnc
  call LimitChange_vector(POSITIVE,rx_any,hxc(:)*qx_any,max_rate_per_step)
  rumn=rum4n+rum3n+rumun+rx_any*ex_sw
!--------------------------------------------------------------------
  if (isnan(rumn(1))) then
     write(LOGUNIT,*) 'rumn is Nan,rum4n,rum3n+rumun,rx_any,ex_sw', &
                rum4n,rum3n,rumun,rx_any,ex_sw
     write(LOGUNIT,*) 'hx,K4n,K3n,hxc,',hx,K4n,cK4n,K3n,hxc
     write(LOGUNIT,*) 'pxlayer2,pxlayer3,D6m,D1m,D2m,Q6c', &
                  pxlayer2,pxlayer3,D6m,D1m,D2m,Q6c,Q6n,Q6p
  endif
!--------------------------------------------------------------------

  qx_any=qpHc-p_qlpc
  ex_sw=insw_vector(qx_any)
  rx_any=rufc*p_qpc
  call LimitChange_vector(POSITIVE,rx_any,hxc(:)*qx_any,max_rate_per_step)
  rump=rum1p+rumup+rx_any*ex_sw

  where (rufc>ZERO)
    rx_any  =   min(  rufc,  ( rumn+ ruQ6n+ ruQ1n )/ p_qlnc)
    rx_any  =   min(  rx_any,( rump +ruQ6p+ ruQ1p )/ p_qlpc)
    reR7c  =  rufc -  max(  rx_any,  ZERO)

    runc  =   rufc- reR7c
  elsewhere
    runc=ZERO
    reR7c=ZERO;rx_any=ZERO
  endwhere
! call LimitChange_vector(POSITIVE,rx_any,cx_any, max_change_per_step)
  jnetHIc(hx,:)=runc

  rupn=max(ZERO,runc*p_qnc-r_xgrazing_c*qnHc)
  cx_any=p_qnc* hxc(:)-hxn(:)
  misn= cx_any*max(rrmc,rugc)/(NZERO+hxc(:)) 
  renn =max(ruQ6n+ruQ1n-rupn-misn,-rumn)

  cx_any=p_qpc* hxc(:)-hxp(:)
  rupp=max(ZERO,runc*p_qpc-r_xgrazing_c*qpHc)
  misp = cx_any*max(rrmc,rugc)/(NZERO+hxc(:)) 
  renp = max(ruQ6p+ruQ1p-rupp-misp,-rump)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Carbon correction: all C which cannot be used for growth due to
  ! lack of nutrients is excreted!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of mortality
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  sm =p_thdo/ ( max(ZERO,iN)+ p_thdo)* p_sd

  r=max(ZERO,hxc/(chm-clm)/p_poro)

  sm  = sm  + p_sd2*r

  rmHc=sm*hxc
  call LimitChange_vector(POSITIVE,rmHc,hxc,max_change_per_step)

  !extra check in case of extreme NC-quotum in BenBac
  rmnun=rmHc*max(ZERO,qnHc-p_qlnc)
  rmnuc=rmnun/p_qNuc
  rmHc=max(ZERO,rmHc-rmnuc)
  rmnun=min(rmnun,rmhc*p_qNuc)

  rmnuc=rmnun/p_qNuc
  rmHc=max(ZERO,rmHc-rmnuc)
  rmnip=rmHc*max(ZERO,qpHc-p_qlpc)
  rmHn=max(ZERO,rmHc*qnHc-rmnun)
  rmHp=max(ZERO,rmHc*qpHc-rmnip)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! add excretion of Detritus and Phosphate due to mortality to collectors vars.
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  jHIQ6c=jHIQ6c+ rmHc*( DONE- p_peZ_R1c)
  jHIQ6n=jHIQ6n+ rmHn*( DONE- p_peZ_R1n)
  jHIQ6p=jHIQ6p+ rmHp*( DONE- p_peZ_R1p)
  jHIKIp=jHIKIp+rmnip

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! fluxes excretion  of LOC and urea due to mortality
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call flux_vector(iiBen,pphxc,ppBenLabileDetritus(p_iQ1(hx),iiC), &
                                                rmHc *p_peZ_R1c+rmnuc)
  call flux_vector(iiBen,pphxn,ppBenLabileDetritus(p_iQ1(hx),iiN), &
                                                rmHn *p_peZ_R1n)
  call flux_vector(iiBen,pphxp,ppBenLabileDetritus(p_iQ1(hx),iiP), &
                                                rmHp *p_peZ_R1p)
  call flux_vector(iiBen,pphxn, ppBenUrea(p_iQ1(hx),iiN),rmnun)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calulate net Fluxes  N
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ex_sw=insw_vector(renn);px_any=DONE
  !-----excretion of nutrients
  if ( hx==iiH2) then
    px_any=min(DONE,pxlayer2*max(ZERO,K14n)/cK4n)
    call flux_vector(iiBen, pphxn,ppK24n, renn*ex_sw*(DONE-px_any))
  endif
  call flux_vector(iiBen, pphxn,ppBenthicAmmonium(p_iK4(hx),iiN),  &
                                                 renn*ex_sw*px_any)
  !-----nutrient uptake
  jNIBIn=-renn*rum4n/(NZERO+ rumn)*(DONE-ex_sw)*insw_vector(cK4n)
  ruK3n =-renn*rum3n/(NZERO+ rumn)*(DONE-ex_sw)
  call flux_vector(iiBen, ppK3n,pphxn,ruK3n)
  ruQun=-renn*rumun/(NZERO+ rumn)*(DONE-ex_sw)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Fluxes  P inlcude P-flux due to mortality
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ex_sw=insw_vector(renp);px_any=DONE
  jHIKIp=jHIKIp+ renp*ex_sw
  !-----excretion of nutrients
  if ( hx==iiH2) then
    px_any=min(ZERO,pxlayer2*max(ZERO,K11p)/(DONE+p_p_ads_K1_ae)/ cK1p)
    call flux_vector(iiBen, pphxp,ppK21p, jHIKIp*(DONE-px_any))
  endif
  call flux_vector(iiBen,pphxp,ppBenthicPhosphate(p_iK1(hx),iiP),jHIKIp*px_any)
  !-----nutrient uptake
  jNIBIp=-renp*rum1p/(NZERO+ rump)*(DONE-ex_sw)
  !organic P
  rx_any=-renp*rumup/(NZERO+ rump)*(DONE-ex_sw)
  ruQ1p=ruQ1p+rx_any

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! fluxes uptake  of LOC (food) and urea (nutrient)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  call flux_vector(iiBen, ppBenLabileDetritus(p_iQ1(hx),iiC),pphxc, & 
                                                       ruQ1c+ruQun/p_qnUc )
  call flux_vector(iiBen, ppBenLabileDetritus(p_iQ1(hx),iiN),pphxn, ruQ1n )
  call flux_vector(iiBen, ppBenLabileDetritus(p_iQ1(hx),iiP),pphxp, ruQ1p )
  call flux_vector(iiBen, ppBenUrea(p_iQ1(hx),iiN),          pphxn, ruQun)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Assigning depends on type of bacteria and layers in which they occur:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  select case (hx)

    case ( iiH1 )
      call sourcesink_flux_vector( iiBen, pphxc,ppG3c,rrtc+ruQun/p_qnUc )
      if (sw_an==0) then
        call flux_vector(iiBen, ppG2o,ppG2o,-( rrtc/ MW_C))
      else
        call flux_vector(iiBen, ppG2o,ppG2o,-(eo* rrtc/ MW_C))
        call flux_vector(iiBen, ppK6r,ppK6r,((DONE-eo)* rrtc/ MW_C*p_qro))
      endif
      call flux_vector(iiBen,ppBenthicAmmonium(p_iK4(hx),iiN), pphxn, jNIBIn)
      call flux_vector(iiBen,ppBenthicPhosphate(p_iK1(hx),iiP),pphxp, jNIBIp )

    case ( iiH2 )
      call sourcesink_flux_vector( iiBen, pphxc,ppG13c,rrtc+ruQun/p_qnUc )

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Respiration in anoxic circumstances produces reduced material
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      if ( .not.combine_anabac) then
!       r=D2m-D1m
!       cmm=log(DONE-0.5*(DONE-exp(-r/DH2m)))*(-DH2m)
!       call GetPenetrationDepthQ1(1,KQ1c,LAYER2,r,Q11c,s,ierr)
!       if ( ierr>0) s=D6m
!       call flux_vector(iiBen, ppDH2m,ppDH2m,&
!                 (s-DH2m(:))*(ruQ1c)/(NZERO+hxc(:)))
!       call flux_vector(iiBen, ppDH2m,ppDH2m,&
!                 (D6m-DH2m(:))*(ruQ6c)/(NZERO+hxc(:)))
      else
        call flux_vector( iiBen, ppK16r,ppK16r, pxlayer2*rrtc/ MW_C* p_qro )
        call flux_vector( iiBen, ppK26r,ppK26r, pxlayer3*rrtc/ MW_C* p_qro )
        ! nutrients for growth are taken of two layers in which H2 is present
        ! according distribution of H2 and nutrients available.
        ! determine the nutrient concentration orginating of anoxic 
        ! layer available for nutrient uptake
!       cx_any=pxlayer3*max(ZERO,K24n)
        px_any=max(ZERO,pxlayer3*K24n/cK4n)
        ex_any=max(ZERO,pxlayer2*K14n/cK4n)
!       cx_any=(px_any+ex_any)
!       px_any=px_any/cx_any
!       ex_any=ex_any/cx_any
        ! calculate uptake which take place in denitrification layer
!       rx_any=jNIBIn*cx_any/(NZERO+cK4n)
        rx_any=jNIBIn*px_any
        ! if necessary limit flux
!       call LimitChange_vector(POSITIVE,rx_any,K24n,DONE/LocalDelta,s)
        !reset s if conc. in anoxic layer is negative (should not be happen)
!       s=s*insw_vector(cK4n-cx_any)
!       call flux_vector(iiBen,ppK24n,pphxn, rx_any*s)
!       call flux_vector(iiBen,ppBenthicAmmonium(p_iK4(hx),iiN),pphxn, &
!                                                         (jNIBIn-rx_any)*s)
        call findnan(rx_any,NO_BOXES_XY,iout)
        if (iout.gt.0) write(LOGUNIT,*) 'BB rx_any' ,rx_any 
        call findnan(jNIBIN,NO_BOXES_XY,iout)
        if (iout.gt.0) write(LOGUNIT,*) 'BB jNIBIn',jNIBIn 
        call flux_vector(iiBen,ppK24n,pphxn, rx_any)
        call flux_vector(iiBen,ppBenthicAmmonium(p_iK4(hx),iiN),pphxn, &
                                                          (jNIBIn-rx_any))

        ! according distribution of H2 and nutrients available.
        ! determine the nutrient conc. orginating of denitrification layer
        ! available for nutrient uptake and take in account different properties
        ! for adsorption in the denitrification layer and the anoxic layer.
        cx_any=max(ZERO,pxlayer2*max(ZERO,K11p)/(DONE+p_p_ads_K1_ae))
        ! calculate uptake which take place in anoxic layer :
        ! total flux -flux in denitrification layer
        ! if K21p negative is the term in the min function may be larger than 1.
        rx_any=jNIBIp*(DONE- min(DONE,cx_any/(NZERO+M1p/to_pm3pw)))
        !limit flux in the layer with lowest adsorption capacity and
        !make flux r 0 if conc. in anoxic layer is <0 (should not be happen)
        call LimitChange_vector(POSITIVE,rx_any,K21p,DONE/LocalDelta,s)
        t=insw_vector(cx_any);s=DONE*(DONE-t)+s*t
        call flux_vector(iiBen,ppK21p,pphxp,rx_any*s)
        call flux_vector(iiBen, ppBenthicPhosphate(p_iK1(hx),iiP),pphxp, &
                           (jNIBIp -rx_any)*s)
      endif

    case ( iiH3 )
      call sourcesink_flux_vector( iiBen, pphxc,ppG23c,rrtc+ruQun/p_qnUc )

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Respiration in anoxic circumstances produces reduced material
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      call flux_vector( iiBen, ppK26r,ppK26r, rrtc/ MW_C* p_qro )

      if ( .not.combine_anabac) then
!       r=p_d_tot-D2m
!       cmm=log(DONE-0.5*(DONE-exp(-r/DH3m)))*(-DH3m)
!       call GetPenetrationDepthQ1(2,KQ1c,LAYER3,r,Q21c,s,ierr)
!       if ( ierr>0) s=D6m
!       call flux_vector(iiBen, ppDH3m,ppDH3m,&
!                 (s-DH3m(:))*(ruQ1c)/(NZERO+hxc(:)))
!       call flux_vector(iiBen, ppDH3m,ppDH3m,&
!                 (D6m-DH3m(:))*(ruQ6c)/(NZERO+hxc(:)))
      endif
  end select

  if (p_sulR6 <p_suhR6 .and. p_cuR6np > DONE )  then
     jHIQ6c=jHIQ6c+reR7c
  else
     call flux_vector( iiBen, pphxc,pphxc, -reR7c )
  endif

  call flux_vector( iiBen, pphxc,ppQ6c, jHIQ6c )
  call flux_vector( iiBen, pphxn,ppQ6n, jHIQ6n )
  call flux_vector( iiBen, pphxp,ppQ6p, jHIQ6p )

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculation of changes due to uptake of detritus in distribution
  ! of state variables (Dx.m is a undetermined source)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  if (hx==iiH2) cmm=FindMidPointExp_vector(D6m,d1_from=D1m,p_to=p_d_tot)
  s=(cmm-abs(D6m)) *(jHIQ6c-ruQ6c)/(NZERO+Q6c(:))
  call LimitChange_vector(NEGATIVE,s,abs(D6m),max_change_per_step)
  call flux_vector(iiBen, ppD6m,ppD6m,s*sign(DONE,D6m))
  if (hx==iiH2) cmm=FindMidPointExp_vector(D7m,d1_from=D1m,p_to=p_d_tot)
  s=(cmm-abs(D7m)) *(jHIQ6n-ruQ6n)/(NZERO+Q6n(:))
  call LimitChange_vector(NEGATIVE,s,abs(D7m),max_change_per_step)
  call flux_vector(iiBen, ppD7m,ppD7m,s*sign(DONE,D7m))
  if (hx==iiH2) cmm=FindMidPointExp_vector(D8m,d1_from=D1m,p_to=p_d_tot)
  s=(cmm-abs(D8m)) *(jHIQ6p-ruQ6p)/(NZERO+Q6p(:))
  call LimitChange_vector(NEGATIVE,s,abs(D8m),max_change_per_step)
  call flux_vector(iiBen, ppD8m,ppD8m,s*sign(DONE,D8m))


  end

!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
