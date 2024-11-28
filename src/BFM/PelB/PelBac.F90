#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: PelBac
!
! DESCRIPTION
!   This process describes the dynamics of bacterioplankton
!
!
!

!   This file is generated directly from OpenSesame model code, using a code
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine PelBacDynamics
!
! !USES:

  ! For the following Pelagic-states fluxes are defined: B1c, R6c, B1n, R6n, &
  ! B1p, R6p, R1c, R1n, R1p, R2c, O2o, N6r, N4n, N1p, N3n, R7c
  ! The following Pelagic 1-d global boxvars are modified : flPTN6r
  ! The following Pelagic 1-d global boxvars are used: ETW, qnB1c, qpB1c, &
  ! eO2mO2, qpR6c, qnR6c
  ! The following 0-d global parameters are used: p_peZ_R1c, p_peZ_R1n, &
  ! p_peZ_R1p, p_qro
  ! The following global constants are used: RLEN
  ! The following constants are used: MW_C

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO,DONE,NZERO

#ifdef NOPOINTERS
  use mem,  ONLY: D3STATE
#else
  use mem,ONLY:B1c,B1n,B1p, R6c,RZc, R1c,R1n,R1p, R2c, O2o, N4n,N3n, N1p
#endif
  use mem, ONLY: ppR6c,ppR6n,ppR6p, ppB1c,ppB1n,ppB1p, ppO2o,ppO3c, &
    iiPel, ppR1c,ppR1n,ppR1p, ppR2c, ppN6r,ppN4n,ppN1p,ppN3n, &
    qnB1c, qpB1c, qpR6c, qnR6c, jnetB1c, flPTN6r,flB1N4n, jPLO3c, &
    sN4N3n,max_change_per_step,max_rate_per_step,Nun,Nup,p_xfree_R2, &
    flux_vector, sourcesink_flux_vector,ETW,HCO3, Depth,NO_BOXES, &
    fr_lim_BN_n,fr_lim_B1_n,fr_lim_B1_p,iiConsumption, &
    ppMicroZooPlankton,iiMicroZooPlankton,iiC,fl_xgrazing_B1c
  use mem,only: ppBac,Bac,sN4N3n,pNaNt,flR1O3c,flR1N4n
  use constants, ONLY: MW_C, p_qnUc,POSITIVE
  use mem_Param, ONLY:p_peZ_R1c,p_peZ_R1n,p_peZ_R1p, p_qro, p_qon_nitri,p_xeff_an
  use mem_PelBac
  use LimitRates, ONLY:LimitChange_vector,DoubleLimitChange_vector
  use global_interface,only:FindNaNInRates
  use SourceFunctions,only:Source_D3_withgroup

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following vector functions are used:eTq_vector, MM_vector, insw_vector
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  use mem_globalfun,   ONLY: eTq_vector, insw_vector

! use mem,  ONLY: Output2d_1,Output2d_2,Output2d_3, Output2d_4

!
!
! !AUTHORS
!   Original version by J.W. Baretta
!    Giovanni Coppini (UNIBO), Hanneke Baretta-Bekker, Marcello Vichi (INGV)
!    Piet Ruardij (NIOZ)
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

  real(RLEN),external  ::GetDelta

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer                         :: iout
  real(RLEN)                      :: xeff,r_xscalar
  real(RLEN),dimension(NO_BOXES)  :: Bhc,Bhn,qnhc
  real(RLEN),dimension(NO_BOXES)  :: rmBc,rmBn,rmBp
  real(RLEN),dimension(NO_BOXES)  :: rupn,rupp
  real(RLEN),dimension(NO_BOXES)  :: pqun3,pqup,pquR1n,et,eo,ea
  real(RLEN),dimension(NO_BOXES)  :: flB1N6r
  real(RLEN),dimension(NO_BOXES)  :: rumc,rrmc
  real(RLEN),dimension(NO_BOXES)  :: rR1c,rR1n,rR1p
  real(RLEN),dimension(NO_BOXES)  :: ruR1c,ruR1n,ruR1p
  real(RLEN),dimension(NO_BOXES)  :: ruR2c
  real(RLEN),dimension(NO_BOXES)  :: ruR6c,ruR6n,ruR6p
                                     !gross nutrient uptake heterotrophic bact.
  real(RLEN),dimension(NO_BOXES)  :: rumHn,rumHp,renHn,runHn,renHp,runHp
  real(RLEN),dimension(NO_BOXES)  :: rugHc,rufHc,runHc
  real(RLEN),dimension(NO_BOXES)  :: rumH4n,rum3n,rumHun,rumH1p,rumHup
  real(RLEN),dimension(NO_BOXES)  :: rumB4n,rumBun,rumB1p,rumBup
  real(RLEN),dimension(NO_BOXES)  :: rugBc,rumtn,rumtp,rrBc
  real(RLEN),dimension(NO_BOXES)  :: misn,misp,luxn,luxp
  real(RLEN),dimension(NO_BOXES)  :: renBn,renBp
  real(RLEN),dimension(NO_BOXES)  :: reR7c
  real(RLEN),dimension(NO_BOXES)  :: rutc
  real(RLEN),dimension(NO_BOXES)  :: runNc
  real(RLEN),dimension(NO_BOXES)  :: puR6
  real(RLEN),dimension(NO_BOXES)  :: puhR1,puhR2
  real(RLEN),dimension(NO_BOXES)  :: iN
  real(RLEN),dimension(NO_BOXES)  :: puhR6,ex_sw,ex_lim1,ex_lim2
  real(RLEN),dimension(NO_BOXES)  :: rx_any,cx_any,qx_any,sx_any,px_any
  real(RLEN),dimension(NO_BOXES)  :: qpR1c,qnR1c
  real(RLEN),dimension(NO_BOXES)  :: rR6c
  real(RLEN),dimension(NO_BOXES)  :: cR1c,cR1n

  real(RLEN),dimension(NO_BOXES)  :: rumNc,sumNac,sumNbc
  real(RLEN),dimension(NO_BOXES)  :: rugNc,sugNac,sugNbc
  real(RLEN),dimension(NO_BOXES)  :: rN4N3n,qnnc
  real(RLEN),dimension(NO_BOXES)  :: uNa,uNb
  real(RLEN),dimension(NO_BOXES)  :: rrNc,rrHc
  real(RLEN),dimension(NO_BOXES)  :: pNun,pquaR1n,pqubR1n
  real(RLEN),dimension(NO_BOXES)  :: sumaN4n,sumaNun
  real(RLEN),dimension(NO_BOXES)  :: sumbN4n,sumbNun
  real(RLEN),dimension(NO_BOXES)  :: sumaN1p,sumbN1p
  real(RLEN),dimension(NO_BOXES)  :: rumN4n,rumNun,rumN1p
  real(RLEN),dimension(NO_BOXES)  :: Ban,Bap,renNn,renNp

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  fl_xgrazing_B1c=Source_D3_withgroup(ppB1c, &
                ppMicroZooPlankton,iiMicroZooPlankton,iiC,iiConsumption)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Temperature effect on pelagic bacteria:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  if (p_version ==5) then
     Bhc=max(ZERO,B1c-Bac)
     cx_any=Bac*p_qnc
     ! distribute nutrients over the B1c and Bac
     ! in case of tool low N content set both on its minimum, which lead
     ! for both to the same mortality. 
     qnhc=max(p_qnlc,(B1n-cx_any)/(NZERO+Bhc))
     qnnc=max(p_qnlc,(B1n-Bhc*qnhc)/(NZERO+Bac))
  else
    Bhc=B1c
    qnhc=qnB1c
  endif

  et  =   eTq_vector(  ETW,  p_q10)

  !All these rates will contain fluxes which originates from both the
  ! heterotrophic and nitrifying bacteria.

  rmBc=ZERO;rmBn=ZERO;rmBp=ZERO

  rR1c=ZERO;rR1n=ZERO;rR1p=ZERO
  rR6c=ZERO

  rrBc=ZERO;rugBc=ZERO

  rumB4n=ZERO;rumBun=ZERO;rumB1p=ZERO;rumBup=ZERO

  if (p_version >=4 )  then
     !correct Labile Organic Carbon for  urea
     cR1c=R1c-Nun/p_qnUc;cR1n=R1n-Nun
  else
     cR1c=R1c;cR1n=R1n
  endif


  !nitrifires must always enough N,
  !calc nutrient-state for nitrifiers based only on P
  !Calculate  correct quotum for heterotrophic bacteria
  iN=max(ZERO,min((qnhc-p_qnlc)/(p_qnc-p_qnlc),(qpB1c-p_qplc)/(p_qpc-p_qplc)))

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !  Mortality:
  !   1. first order mortality: old definition
  !   2. density dependent mortality due to virus infection
  !
  !   It is assumed the mortality is distributed in the same way over
  !   LOC (R1) and detritus (R6) s for phytoplankton and microzooplankton.
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if (p_version ==5 )  then

     rx_any= (  p_thdo/ ( iN+ p_thdo)* p_sd+( p_sd2* Bhc))*Bhc
     call LimitChange_vector(POSITIVE,rx_any,Bhc,max_change_per_step)
     rmBc= rmBc+rx_any
     rmBn= rmBn+rx_any*qnhc
     rmBp=rmBp +rx_any*qpB1c
  else
     rmBc  =  ( p_sd* et+( p_sd2* B1c)* B1c)
     rmBn=     rmBc*qnB1c
     rmBp=     rmBp*qpB1c
  endif

  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Calculate quota in R1c
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    !labile nutrient carbon is fully copied to nutrient availability..
    !a measure  for uptake of R1c is the quotum avail nutr's (dissolved +LABILE)
    !to R1C
    qnR1c  =   cR1n/ (NZERO + cR1c)
    qpR1c  =   R1p/ (NZERO + cR1c)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Potential available nutrients for uptake by bacteria
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  r_xscalar=NZERO+p_lN3N4n
  pqun3 = r_xscalar/(r_xscalar+N4n)
  r_xscalar=NZERO+p_lureaN4n
  pquR1n = r_xscalar/(r_xscalar+N4n)
  r_xscalar=NZERO+p_lN1
  pqup  =  r_xscalar/(r_xscalar+ N1p)

   px_any=fr_lim_B1_n
   rumH4n  =  et*Bhc* p_qun* N4n
   call DoubleLimitChange_vector(POSITIVE,rumH4n,N4n,px_any,max_change_per_step)
   rum3n  =  et*Bhc* p_qun* N3n* pqun3
   call DoubleLimitChange_vector(POSITIVE,rum3n,N3n,px_any,max_change_per_step)
   rumHun  =  et*Bhc* p_qun*Nun* pquR1n
   call DoubleLimitChange_vector(POSITIVE,rumHun,Nun,px_any,max_change_per_step)
   px_any=fr_lim_B1_p
   rumH1p  =  max(ZERO,et*p_qup* N1p* B1c )
   call DoubleLimitChange_vector(POSITIVE,rumH1p,N1p,px_any,max_change_per_step)
   rumHup  =  max(ZERO,et*p_qup* Nup* B1c*pqup )
   call DoubleLimitChange_vector(POSITIVE,rumHup,Nup,px_any,max_change_per_step)

   rumB4n= rumB4n +rumH4n
   rumBun= rumBun +rumHun
   rumB1p= rumB1p +rumH1p
   rumBup= rumBup +rumHup

 ! call findnega(rumB4n,NO_BOXES,iout)
 ! if (iout.gt.0) then
 !   write(LOGUNIT,*) '1: rumB4n, rumH4n:', &
 !     rumB4n(iout),rumH4n(iout)
 ! endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Substrate availability
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ! correction of food avilabilities dependent on internal quota
    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    xeff=DONE-p_pu_ra
    if (p_sulR1 > ZERO ) then
      if (p_version==5) then
        puhR1=max(ZERO,DONE-min(abs(qnR1c-xeff*p_qnc)/(xeff*p_qnc), &
                          abs(qpR1c-xeff*p_qpc)/(xeff*p_qpc)), iN)
        puhR2=max(ZERO,iN)
      else
        puhR1=max(ZERO,DONE-min(abs(qnR1c-xeff*p_qnc)/(xeff*p_qnc), &
                          abs(qpR1c-xeff*p_qpc)/(xeff*p_qpc)))
      endif
    else
       puhR1=DONE
    endif

    puR6  =   DONE;puhR6=DONE
    if (p_sulR6 <p_suR6 .and. p_cuR6np > DONE ) then
      if (p_version==5) then
        puhR6=max(ZERO,DONE-max(abs(qnR6c*p_cuR6np-xeff*p_qnc)/(xeff*p_qnc), &
                          abs(qpR6c*p_cuR6np-xeff*p_qpc)/(xeff*p_qpc)))
      else
        puhR6=max(ZERO,DONE-min(abs(qnR6c*p_cuR6np-xeff*p_qnc)/(xeff*p_qnc), &
                          abs(qpR6c*p_cuR6np-xeff*p_qpc)/(xeff*p_qpc)))
      endif
      puR6=(p_suR6*puhR6 +p_sulR6*(DONE-puhR6))/p_suR6
    endif

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ! oxygen environment:
    ! To provide a faster switching between the two metabolic pathways the
    ! oxygen dependence eo has been changed from the standard
    !     eo = MM(O2.o, p_chdo) to the cubic one written below.
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    eo = DONE-exp(- max(NZERO,O2o(:)) / p_chdo)

    ! rest respiration
    rrmc  =   p_srs* Bhc* et

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Calculate amount for R1, R6, and R2 and total amount of substrate avilable
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    !effect of anaerobic circumstances on uptake rate of food
    ea=p_xeff_an+(DONE-p_xeff_an)*eo

    ruR1c  =   ea*et*(p_suhR1* puhR1 + p_sulR1*(DONE-puhR1)) * cR1c
    !Assume the faecelpellets in R6 are not eaten
    ruR6c  =   ea*et*p_suR6* puR6* max(ZERO,R6c-RZc)
    ruR2c  =   ea*et*p_suR2*puhR2*p_xfree_R2* R2c

    rutc  =   NZERO + ruR6c+ ruR2c+ ruR1c

  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Potential uptake by bacteria
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    rumc  =   eo* p_sum* et*  Bhc
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Actual uptake by bacteria
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    rugHc  =   min(  rumc,  rutc)
    rufHc=rugHc*(DONE-p_pu_ra)-rrmc
    ! if net-gain is less than the costs for the rest-respration the bacteria
    ! stop all activity.
    ex_sw=insw_vector(rufHc)
    rufHc=rufHc*ex_sw
    rugHc=rugHc*ex_sw

    rrHc  =  rrmc +  p_pu_ra* rugHc
    rufHc=rugHc-rrHc

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Carbon fluxes into bacteria
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ruR6c  =   rugHc* ruR6c/ rutc
    ruR2c  =   rugHc* ruR2c/ rutc
    ruR1c  =   rugHc* ruR1c/ rutc
    call LimitChange_vector(POSITIVE,ruR1c,cR1c,max_change_per_step)
    rugHc= ruR6c+ruR2c+ ruR1c

    call flux_vector( iiPel, ppR6c,ppB1c, ruR6c )
    call flux_vector( iiPel, ppR2c,ppB1c, ruR2c )

    !Collect...
    rR1c=rR1c+ruR1c
    rugBc=rugBc+rugHc

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Organic Nitrogen and Phosphrous uptake
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    ruR6n  =   qnR6c* ruR6c
    ruR6p  =   qpR6c* ruR6c
    if (p_sulR6 <p_suR6 .and. p_cuR6np > DONE ) then
       ruR6n=ruR6n*p_cuR6np
       ruR6p=ruR6p*p_cuR6np
    endif

    call flux_vector( iiPel, ppR6n,ppB1n, ruR6n )
    call flux_vector( iiPel, ppR6p,ppB1p, ruR6p )

    ruR1n  =   qnR1c* ruR1c
    ruR1p  =   qpR1c* ruR1c

    !Collect...
    rR1n=rR1n+ruR1n
    rR1p=rR1p+ruR1p

    rrBc=rrBc+rrHc

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !        Calculation of net growth corrected for nutrient limitation
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    !N cyclus 
    rx_any=max(ZERO,rufHc*(qnhc(:)-p_qnlc))
    !a too high nutrient content of food is related to the food upake  (r< 0)
    !a too low nutrient content of food is related to net growth       (r> 0)
    call LimitChange_vector(POSITIVE,rx_any,Bhn,max_rate_per_step)
    rumHn=rum3n+rumH4n+rumHun+rx_any

    !P cyclus  
    rx_any=max(ZERO,rufHc*(qpB1c(:)-p_qplc))
    call LimitChange_vector(POSITIVE,rx_any,Bhc*qpB1c,max_rate_per_step)
    rumHp=rumH1p +rumHup+rx_any

    rx_any  =   min(  rufHc,( rumHn+ ruR6n+ ruR1n )/ p_qnlc)
    rx_any  =   min( rx_any,( rumHp +ruR6p+ ruR1p )/ p_qplc)
    reR7c =max(ZERO, rufHc - rx_any)
    runHc  =max(ZERO, rufHc- reR7c)

  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Inorganic Nutrient uptake by heterotrophic bacteria
  !  a too high nutrient uptake rate is related to the food upake (r<0)
  !  a too low nutrient uptake rate  is related to nutreint limitation  (r>0)
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    !N cyclus 
    !estimate N needed for net.growth net prod-grazing-main.respiration
    rupn=max(ZERO,runHc*p_qnc-fl_xgrazing_B1c*qnhc)
    ! intracellular missing amount of N
    rx_any=max(ZERO,(max(rrHc,rugHc))*(p_qnc-qnhc))
    misn= max(ZERO, rx_any)
    luxn= max(ZERO,-rx_any)
    runHn=min(rumH4n+rumHun-luxn,rupn+ misn)
    renHn=max(ruR6n+ruR1n-rupn,-runHn)

    !P cyclus  for all types of bacteria
    rupp=max(ZERO,runHc*p_qpc-fl_xgrazing_B1c*qpB1c)
    ! intracellular missing amount of N
    rx_any=max(ZERO,(max(rrHc,rugHc))*(p_qpc-qpB1c))
    misp= max(ZERO, rx_any)
    luxp= max(ZERO,-rx_any)
    runHp=min(rumH1p+rumHup-luxp,rupp+ misp)
    renHp=max(ruR6p+ruR1p-rupp,-runHp)

  if ( p_version>=4) then
     px_any= &
       max(ZERO,min((qnnc-p_qnlc)/(p_qnc-p_qnlc),(qpB1c-p_qplc)/(p_qpc-p_qplc)))
!    px_any=max(ZERO,(qpB1c-p_qplc)/(p_qpc-p_qplc))
     rx_any= (  p_thdo/ ( px_any+ p_thdo)* p_sd+( p_sd2* Bac))*Bac
     call flux_vector( iiPel, ppBac,ppBac, -rx_any)
     rmBc=rmBc +rx_any
     rmBn=rmBn +rx_any*p_qnc
     rmBp=rmBp +rx_any*qpB1c

     !max growth rate
     uNa=pNaNt
     uNb=(DONE-pNaNt)
     sumNac=eo* et *p_suNa*pNaNt
     sumNbc=eo* et *p_suNb*(DONE-pNaNt)
     r_xscalar=NZERO+p_lureaNbN4n
     pqubR1n = r_xscalar/(r_xscalar+N4n)
     r_xscalar=NZERO+p_lureaNaN4n
     pquaR1n = r_xscalar/(r_xscalar+N4n)

     !max nutrient uptake rate
     sumaNun=et*p_qunNa*Nun* pquaR1n*uNa
     sumbNun=et*p_qunNb*Nun* pqubR1n*uNb
     sumaN4n=et*p_qunNa*N4n*uNa
     sumbN4n=et*p_qunNb*N4n*uNb
     sumaN1p=et*p_qup*N1p*uNa
     sumbN1p=et*p_qup*N1p*uNb

     !limitation nutrient uptake ammonium
     rx_any=sumaN4n*Bac*uNa
     call DoubleLimitChange_vector(POSITIVE,rx_any,N4n,fr_lim_BN_n, &
                                        max_change_per_step,ex_lim1)
     sumaN4n=sumaN4n*ex_lim1

     rx_any=sumaNun*Bac*uNa
     call DoubleLimitChange_vector(POSITIVE,rx_any,Nun,fr_lim_BN_n, &
                                        max_change_per_step,ex_lim1)
     sumaNun=sumaNun*ex_lim1

     rx_any=sumaN1p*Bac*uNa
     call DoubleLimitChange_vector(POSITIVE,rx_any,N1p,fr_lim_B1_p, &
                                        max_change_per_step,ex_lim1)
     sumaN1p=sumaN1p*ex_lim1

     ! total gross specific nutrient uptake archea
     sx_any=sumaN4n+sumaNun
     ! determine limitation between C and N uptake
     px_any=min(sx_any,sumNac*p_qnNc)/(NZERO+sx_any)
     !correct ammonum and urea uptake
     sumaN4n=sumaN4n*px_any
     sumaNun=sumaNun*px_any

     rx_any=sumbN4n*Bac*uNb
     call DoubleLimitChange_vector(POSITIVE,rx_any,N4n,fr_lim_BN_n, &
                                        max_change_per_step,ex_lim1)
     sumbN4n=sumbN4n*ex_lim1

     rx_any=sumbNun*Bac*uNb
     call DoubleLimitChange_vector(POSITIVE,rx_any,Nun,fr_lim_BN_n, &
                                        max_change_per_step,ex_lim1)
     sumbNun=sumbNun*ex_lim1

     rx_any=sumbN1p*Bac*uNb
     call DoubleLimitChange_vector(POSITIVE,rx_any,N1p,fr_lim_B1_p, &
                                        max_change_per_step,ex_lim1)
     sumbN1p=sumbN1p*ex_lim1

     ! total gross specific nutrient uptake bacteria
     sx_any=sumbN4n+sumaNun
     ! determine limitation between C and N uptake
     px_any=min(sx_any,sumNbc*p_qnNc)/(NZERO+sx_any)
     !correct ammonum and urea uptake
     sumbN4n=sumbN4n*px_any
     sumbNun=sumbNun*px_any

     ! uptake rates for ammonium and urea
     rumN4n=(sumbN4n*uNb+sumaN4n*uNa)*Bac
     rumNun=(sumbNun*uNb+sumaNun*uNa)*Bac
     rumN1p=(sumbN1p*uNb+sumaN1p*uNa)*Bac

     !r_xscalar is a quotum...
     ! correct for the fact that the nitrifiers need some nitrogen to grow
     ! which cannot be used for nitrification
     r_xscalar=DONE/(p_qnNc*(DONE+p_qnlc/p_qnNc))
     !realized uptake
     sugNac=min(sumNac,(sumaN4n+sumaNun)*r_xscalar)
     sugNbc=min(sumNbc,(sumbN4n+sumbNun)*r_xscalar)


     !nitrifier which has no netto production stops...
     uNa=uNa*insw_vector(sugNac*uNa*(DONE-p_pu_raN)-et*p_srsN)
     uNb=uNb*insw_vector(sugNbc*uNb*(DONE-p_pu_raN)-et*p_srsN)
     rugNc=Bac*(sugNac*uNa+sugNbc*uNb)

     runNc=rugNc-et*p_srsN *(uNa+uNb)*Bac
     qx_any=(qpB1c(:)-p_qplc)
     ex_sw=insw_vector(qx_any)
     px_any  = min( runNc,(rumN1p +runNc*ex_sw*qpB1c)/ p_qplc)/runNc

     !nitrifier is too limited by phosphate stops...
     uNa=uNa*insw_vector(px_any*sugNac*uNa*(DONE-p_pu_raN)-et*p_srsN)
     uNb=uNb*insw_vector(px_any*sugNbc*uNb*(DONE-p_pu_raN)-et*p_srsN)
     rugNc=Bac*(sugNac*uNa+sugNbc*uNb)
     rumNc=Bac*(sumNac*uNa+sumNbc*uNb)

     sx_any=NZERO+uNa*(sumaN4n+sumaNun)+uNb*(sumbN4n+sumbNun)

     !part of uptake od dissolved which orginates from urea
     pNun=(uNa*sumaNun+uNb*sumbNun)/ sx_any

     rN4N3n=rugNc*p_qnNc

#ifdef INCLUDE_PELCO2
      !limitation
       rx_any=rugNc/MW_C
       call LimitChange_vector(POSITIVE,rx_any,HCO3,max_change_per_step,ex_lim2)
#endif
      rrNc =rugNc *p_pu_raN + et * p_srsN * Bac*(uNa+uNb)
      rx_any=rN4N3n*p_qon_nitri+rrNc/MW_C
      call LimitChange_vector(POSITIVE,rx_any,O2o,max_change_per_step,ex_lim1)
      rN4N3n=rN4N3n*min(ex_lim1,ex_lim2)
      rugNc=rugNc*min(ex_lim1,ex_lim2)
      rumNc=rumNc*min(ex_lim1,ex_lim2)
      !In case of negative growth: we assume direct cystis formation
      rrNc =rugNc *p_pu_raN + et * p_srsN * Bac* (uNa+uNb)
      !Primary porduction:uptake od carbon by nitrifiers
      runNc=max(ZERO,rugNc-rrNc)
      ex_sw=insw_vector(runNc)
      rugNc=rugNc*ex_sw
      rrNc =rugNc *p_pu_raN + et * p_srsN * Bac*(uNa+uNb)
      call flux_vector(iiPel, ppBac,ppBac,  rugNc)
      call flux_vector(iiPel, ppO3c,ppB1c,  rugNc)
      call flux_vector( iiPel,ppO2o,ppO2o, rugNc/MW_C)
      call flux_vector(iiPel, ppBac,ppBac,-rrNc)
      rrBc=rrBc+ rrNc
      rugBc=rugBc+ rugNc

      rupn=p_qnc*runNc
      rupp=p_qpc*runNc

      sN4N3n=rN4N3n/(NZERO+ N4n)
      rx_any=pNun*rugNc*p_qnNc
      !Collect rates: (See PelGlobal and PelChem)
      flR1N4n=flR1N4n+ rx_any
      flR1O3c=flR1O3c+ rx_any/p_qnUc

      renNn=max(-rupn,-(rumN4n+rumNun-rN4N3n))
      renNp=max(-rupp,-rumN1p)
      !subtract N used for nitrification.
      ! The left over can be used for filling the quota.
      rN4N3n=min(rumN4n+rumNun,rN4N3n)
      rumB4n=rumB4n + max(ZERO,rumN4n-rN4N3n*(DONE-pNun))
      rumBun=rumBun + max(ZERO,rumNun-rN4N3n*pNun)
      rumB1p=rumB1p +rumN1p
  else  
    rugNc=ZERO;rrNc=ZERO;renNp=ZERO;renNn=ZERO
  endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Respiration calculation + flux
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  call sourcesink_flux_vector( iiPel, ppB1c,ppO3c, rrBc )
  jPLO3c(1)=jPLO3c(1)+sum(rrBc*Depth)

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Pelagic bacteria are a wide functional group comprising both aerobic and
  ! anaerobic bacteria. At (very) low Oxygen concentrations bacteria use
  ! N6.r as electron acceptor in the respiration process. However, if N3.n
  ! is still present, the rate of consumption of N6.r is converted in N3.n
  ! consumption (see ChemicalProcesses.p).
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  call flux_vector( iiPel, ppO2o,ppO2o,-rrNc-( eo* rrHc/ MW_C))
  flB1N6r  =  ( DONE- eo)* rrHc/ MW_C* p_qro
  call flux_vector( iiPel, ppN6r,ppN6r, flB1N6r )
  flPTN6r  =   flPTN6r+ flB1N6r

  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Dissolved Nitrogen dynamics
  ! insw: No excretion if net. growth <0
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! add ren and  for nitrifiers
  renBn=renHn+renNn
  renBp=renHp+renNp

  rumtn=rumB4n+rumBun+rum3n
  rumtp=rumB1p+rumB1p

  ! excess of nutrients : renBn > 0
  ex_sw=insw_vector(renBn)
  call flux_vector( iiPel, ppB1n,ppN4n,  renBn*ex_sw )
  flB1N4n=flB1N4n+renBn*ex_sw

  ! shortage of nutrients : renBn < 0 --> Nutrient uptake
  rx_any=-renBn*(DONE-ex_sw)

  !ammonium  ------------------
  px_any=rumB4n/(NZERO+ rumtn)
  call findnega(px_any,NO_BOXES,iout)
  if (iout.gt.0) then
    write(LOGUNIT,*) 'rx_any,renBn,rumHn,rumtn,rumH4n,rumN4n,uN:', &
      rx_any(iout),renBn(iout),rumHn(iout),rumtn(iout), &
       rumH4n(iout),rumN4n(iout),rN4N3n(iout)*(DONE-pNun(iout)),uNa(iout)
    write(LOGUNIT,*) 'Nun,N4n,Bac,B1c:',N4n(iout),Nun(iout),Bac(iout),B1c(iout)
    write(LOGUNIT,*) 'rumB4n',rumB4n
  endif
  call flux_vector(iiPel, ppN4n,ppB1n, rx_any* px_any )
  !nitrate  ------------------
  px_any= rum3n/(NZERO+ rumtn)
  call flux_vector(iiPel, ppN3n,ppB1n, rx_any*px_any )
  !oxygen which freed by bacteria when they take up nitrate
  call flux_vector(iiPel, ppO2o,ppO2o,rx_any*px_any *p_qon_nitri)
  ! uptake of urea by bacteria and nitrifiers
  !Collect all R1'rates
  px_any=rumBun/(NZERO +rumtn)
  rR1n=rR1n              +rx_any *px_any
  ! CO2 produced by heterotrophic bacteria by urease
  flR1O3c=flR1O3c        +rx_any *px_any/p_qnUc
  call findnega(rR1n,NO_BOXES,iout)

  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Dissolved Phosphorus dynamics
  ! insw: No excretion if net. growth <0
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ! excess of nutrients : renBp > 0
  ex_sw=insw_vector(renBp)
  call flux_vector( iiPel, ppB1p,ppN1p, renBp* ex_sw )

  ! shortage of nutrients : renBp < 0 --> Nutrient uptake
  rx_any= -renBp*(DONE-ex_sw)
  call flux_vector( iiPel, ppN1p,ppB1p, rx_any *rumB1p/(NZERO +rumtp) )
  rR1p=rR1p+                            rx_any *rumBup/(NZERO +rumtp)

  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Uptake fluxes of R1c
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   call flux_vector( iiPel, ppR1c,ppB1c, rR1c )
   call flux_vector( iiPel, ppR1n,ppB1n, rR1n )
   call flux_vector( iiPel, ppR1p,ppB1p, rR1p )
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Excretion fluxes + correction net prod.:
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  if (p_sulR6 <p_suR6 .and. p_cuR6np > DONE )  then
     rR6c=rR6c+reR7c
  else
!    old option:
!    call flux_vector( iiPel, ppB1c,ppR7c, reR7c )
  endif
  rR6c=rR6c+ rmBc*( DONE- p_peZ_R1c)
  call flux_vector( iiPel, ppB1c,ppR6c, rR6c )
  call flux_vector( iiPel, ppB1n,ppR6n, rmBn*( DONE- p_peZ_R1n) )
  call flux_vector( iiPel, ppB1p,ppR6p, rmBp*( DONE- p_peZ_R1p) )

  call flux_vector( iiPel, ppB1c,ppR1c, rmBc* p_peZ_R1c )
  call flux_vector( iiPel, ppB1n,ppR1n, rmBn* p_peZ_R1n )
  call flux_vector( iiPel, ppB1p,ppR1p, rmBp* p_peZ_R1p )

  rx_any=rugBc -rrBc+rrmc-reR7c
  jnetB1c(1)=jnetB1c(1)+sum(Depth*rx_any)

  end

!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
