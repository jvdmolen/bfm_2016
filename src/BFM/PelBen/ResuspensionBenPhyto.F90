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
  subroutine ResuspensionBenPhyto
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      use global_mem, ONLY:RLEN,ZERO,LOGUNIT,NZERO,DONE
      use mem, ONLY: BoxNumberXY,  &
        NO_BOXES,NO_BOXES_Z,NO_BOXES_XY,Dfm,Dcm,ppDfm,ppDcm,Depth,  &
        flux,flux_vector,iiBenPhyto,PelBoxAbove,jbotQ2c,jbotQ2n, &
        iiBen,iiPel,iiC,iiN,iiP,iiL,iiS,ppBenPhyto,BenPhyto,ppPhytoPlankton,&
        iiPhytoPlankton, &
        CoupledtoBDc,jbotBPc,Q2c,Q2n,ppR2c,ppR2n,ppQ2c,ppQ2n, &
        max_change_per_step,LocalDelta
      use mem_Param,ONLY:p_d_tot
      use mem_BenPhyto,ONLY:p_useparams,p_xresus,&
                                             CalculateBenPhyto_vector

      use mem_globalfun,only:insw_vector
      use constants, ONLY: INTEGRAL,AVERAGE,ANY,POSITIVE
      use LimitRates, ONLY:LimitChange_vector
      use mem_pelGlobal,ONLY:qlR2P1
!     use mem,  ONLY: Output2d_1,Output2d_2,Output2d_3, Output2d_4
!
!
! !AUTHORS
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
       implicit none

       integer                                        :: bp
       integer                                        :: layer,ll
       integer                                        :: nrphyto
       REAL(RLEN),dimension(NO_BOXES_XY)              :: thlayer
       REAL(RLEN),dimension(NO_BOXES_XY)              :: dummy
       REAL(RLEN),dimension(NO_BOXES_XY)              :: r,rx_any,cx_any,c_m
       REAL(RLEN),dimension(NO_BOXES_XY)              :: limit_1,limit_2
       REAL(RLEN),dimension(NO_BOXES_XY)              :: pc,pl
       REAL(RLEN),dimension(NO_BOXES_XY)              :: bpc
       REAL(RLEN),dimension(NO_BOXES_XY)              :: bpn
       REAL(RLEN),dimension(NO_BOXES_XY)              :: bpp
       REAL(RLEN),dimension(NO_BOXES_XY)              :: bpl
       REAL(RLEN),dimension(NO_BOXES_XY)              :: bps
       REAL(RLEN)                 :: sc,sn,sp,sl
       REAL(RLEN)                 :: sqc,sqn
       REAL(RLEN)                 :: ss
       REAL(RLEN)                                     :: max_change

       !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
       ! user defined external functions
       !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
       integer, external  :: D3toD1
       integer, external  :: D2toD1
       !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

       max_change=0.05D+00/LocalDelta
       dummy=ZERO
       do bp=1,iiBenPhyto
         nrphyto=p_useparams(bp)
         jbotBPc(bp,:)=ZERO
         ll=CoupledtoBDc(nrphyto)
         ! only resuspension in case of more than a low biomass of BenPhyto.
         ! See further at ModuleBenPhyto.F90
         if (ll>0.and.ll<iiPhytoPlankton.and.p_xresus(bp)>ZERO) then
           ! rrESS is reuspension/d-->
           !               in one timestep rrESS*LocalDelta is resuspended.
           ! mg Silt/m2/d*d/(mg silt/m2)*1m
           call CouplingBioSilt(2,thlayer)
           bpc=max(ZERO,BenPhyto(bp,iiC))
           thlayer=max(ZERO,thlayer*p_xresus(bp)*LocalDelta) &
                                                    *bpc/(bpc+1000.0D+00)
           thlayer=thlayer*insw_vector(thlayer-0.00004)
           thlayer=max(ZERO,thlayer)      !JM added
           bpn=max(ZERO,BenPhyto(bp,iiN))
           bpp=max(ZERO,BenPhyto(bp,iiP))
           bpl=max(ZERO,BenPhyto(bp,iiL))
           bps=max(ZERO,BenPhyto(bp,iiS))
           c_m=ZERO
           ! calculate fraction of silt which is resuspended...
           pc=  CalculateBenPhyto_vector(iiC,INTEGRAL,c_m,thlayer)
           pl=  CalculateBenPhyto_vector(iiL,INTEGRAL,c_m,thlayer)
           !Calculate the relative concentration per m3:
!          ! correct for low concentrations  unit: 1/d
           ! reset rates at negative concentrations
           limit_1=insw_vector(bpc)*insw_vector(min(bpc,bpl))
           pc=pc/LocalDelta*limit_1
           pl=pl/LocalDelta*limit_1
           cx_any=DONE
           ! resuspension is limited at low concentrations
!JM           call LimitChange_vector(POSITIVE,pc,cx_any,max_change,limit_1)
!JM           call LimitChange_vector(POSITIVE,pl,cx_any,max_change,limit_2)
           call LimitChange_vector(ANY,pc,cx_any,max_change,limit_1)
           call LimitChange_vector(ANY,pl,cx_any,max_change,limit_2)
           limit_1=min(limit_1,limit_2);pl=pl*limit_1;pc=pc*limit_1
           DO BoxNumberXY=1,NO_BOXES_XY
             layer=PelBoxAbove(BoxNumberXY)
             if (thlayer(BoxNumberXY).gt.NZERO.and. &
                                     bpc(BoxNumberXY).gt.NZERO) then
               ! unit: 1/d mgSilt/m3/d * ( m2 *d /mgSilt) = 1/(m*d)
               ! unit: 1/(m *d) * mgC/m2 =mgC/m3/d
               sc=pc(BoxNumberXY)*bpc(BoxNumberXY)/Depth(layer)
               ll=ppPhytoPlankton(nrphyto,iiC);call flux(layer,iiPel,ll,ll,sc)
               sn=pl(BoxNumberXY)*bpn(BoxNumberXY)/Depth(layer)
               ll=ppPhytoPlankton(nrphyto,iiN);call flux(layer,iiPel,ll,ll,sn)
               sp=pl(BoxNumberXY)*bpp(BoxNumberXY)/Depth(layer)
               ll=ppPhytoPlankton(nrphyto,iiP);call flux(layer,iiPel,ll,ll,sp)
               sl=pl(BoxNumberXY)*bpl(BoxNumberXY)/Depth(layer)
               ll=ppPhytoPlankton(nrphyto,iiL);call flux(layer,iiPel,ll,ll,sl)
               ll=ppPhytoPlankton(nrphyto,iiS)
               if (ll.gt.0)  then
                 ss=pl(BoxNumberXY)*bps(BoxNumberXY)/Depth(layer)
                 call flux(layer,iiPel,ll,ll,ss)
               endif
               sqc=pc(BoxNumberXY)* min(Q2c(BoxNumberXY), &
                             bpc(BoxNumberXY)*qlR2P1)/Depth(layer)
               call flux(layer,iiPel,ppR2c,ppR2c,sqc)
               sqn=pl(BoxNumberXY)* min(Q2n(BoxNumberXY), &
                             bpc(BoxNumberXY)*qlR2P1)/Depth(layer)
               call flux(layer,iiPel,ppR2n,ppR2n,sqn)
             endif
           enddo
           if (sum(pc).gt.NZERO) then
           ! just as with the silt model  in case of a resuspension event
           ! the diatoms are directly distributed over the whole water column
             jbotBPc(bp,:)=pc*bpc
             ll=ppBenPhyto(bp,iiC);       call flux_vector(iiBen,ll,ll,-pc*bpc)
             ll=ppBenPhyto(bp,iiN);       call flux_vector(iiBen,ll,ll,-pl*bpn)
             ll=ppBenPhyto(bp,iiP);       call flux_vector(iiBen,ll,ll,-pl*bpp)
             ll=ppBenPhyto(bp,iiL);       call flux_vector(iiBen,ll,ll,-pl*bpl)
             ll=ppBenPhyto(bp,iiS)
             if (ll >0) call flux_vector(iiBen,ll,ll,-pl*bps)
             ! TEP produce by benthic diatoms is resuspendend too!
             cx_any=min(Q2c,bpc*qlR2P1)
             call flux_vector(iiBen,ppQ2c,ppQ2c,-pc*cx_any)
             jbotQ2c=pc*cx_any
             cx_any=min(Q2n,bpn*qlR2P1)
             call flux_vector(iiBen,ppQ2n,ppQ2n,-pl*cx_any)
             jbotQ2n=pc*cx_any
             c_m=p_d_tot
             !calculate depth where the median of the concentration is in
             ! the layer between from thlayer to Inf
             ! because a layer with thickness 'thelayer' is removed. The average
             ! will shifted in upwards direction.
             r=CalculateBenPhyto_vector(iiC,AVERAGE,thlayer,c_m) -thlayer
             rx_any=pc* (r-Dfm)
             call LimitChange_vector(ANY,rx_any,Dfm,max_change)
             call flux_vector(iiBen, ppDfm,ppDfm,rx_any)
             r=CalculateBenPhyto_vector(iiL,AVERAGE,thlayer,c_m) -thlayer
             rx_any=pl* (r-Dcm)
             call LimitChange_vector(ANY,r,Dcm,max_change)
             call flux_vector(iiBen, ppDcm,ppDcm,r)
           else
             jbotQ2c=ZERO
             jbotQ2n=ZERO
           endif
         endif
       enddo
       end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
