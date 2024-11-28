#include "DEBUG.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Settling
!
! DESCRIPTION
!   This process describes the dynamics of sedimentation and
!    deposition of phytoplankton (P1, P2, P3, P4, ...) and detritus (R6)
!    in the benthic system.
!    A burial velocity is defined, which controls the magnitude
!    of the inflow rate of detritus from the pelagic form R6 to
!    the benthic form Q6.
!    The processes described here taks only place in the lowest
!    boxes in the Z-direction.
!
!

!   This file is generated directly from OpenSesame model code, using a code
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine SettlingDynamics
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,ZERO,NZERO,DONE
  use mem,ONLY: Dfm,Dcm,R2c,R2n,R6c,R6n,R6p,R6s,R9x, RZc, D3STATE
  use mem,ONLY: Q6c,Q6n,Q6p,Q6s, D1m, D6m,D7m,D8m,D9m
  use mem,ONLY: ppD6m,ppD7m,ppD8m,ppD9m
  use mem,ONLY: jbotQ6c,jbotQ6n,jbotQ6p,jbotQ6s ,jrESS
  use mem,ONLY: jbotR6c,jbotR6n,jbotR6p,jbotR6s
  use mem,ONLY: ppR2c,ppR2n,ppR6c,ppR6n,ppR6p,ppR6s,ppRZc,iiC,iiN,iiP,iiL,iiS, &
    ppDfm,ppDcm, ppPhytoPlankton, iiPhytoPlankton,PhytoPlankton, &
    ppBenPhyto,BenPhyto, &
    iiPel,iiBen,iiR2, ppQ1c,ppQ1n,ppQ1p, ppQ2c,ppQ2n,ppQ6c,ppQ6s,ppQ6n,ppQ6p,  &
    NO_BOXES_XY,BoxNumber,BoxNumberXY, Depth,OCDepth, dry_z,  &
    sediR6,sediRZ,sediPI,sediR9, jRIQIc,jRIQIs, &
    PelBoxAbove,CoupledtoBDc,LocalDelta, flux,max_change_per_step
  use mem_Param,ONLY:CalcPhytoPlankton,p_dry_ben, p_poro,&
                                     p_pe_R1c,p_pe_R1n,p_pe_R1p
  use mem_Settling
  use constants,only:ANY,NEGATIVE,POSITIVE,INTEGRAL,AVERAGE
  use botflux,only:addbotflux,openbotflux
  use mem_BenPhyto,only:CalculateBenPhyto,p_hBPc
  use mem_Phyto,only:p_iRi,p_qnRc,p_qpRc,p_qsRc,p_xqn,p_xqp,p_xqs,p_qchlc,p_qnR2c
  use mem_PelGlobal,only:p_rrPIm,p_rrR6m
  use mem_Param,  ONLY:p_d_tot,p_qR6cQ9x
  use Track, ONLY:define_track_bot_rate
  use LimitRates, ONLY:LimitChange
  use mem_globalfun,   ONLY: IntegralExpDist,SetExpDist
  use mem_PelGlobal,ONLY:qlR2P1
  use global_interface,only:RecalcPenetrationDepth,CorrectConcNearBed,FindMidPointExp
  use mem_globalfun,   ONLY: IntegralExpDist,SetExpDist
! use mem,  ONLY: Output2d_1,Output2d_2,Output2d_3, Output2d_4
!
!
! !AUTHORS
!   ERSEM-team
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
  ! Local Vectors used  of group vectors
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN), dimension(:), pointer  ::lcl_PhytoPlankton
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer     :: i,j,k,l
  real(RLEN)  :: thlayer,alpha
  real(RLEN)  :: max_change
  real(RLEN)  :: corr
  real(RLEN)  :: sedi,Pc,Px_m2_c
  real(RLEN)  :: ruQIc,ruQIn,ruQIp,ruQIl
  real(RLEN)  :: ruQ2c,ruQ2n,ruQ6c,ruQ6n,ruQ6p,ruQ6s
  real(RLEN)  :: ruQ6x,ruQ1x,ruQIx
  real(RLEN)  :: ex_dry,aRc,aRx,s_m,qx_1tol,s_R6m
  real(RLEN)  :: px_lim,rin
  real(RLEN)  :: qx_any,sx_any,mx_any,px_any !(any quotum,specific rate,depth,fraction)
  real(RLEN)  :: cx_any,rx_any !(any concentration,rate)
  real(RLEN)  :: cx_attop !Concentration of Q6* at depth 0 =( sediment-eaterinterface)
  real(RLEN)  :: height
  real(RLEN)  :: newDm,rout
  real(RLEN)  :: cx_R2producers_c,R6_insilt,px_R6minusRZ
  real(RLEN),dimension(NO_BOXES_XY) :: thben

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! user defined external functions
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN), external  :: GetDelta
  external              :: GetLayer
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! user defined external functions
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  Call CouplingBioSilt(2,thben);thben=-min(ZERO,thben);thben=thben*GetDelta()
  DO BoxNumberXY=1,NO_BOXES_XY
    BoxNumber=PelBoxAbove(BoxNumberXY)
       cx_R2producers_c=NZERO
      ! No sedimentation when there is nearly no water onn tidal flat
      ex_dry=DONE
      if (p_dry_ben) ex_dry=min(DONE,Depth(BoxNumber)*dry_z(1)/p_height)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! sinking of pelagic phytoplankton
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      max_change=0.05D+00/LocalDelta
      do i= 1,iiPhytoPlankton
      !determine total amount of diatoms
      ! (all phytoplankton which porudce TEP)
        lcl_PhytoPlankton => PhytoPlankton(i,iiC)
        if (p_iRi(i)==iiR2)  cx_R2producers_c= &
                   cx_R2producers_c+lcl_PhytoPlankton(BoxNumber)
      enddo
      do i= 1,iiPhytoPlankton
        sedi = sediPI(i,BoxNumber)
        ! In this do-loop only the phytplankton types are checked on sedimentation
        ! which are not coupled to a benthic phytoplankton.
        if ( CoupledtoBDc(i)==0 .and. sedi> ZERO.and.   &
              CalcPhytoPlankton(i).and.p_burvel_PI(i) > ZERO  ) then
          Px_m2_c=sum(lcl_PhytoPlankton*Depth)
          s_m=p_burvel_PI(i)*ex_dry*min(DONE,sedi/(NZERO+p_rrPIm(i)))
          call CorrectConcNearBed(Depth(BoxNumber),s_m,p_height,1.0D+20,corr)
          !C:-------------------------
          j=ppPhytoPlankton(i,iiC)
          Pc=lcl_PhytoPlankton(BoxNumber)
          ruQIc = corr* s_m* Pc
          call LimitChange(POSITIVE,ruQIc,Px_m2_c,max_change)
          ruQ1x = p_pe_R1c* ruQIc;ruQ6x = ruQIc- ruQ1x
          call addbotflux(POSITIVE,BoxNumberXY,iiPel,j,iiBen,ppQ1c,ruQ1x)
          call addbotflux(POSITIVE,BoxNumberXY,iiPel,j,iiBen,ppQ6c,ruQ6x)
          if (p_iRi(i)==iiR2) then
            ruQ2c=ruQIc*min(R2c(BoxNumber)/(NZERO+cx_R2producers_c),qlR2P1)
            call addbotflux(POSITIVE,BoxNumberXY,iiPel,ppR2c,iiBen,ppQ2c,ruQ2c)
            if (p_qnR2c > NZERO) then
              ruQ2n=ruQ2c* R2n(Boxnumber)/(NZERO+R2c(Boxnumber))
              call addbotflux(POSITIVE,BoxNumberXY,iiPel,ppR2n,iiBen,ppQ2n,ruQ2c)
            endif
          endif
          !N:-------------------------
          j=ppPhytoPlankton(i,iiN)
          lcl_PhytoPlankton => PhytoPlankton(i,iiN)
          qx_any=lcl_PhytoPlankton(BoxNumber)/(NZERO+Pc)
          qx_any=min(qx_any,p_xqn(i)*p_qnRc(i))
          ruQIn = ruQIc* qx_any
          ruQ1x = p_pe_R1n* ruQIn;ruQ6x = ruQIn- ruQ1x
          call addbotflux(POSITIVE,BoxNumberXY,iiPel,j,iiBen,ppQ1n,ruQ1x)
          call addbotflux(POSITIVE,BoxNumberXY,iiPel,j,iiBen,ppQ6n,ruQ6x)
          call define_track_bot_rate(iiPel,j,BoxNumber,BoxNumberXY,0,0,-ruQIn)
          !P:-------------------------
          j=ppPhytoPlankton(i,iiP)
          lcl_PhytoPlankton => PhytoPlankton(i,iiP)
          qx_any=lcl_PhytoPlankton(BoxNumber)/(NZERO+Pc)
          qx_any=min(qx_any,p_xqp(i)*p_qpRc(i))
          ruQIp = ruQIc* qx_any
          ruQ1x = p_pe_R1p* ruQIp;ruQ6x = ruQIp- ruQ1x
          call addbotflux(POSITIVE,BoxNumberXY,iiPel,j,iiBen,ppQ1p,ruQ1x)
          call addbotflux(POSITIVE,BoxNumberXY,iiPel,j,iiBen,ppQ6p,ruQ6x)
          call define_track_bot_rate(iiPel,j,BoxNumber,BoxNumberXY,0,0,-ruQIp)
          !L:-------------------------
          j=ppPhytoPlankton(i,iiL)
          lcl_PhytoPlankton => PhytoPlankton(i,iiL)
          qx_any=lcl_PhytoPlankton(BoxNumber)/(NZERO+Pc)
          ruQIl = ruQIc* qx_any
          call openbotflux(NEGATIVE,BoxNumberXY,iiPel,j,-ruQIl)
          !Si:-------------------------
          j=ppPhytoPlankton(i,iiS)
          if ( j> 0) then
           lcl_PhytoPlankton => PhytoPlankton(i,iiS)
           qx_any=lcl_PhytoPlankton(BoxNumber)/(NZERO+Pc)
           qx_any=min(qx_any,p_xqs(i)*p_qsRc(i))
           ruQ6s = ruQIc* qx_any
           call addbotflux(POSITIVE,BoxNumberXY,iiPel,j,iiBen,ppQ6s,ruQ6s)
           call define_track_bot_rate(iiPel,j,BoxNumber,BoxNumberXY,0,0,-ruQ6s)
          endif
        end if
      end do

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! settling of resuspended diatoms in the pelagic
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      do i= 1,iiPhytoPlankton
        sedi= sediPI(i,BoxNumber)
        k= CoupledtoBDc(i)
        if ( k.ne.0 .and. sedi.gt.ZERO.and. CalcPhytoPlankton(i) ) then
        ! always settling pelagic phytoplankton is coupled to a BenPhyto type.
        ! See further at ModuleBenPhyto.F90
          k=mod(k,iiPhytoPlankton)
          s_m=p_burvel_PI(i)*min(DONE,sedi/(NZERO+p_rrPIm(i)))
!         height=s*GetDelta()
          height=p_height
          call CorrectConcNearBed(Depth(BoxNumber)*dry_z(BoxNumberXY), &
                                                  sedi,height,1.0D+20,corr)
          !C:-----------------------------------------
          j=ppPhytoPlankton(i,iiC);l=ppBenPhyto(k,iiC)
          lcl_PhytoPlankton => PhytoPlankton(i,iiC)
          Px_m2_c=sum(lcl_PhytoPlankton*Depth)
          Pc = lcl_PhytoPlankton(BoxNumber)
          ! netto sedimentation rate
          ruQIc = max(ZERO,corr* s_m*Pc)
          call LimitChange(POSITIVE,ruQIc,Px_m2_c,max_change)
          if ( ruQIc.gt.1.0D-5 ) then
            call addbotflux(POSITIVE,BoxNumberXY,iiPel,j,iiBen,l,ruQIc)
            !--R2->Q2 dynmaics is coupled to sinking of Benthis diatom
            ! sinking=sinking benthic diatoms/(total diatoms*R2c
            if (p_iRi(i)==iiR2) then
             ruQ2c=ruQIc &
                 *min(R2c(BoxNumber)/(NZERO+cx_R2producers_c),qlR2P1)
             call addbotflux(POSITIVE,BoxNumberXY,iiPel,ppR2c,iiBen,ppQ2c,ruQ2c)
             if (p_qnR2c > NZERO) then
               ruQ2n=ruQ2c* R2n(Boxnumber)/(NZERO+R2c(Boxnumber))
               call addbotflux(POSITIVE,BoxNumberXY,iiPel,ppR2n,iiBen,ppQ2n,ruQ2n)
             endif
            endif
            !N:-----------------------------------------
            lcl_PhytoPlankton => PhytoPlankton(i,iiN)
            qx_any=(p_qnRc(i)*NZERO+lcl_PhytoPlankton(BoxNumber))/(NZERO+Pc)
            qx_any=min(qx_any,p_xqn(i)*p_qnRc(i))
            j=ppPhytoPlankton(i,iiN);l=ppBenPhyto(k,iiN)
            ruQIx= ruQIc * qx_any
            call addbotflux(POSITIVE,BoxNumberXY,iiPel,j,iiBen,l,ruQIx)
            !P:-----------------------------------------
            lcl_PhytoPlankton => PhytoPlankton(i,iiP)
            qx_any=(p_qpRc(i)*NZERO+lcl_PhytoPlankton(BoxNumber))/(NZERO+Pc)
            qx_any=min(qx_any,p_xqp(i)*p_qpRc(i))
            j=ppPhytoPlankton(i,iiP);l=ppBenPhyto(k,iiP)
            ruQIx= ruQIc * qx_any
            call addbotflux(POSITIVE,BoxNumberXY,iiPel,j,iiBen,l,ruQIx)
            !Chla:-----------------------------------------
            lcl_PhytoPlankton => PhytoPlankton(i,iiL)
            j=ppPhytoPlankton(i,iiL);l=ppBenPhyto(k,iiL)
            qx_any=(p_qchlc(i)*NZERO+lcl_PhytoPlankton(BoxNumber))/(NZERO+Pc)
            qx_any=min(qx_any,p_xqn(i)*p_qchlc(i))
            ruQIx= ruQIc * qx_any
            call addbotflux(POSITIVE,BoxNumberXY,iiPel,j,iiBen,l,ruQIx)
            !Si:-----------------------------------------
            j=ppPhytoPlankton(i,iiS)
            if ( j>0 ) then
              l=ppBenPhyto(k,iiS);lcl_PhytoPlankton => PhytoPlankton(i,iiS)
              qx_any=lcl_PhytoPlankton(BoxNumber)/(NZERO+Pc)
              qx_any=(p_qsRc(i)*NZERO+lcl_PhytoPlankton(BoxNumber))/(NZERO+Pc)
              qx_any=min(qx_any,p_xqs(i)*p_qsRc(i))
              ruQIx= ruQIc * qx_any
              call addbotflux(POSITIVE,BoxNumberXY,iiPel,j,iiBen,l,ruQIx)
            endif
            !calculate the thickness in whichthe amount ruQIc should fit
            !for the actual Dfm
            lcl_PhytoPlankton => BenPhyto(k,iiC)
            Px_m2_c=lcl_PhytoPlankton(BoxNumberXY)
            !specific rate
            sx_any=min(DONE,ruQIc/(NZERO+Px_m2_c))
            call LimitChange(ANY,sx_any,DONE,max_change)
            mx_any=max(Dfm(BoxNumberXY),thben(BoxNumberXY))
            !Calculate the layer thickness in which the sedimentated fraction fits.
            mx_any=CalculateBenPhyto(iiC,INTEGRAL,BoxNumberXY, &
              ZERO,mx_any,inverse=sx_any*GetDelta())
            !benthic diatoms will distributed over the thickneess in the
            ! in a layer with ththickness of the layer of the sedimeted silt.
            ! Or in case of high concentration diatoms the layer will be thicker
            ! than the silt layer.
            mx_any=NZERO+max(mx_any,thben(BoxNumberXY))

            !fraction
            !Check if the concentration is not higher that a maximum parameter
            px_any=max(DONE,ruQIc*GetDelta()/mx_any / &
                                    p_poro(BoxNumberXY)/ p_hBPc(k) )

            mx_any= mx_any +&
            CalculateBenPhyto(iiC,AVERAGE,BoxNumberXY,mx_any,p_d_tot)

            !Calculate in change in Dfm (distribution of C in bentphyto)
            rx_any=sx_any*(mx_any*px_any -Dfm(BoxNumberXY))
            call LimitChange(ANY,rx_any,Dfm(BoxNumberXY),max_change)
            call flux(BoxNumberXY,iiBen,ppDfm,ppDfm, rx_any)

            !Calculate in change in Dfm (distribution of C in bentphyto)
            lcl_PhytoPlankton => PhytoPlankton(i,iiL)
            ruQix=ruQIc*lcl_PhytoPlankton(BoxNumber)/(NZERO+Pc)
            lcl_PhytoPlankton => BenPhyto(k,iiL)
            sx_any=ruQix/(NZERO+lcl_PhytoPlankton(BoxNumberXY))
            rx_any=sx_any*(mx_any*px_any -Dcm(BoxNumberXY))
            call LimitChange(ANY,rx_any,Dcm(BoxNumberXY),max_change)
            call flux(BoxNumberXY,iiBen,ppDcm,ppDcm, rx_any)
          endif
        endif
      enddo

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! settling of faecel-pellets detritus (RZ)
      ! Near the sediment fecalpellets fall apart with rate of p_smRZ
      ! (due physical stress near the bed)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      call flux(BoxNumber, iiPel, ppRZc, ppRZc, -p_smRZ*RZc(BoxNumber))

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      !settling of detritus
      !               -passive sedimentation
      !               -coupled to silt (R9x) using parameter p_qR6cR9x
      !Calculate sinking rate corrected for RZ
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! determine the number of layers from which is silt sedimentated
      ! the layers 1:l are used to average the quaota of R6* (cnp,s)
      mx_any=sum(sediR9*Depth)/OCDepth(1)*GetDelta()
      call getlayer(mx_any,l)
      mx_any=NZERO+sum(Depth(1:l))

      ! determination of the passive sedimentation rate
      s_R6m=ex_dry*max(ZERO,min(sediR6(BoxNumber), &
        (sediR6(BoxNumber)*R6c(BoxNumber)-sediRZ(BoxNumber)*RZc(BoxNumber) &
          -sediR9(BoxNumberXY)*R9x(BoxNumberXY)*p_qR6cQ9x)/ &
          (NZERO+max(ZERO,R6c(BoxNumber) &
             -RZc(BoxNumber)-R9x(BoxNumberXY)*p_qR6cQ9x))))
      s_m=p_burvel_R6cnp*s_R6m/(NZERO+p_rrR6m)
      !fecellpellets must first fall apart before it can be "catched" by
      !the sediment, hence only "slow" detritus  is catched
      px_R6minusRZ=(R6c(BoxNumber)-RZc(BoxNumber))/(NZERO+R6c(BoxNumber))

      ! all detritus from 0:p_height will sedimentates
      height=p_height

      !correction  for power profile near sediment
      call CorrectConcNearBed(Depth(BoxNumber) &
                             *dry_z(BoxNumberXY),s_m,height,1.0D+20,corr)


      !Calculate sedimentation of C -detritus
      aRc=sum(Depth(1:l)*px_R6minusRZ*R6c(1:l))/mx_any
      aRx=sum(Depth(1:l)*R9x(1:l))/mx_any
      !determine part of R6 whic is coupled to silt
      R6_insilt=min(aRc,aRx*p_qR6cQ9x)
      !derive from silt sedimentation the sedimentation rate of fraction of R6
      rx_any=-min(ZERO,jrESS(BoxNumberXY)*p_qR6cQ9x)
      ruQ6c = max(jRIQIc(BoxNumberXY),corr * s_m *max(ZERO,aRc-R6_insilt)) &
                                                                     +rx_any

      !Precorrection to limit flux to the  netto sedimentation
      !( -resuspension of detritus)
      cx_any=sum(R6c*Depth)+jbotQ6c(BoxNumberXY)*GetDelta()
      call LimitChange(POSITIVE,ruQ6c,cx_any,max_change,px_lim)
      !Calculate fraction of precorrected netto sedimentation h
      !corrected for resuspension and (already added) pseudo-faeces flux
      px_any=min(DONE,max(ZERO, &
            (ruQ6c*px_lim-jRIQIc(BoxNumberXY))/(NZERO+ruQ6c)))

      if (aRc.gt.NZERO) then
        ruQ6c= px_any* ruQ6c
        !Calculate sedimentation of N -detritus
        cx_any=sum(Depth(1:l)*px_R6minusRZ*R6n(1:l))/mx_any
        qx_1tol=cx_any/aRc
        qx_any=R6n(BoxNumberXY)/(NZERO+R6c(BoxNumberXY))
        !units: mmolN/m3 =moln/mgC * (-,(mg Silt/m3)/(mgC/m-3) *(grC/gr Silt)
        R6_insilt=cx_any*min(DONE,aRx/aRc*p_qR6cQ9x)
        rx_any=-min(ZERO,jrESS(BoxNumberXY)*p_qR6cQ9x*qx_1tol)
        !units   (-) * (-) * (m/d) * (-) *(-,molN/m3-molN/m3) +mmol N/m2/d
        ruQ6n = corr*s_m*max(ZERO,cx_any-R6_insilt)+rx_any
        !Check
        cx_any=sum(R6n*Depth)+jbotQ6n(BoxNumberXY)*GetDelta()
        call LimitChange(ANY,ruQ6n,cx_any,max_change,px_lim)
        px_any=min(DONE,max(ZERO,  &
                 (ruQ6n*px_lim-jRIQIc(BoxNumberXY)*qx_any)/(NZERO+ruQ6n)))
        ruQ6n = ruQ6n *px_any

        !Calculate sedimentation of P -detritus
        cx_any=sum(Depth(1:l)*px_R6minusRZ*R6p(1:l))/mx_any
        qx_1tol=cx_any/aRc
        qx_any=R6p(BoxNumberXY)/(NZERO+R6c(BoxNumberXY))
        !units:mmol P/m3 =mol P/mg C *(-,(mg Silt/m3)/(mgC/m-3) *(grC/gr Silt)
        R6_insilt=cx_any*min(DONE,aRx/aRc*p_qR6cQ9x)
        rx_any=-min(ZERO,jrESS(BoxNumberXY)*p_qR6cQ9x*qx_1tol)
        ruQ6p = corr *s_m* max(ZERO,cx_any-R6_insilt)+rx_any
        !Check
        cx_any=sum(R6p*Depth)+jbotQ6p(BoxNumberXY)*GetDelta()
        call LimitChange(ANY,ruQ6p,cx_any,max_change,px_lim)
        px_any=min(DONE,max(ZERO, &
                  (ruQ6p*px_lim-jRIQIc(BoxNumberXY)*qx_any)/(NZERO+ruQ6p)))
        ruQ6p = ruQ6p *px_any

        !Calculate sedimentation of Si -detritus
        s_m=p_burvel_R6s*s_R6m/(NZERO+p_rrR6m)
        cx_any=sum(Depth(1:l)*px_R6minusRZ*R6s(1:l))/mx_any
        qx_1tol=cx_any/aRc
        rx_any=-min(ZERO,jrESS(BoxNumberXY)*p_qR6cQ9x*qx_1tol)
        R6_insilt=cx_any*min(DONE,aRx/aRc*p_qR6cQ9x)
        ruQ6s = max(jRIQIs(BoxNumberXY),corr *s_m*max(ZERO,cx_any-R6_insilt)) &
                                                                +rx_any

        !Check
        cx_any=sum(R6s*Depth)+jbotQ6s(BoxNumberXY)*GetDelta()
        call LimitChange(ANY,ruQ6s,cx_any,max_change,px_lim)
        px_any=min( &
             DONE,max(ZERO, (ruQ6s*px_lim-jRIQIs(BoxNumberXY))/(NZERO+ruQ6s)))
        ruQ6s = ruQ6s *px_any

!       Output2d_1=ruQ6c
!       Output2d_2=ruQ6n
!       Output2d_3=ruQ6p
!       Output2d_4=ruQ6s
        call addbotflux(POSITIVE,BoxNumberXY,iiPel,ppR6c,iiBen,ppQ6c,ruQ6c)
        call addbotflux(POSITIVE,BoxNumberXY,iiPel,ppR6n,iiBen,ppQ6n,ruQ6n)
        call addbotflux(POSITIVE,BoxNumberXY,iiPel,ppR6p,iiBen,ppQ6p,ruQ6p)
        call addbotflux(POSITIVE,BoxNumberXY,iiPel,ppR6s,iiBen,ppQ6s,ruQ6s)
      endif
      !prevent that RZc is larger then R6c in bottom layer..
      cx_any=R6c(BoxNumberXY)+ &
                    (jbotR6c(BoxNumberXY)+jbotQ6c(BoxNumberXY)) *GetDelta()
      rx_any=min(ZERO,cx_any-RZc(BoxNumberXY))/GetDelta()
      call openbotflux(NEGATIVE,BoxNumberXY,iiPel,ppRZc,rx_any)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      !Calculate redistribution of detritus
      !Redistribution calculation keep total detritus concentration below
      ! oxygen penetration depth (D1m) keep constant
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
       j=BoxNumberXY
       ! Als jrEss negtive : sedimentation of detritus.....
       rin= -(jbotR6c(j)+jbotQ6c(j))
       cx_any=rin*LocalDelta
       if ( rin.gt.NZERO) then
         rout= ( -abs(D6m(j))) *(rin)/(NZERO+Q6c(j))
         call RecalcPenetrationDepth(D1m(j),p_d_tot, D6m(j),p_hDQ6m, &
                                                         cx_any,Q6c(j),newDm)
         rx_any=(abs(newDm)- abs(D6m(j)))*sign(DONE,newDm)/LocalDelta
         ! r and rout both negative: the maximum give the smallest change
         rout=max(rout,rx_any)
       else
         alpha=SetExpDist(D6m(j),Dmm=p_d_tot)
         cx_attop=Q6c(j)/IntegralExpDist( -alpha,p_d_tot)
         thlayer=log(DONE+cx_any*alpha/(NZERO+cx_attop))/(-alpha)
         mx_any=FindMidPointExp(D6m(j),thlayer,p_d_tot)-thlayer
         rout= (mx_any -abs(D6m(j))) +max(ZERO,D6m(j)-p_d_tot)
         rout= rout*(rin)/(NZERO+Q6c(j))
       endif
       call flux(j,iiBen,ppD6m,ppD6m,rout*sign(DONE,D6m(j)))

       rin= -(jbotR6n(j)+jbotQ6n(j))
       cx_any=rin*LocalDelta
       if ( rin.gt.NZERO) then
         rout= ( -abs(D7m(j))) *(rin)/(NZERO+Q6n(j))
         call RecalcPenetrationDepth(D1m(j),p_d_tot, D7m(j),p_hDQ6m, &
                                                         cx_any,Q6n(j),newDm)
         rx_any=(abs(newDm)- abs(D7m(j)))*sign(DONE,newDm)/LocalDelta
         rout=max(rout,rx_any)
       else
         alpha=SetExpDist(D7m(j),Dmm=p_d_tot)
         cx_attop=Q6n(j)/IntegralExpDist( -alpha,p_d_tot)
         thlayer=log(DONE+cx_any*alpha/(NZERO+cx_attop))/(-alpha)
         mx_any=FindMidPointExp(D7m(j),thlayer,p_d_tot)-thlayer
         rout= (mx_any -abs(D7m(j))) +max(ZERO,D7m(j)-p_d_tot)
         rout= rout*(rin)/(NZERO+Q6n(j))
       endif
       call flux(j,iiBen,ppD7m,ppD7m,rout*sign(DONE,D7m(j)))

       rin= -(jbotR6p(j)+jbotQ6p(j))
       cx_any=rin*LocalDelta
       if ( rin.gt.NZERO) then
         rout= (  -abs(D8m(j))) *(rin)/(NZERO+Q6p(j))
         call RecalcPenetrationDepth(D1m(j),p_d_tot, D8m(j),p_hDQ6m, &
                                                         cx_any,Q6p(j),newDm)
         rx_any=(abs(newDm)- abs(D8m(j)))*sign(DONE,newDm)/LocalDelta
         rout=max(rout,rx_any)
       else
         alpha=SetExpDist(D8m(j),Dmm=p_d_tot)
         cx_attop=Q6p(j)/IntegralExpDist( -alpha,p_d_tot)
         thlayer=log(DONE+cx_any*alpha/(NZERO+cx_attop))/(-alpha)
         mx_any=FindMidPointExp(D8m(j),thlayer,p_d_tot)-thlayer
         rout= (mx_any -abs(D8m(j))) +max(ZERO,D8m(j)-p_d_tot)
         rout= rout*(rin)/(NZERO+Q6p(j))
       endif
       call flux(j,iiBen,ppD8m,ppD8m,rout*sign(DONE,D8m(j)))

       rin= -(jbotR6s(j)+jbotQ6s(j))
       cx_any=rin*LocalDelta
       if ( rin.gt.NZERO) then
         rout= ( -abs(D9m(j))) *(rin)/(NZERO+Q6s(j))
         call RecalcPenetrationDepth(D1m(j),p_d_tot, D9m(j),p_hDQ6m, &
                                                         cx_any,Q6s(j),newDm)
         rx_any=(abs(newDm)- abs(D9m(j)))*sign(DONE,newDm)/LocalDelta
         rout=max(rout,rx_any)
       else
         alpha=SetExpDist(D9m(j),Dmm=p_d_tot)
         cx_attop=Q6s(j)/IntegralExpDist( -alpha,p_d_tot)
         thlayer=log(DONE+cx_any*alpha/(NZERO+cx_attop))/(-alpha)
         mx_any=FindMidPointExp(D9m(j),thlayer,p_d_tot)-thlayer
         rout= (mx_any -abs(D9m(j))) +max(ZERO,D9m(j)-p_d_tot)
         rout= rout*(rin)/(NZERO+Q6s(j))
       endif
       call flux(j,iiBen,ppD9m,ppD9m,rout*sign(DONE,D9m(j)))
  enddo

  end
  subroutine getlayer(depth,layer)
  use global_mem, ONLY:RLEN,ZERO
  use mem,only:OCDepth,NO_BOXES
  implicit none
  REAL(RLEN),intent(IN)   :: depth
  integer,intent(OUT)     :: layer
  REAL(RLEN)              ::r
  integer                 :: nlist(NO_BOXES)

   r=max(ZERO,OCDepth(1)-depth)
   nlist=minloc(abs(OCdepth-r))
   layer=min(NO_BOXES,max(1,nlist(1)))
   return
  end

!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
