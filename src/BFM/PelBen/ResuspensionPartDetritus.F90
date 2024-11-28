#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: ResuspensionPartDetritus
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
  subroutine ResuspensionPartDetritus
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

       use global_mem, ONLY:RLEN,ZERO,LOGUNIT,NZERO,DONE
       use mem, ONLY: BoxNumberXY, NO_BOXES_XY,PelBoxAbove,  &
         Depth,flux,flux_vector,jrESS,max_change_per_step, &
         iiBen,iiPel,ppR6c,ppR6n,ppR6p,ppR6s, ppQ6c,ppQ6n,ppQ6p,ppQ6s, &
         Q6c,Q6n,Q6p,Q6s, iiBenPhyto,BenPhyto,iiC,CoupledtoBDc, &
         LocalDelta, D6m,D7m,D8m,D9m,jbotQ6c,jbotQ6n,jbotQ6p,jbotQ6s
       use mem_Param,ONLY:p_d_tot_2,p_poro
       use mem_Settling,ONLY: p_xresusQ6
       use mem_BenPhyto,ONLY:p_useparams,CalculateBenPhyto_vector

       use constants,ONLY:INTEGRAL
       use mem_globalfun,   ONLY: IntegralExpDist_vector,SetExpDist_vector
       use LimitRates, ONLY:LimitChange_vector
!      use mem,  ONLY: Output2d_1,Output2d_2,Output2d_3, Output2d_4
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
       
       integer                                        :: layer,l,nrphyto,iout
       REAL(RLEN),dimension(NO_BOXES_XY)              :: thlayer,rlim
       REAL(RLEN),dimension(NO_BOXES_XY)              :: alpha,cmm,pxperstep,h
       REAL(RLEN),dimension(NO_BOXES_XY)              :: pxc
       REAL(RLEN),dimension(NO_BOXES_XY)              :: pxn
       REAL(RLEN),dimension(NO_BOXES_XY)              :: pxp
       REAL(RLEN),dimension(NO_BOXES_XY)              :: pxs
!      REAL(RLEN),dimension(NO_BOXES_XY)              :: psilt
!      REAL(RLEN),dimension(NO_BOXES_XY)              :: silt
       REAL(RLEN),dimension(NO_BOXES_XY)              :: bpc,tbpc,lim_r
       REAL(RLEN)                                     :: srESS3
       REAL(RLEN)                                     :: r_xany !any rate /m3
       REAL(RLEN)                                     :: scalar
       
       !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
       ! user defined external functions
       !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
       integer, external  :: D3toD1
       integer, external  :: D2toD1
       !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

       if (p_xresusQ6.lt.NZERO) then
         jbotQ6c=ZERO; jbotQ6n=ZERO; jbotQ6p=ZERO; jbotQ6s=ZERO;
         return
       endif
       !jrESS is reuspension/d->in one timestep jrESS*LocalDelta is resuspended.
       ! mg Silt/m2/d*d/(mg silt/m2)*1m
       r_xany=ZERO; tbpc=ZERO;
       call CouplingBioSilt(2,thlayer)
       thlayer=max(ZERO,thlayer*LocalDelta)
       do l=1,iiBenPhyto
         nrphyto=p_useparams(l)
         if (CoupledtoBDc(nrphyto)>0) then
           h=ZERO;bpc=max(ZERO,BenPhyto(l,iiC))
           scalar=maxval(bpc)
           ! calculate fraction of silt which is resuspended...
           if (scalar >1.0D-5)  &
           tbpc=tbpc+  &
               bpc*CalculateBenPhyto_vector(iiC,INTEGRAL,h,thlayer)
         endif
       enddo
       lim_r=1000.0/(tbpc+1000.0)
       where ( jrESS > NZERO) 
         cmm=p_d_tot_2; pxperstep=p_xresusQ6*lim_r/LocalDelta
         alpha=SetExpDist_vector(D6m)
         pxc= pxperstep/ IntegralExpDist_vector(-alpha, cmm)* &
              IntegralExpDist_vector( -alpha,thlayer) 
         alpha=SetExpDist_vector(D7m)
         pxn= pxperstep/ IntegralExpDist_vector(-alpha, cmm)* &
               IntegralExpDist_vector( -alpha,thlayer) 
         alpha=SetExpDist_vector(D8m)
         pxp= pxperstep/ IntegralExpDist_vector(-alpha, cmm)* &
               IntegralExpDist_vector( -alpha,thlayer) 
         alpha=SetExpDist_vector(D9m)
         pxs= pxperstep/ IntegralExpDist_vector(-alpha, cmm)* &
               IntegralExpDist_vector( -alpha,thlayer) 
       elsewhere
         pxc=ZERO; pxn=ZERO; pxp=ZERO; pxs=ZERO;
       endwhere
         if (isnan(pxc(1))) then
           write(LOGUNIT,*) 'ResusPartDetr: pxc,cmm,alpha,D6m,thlayer:', &
                                               pxc,cmm,alpha,D6m,thlayer
         endif
       jbotQ6c=pxc*Q6c
       jbotQ6n=pxn*Q6n
       jbotQ6p=pxp*Q6p
       jbotQ6s=pxs*Q6s
       call findnan(pxp,NO_BOXES_XY,iout)
       call LimitChange_vector(1,jbotQ6c,Q6c,max_change_per_step,rlim)
       pxc=pxc*rlim
       call LimitChange_vector(1,jbotQ6n,Q6n,max_change_per_step,rlim)
       pxn=pxn*rlim
       call LimitChange_vector(1,jbotQ6p,Q6p,max_change_per_step,rlim)
       pxp=pxp*rlim
       call LimitChange_vector(1,jbotQ6s,Q6s,max_change_per_step,rlim)
       pxs=pxs*rlim
       DO BoxNumberXY=1,NO_BOXES_XY
         layer=PelBoxAbove(BoxNumberXY)
         scalar= jrESS(BoxNumberXY)
         if (scalar> NZERO) then
           srESS3=DONE/Depth(layer)
           r_xany=Q6c(BoxNumberXY)*pxc(BoxNumberXY)*srESS3       
           call flux(layer,iiPel,ppR6c,ppR6c,r_xany)
           r_xany=Q6n(BoxNumberXY)*pxn(BoxNumberXY)*srESS3
           call flux(layer,iiPel,ppR6n,ppR6n,r_xany)
           r_xany=Q6p(BoxNumberXY)*pxp(BoxNumberXY)*srESS3
           call flux(layer,iiPel,ppR6p,ppR6p,r_xany)
           r_xany=Q6s(BoxNumberXY)*pxs(BoxNumberXY)*srESS3
           call flux(layer,iiPel,ppR6s,ppR6s,r_xany)
         endif
       enddo
       ! g/m2  =gm/m2  
       jbotQ6c=pxc*Q6c
       jbotQ6n=pxn*Q6n
       jbotQ6p=pxp*Q6p
       jbotQ6s=pxs*Q6s
       call flux_vector(iiBen,ppQ6c,ppQ6c,-jbotQ6c)
       call flux_vector(iiBen,ppQ6n,ppQ6n,-jbotQ6n)
       call flux_vector(iiBen,ppQ6p,ppQ6p,-jbotQ6p)
       call flux_vector(iiBen,ppQ6s,ppQ6s,-jbotQ6s)
       !detritus is dirct redistributed over the layeras according
       !the distribution of the resuspendend silt.
      end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
