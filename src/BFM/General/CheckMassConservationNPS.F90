#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: CheckMassConservationNPS
!
! DESCRIPTION
!   !

!   This file is generated directly from OpenSesame model code, using a code
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  subroutine CheckMassConservationNPSDynamics
!
! !USES:
  ! The following Pelagic-states are used (NOT in fluxes): &
  ! B1p, N1p, N3n, N4n, P1s, R6s, N5s
  ! The following Benthic-states are used (NOT in fluxes): Q1p, Q6p, Q1n, Q6n, &
  ! Q6s, Q11p, K1p, K11p, Q11n, K4n, K14n, K21p, G4n, K3n, K24n, K5s
  ! The following Benthic 1-d global boxvars got a value: totBENp, totBENn, &
  ! totBENs
  ! The following 0-d global parameters are used: CalcBenthicFlag
  ! The following global constants are used: RLEN
  ! The following constants are used: BENTHIC_RETURN, BENTHIC_BIO, BENTHIC_FULL

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN
#ifdef NOPOINTERS
  use mem,  ONLY: D2STATE,D3STATE
#else
  use mem, ONLY: B1p,  N1p, &
    B1n, N3n, N4n, &
    N5s, D2STATE
  use mem, ONLY: Q1p, Q6p, Q1n, Q6n, Q6s, &
    K1p, K11p, K4n, K14n, &
    K3n, K5s, K15s, D3STATE
#endif
  use mem,ONLY: ppB1p,ppN1p, ppB1n, ppN3n, ppN4n, ppN5s, Depth
  use mem, ONLY: ppQ6p, ppQ6n, ppQ6s, ppK1p, ppK11p, ppK4n, ppK14n, &
    ppK3n, ppK5s,ppK15s, totBENp, totBENn, totBENs, &
    totPELp,totPELn,totPELs, totSYSp,totSYSn,totSYSs, NO_BOXES_XY, &
    NO_BOXES_Z, NO_BOXES_X,NO_BOXES_Y,NO_BOXES,  &
    BoxNumberX,BoxNumberY,BoxNumberZ,BoxNumberXY,  &
    ppMicroZooplankton,ppMesoZooPlankton,MicroZooplankton,MesoZooplankton, &
    iiMicroZooplankton,iiMesoZooplankton,iiC,iiN,iiP,iiS, &
    iiPhytoPlankton,ppPhytoPlankton,PhytoPlankton, &
    iiBenBacteria,ppBenBacteria,BenBacteria, &
    iiBenLabileDetritus,ppBenLabileDetritus,BenLabileDetritus,&
    iiBenUrea,ppBenUrea,BenUrea,&
    iiPelDetritus,ppPelDetritus,PelDetritus, &
    iiBenthicAmmonium,ppBenthicAmmonium, BenthicAmmonium, &
    iiBenthicPhosphate,ppBenthicPhosphate,BenthicPhosphate, &
    iiBenOrganisms,ppBenOrganisms,BenOrganisms,iiBenPhyto,ppBenPhyto,BenPhyto, &
    iiSuspensionFeeders,ppSuspensionFeeders,SuspensionFeeders, &
    jtotBENPELn,jtotBENPELp,jtotBENPELs,PELBOTTOM,jtPELn,jtPELp,jtPELs, &
    jtBENn,jtBENp,jtBENs,Source_D3_vector,Source_D2_vector

  use constants,  ONLY: ZERO,BENTHIC_RETURN, BENTHIC_BIO, BENTHIC_FULL
  use mem_Param,  ONLY: CalcBenthicFlag
  use mem_MesoZoo, ONLY: p_qnMec=>p_qnc,p_qpMec=>p_qpc
  use mem_MicroZoo, ONLY: p_qnMic=>p_qnc,p_qpMic=>p_qpc
!
!
! !AUTHORS
!   Piet Ruardij
!
!
! !REVISION_HISTORY
!   Created at Mon Nov 21 09:44:23 CET 2005
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
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN),dimension(:),pointer  :: sc
  real(RLEN),dimension(NO_BOXES_Z) :: s
  real(RLEN),dimension(NO_BOXES_Z) :: d
  real(RLEN),dimension(:),pointer  :: rc
  real(RLEN),dimension(NO_BOXES_XY):: r,fl
  integer                           ::i,j,k,l
  integer                           ::f,t

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! user defined external functions
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer, external  :: D3toD1
  integer, external  :: D2toD1
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  BoxNumberZ = NO_BOXES_Z
  DO BoxNumberY=1,NO_BOXES_Y
    DO BoxNumberX=1,NO_BOXES_X
      f=D3toD1(BoxNumberX,BoxNumberY,1)
      t=D3toD1(BoxNumberX,BoxNumberY,BoxNumberZ)
      BoxNumberXY=D2toD1(BoxNumberX,BoxNumberY)

      totPELp(BoxNumberXY)=ZERO
      totPELn(BoxNumberXY)=ZERO
      totPELs(BoxNumberXY)=ZERO

      jtotBENPELp(BoxNumberXY)=ZERO
      jtotBENPELn(BoxNumberXY)=ZERO
      jtotBENPELs(BoxNumberXY)=ZERO

      jtPELp(BoxNumberXY)=ZERO
      jtPELn(BoxNumberXY)=ZERO
      jtPELs(BoxNumberXY)=ZERO

      d=Depth(f:t)
      do i=1, iiMicroZooplankton
         l= ppMicroZooplankton(i,iiC)
         sc=> MicroZooplankton(i,iiC)
         k=ppMicroZooPlankton(i,iiN)
         if ( k==0) then
           s=sc(f:t);s=s*p_qnMic(i)
           fl=PELBOTTOM(:,l)*p_qnMic(i)
         else
           rc=> MicroZooplankton(i,iiN);s=rc(f:t)
           fl=PELBOTTOM(:,k)
         endif
         totPELn(BoxNumberXY)=totPELn(BoxNumberXY) +sum(s*d)
         jtotBENPELn(BoxNumberXY)=jtotBENPELn(BoxNumberXY)+fl(BoxNumberXY)
         if ( k==0) then
           s=Source_D3_vector(l,0)*p_qnMic(i)
         else
           s=Source_D3_vector(k,0)
         endif
         jtPELn(BoxNumberXY)=jtPELn(BoxNumberXY)+sum(s*d)
         k=ppMicroZooPlankton(i,iiP)
         if ( k==0) then
           s=sc(f:t);s=s*p_qpMic(i)
           fl=PELBOTTOM(:,l)*p_qpMic(i)
         else
           rc=> MicroZooplankton(i,iiP);s=rc(f:t)
           fl=PELBOTTOM(:,k)
         endif
         totPELp(BoxNumberXY)=totPELp(BoxNumberXY) +sum(s*d)
         jtotBENPELp(BoxNumberXY)=jtotBENPELp(BoxNumberXY)+fl(BoxNumberXY)
         if ( k==0) then
           s=Source_D3_vector(l,0)*p_qpMic(i)
         else
           s=Source_D3_vector(k,0)
         endif
         jtPELp(BoxNumberXY)=jtPELp(BoxNumberXY)+sum(s*d)
      enddo
      do i=1, iiMesoZooplankton
         l= ppMesoZooplankton(i,iiC)
         sc=> MesoZooplankton(i,iiC)
         k=ppMesoZooPlankton(i,iiN)
         if ( k==0) then
           s=sc(f:t);s=s*p_qnMec(i)
           fl=PELBOTTOM(:,l)*p_qnMec(i)
         else
           rc=> MesoZooplankton(i,iiN);s=rc(f:t)
           fl=PELBOTTOM(:,k)
         endif
         totPELn(BoxNumberXY)=totPELn(BoxNumberXY) +sum(s*d)
         jtotBENPELn(BoxNumberXY)=jtotBENPELn(BoxNumberXY)+fl(BoxNumberXY)
         if ( k==0) then
           s=Source_D3_vector(l,0)*p_qnMec(i)
         else
           s=Source_D3_vector(k,0)
         endif
         jtPELn(BoxNumberXY)=jtPELn(BoxNumberXY)+sum(s*d)
         k=ppMesoZooPlankton(i,iiP)
         if ( k==0) then
           s=sc(f:t);s=s*p_qpMec(i)
           fl=PELBOTTOM(:,l)*p_qPMec(i)
         else
           rc=> MesoZooplankton(i,iiP);s=rc(f:t)
           fl=PELBOTTOM(:,k)
         endif
         totPELp(BoxNumberXY)=totPELp(BoxNumberXY) +sum(s*d)
         jtotBENPELp(BoxNumberXY)=jtotBENPELp(BoxNumberXY)+fl(BoxNumberXY)
         if ( k==0) then
           s=Source_D3_vector(l,0)*p_qpMec(i)
         else
           s=Source_D3_vector(k,0)
         endif
         jtPELp(BoxNumberXY)=jtPELp(BoxNumberXY)+sum(s*d)
      enddo
      do i=1, iiPhytoPlankton
         sc=> PhytoPlankton(i,iiN); s=sc(f:t)
         totPELn(BoxNumberXY)=totPELn(BoxNumberXY) +sum(s*d)
         k= ppPhytoPlankton(i,iiN)
         s=Source_D3_vector(k,0)
         jtPELn(BoxNumberXY)=jtPELn(BoxNumberXY)+sum(s*d)
         fl=PELBOTTOM(:,k)
         jtotBENPELn(BoxNumberXY)=jtotBENPELn(BoxNumberXY)+fl(BoxNumberXY)
         sc=> PhytoPlankton(i,iiP); s=sc(f:t)
         totPELp(BoxNumberXY)=totPELp(BoxNumberXY) +sum(s*d)
         k= ppPhytoPlankton(i,iiP)
         s=Source_D3_vector(k,0)
         jtPELp(BoxNumberXY)=jtPELp(BoxNumberXY)+sum(s*d)
         fl=PELBOTTOM(:,k)
         jtotBENPELp(BoxNumberXY)=jtotBENPELp(BoxNumberXY)+fl(BoxNumberXY)
         k=ppPhytoPlankton(i,iiS)
         if (k>0) then
           sc=> PhytoPlankton(i,iiS);s=sc(f:t)
           totPELs(BoxNumberXY)=totPELs(BoxNumberXY) +sum(s*d)
           fl=PELBOTTOM(:,k)
           s=Source_D3_vector(k,0)
           jtPELs(BoxNumberXY)=jtPELs(BoxNumberXY)+sum(s*d)
           jtotBENPELs(BoxNumberXY)=jtotBENPELs(BoxNumberXY)+fl(BoxNumberXY)
         endif
      enddo
      do i=1, iiPelDetritus
         k=ppPelDetritus(i,iiN)
         if ( k>0) then
           sc=> PelDetritus(i,iiN);s=sc(f:t)
           totPELn(BoxNumberXY)=totPELn(BoxNumberXY) +sum(s*d)
           s=Source_D3_vector(k,0)
           jtPELn(BoxNumberXY)=jtPELn(BoxNumberXY)+sum(s*d)
           fl=PELBOTTOM(:,k)
           jtotBENPELn(BoxNumberXY)=jtotBENPELn(BoxNumberXY)+fl(BoxNumberXY)
         endif
         k=ppPelDetritus(i,iiP)
         if ( k>0) then
           sc=> PelDetritus(i,iiP);s=sc(f:t)
           totPELp(BoxNumberXY)=totPELp(BoxNumberXY) +sum(s*d)
           s=Source_D3_vector(k,0)
           jtPELp(BoxNumberXY)=jtPELp(BoxNumberXY)+sum(s*d)
           fl=PELBOTTOM(:,k)
           jtotBENPELp(BoxNumberXY)=jtotBENPELp(BoxNumberXY)+fl(BoxNumberXY)
         endif
         k=ppPelDetritus(i,iiS)
         if ( k>0) then
           sc=> PelDetritus(i,iiS);s=sc(f:t)
           totPELs(BoxNumberXY)=totPELs(BoxNumberXY) +sum(s*d)
           s=Source_D3_vector(k,0)
           jtPELs(BoxNumberXY)=jtPELs(BoxNumberXY)+sum(s*d)
           fl=PELBOTTOM(:,k)
           jtotBENPELs(BoxNumberXY)=jtotBENPELs(BoxNumberXY)+fl(BoxNumberXY)
         endif
      enddo

      totPELp(BoxNumberXY)=totPELp(BoxNumberXY)+ sum((B1p(f:t)+ N1p(f:t))* d)
      totPELn(BoxNumberXY)=totPELn(BoxNumberXY)+ &
                                       sum((B1n(f:t)+ N3n(f:t)+ N4n(f:t))* d)
      totPELs(BoxNumberXY)=totPELs(BoxNumberXY)+ sum(N5s(f:t)* d)

      s=Source_D3_vector(ppB1n,0)
      jtPELn(BoxNumberXY)=jtPELn(BoxNumberXY)+sum(s*d)
      s=Source_D3_vector(ppN3n,0)
      jtPELn(BoxNumberXY)=jtPELn(BoxNumberXY)+sum(s*d)
      s=Source_D3_vector(ppN4n,0)
      jtPELn(BoxNumberXY)=jtPELn(BoxNumberXY)+sum(s*d)
      s=Source_D3_vector(ppB1p,0)
       jtPELp(BoxNumberXY)=jtPELp(BoxNumberXY)+sum(s*d)
      s=Source_D3_vector(ppN1p,0)
       jtPELp(BoxNumberXY)=jtPELp(BoxNumberXY)+sum(s*d)
      s=Source_D3_vector(ppN5s,0)
       jtPELs(BoxNumberXY)=jtPELs(BoxNumberXY)+sum(s*d)

      jtotBENPELp(BoxNumberXY)=jtotBENPELp(BoxNumberXY) &
        +PELBOTTOM(BoxNumberXY,ppB1p)+PELBOTTOM(BoxNumberXY,ppN1p)
      jtotBENPELn(BoxNumberXY)=jtotBENPELn(BoxNumberXY) &
        +PELBOTTOM(BoxNumberXY,ppB1n)+PELBOTTOM(BoxNumberXY,ppN3n)  &
        +PELBOTTOM(BoxNumberXY,ppN4n)
      jtotBENPELs(BoxNumberXY)=jtotBENPELs(BoxNumberXY) &
        +PELBOTTOM(BoxNumberXY,ppN5s)

      jtPELn(BoxNumberXY)=jtPELn(BoxNumberXY)-jtotBENPELn(BoxNumberXY)
      jtPELp(BoxNumberXY)=jtPELp(BoxNumberXY)-jtotBENPELp(BoxNumberXY)
      jtPELs(BoxNumberXY)=jtPELs(BoxNumberXY)-jtotBENPELs(BoxNumberXY)


      totBENn(BoxNumberXY) =   ZERO
      totBENp(BoxNumberXY) =   ZERO
      totBENs(BoxNumberXY) =   ZERO

      jtBENn(BoxNumberXY) =   ZERO
      jtBENp(BoxNumberXY) =   ZERO
      jtBENs(BoxNumberXY) =   ZERO

      select case ( CalcBenthicFlag)
        case ( 0 )
        case ( BENTHIC_RETURN )  ! Simple benthic return
          ! Mass conservation variables
          totBENp(BoxNumberXY)  =  ( Q1p(BoxNumberXY)+ Q6p(BoxNumberXY))
          totBENn(BoxNumberXY)  =  ( Q1n(BoxNumberXY)+ Q6n(BoxNumberXY))
          totBENs(BoxNumberXY)  =  ( Q6s(BoxNumberXY))
        case ( BENTHIC_BIO )  ! Intermediate benthic return
         do i=1, iiBenPhyto
           rc => BenPhyto(i,iiN)
           totBENn(BoxNumberXY)=totBENn(BoxNumberXY)+rc(BoxNumberXY)
           r=Source_D2_vector(ppBenPhyto(i,iiN),0)
           jtBenn(BoxNumberXY)=jtBenn(BoxNumberXY)+r(BoxNUmberXY)
           rc => BenPhyto(i,iiP)
           totBENp(BoxNumberXY)=totBENp(BoxNumberXY)+rc(BoxNumberXY)
           r=Source_D2_vector(ppBenPhyto(i,iiP),0)
           jtBenp(BoxNumberXY)=jtBenp(BoxNumberXY)+r(BoxNUmberXY)
           j=ppBenPhyto(i,iiS)
           if ( j>0) then
             rc => BenPhyto(i,iiS)
             totBENs(BoxNumberXY)=totBENs(BoxNumberXY)+rc(BoxNumberXY)
             r=Source_D2_vector(j,0)
             jtBens(BoxNumberXY)=jtBens(BoxNumberXY)+r(BoxNUmberXY)
           endif
         enddo
         do i=1, iiBenOrganisms
           rc => BenOrganisms(i,iiN)
           totBENn(BoxNumberXY)=totBENn(BoxNumberXY)+rc(BoxNumberXY)
           r=Source_D2_vector(ppBenOrganisms(i,iiN),0)
           jtBenn(BoxNumberXY)=jtBenn(BoxNumberXY)+r(BoxNUmberXY)
           rc => BenOrganisms(i,iiP)
           totBENp(BoxNumberXY)=totBENp(BoxNumberXY)+rc(BoxNumberXY)
           r=Source_D2_vector(ppBenOrganisms(i,iiP),0)
           jtBenp(BoxNumberXY)=jtBenp(BoxNumberXY)+r(BoxNUmberXY)
         enddo
         do i=1, iiSuspensionFeeders
           rc => SuspensionFeeders(i,iiN)
           totBENn(BoxNumberXY)=totBENn(BoxNumberXY)+rc(BoxNumberXY)
           r=Source_D2_vector(ppSuspensionFeeders(i,iiN),0)
           jtBenn(BoxNumberXY)=jtBenn(BoxNumberXY)+r(BoxNUmberXY)
           rc => SuspensionFeeders(i,iiP)
           totBENp(BoxNumberXY)=totBENp(BoxNumberXY)+rc(BoxNumberXY)
           r=Source_D2_vector(ppSuspensionFeeders(i,iiP),0)
           jtBenp(BoxNumberXY)=jtBenp(BoxNumberXY)+r(BoxNUmberXY)
         enddo
         do i=1, iiBenLabileDetritus
           rc => BenLabileDetritus(i,iiN)
           totBENn(BoxNumberXY)=totBENn(BoxNumberXY)+rc(BoxNumberXY)
           r=Source_D2_vector(ppBenLabileDetritus(i,iiN),0)
           jtBenn(BoxNumberXY)=jtBenn(BoxNumberXY)+r(BoxNUmberXY)
           rc => BenLabileDetritus(i,iiP)
           totBENp(BoxNumberXY)=totBENp(BoxNumberXY)+rc(BoxNumberXY)
           r=Source_D2_vector(ppBenLabileDetritus(i,iiP),0)
           jtBenp(BoxNumberXY)=jtBenp(BoxNumberXY)+r(BoxNUmberXY)
         enddo
         do i=1, iiBenUrea
           rc => BenUrea(i,iiN)
           totBENn(BoxNumberXY)=totBENn(BoxNumberXY)+rc(BoxNumberXY)
           r=Source_D2_vector(ppBenUrea(i,iiN),0)
           jtBenn(BoxNumberXY)=jtBenn(BoxNumberXY)+r(BoxNUmberXY)
         enddo
         do i=1, iiBenBacteria
           rc => BenBacteria(i,iiN)
           totBENn(BoxNumberXY)=totBENn(BoxNumberXY)+rc(BoxNumberXY)
           r=Source_D2_vector(ppBenBacteria(i,iiN),0)
           jtBenn(BoxNumberXY)=jtBenn(BoxNumberXY)+r(BoxNUmberXY)
           rc => BenBacteria(i,iiP)
           totBENp(BoxNumberXY)=totBENp(BoxNumberXY)+rc(BoxNumberXY)
           r=Source_D2_vector(ppBenBacteria(i,iiP),0)
           jtBenp(BoxNumberXY)=jtBenp(BoxNumberXY)+r(BoxNUmberXY)
         enddo
          ! Mass conservation variables
          totBENp(BoxNumberXY) = totBENp(BoxNumberXY) &
              +Q6p(BoxNumberXY)+ K1p(BoxNumberXY)+ K11p(BoxNumberXY)
          r=Source_D2_vector(ppQ6p,0)
          jtBenp(BoxNumberXY)=jtBenp(BoxNumberXY)+r(BoxNUmberXY)
          r=Source_D2_vector(ppK1p,0)
          jtBenp(BoxNumberXY)=jtBenp(BoxNumberXY)+r(BoxNUmberXY)
          r=Source_D2_vector(ppK11p,0)
          jtBenp(BoxNumberXY)=jtBenp(BoxNumberXY)+r(BoxNUmberXY)
          totBENn(BoxNumberXY) = totBENn(BoxNumberXY) &
              +Q6n(BoxNumberXY)+ K4n(BoxNumberXY)+ K14n(BoxNumberXY)
          r=Source_D2_vector(ppQ6n,0)
          jtBenn(BoxNumberXY)=jtBenn(BoxNumberXY)+r(BoxNUmberXY)
          r=Source_D2_vector(ppK4n,0)
          jtBenn(BoxNumberXY)=jtBenn(BoxNumberXY)+r(BoxNUmberXY)
          r=Source_D2_vector(ppK14n,0)
          jtBenn(BoxNumberXY)=jtBenn(BoxNumberXY)+r(BoxNUmberXY)
          totBENs(BoxNumberXY)  =  ( Q6s(BoxNumberXY))
          r=Source_D2_vector(ppQ6s,0)
          jtBens(BoxNumberXY)=jtBens(BoxNumberXY)+r(BoxNUmberXY)

          jtBENn(BoxNumberXY)=jtBENn(BoxNumberXY)-jtotBENPELn(BoxNumberXY)
          jtBENp(BoxNumberXY)=jtBENp(BoxNumberXY)-jtotBENPELp(BoxNumberXY)
          jtBENs(BoxNumberXY)=jtBENs(BoxNumberXY)-jtotBENPELs(BoxNumberXY)
       case ( BENTHIC_FULL )  ! Full benthic nutrients
        do i=1, iiBenPhyto
           rc => BenPhyto(i,iiN)
           totBENn(BoxNumberXY)=totBENn(BoxNumberXY)+rc(BoxNumberXY)
           r=Source_D2_vector(ppBenPhyto(i,iiN),0)
           jtBenn(BoxNumberXY)=jtBenn(BoxNumberXY)+r(BoxNUmberXY)
           rc => BenPhyto(i,iiP)
           totBENp(BoxNumberXY)=totBENp(BoxNumberXY)+rc(BoxNumberXY)
           r=Source_D2_vector(ppBenPhyto(i,iiP),0)
           jtBenp(BoxNumberXY)=jtBenp(BoxNumberXY)+r(BoxNUmberXY)
           j=ppBenPhyto(i,iiS)
           if ( j>0) then
             rc => BenPhyto(i,iiS)
             totBENs(BoxNumberXY)=totBENs(BoxNumberXY)+rc(BoxNumberXY)
             r=Source_D2_vector(j,0)
             jtBens(BoxNumberXY)=jtBens(BoxNumberXY)+r(BoxNUmberXY)
           endif
         enddo
         do i=1, iiBenOrganisms
           rc => BenOrganisms(i,iiN)
           totBENn(BoxNumberXY)=totBENn(BoxNumberXY)+rc(BoxNumberXY)
           r=Source_D2_vector(ppBenOrganisms(i,iiN),0)
           jtBenn(BoxNumberXY)=jtBenn(BoxNumberXY)+r(BoxNUmberXY)
           rc => BenOrganisms(i,iiP)
           totBENp(BoxNumberXY)=totBENp(BoxNumberXY)+rc(BoxNumberXY)
           r=Source_D2_vector(ppBenOrganisms(i,iiP),0)
           jtBenp(BoxNumberXY)=jtBenp(BoxNumberXY)+r(BoxNUmberXY)
         enddo
         do i=1, iiSuspensionFeeders
           rc => SuspensionFeeders(i,iiN)
           totBENn(BoxNumberXY)=totBENn(BoxNumberXY)+rc(BoxNumberXY)
           r=Source_D2_vector(ppSusPensionFeeders(i,iiN),0)
           jtBenn(BoxNumberXY)=jtBenn(BoxNumberXY)+r(BoxNUmberXY)
           rc => SuspensionFeeders(i,iiP)
           totBENp(BoxNumberXY)=totBENp(BoxNumberXY)+rc(BoxNumberXY)
           r=Source_D2_vector(ppSusPensionFeeders(i,iiP),0)
           jtBenp(BoxNumberXY)=jtBenp(BoxNumberXY)+r(BoxNUmberXY)
         enddo
         do i=1, iiBenLabileDetritus
           rc => BenLabileDetritus(i,iiN)
           totBENn(BoxNumberXY)=totBENn(BoxNumberXY)+rc(BoxNumberXY)
           r=Source_D2_vector(ppBenLabileDetritus(i,iiN),0)
           jtBenn(BoxNumberXY)=jtBenn(BoxNumberXY)+r(BoxNUmberXY)
           rc => BenLabileDetritus(i,iiP)
           totBENp(BoxNumberXY)=totBENp(BoxNumberXY)+rc(BoxNumberXY)
           r=Source_D2_vector(ppBenLabileDetritus(i,iiP),0)
           jtBenp(BoxNumberXY)=jtBenp(BoxNumberXY)+r(BoxNUmberXY)
         enddo
         do i=1, iiBenUrea
           rc => BenUrea(i,iiN)
           totBENn(BoxNumberXY)=totBENn(BoxNumberXY)+rc(BoxNumberXY)
           r=Source_D2_vector(ppBenUrea(i,iiN),0)
           jtBenn(BoxNumberXY)=jtBenn(BoxNumberXY)+r(BoxNUmberXY)
         enddo
         do i=1, iiBenBacteria
           rc => BenBacteria(i,iiN)
           totBENn(BoxNumberXY)=totBENn(BoxNumberXY)+rc(BoxNumberXY)
           r=Source_D2_vector(ppBenBacteria(i,iiN),0)
           jtBenn(BoxNumberXY)=jtBenn(BoxNumberXY)+r(BoxNUmberXY)
           rc => BenBacteria(i,iiP)
           totBENp(BoxNumberXY)=totBENp(BoxNumberXY)+rc(BoxNumberXY)
           r=Source_D2_vector(ppBenBacteria(i,iiP),0)
           jtBenp(BoxNumberXY)=jtBenp(BoxNumberXY)+r(BoxNUmberXY)
         enddo
         do i=1,iiBenthicPhosphate
           rc => BenthicPhosphate(i,iiP)
           totBENp(BoxNumberXY)=totBENp(BoxNumberXY)+rc(BoxNumberXY)
           r=Source_D2_vector(ppBenthicPhosphate(i,iiP),0)
           jtBenp(BoxNumberXY)=jtBenp(BoxNumberXY)+r(BoxNUmberXY)
         enddo
         do i=1,iiBenthicAmmonium
           rc => BenthicAmmonium(i,iiN)
           totBENn(BoxNumberXY)=totBENn(BoxNumberXY)+rc(BoxNumberXY)
           r=Source_D2_vector(ppBenthicAmmonium(i,iiN),0)
           jtBenn(BoxNumberXY)=jtBenn(BoxNumberXY)+r(BoxNUmberXY)
         enddo
         totBENn(BoxNumberXY) =totBENn(BoxNumberXY)  &
                                          + Q6n(BoxNumberXY)+K3n(BoxNumberXY)
         r=Source_D2_vector(ppQ6n,0); 
         jtBenn(BoxNumberXY)=jtBenn(BoxNumberXY)+r(BoxNUmberXY)
         r=Source_D2_vector(ppK3n,0); 
         jtBenn(BoxNumberXY)=jtBenn(BoxNumberXY)+r(BoxNUmberXY)
         totBENp(BoxNumberXY) =totBENp(BoxNumberXY) + Q6p(BoxNumberXY)
         r=Source_D2_vector(ppQ6p,0); 
         jtBenp(BoxNumberXY)=jtBenp(BoxNumberXY)+r(BoxNUmberXY)
         totBENs(BoxNumberXY) =totBENs(BoxNumberXY)  &
                        + K5s(BoxNumberXY)+K15s(BoxNumberXY)+Q6s(BoxNumberXY)
         r=Source_D2_vector(ppK5s,0);  
         jtBens(BoxNumberXY)=jtBens(BoxNumberXY)+r(BoxNUmberXY)
         r=Source_D2_vector(ppK15s,0); 
         jtBens(BoxNumberXY)=jtBens(BoxNumberXY)+r(BoxNUmberXY)
         r=Source_D2_vector(ppQ6s,0); 
         jtBens(BoxNumberXY)=jtBens(BoxNumberXY)+r(BoxNUmberXY)

          jtBENn(BoxNumberXY)=jtBENn(BoxNumberXY)+jtotBENPELn(BoxNumberXY)
          jtBENp(BoxNumberXY)=jtBENp(BoxNumberXY)+jtotBENPELp(BoxNumberXY)
          jtBENs(BoxNumberXY)=jtBENs(BoxNumberXY)+jtotBENPELs(BoxNumberXY)
      end select

      totSYSn(BoxNumberXY)=totPELn(BoxNumberXY)+totBENn(BoxNumberXY)
      totSYSp(BoxNumberXY)=totPELp(BoxNumberXY)+totBENp(BoxNumberXY)
      totSYSs(BoxNumberXY)=totPELs(BoxNumberXY)+totBENs(BoxNumberXY)
    enddo
  enddo
  end
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
