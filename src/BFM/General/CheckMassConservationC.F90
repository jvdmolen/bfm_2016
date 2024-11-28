#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: CheckMassConservationC
!
! DESCRIPTION
!   ! Calculation of the budgets/massconservation  for inorganic C, organic C and Alkalinity (H+)
! AUTHOR 
!        Piet Ruardij
! !INTERFACE
      subroutine CheckMassConservationCDynamics
!
      use constants, ONLY:RLEN,ZERO,SEC_PER_DAY
      use mem,only:D3SINK,D3SOURCE,D2SINK,D2SOURCE,D2STATE,D3STATE
      use mem,only:jupPELc,jminPELc,jupBENc,jminBENc,jY3O3c,jtotbenpelc
      use mem,only:B1c,Q1c,Q6c,totPelc,totBENc,totPELInc,totBENInc
      use mem,only:totPELh,totBENh,N3n,N6r,PELBOTTOM,Depth
      use mem,only:NO_BOXES_X,NO_BOXES_Y,NO_BOXES_Z,NO_BOXES
      use mem,only:BoxNumberX,BoxNumberY,BoxNumberXY,BoxNumberZ
      use mem,only:iiBen,iiPel,ppB1c,ppO3c,ppO3h,ppG3c,ppG13c,ppG23c, &
         ppG3h,ppG13h,ppG23h,iiC,iiH2,iiH3,iiHN,ppK6r,ppK16r,ppK26r,ppK3n
      use mem,only:ppMicroZooplankton,ppMesoZooPlankton,&
        iiMicroZooplankton,iiMesoZooPlankton,&
        iiPhytoPlankton,ppPhytoPlankton,iiBenBacteria,ppBenBacteria,&
        iiBenOrganisms,ppBenOrganisms,iiSuspensionFeeders,ppSuspensionFeeders
      use mem,only:MicroZooplankton,MesoZooPlankton,PhytoPlankton, &
        BenBacteria, BenOrganisms,SuspensionFeeders,ppPelDetritus, &
        iiPelDetritus,PelDetritus,iiBenLabileDetritus,BenLabileDetritus, &
        BenPhyto,iiBenPhyto,ppBenthicCO2
      use mem_Param,only:CalcBenthicFlag,p_qro,combine_anabac
      use mem_CO2, ONLY: p_qhK4K3n,p_qhATo
      use constants,  ONLY: BENTHIC_RETURN, BENTHIC_BIO, BENTHIC_FULL

      implicit none
      real(RLEN),dimension(:),pointer  :: sc
      real(RLEN),dimension(NO_BOXES_Z) :: d
      integer                           ::i,j,nr
      integer                           ::f,t
    
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! user defined external functions
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
       integer, external  :: D3toD1
       integer, external  :: D2toD1
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


      jtotbenpelc=ZERO
      jupPELc=ZERO  ; jupBENc=ZERO
      jminPELc=ZERO ; jminBENc=ZERO
      totPELc=ZERO  ; totBENc = ZERO
      totPELInc=ZERO; totBENInc=ZERO
      totPELh=ZERO  ; totBENh=ZERO

      BoxNumberZ = NO_BOXES_Z
      DO BoxNumberY=1,NO_BOXES_Y
        DO BoxNumberX=1,NO_BOXES_X
         f=D3toD1(BoxNumberX,BoxNumberY,1)
         t=D3toD1(BoxNumberX,BoxNumberY,BoxNumberZ)
         BoxNumberXY=D2toD1(BoxNumberX,BoxNumberY)

         d=Depth(f:t)

         nr=ppO3c
         do i=1,iiPhytoPlankton
            j=ppPhytoPlankton(i,iiC)
            if (ppO3c.eq.0) nr=j
            jupPELc(BoxNumberXY)=jupPELc(BoxNumberXY) &
                    +sum(D3SOURCE(j,nr,f:t)*d)*SEC_PER_DAY 
            if (ppO3c.gt.0)  &
              jminPELc(BoxNumberXY)=jminPELc(BoxNumberXY) &
                    +sum(D3SINK(j,nr,f:t)*d)*SEC_PER_DAY
            sc=> PhytoPlankton(i,iiC); 
            totPELc(BoxNumberXY)=totPELc(BoxNumberXY) +sum(sc(f:t)*d)
            jtotbenpelc(BoxNumberXY)=jtotbenpelc(BoxNumberXY) &
                                         +PELBOTTOM(j,BoxNumberXY)
         enddo
         do i=1,iiMesoZooPlankton
            j=ppMesoZooPlankton(i,iiC)
            if (ppO3c.eq.0) nr=j
            jminPELc(BoxNumberXY)=jminPELc(BoxNumberXY) &
                    +sum(D3SINK(j,nr,f:t)*d) *SEC_PER_DAY
            sc=> MesoZooplankton(i,iiC)
            totPELc(BoxNumberXY)=totPELc(BoxNumberXY) +sum(sc(f:t)*d)
            jtotbenpelc(BoxNumberXY)=jtotbenpelc(BoxNumberXY)&
                                 +PELBOTTOM(j,BoxNumberXY)
         enddo
         do i=1,iiMicroZooPlankton
            j=ppMicroZooPlankton(i,iiC)
            if (ppO3c.eq.0) nr=j
            jminPELc(BoxNumberXY)=jminPELc(BoxNumberXY) &
                    +sum(D3SINK(j,nr,f:t)*d) *SEC_PER_DAY
            sc=> MicroZooplankton(i,iiC)
            totPELc(BoxNumberXY)=totPELc(BoxNumberXY) +sum(sc(f:t)*d)
            jtotbenpelc(BoxNumberXY)=jtotbenpelc(BoxNumberXY)&
                                      +PELBOTTOM(j,BoxNumberXY)
         enddo
         do i=1, iiPelDetritus
           j= ppPelDetritus(i,iiC)
           sc=> PelDetritus(i,iiC)
           totPELc(BoxNumberXY)=totPELc(BoxNumberXY) +sum(sc(f:t)*d)
           jtotbenpelc(BoxNumberXY)=jtotbenpelc(BoxNumberXY)&
                                   +PELBOTTOM(j,BoxNumberXY)
         enddo

         totPELc(BoxNumberXY)=totPELc(BoxNumberXY)+ sum(B1c(f:t)* d)
         !jtotbenpelc and B1c : no flux to sediment of B1
         if ( ppO3c.eq.0) nr=ppB1c
         jminPELc(BoxNumberXY)=jminPELc(BoxNumberXY) &
                    +sum(D3SINK(ppB1c,nr,f:t)*d)*SEC_PER_DAY

         if ( ppO3c.gt.0) & 
         totPELInc(BoxNumberXY)=sum(D3STATE(ppO3c,f:t)*d)
         if (ppO3h.gt.0) totPELh(BoxNumberXY)=sum((D3STATE(ppO3h,f:t) &
                         +N3n(f:t)*p_qhK4K3n-N6r(f:t)*p_qhAto/p_qro)*d) 
         
         !respiration of FIlterFeeders(Y3), which take partly place in pelagic
         jminPELc(BoxNumberXY)=jminPELc(BoxNumberXY)+jY3O3c(BoxNumberXY)

         select case ( CalcBenthicFlag)
          case ( 0 )
          case ( BENTHIC_RETURN )  ! Simple benthic return
            ! Mass conservation variables
            totBENc  =  ( Q1c+ Q6c)
          case ( BENTHIC_BIO,BENTHIC_FULL )  ! Intermediate benthic return
            nr=ppG3c
            do i=1, iiBenOrganisms
              j=ppBenOrganisms(i,iiC)
              if ( ppG3c.eq.0) nr=j
              jminBENc(BoxNumberXY)=jminBENc(BoxNumberXY) &
                    +D2SINK(j,nr,BoxNumberXY) *SEC_PER_DAY
              sc => BenOrganisms(i,iiC)
              totBENc(BoxNumberXY)=totBENc(BoxNumberXY)+sc(BoxNumberXY)
            enddo
            do i=1, iiSuspensionFeeders
              j=ppSusPensionFeeders(i,iiC)
              if ( ppG3c.eq.0) nr=j
              jminBENc(BoxNumberXY)=jminBENc(BoxNumberXY) &
                    +D2SINK(j,nr,BoxNumberXY) *SEC_PER_DAY
              sc => SusPensionFeeders(i,iiC)
              totBENc(BoxNumberXY)=totBENc(BoxNumberXY)+sc(BoxNumberXY)
            enddo
            do i=1, iiBenLabileDetritus
              sc => BenLabileDetritus(i,iiC)
              totBENc(BoxNumberXY)=totBENc(BoxNumberXY)+sc(BoxNumberXY)
            enddo
            do i=1, iiBenPhyto
              sc => BenPhyto(i,iiC)
              totBENc(BoxNumberXY)=totBENc(BoxNumberXY)+sc(BoxNumberXY)
              nr=ppG3c;if ( ppG3c.eq.0) nr=j
              jupBENc(BoxNumberXY)=jupBENc(BoxNumberXY)&
                            +D2SOURCE(nr,j,BoxNumberXY) *SEC_PER_DAY
            enddo
            do i=1, iiBenBacteria
              j=ppBenBacteria(i,iiC)
              if (i.ne.iiH2) then
                nr=ppG3c;if ( ppG3c.eq.0) nr=j
                if (i.eq.iiHN)jupBENc(BoxNumberXY)=jupBENc(BoxNumberXY)&
                            +D2SOURCE(nr,j,BoxNumberXY) *SEC_PER_DAY
                !Aerobic bacteria+ nitrifiers
                jminBENc(BoxNumberXY)=jminBENc(BoxNumberXY) &
                    +D2SINK(j,nr,BoxNumberXY) *SEC_PER_DAY
              elseif (combine_anabac.and.i==iiH2) then
                !Anaerobic bacteria
                nr=ppG13c;if (ppG13c.eq.0) nr=j
                jminBENc(BoxNumberXY)=jminBENc(BoxNumberXY) &
                    +D2SINK(j,nr,BoxNumberXY) *SEC_PER_DAY
                nr=ppG23c;if (ppG23c.gt.0) & 
                   jminBENc(BoxNumberXY)=jminBENc(BoxNumberXY) &
                      +D2SINK(j,nr,BoxNumberXY) *SEC_PER_DAY
              else
                !Anaerobic bacteria
                nr=ppBenthicCO2(i,iiC);if (nr.eq.0) nr=j
                jminBENc(BoxNumberXY)=jminBENc(BoxNumberXY) &
                      +D2SINK(j,nr,BoxNumberXY) *SEC_PER_DAY
              endif
              sc => BenBacteria(i,iiC)
              totBENc(BoxNumberXY)=totBENc(BoxNumberXY)+sc(BoxNumberXY)
            enddo
            ! Mass conservation variables
            totBENc(BoxNumberXY) = totBENc(BoxNumberXY)+Q6c(BoxNumberXY)
            if ( ppG3c.gt.0) & 
               totBENInc(BoxNumberXY)=D2STATE(ppG3c,BoxNumberXY) &
               +D2STATE(ppG13c,BoxNumberXY)+D2STATE(ppG23c,BoxNumberXY)
            if ( ppG3h.gt.0)  totBENh(BoxNumberXY)= &
                D2STATE(ppG3h,BoxNumberXY) &
               +D2STATE(ppG13h,BoxNumberXY)+D2STATE(ppG23h,BoxNumberXY) &
               -p_qhAto/p_qro *(D2STATE(ppK6r,BoxNumberXY) &
               +D2STATE(ppK16r,BoxNumberXY)+D2STATE(ppK26r,BoxNumberXY)) &
               +D2STATE(ppK3n,BoxNumberXY)*p_qhK4K3n
         end select
        enddo
      enddo

      return
      end
