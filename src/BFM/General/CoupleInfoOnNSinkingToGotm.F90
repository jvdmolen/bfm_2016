#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: General/CoupleInfoOnNSinkingToGotm.F90
!
! DESCRIPTION
!   !
!
! !INTERFACE
  subroutine CoupleInfoOnNSinkingToGotm

! use global_mem, ONLY:RLEN,LOGUNIT,ZERO,NZERO,DONE
  use mem, ONLY: ppR2c,ppR2n,ppR6p, ppR6c, ppR6n, ppR6s, ppRZc,ppB1c, &
    ppB1c,ppB1n,ppB1p,ppBac, ppP6c,ppR3c,ppPcc, &
    ppMicroZooPlankton, ppMesoZooPlankton, ppPhytoPlankton,iiPhytoPlankton,  &
    iiMicroZooPlankton,iiMesoZooPlankton, &
    iiP,iiC,iiN,iiL,iiS, &
    iiPELSINKREF, &
    ppsediB1,ppsediR2,ppsediRZ,ppsediR6,ppsediR6s,ppsediPI,ppsediMiZ,ppsediMeZ
!     sediB1,  sediR2,  sediRZ, sediR6,   sediR6s, sediPI,   sediMeZ,  sediMiZ

    implicit none
    integer                   ::i,j,j0
 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Sedimentation scheme
  ! A diagnotic variable in which the sedimentation rate is coupled in
  ! this way to the state variables:
  ! For example for all R6 state varables in bio_bfm.F90 the sinking rates
  ! as defined in array sediR6 are used. In bio_bfm.F90 a recalculation take
  ! place from rates per day to rates per second and if necessay a limitation
  ! of the sinking rate for shallow grid points. 
  ! The sequence num ber in which teh sinking rate is defined
  ! is given to the an array element of iiPelSINKREF withe sequence number 
  ! of the the C-concistuent of the state variable for which the sinking
  ! rate is defined 
  ! For the other constituents in the values  are indiect copied to the
  ! the array elemnt of iiPelSinkREF througingn them a neative value of the
  ! sequence number of the state variable
  ! the values for R6c (iiPelSINKREF<0 )
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   ! These assigments are done only once:the first time this routine is called
     !Detritus 
     iiPELSINKREF(ppR2c)=  ppsediR2
     iiPELSINKREF(ppR2n)=  -ppR2c
     iiPELSINKREF(ppRZc)=  ppsediRZ
     iiPELSINKREF(ppR6c)=  ppsediR6
     iiPELSINKREF(ppR6n)=-ppR6c
     iiPELSINKREF(ppR6p)=-ppR6c
     iiPELSINKREF(ppR6s)=ppsediR6s
     !Bacteria
     iiPELSINKREF(ppB1c)=ppsediB1
     iiPELSINKREF(ppB1n)=-ppB1c
     iiPELSINKREF(ppB1p)=-ppB1c
     iiPELSINKREF(ppBac)=-ppB1c
     !PhytoPLankton
     do i = 1,iiPhytoPlankton
        j0=ppPhytoPlankton(i,iiC)
        iiPELSINKREF(j0)=ppsediPI(i)
        j=ppPhytoPlankton(i,iiN)
        iiPELSINKREF(j)=-j0
        j=ppPhytoPlankton(i,iiP)
        iiPELSINKREF(j)=-j0
        j=ppPhytoPlankton(i,iiL)
        iiPELSINKREF(j)=-j0
        j=ppPhytoPlankton(i,iiS)
        if (j>0) iiPELSINKREF(j)=-j0
     enddo
     !Extra state var. for Coupled to  to sink rate of Phaeocystis C
     iiPELSINKREF(ppR3c)=-ppP6c
     iiPELSINKREF(ppPcc)=-ppP6c
     !MesoZooPlankton
     do i = 1 , iiMesoZooPlankton
        j0=ppMesoZooPlankton(i,iiC)
        iiPELSINKREF(j0)=ppsediMeZ(i)
        j=ppMesoZooPlankton(i,iiN)
        if (j>0) iiPELSINKREF(j)=-j0
        j=ppMesoZooPlankton(i,iiP)
        if (j>0) iiPELSINKREF(j)=-j0
     enddo
     !MicroZooPlankton
     do i = 1 , iiMicroZooPlankton
        j0=ppMicroZooPlankton(i,iiC)
        iiPELSINKREF(j0)=ppsediMiZ(i)
        j=ppMicroZooPlankton(i,iiN)
        if (j>0) iiPELSINKREF(j)=-j0
        j=ppMicroZooPlankton(i,iiP)
        if (j>0) iiPELSINKREF(j)=-j0
     enddo

     end
