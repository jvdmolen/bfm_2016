#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !MODULE: SourceFunctions
!
! DESCRIPTION
!   !

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
      module SourceFunctions
! !USES:
          use global_mem, only: RLEN, ZERO,DONE
          use mem,only:D3SINK,D3SOURCE,D2SINK,D2SOURCE
          use mem,only:iiProduction,iiConsumption,iiTotal
          use constants,only:SEC_PER_DAY

!
! !AUTHORS
!   Piet Ruardij   
!
!
! CHANGE_LOG
!   
!
! COPYING
!   
!   Copyright (C) 2004 P. Ruardij, the mfstep group, the ERSEM team 
!   (rua@nioz.nl, vichi@bo.ingv.it)
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation;
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!  functions 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

          contains

!EOP
!-------------------------------------------------------------------------!
!BOP

!
! FUNCTION
!   Source_D3_withgroup
!
! DESCRIPTION
!   
!  
! !INTERFACE

          function Source_D3_withgroup(iistate,ppgroup,lgroup,type,mode)
          implicit none

          integer, intent(IN)             ::iistate
          integer, intent(IN)             ::lgroup
          integer, intent(IN)             ::type
          integer, intent(IN)             ::mode
          interface                         ! Specification
            integer function ppgroup(n,iiC) ! Specification
            integer,intent(IN)   ::n,iiC    ! Specification
            end function                    ! Specification
          end  interface   
!JM          real(RLEN) :: Source_D3_withgroup(size(D3SOURCE,DIM=3))
!JM          real(RLEN) :: fill(size(D3SOURCE,DIM=3))
          real(RLEN) :: Source_D3_withgroup(size(D3SOURCE,DIM=1))
          real(RLEN) :: fill(size(D3SOURCE,DIM=1))
          integer    :: i,j
          real(RLEN) :: l

          fill=ZERO;l=DONE;if ( mode ==iiTotal ) l=-l
          ! Array in sum is by sum seen as 2D-array: DIM=1 and NOT 2
          do i=1,lgroup
            j=ppgroup(i,type)
            if (mode.eq.iiConsumption.or.mode.eq.iiToTal) &
               fill=fill + l* D3SINK(:,iistate,j)
            if (mode.eq.iiProduction.or.mode.eq.iiToTal) &
               fill=fill + D3SOURCE(:,iistate,j)
          enddo
          Source_D3_withgroup=fill*SEC_PER_DAY
        end function Source_D3_withgroup

          function Source_D3_withstate(iistate,jjstate,mode)
          implicit none

          integer, intent(IN)             ::iistate
          integer, intent(IN)             ::jjstate
          integer, intent(IN)             ::mode
!JM          real(RLEN) :: Source_D3_withstate(size(D3SOURCE,DIM=3))
!JM          real(RLEN) :: fill(size(D3SOURCE,DIM=3))
          real(RLEN) :: Source_D3_withstate(size(D3SOURCE,DIM=1))
          real(RLEN) :: fill(size(D3SOURCE,DIM=1))
          real(RLEN) :: l

          fill=ZERO;l=DONE;if ( mode ==iiTotal ) l=-l
          ! Array in sum is by sum seen as 2D-array: DIM=1 and NOT 2
          if (mode.eq.iiConsumption.or.mode.eq.iiToTal) &
               fill=fill + l*D3SINK(:,iistate,jjstate)
          if (mode.eq.iiProduction.or.mode.eq.iiToTal) &
               fill=fill + D3SOURCE(:,iistate,jjstate)
          Source_D3_withstate=fill*SEC_PER_DAY
        end function Source_D3_withstate

          function Source_D2_withgroup(iistate,ppgroup,lgroup,type,mode)

          implicit none

          integer, intent(IN)             ::iistate
          integer, intent(IN)             ::lgroup
          integer, intent(IN)             ::type
          integer, intent(IN)             ::mode
          interface                         ! Specification
            integer function ppgroup(n,iiC) ! Specification
            integer,intent(IN)   ::n,iiC    ! Specification
            end function                    ! Specification
          end  interface   
!JM          real(RLEN) :: Source_D2_withgroup(size(D2SOURCE,DIM=3))
!JM          real(RLEN) :: fill(size(D2SOURCE,DIM=3))
          real(RLEN) :: Source_D2_withgroup(size(D2SOURCE,DIM=1))
          real(RLEN) :: fill(size(D2SOURCE,DIM=1))
          integer    :: i,j
          real(RLEN) :: l

          fill=ZERO;l=DONE;if ( mode ==iiTotal ) l=-l
          ! Array in sum is by sum seen as 2D-array: DIM=1 and NOT 2
          do i=1,lgroup
            j=ppgroup(i,type)
            if (mode.eq.iiConsumption.or.mode.eq.iiToTal) &
               fill=fill + l* D2SINK(:,iistate,j)
            if (mode.eq.iiProduction.or.mode.eq.iiToTal) &
               fill=fill + D2SOURCE(:,iistate,j)
          enddo
          Source_D2_withgroup=fill*SEC_PER_DAY
        end function Source_D2_withgroup

          function Source_D2_withstate(iistate,jjstate,mode)
          implicit none

          integer, intent(IN)             ::iistate
          integer, intent(IN)             ::jjstate
          integer, intent(IN)             ::mode
!JM          real(RLEN) :: Source_D2_withstate(size(D2SOURCE,DIM=3))
!JM          real(RLEN) :: fill(size(D2SOURCE,DIM=3))
          real(RLEN) :: Source_D2_withstate(size(D2SOURCE,DIM=1))
          real(RLEN) :: fill(size(D2SOURCE,DIM=1))
          real(RLEN) :: l

          fill=ZERO;l=DONE;if ( mode ==iiTotal ) l=-l
          ! Array in sum is by sum seen as 2D-array: DIM=1 and NOT 2
          if (mode.eq.iiConsumption.or.mode.eq.iiToTal) &
               fill=fill + l*D2SINK(:,iistate,jjstate)
          if (mode.eq.iiProduction.or.mode.eq.iiToTal) &
               fill=fill + D2SOURCE(:,iistate,jjstate)
          Source_D2_withstate=fill*SEC_PER_DAY
        end function Source_D2_withstate

      end module SourceFunctions
