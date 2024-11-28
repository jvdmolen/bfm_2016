#include "DEBUG.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: FindNaNInRates
!
! DESCRIPTION
!
! !INTERFACE
  subroutine FindExtremesInRates(mode,iiSys,ppState,message)
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,LOGUNIT
  use mem,    ONLY: iiBen, iiPel,Source_D2_vector,Source_D3_vector, &
              NO_BOXES_XY,NO_BOXES,D2STATE,D3STATE
#ifdef BFM_GOTM
  use bio_var,ONLY: var_names, stPelStates,stBenStates
#else
  use api_bfm,only: var_names, stPelStates,stBenStates
#endif

  use mem_Param,only:nan_check
  use BFM_ERROR_MSG,only:set_warning_for_getm
!
!
! !AUTHORS
!   ERSEM team	

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
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

  integer,intent(IN)                   ::mode
  integer,intent(IN)                   ::iiSys
  integer,intent(IN)                   ::ppState
  character(len=*)                     ::message

  real(RLEN),dimension(NO_BOXES)    :: r3,s3
  real(RLEN),dimension(NO_BOXES_XY) :: r2,s2

  integer                           :: jout
  real(RLEN), external  :: GetDelta

  jout=0
  if ( .not.nan_check) return
  select case (iiSys)
    case (iiPel)
      r3=Source_D3_vector(ppState,0)
      if (mode==3) s3=D3STATE(ppSTATE,:)+GetDelta()*r3
      select case  (mode)
        case(1) ;call findnan(r3,NO_BOXES,jout)
        case(2) ;call findlarge(r3,NO_BOXES,1.0D+02,jout)
        case(3) ;call findnega(s3,NO_BOXES,jout)
      end select
      if ( jout>0) then
         write(logunit,'(A)') message
         write(logunit,'(A,A,A,I2)') &
                'in rate ',trim(var_names(stPelStateS+ppState-1)),' layer:',jout
         if (mode.eq.2)write(logunit,'(''Rate of '',A,''('',I2,'')='',G13.6)') &
              trim(var_names(stPelStateS+ppState-1)),jout,r3(jout)
         write(logunit,'(A,''('',I2,'')='',G13.6)') &
              trim(var_names(stPelStateS+ppState-1)),jout,D3STATE(ppSTATE,jout)
         call set_warning_for_getm
      endif
    case (iiBen)
      r2=Source_D2_vector(ppState,0)
      if (mode==3) s2=D2STATE(ppSTATE,:)+GetDelta()*r2
      select case  (mode)
        case(1) ;call findnan(r2,NO_BOXES_XY,jout)
        case(2) ;call findlarge(r2,NO_BOXES_XY,1.0D+04,jout)
        case(3) ;call findnega(s2,NO_BOXES_XY,jout)
      end select
      if ( jout>0) then
         write(logunit,'(A)') message
         write(logunit,'(A,A,A,I2)') &
               'in rate ',trim(var_names(stBenStateS+ppState-1)),' layer:',jout
         if (mode.eq.2)write(logunit,'(''Rate of '',A,''('',I2,'')='',G13.6)') &
              trim(var_names(stBenStateS+ppState-1)),jout,r2(jout)
         write(logunit,'(A,''('',I2,'')='',G13.6)') &
              trim(var_names(stBenStateS+ppState-1)),jout,D2STATE(ppSTATE,jout)
         call set_warning_for_getm
      endif
   end select
  end subroutine FindExtremesInRates

  subroutine FindLargeInRates(iiSys,ppState,message)
  IMPLICIT NONE
  integer,intent(IN)                   ::iiSys
  integer,intent(IN)                   ::ppState
  character(len=*)                     ::message
     call FindExtremesInRates(2,iiSys,ppState,message)
  end subroutine FindLargeInRates
  subroutine FindNegaInRates(iiSys,ppState,message)
  IMPLICIT NONE
  integer,intent(IN)                   ::iiSys
  integer,intent(IN)                   ::ppState
  character(len=*)                     ::message
     call FindExtremesInRates(3,iiSys,ppState,message)
  end subroutine FindNegaInRates
  subroutine FindNaNInRates(iiSys,ppState,message)
  IMPLICIT NONE
  integer,intent(IN)                   ::iiSys
  integer,intent(IN)                   ::ppState
  character(len=*)                     ::message
     call FindExtremesInRates(1,iiSys,ppState,message)
  end subroutine FindNaNInRates
 subroutine findnan( vector,n,iout)
      use global_mem, only:RLEN,LOGUNIT
      use BFM_ERROR_MSG,only:set_warning_for_getm

      implicit none
      integer,intent(IN)               :: n
      REAL(RLEN),intent(IN)            :: vector(n)
      integer,intent(OUT)              :: iout

      integer      ::i

          do i=1,n
            if (isnan(vector(i))) then
                write(LOGUNIT,*) '-------------case:NaN----------'
                call set_warning_for_getm
                iout=i
                return
            endif
          enddo

      iout=0
      return
   end subroutine findnan

subroutine findsmall( vector,n,small,iout)
      use global_mem, only:RLEN,LOGUNIT
      use BFM_ERROR_MSG,only:set_warning_for_getm

      implicit none
      integer,intent(IN)               :: n
      REAL(RLEN),intent(IN)            :: vector(n)
      REAL(RLEN),intent(IN)            :: small
      integer,intent(OUT)              :: iout

      integer      ::i

          do i=1,n
            if (abs(vector(i))<small) then
                write(LOGUNIT,*) '-------------case:< 1.0-e8----------'
                call set_warning_for_getm
                iout=i
                return
            endif
          enddo

      iout=0
      return
   end subroutine findsmall


 subroutine findlarge( vector,n,large,iout)
      use global_mem, only:RLEN,LOGUNIT
      use bennut_interface,only:imaxloc
      use BFM_ERROR_MSG,only:set_warning_for_getm

      implicit none
      integer,intent(IN)               :: n
      REAL(RLEN),intent(IN)            :: vector(n)
      REAL(RLEN),intent(IN)            :: large
      integer,intent(OUT)              :: iout

      integer      ::i
      character(len=10)  ::msg


        i=imaxloc(n,abs(vector))
        if (abs(vector(i))>large) then
            write(msg,'(G10.2)') large
            write(LOGUNIT,*) '-------------case:>',msg,'----------'
            call set_warning_for_getm
            iout=i
            return
        endif


      iout=0
      return
   end subroutine findlarge

 subroutine findnega(vector,n,iout)
      use global_mem, only:RLEN,LOGUNIT,ZERO
      use BFM_ERROR_MSG,only:set_warning_for_getm

      implicit none
      integer,intent(IN)               :: n
      REAL(RLEN),intent(IN)            :: vector(n)
      integer,intent(OUT)              :: iout

      integer      ::i

          do i=1,n
            if (vector(i)<ZERO) then
                write(LOGUNIT,*) '-------------case:<0 ----------'
                call set_warning_for_getm
                iout=i
                return
            endif
          enddo

      iout=0
      return
   end subroutine findnega

   subroutine isinf(x,iout)
      use global_mem, only:RLEN,LOGUNIT
      use string_functions, ONLY: getseq_number,empty
      implicit none
      REAL(RLEN),intent(IN)            :: x
      integer,intent(OUT)              :: iout
      character(len=80)                :: mess

      iout=0
      write(mess,'(''x='',G13.6)') X
      if (index('Infinity',mess(1:len_trim(mess)))>0) iout=1
!     write(LOGUNIT,*) mess(1:len_trim(mess)),iout
   end subroutine isinf
!----------------------------------------------------------------------------
! !ROUTINE: FindInfInRates
!
! DESCRIPTION
!   !	This submodel calls all other submodels
!
!

!
! !INTERFACE
  subroutine FindInfInRates(iiSys,ppState,message)
!
! !USES:

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  use global_mem, ONLY:RLEN,LOGUNIT
  use mem,    ONLY: iiBen, iiPel,Source_D2_vector,Source_D3_vector, &
              NO_BOXES_XY,NO_BOXES,D2STATE,D3STATE
#ifdef BFM_GOTM
  use bio_var,ONLY: var_names, stPelStates,stBenStates
#else
  use api_bfm,only: var_names, stPelStates,stBenStates
#endif

  use BFM_ERROR_MSG,only:set_warning_for_getm
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following global functions are used:ResetTotMassVar
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
!
! !AUTHORS
!   ERSEM team	
!
!
!
! !REVISION_HISTORY
!   !
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

  integer,intent(IN)                   ::iiSys
  integer,intent(IN)                   ::ppState
  character(len=*)                     ::message

  real(RLEN),dimension(NO_BOXES)    :: r3
  real(RLEN),dimension(NO_BOXES_XY) :: r2
  integer                           :: jout

  jout=0
  select case (iiSys)
    case (iiPel)
      r3=Source_D3_vector(ppState,0)
      call findinf(r3,NO_BOXES,jout)
      if (jout>0) then
         write(logunit,'(A)') message
         write(logunit,'(A,A,A,I2)') &
                'Inf in rate ',trim(var_names(stPelStateS+ppState-1)),' layer:',jout
         write(logunit,'(A,''('',I2,'')='',G13.6)') &
              trim(var_names(stPelStateS+ppState-1)),jout,D3STATE(ppSTATE,jout)
         call set_warning_for_getm
      endif
    case (iiBen)
      r2=Source_D2_vector(ppState,0)
      call findinf(r2,NO_BOXES_XY,jout)
      if ( jout>0) then
         write(logunit,'(A)') message
         write(logunit,'(A,A,A,I2)') &
               'Inf in rate ',trim(var_names(stBenStateS+ppState-1)),' layer:',jout
         write(logunit,'(A,''('',I2,'')='',G13.6)') &
              trim(var_names(stBenStateS+ppState-1)),jout,D2STATE(ppSTATE,jout)
         call set_warning_for_getm
      endif
   end select
  end subroutine FindInfInRates

 subroutine findinf( vector,n,iout )
      use global_mem, only:RLEN,LOGUNIT
      use BFM_ERROR_MSG,only:set_warning_for_getm

      implicit none
      integer,intent(IN)               :: n
      REAL(RLEN),intent(IN)            :: vector(n)
      integer,intent(OUT)              :: iout

      integer      ::i,jout

      do i=1,n
        call isinf(vector(i),jout)
        if (iout>0) then
            write(LOGUNIT,*) '-------------case:Inf----------'
            call set_warning_for_getm
            iout=jout
            return
        endif
      enddo

      iout=0
      return
   end subroutine findinf



!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
