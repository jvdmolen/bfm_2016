#include "DEBUG.h"
#include "INCLUDE.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !MODULE: botflux
!
! DESCRIPTION
!   !

!   This file is geneRAted directly from OpenSesame model code, using a code
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
      module botflux
! !USES:
        USE global_mem,      ONLY:RLEN,ZERO,DONE,NZERO,LOGUNIT
        use mem,only:iiBen,iiPel,D3STATE, &
          NO_BOXES,NO_BOXES_Y,NO_BOXES_X,NO_BOXES_Z,NO_BOXES_XY, &
          BoxNumberZ,BoxNumberX,BoxNumberY,BoxNumberXY,NO_D3_BOX_STATES, &
          PELBOTTOM,flux_vector,flux,PelBoxAbove
        USE BFM_ERROR_MSG, ONLY: set_warning_for_getm,BFM_ERROR
        use bio_var, only: var_names,stPelStates,stBenStateS
        use constants,only: ANY,POSITIVE,INITIALIZE,ADD,NEGATIVE

      real(RLEN),private                       ::r
      integer,private                          ::i,j,k,l
      character(len=20),private                ::t,tsource,tgoal
      character(len=80),private                ::msg
      real(RLEN),public,pointer                ::BENBOUND(:,:)
     
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!  functions
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

     contains

      subroutine initbotflux
        implicit none
        integer,external               ::D2toD1,D3toD1
        BoxNumberZ = 1
        DO BoxNumberY=1,NO_BOXES_Y
          DO BoxNumberX=1,NO_BOXES_X
            BoxNumberXY=D2toD1(BoxNumberX,BoxNumberY)
            k=D3toD1(BoxNumberX,BoxNumberY,BoxNumberZ)
            !link to layer at the sediment water interface
            PelBoxAbove(BoxNumberXY)=k
          enddo
        enddo
!JM        BENBOUND=>D3STATE(:,1:NO_BOXES:NO_BOXES_Z)
        BENBOUND=>D3STATE(1:NO_BOXES_Z,:)
        write(logunit,*) 'PelBoxAbove;',PelBoxAbove
      end subroutine initbotflux

      subroutine setbotflux(mode,n,Depth3D)
      implicit none
      integer,   intent(IN)          ::mode
      integer,   intent(IN)          ::n
      real(RLEN),intent(IN),optional ::Depth3D(n)

      select case (mode)
        case (INITIALIZE)
          PELBOTTOM=ZERO
          !only the first time.....
        case (ADD)
          do j=1,NO_BOXES_XY
            k=PelBoxAbove(j)
            if (k.le.n) then
              r=Depth3D(k)
              if ( r.gt.ZERO) then
                do i=1,NO_D3_BOX_STATES
                  call flux(k,iiPel,i,i, PELBOTTOM(j,i)/r)
                enddo
              endif
            endif
          enddo
      end select
      end subroutine setbotflux

      subroutine checkbotflux(mode,n,s,Depth3D)
      implicit none
      integer,   intent(IN)          ::mode
      integer,   intent(IN)          ::n
      real(RLEN),intent(IN)          ::s
      real(RLEN),intent(IN)          ::Depth3D(n)
      do j=1,NO_BOXES_XY
        k=PelBoxAbove(j)
        if (k.le.n) then
          do i=1,NO_D3_BOX_STATES
            select case (n) 
              case(1) ;r=D3STATE(k,i)*Depth3D(k)
              case(2) 
                l=k+NO_BOXES_Z-1 
                r=NZERO+sum(D3STATE(k:l,i)*Depth3D(k:l))
            end select
            r=abs(PELBOTTOM(j,i)/r)
            if ( r>s ) then     
              tgoal=var_names(i)
              write(LOGUNIT,*)  &
               'Rate of ',trim(tgoal),' at water/sediment interface >',s
              call set_warning_for_getm
            endif
          enddo
        endif
      enddo
      end subroutine checkbotflux


      subroutine openbotflux_vector(mode,Sys,State,Rate_vector,Depth_vector)
      implicit none
      integer,intent(IN)            ::mode
      integer,intent(IN)            ::Sys
      integer,intent(IN)            ::State
      real(RLEN),intent(IN)         ::Rate_vector(NO_BOXES_XY)
      real(RLEN),intent(IN),optional::Depth_vector(NO_BOXES_XY)

      logical                       ::warning=.false.
      if (minval(Rate_vector)<ZERO.and.mode==POSITIVE) then
        warning=.true.;t='negative'
        do i=1,NO_BOXES_XY
           if (Rate_vector(i)<ZERO) then
             write(LOGUNIT,*)  &
                'Error: In ',i,' a ',trim(t),' Rate is found'
           endif
        enddo
      elseif (maxval(Rate_vector)>ZERO.and.mode==NEGATIVE) then
        warning=.true.;t='positive'
        do i=1,NO_BOXES_XY
           if (Rate_vector(i)>ZERO) then
             write(LOGUNIT,*)  &
                'Error: In ',i,' a ',trim(t),' Rate is found'
           endif
        enddo
      endif
      if (warning) then
        tsource=var_names(State)
        write(msg,*) 'SinkFlux from ',trim(tsource)
        call BFM_ERROR("openbotflux_vector",trim(msg))
      endif
      if (Sys.ne.iiPel) &
        call BFM_ERROR("openbotflux", &
              "open flux can only be defined for a pelagic state!")
      if (present(Depth_vector)) then
        do i=1,NO_BOXES_XY
          j=PelBoxAbove(i)
          PELBOTTOM(j,State)=PELBOTTOM(j,State) +&
                              Rate_vector(i)/Depth_vector(i)
        enddo
      else
        do i=1,NO_BOXES_XY
          j=PelBoxAbove(i);
          PELBOTTOM(j,State)=PELBOTTOM(j,State)+Rate_vector(i)
        enddo
      endif
      end subroutine openbotflux_vector

      subroutine addbotflux_vector(mode,SourceSys,Source,GoalSyS,Goal, &
                             Rate_vector,Depth_vector)
      implicit none
      integer,intent(IN)            ::mode
      integer,intent(IN)            ::SourceSys
      integer,intent(IN)            ::Source
      integer,intent(IN)            ::GoalSys
      integer,intent(IN)            ::Goal
      real(RLEN),intent(IN)         ::Rate_vector(NO_BOXES_XY)
      real(RLEN),intent(IN),optional::Depth_vector(NO_BOXES_XY)


      select case (mode)
        case (ANY)       ;l=ZERO;r=ZERO
        case (POSITIVE)  ;l= DONE;r=minval(Rate_vector);t='negative'
      end select
      if (l*r<ZERO )  then
        do i=1,NO_BOXES_XY
          if (l*Rate_vector(i)<ZERO) then
write(LOGUNIT,*) 'Error: In ',i,' a ',trim(t),' Rate is found'
          write(msg,*) 'Error: In ',i,' a ',trim(t),' Rate is found'
          call BFM_ERROR("addbotflux_vector",trim(msg))
          endif
        enddo
        tsource=var_names(Source+(SourceSys/iiBen)*(stBenStateS-1))
!JM        tgoal  =var_names(Goal+(iiBen==GoalSys)*(stBenStateS-1))
        if (iiBen==GoalSys) then
          tgoal  =var_names(Goal+(stBenStateS-1))
        else
          tgoal  =var_names(Goal)
        endif
write(LOGUNIT,*) 'Flux from ',trim(tsource),' to ',trim(tgoal)
        write(msg,*) 'Flux from ',trim(tsource),' to ',trim(tgoal)
        call BFM_ERROR("addbotflux_vector",trim(msg))
      endif
if (SourceSys.eq.GoalSys) &
write(LOGUNIT,*) 'wrong use of routine:SourceSys and GoalSys equal'
      if (SourceSys.eq.GoalSys) &
        call BFM_ERROR( "addbotflux_vector", &
                 'wrong use of routine:SourceSys and GoalSys equal')
      select case (SourceSys)
        case(iiBen)
          call flux_vector(iiBen,Source,Source,-Rate_vector)
          do i=1,NO_BOXES_XY
             j=PelBoxAbove(i)
             PELBOTTOM(j,Goal)=PELBOTTOM(j,Goal)+ Rate_vector(i)
          enddo
        case(iiPel)
          if (present(Depth_vector)) then
            do i=1,NO_BOXES_XY
              j=PelBoxAbove(i)
              PELBOTTOM(j,Source)=PELBOTTOM(j,Source) &
                             - Rate_vector(i)/Depth_vector(i)
            enddo
            call flux_vector(iiBen,Goal,Goal,+Rate_vector/Depth_vector)
          else
            do i=1,NO_BOXES_XY
              j=PelBoxAbove(i)
              PELBOTTOM(j,Source)=PELBOTTOM(j,Source)- Rate_vector(i)
            enddo
            call flux_vector(iiBen,Goal,Goal,+Rate_vector)
          endif
      end select
      end subroutine addbotflux_vector

      subroutine addbotflux(mode,BenBox,SourceSys,Source,&
                          GoalSyS,Goal,Rate,Depth)
      implicit none
      integer,intent(IN)            ::mode
      integer,intent(IN)            ::BenBox
      integer,intent(IN)            ::SourceSys
      integer,intent(IN)            ::Source
      integer,intent(IN)            ::GoalSys
      integer,intent(IN)            ::Goal
      real(RLEN),intent(IN)         ::Rate
      real(RLEN),intent(IN),optional::Depth

      select case (mode)
        case (ANY)       ;l=ZERO;r=ZERO
        case (POSITIVE)  ;l= DONE;t='negative'
      end select
      if ( l*Rate <ZERO )  then
        write(msg,*) 'Error: a ',trim(t),' Rate is found'
        tsource=var_names(Source+(SourceSys/iiBen)*(stBenStateS-1))
        tgoal  =var_names(Goal+(GoalSys/iiBen)*(stBenStateS-1))
        write(LOGUNIT,*) 'Flux from ',trim(tsource),' to ',trim(tgoal)
        call BFM_ERROR('addbotflux',trim(msg));
      endif
      if (SourceSys.eq.GoalSys) &
        call BFM_ERROR( "addbotflux", &
               "wrong use of routine:SourceSys and GoalSys equal")
      select case (SourceSys)
        case(iiBen)
          call flux(BenBox,iiBen,Source,Source,-Rate)
          j=PelBoxAbove(BenBox)
          PELBOTTOM(j,Goal)=PELBOTTOM(j,Goal)+ Rate
        case(iiPel)
          j=PelBoxAbove(BenBox)
          if (present(Depth)) then
            PELBOTTOM(j,Source)=PELBOTTOM(j,Source) - Rate/Depth
            call flux(BenBox,iiBen,Goal,Goal,+Rate/Depth)
          else
            PELBOTTOM(j,Source)=PELBOTTOM(j,Source) - Rate
            call flux(BenBox,iiBen,Goal,Goal,+Rate)
          endif
      end select
      end subroutine addbotflux

      subroutine openbotflux(mode,BenBox,Sys,State,Rate,Depth)
      implicit none
      integer,intent(IN)            ::mode
      integer,intent(IN)            ::BenBox
      integer,intent(IN)            ::Sys
      integer,intent(IN)            ::State
      real(RLEN),intent(IN)         ::Rate
      real(RLEN),intent(IN),optional::Depth
      logical                       ::warning=.false.
      if (Sys==iiBen) then
        call BFM_ERROR( "openbotFlux", &
               "openflux can only be defined for a pelagic state!")
        warning=.true.
      elseif (Rate <ZERO .and. mode==POSITIVE )  then
        warning=.true.
      elseif (Rate >ZERO .and. mode==NEGATIVE )  then
        warning=.true.
      endif
      if ( warning ) then
        write(msg,*)'Error: sign of value is opposite parameter in call'
        tsource=var_names(State)
        write(LOGUNIT,*) 'Sink of ',trim(tsource)
        call BFM_ERROR('openbotflux',trim(msg));
      endif
      j=PelBoxAbove(BenBox)
      if (present(Depth)) then
         PELBOTTOM(j,State)=PELBOTTOM(j,State) + Rate/Depth
      else
         PELBOTTOM(j,State)=PELBOTTOM(j,State) + Rate
      endif
      end subroutine openbotflux

      function getbotflux_3D(mode,Depth3D)
      implicit none
      integer,   intent(IN) ::mode
      real(RLEN),intent(IN) ::Depth3D(NO_BOXES)
      real(RLEN)            ::getbotflux_3D(NO_BOXES)

      getbotflux_3D=ZERO
      do j=1,NO_BOXES_XY
        k=PelBoxAbove(j)
        r=Depth3D(k)
        if (r.gt.ZERO) getbotflux_3D(k)=PELBOTTOM(j,mode)/r
      enddo 
      end function getbotflux_3D
!EOP
!------------------------------------------------------------------------
      end module botflux
