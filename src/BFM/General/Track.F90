#include "DEBUG.h"


!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !MODULE: Track
!
! DESCRIPTION
!   Here are the values for the benthic pelagic coupling.
!
!

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
     module Track
!
! !USES:

     use bfm_output, only: var_names
     use bfm_output, only: stPelStateS,stBenStateS
     use BFM_ERROR_MSG, ONLY: BFM_ERROR
     use global_mem, only: LOGUNIT,ALLTRANSPORT,ZERO
     use constants, only: RLEN, ZERO, SEC_PER_DAY
     use mem_Param,  ONLY: p_small,p_check_track
     use mem, ONLY: flux,iiTrack,ii3dTrack,nr_3d_track,ii3dptTrack, &  
          flag_3d_track_bot, fix_3d_track_bot, check_3d_track_bot, &
          ii2dTrack, nr_2d_track,fix_2d_track_bot,ii2dptTrack, & 
          ii3daptTrack

      character(len=80)                :: message
      integer,dimension(:),allocatable :: check_ben_to_pel
!  
!
! !AUTHORS
!   Original version by  P. Ruardij
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
!   the Free Software Foundation;
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
  ! SHARED PUBLIC FUNCTIONS (must be explicited below "contains")

  public define_track_rate,set_exchange_benpel,calc_all_track_rates, &
         benpelcoup_track,set_track_rates,calculate_labeled_fraction,&
         calculate_average_labeled_fraction,fill_sfl_if_tracking,&
         inittransportstatetypes_track,set_all_exchanges_for_tracking_for_your_model,&
         check_track_states,check_benthic_states_on_negative_values
  contains
!----------------------------------------------------------------------------
     
      subroutine calc_all_track_rates(mode)
        use mem,ONLY:  D3SOURCE,D3SINK,D3STATE,NO_D3_BOX_STATES, &
                       D2SOURCE,D2SINK,D2STATE,NO_D2_BOX_STATES, PELBOTTOM, &
                       D3DIAGNOS,D2DIAGNOS ,&
                       D3STATETYPE,D2STATETYPE
        use mem,ONLY:  Depth,NO_BOXES_Y,NO_BOXES_X,NO_BOXES,NO_BOXES_XY ,&
                       NO_D3_BOX_DIAGNOSS,NO_D2_BOX_DIAGNOSS
        use mem, ONLY: iiPel,iiBen

        implicit none
        integer,intent(IN)             ::mode

        integer     :: BoxNumber
        integer     :: BoxNumberZ
        integer     :: BoxNumberY
        integer     :: BoxNumberX
        integer     :: BoxNumberXY

        real(RLEN),dimension(NO_BOXES) :: r,s
        real(RLEN),dimension(NO_D3_BOX_STATES) :: ld3state
        real(RLEN),dimension(NO_D2_BOX_STATES) :: ld2state
        real(RLEN),dimension(NO_D3_BOX_STATES, NO_D3_BOX_STATES) :: ld3sink
        real(RLEN),dimension(NO_D2_BOX_STATES, NO_D2_BOX_STATES) :: ld2sink

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! user defined external functions
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        integer, external  :: D3toD1
        integer, external  :: D2toD1
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        !iiTrack==1 : there is tracking defined however there is no initialization
        !iiTrack==2 : ther is tracking and initialization
        !mode ==0 forced intializtion of the states.
        !mode > 0 no initialization

        if ( iiTrack==1.or.mode ==0) then
          if (iiTrack==0) return;
          write(LOGUNIT,*) 'D3_init_track_states'
          call inittransportstatetypes_track(3,D3STATE,D3STATETYPE, NO_D3_BOX_STATES,  &
                                             NO_BOXES,nr_3d_track,ii3dTrack)
           write(LOGUNIT,*) 'D2_init_track_states'
          call inittransportstatetypes_track(2,D2STATE,D2STATETYPE, NO_D2_BOX_STATES,  &
                                             NO_BOXES_XY,nr_2d_track,ii2dTrack)
          iiTrack=2
        endif
        if ( mode.eq.0) return
        ! For initialization of for set_sedimentation_scheme_tracking first some
        ! variables mut hvae get their values. This is done the first time in PelGlobal.F90
        ! Therfor the initialization for tracking take only place after 
        ! the calculation of the first time of the rates.
        ! Thies means that the tracking and calculation of the rate of cahange in the tracking 
        ! start only at step 2.
        if ( iiTrack==1.or.iiTrack==2) then
           write(LOGUNIT,*) 'set all exchanges for trackings'
          call set_all_exchanges_for_tracking_for_your_model
           write(LOGUNIT,*) 'set sedimentation scheme for tracking'
          call  set_sedimentation_scheme_tracking

          iiTrack=3 ; !if this routine is preformed all initalizations is done.
        elseif ( iiTrack ==3 ) then

          call check_track_states(10,'at start of calculations for tracking')

          !write(LOGUNIT,*) 'D3_track_frac'
          call calculate_labeled_fraction(D3DIAGNOS, D3STATE, NO_D3_BOX_STATES,&
                                       NO_D3_BOX_DIAGNOSS,NO_BOXES,nr_3d_track,ii3dTrack,ii3dptTrack)
          !write(LOGUNIT,*) 'D2_track_frac'
          call calculate_labeled_fraction(D2DIAGNOS, D2STATE, NO_D2_BOX_STATES,&
                                       NO_D2_BOX_DIAGNOSS,NO_BOXES_XY,nr_2d_track,ii2dTrack,ii2dptTrack)
      
!         if ( track_error ==1 ) then
!            write(LOGUNIT,*) 'Before set_track_rates'
!            r=Source_D3_Vector(ppN3n,0)/(1.0D-80+N3n)
!            s=Source_D3_Vector(pptrN3n,0)/(1.0D-80+trN3n)
!            write(LOGUNIT,*) "Compare relrates N3n/trN3n:",r(25),s(25)
!         endif

          !write(LOGUNIT,*) 'D3_track_rates'
          call set_track_rates(D3SOURCE,D3SINK,D3STATE,NO_D3_BOX_STATES,NO_BOXES,&
                                              nr_3d_track,ii3dTrack)
          !write(LOGUNIT,*) 'D2_track_rates'
          call set_track_rates(D2SOURCE,D2SINK,D2STATE,NO_D2_BOX_STATES,NO_BOXES_XY,&
                                              nr_2d_track,ii2dTrack)
!         if ( track_error ==1 ) then
!            write(LOGUNIT,*) 'After set_track_rates'
!            r=Source_D3_Vector(ppN3n,0)/(1.0D-80+N3n)
!            s=Source_D3_Vector(pptrN3n,0)/(1.0D-80+trN3n)
!            write(LOGUNIT,*) "Compare relrates N3n/trN3n:",r(25),s(25)
!         endif

          BoxNumberZ = 1
          DO BoxNumberY=1,NO_BOXES_Y
            DO BoxNumberX=1,NO_BOXES_X
              BoxNumber=D3toD1(BoxNumberX,BoxNumberY,BoxNumberZ)
              BoxNumberXY=D2toD1(BoxNumberX,BoxNumberY)

              ld3state= D3STATE(:,BoxNumber)
              ld2state =D2STATE(:,BoxNumberXY)
              ld3sink =D3SINK(:,:,BoxNumber)
              ld2sink =D2SINK(:,:,BoxNumberXY)
              !write(LOGUNIT,*) 'D2_track_bencoup'
              call benpelcoup_track( BoxNumber,BoxNumberXY,ld3state,NO_D3_BOX_STATES, &
                          ld2state,NO_D2_BOX_STATES, &
                          PELBOTTOM(1:NO_D3_BOX_STATES,BoxNumber),Depth(BoxNumber))

              call definitive_loss_track(iiPel,BoxNumber,ld3state,ld3sink,NO_D3_BOX_STATES,NO_D3_BOX_STATES, &
                          PELBOTTOM(1:NO_D3_BOX_STATES,BoxNumber),Depth(BoxNumber))
              call definitive_loss_track(iiBen,BoxNumberXY,ld2state,ld2sink,NO_D2_BOX_STATES,NO_D3_BOX_STATES, &
                          PELBOTTOM(1:NO_D3_BOX_STATES,BoxNumber),Depth(BoxNumber))

              call calculate_average_labeled_fraction(BoxNumberXY,D2DIAGNOS, D3STATE, NO_D3_BOX_STATES,&
                                 NO_D3_BOX_DIAGNOSS,NO_BOXES, Depth,nr_3d_track,ii3dTrack,ii3daptTrack)

            enddo
          enddo
          call check_track_states(100,'at end of calulations for tracking')
          call calc_track_rates_for_your_model

        endif
      end subroutine calc_all_track_rates
!----------------------------------------------------------------------------
      subroutine set_sedimentation_scheme_tracking
         use mem,only: iiPelSINKREF,NO_D3_BOX_STATES

         implicit none
         integer          ::ml
         integer          ::mt
         integer          ::j

         do mt=1,NO_D3_BOX_STATES-ii3dTrack
               ml=nr_3d_track(mt)
               if ( ml > 0 ) then
                  j=iiPELSINKREF(mt)
                  if ( j> 0 ) then
                     iiPELSINKREF(ml)=-mt
                  elseif ( j<0) then
                     iiPELSINKREF(ml)=j
                  endif
               endif
         enddo
      end subroutine set_sedimentation_scheme_tracking

!----------------------------------------------------------------------------
      subroutine define_track_bot_rate(iiSys,iiState,BoxNumber,BoxNumberXY,jjSys,jjState,rate)
        use mem,ONLY:  D3SOURCE,D3SINK,D3STATE,NO_D3_BOX_STATES, &
                       D2SOURCE,D2SINK,D2STATE,NO_D2_BOX_STATES
        use mem,ONLY:  Depth,NO_BOXES_Y,NO_BOXES_X,NO_STATES

        integer,intent(IN)            ::iiSys
        integer,intent(IN)            ::jjSys
        integer,intent(IN)            ::iiState
        integer,intent(IN)            ::jjState
        real(RLEN),intent(IN)         ::rate
        integer,intent(IN)            ::BoxNumber
        integer,intent(IN)            ::BoxNumberXY

        real(RLEN),dimension(NO_D3_BOX_STATES) :: ld3state
        real(RLEN),dimension(NO_D2_BOX_STATES) :: ld2state

!JM        ld3state= D3STATE(:,BoxNumber)
!JM        ld2state =D2STATE(:,BoxNumberXY)
        ld3state= D3STATE(BoxNumber,:)
        ld2state =D2STATE(BoxNumberXY,:)
        call define_track_bot_rate_interface(BoxNumber,BoxNumberXY, &
                     ld3state,NO_D3_BOX_STATES,ld2state,NO_D2_BOX_STATES, &
                     Depth(BoxNumber),iiSys,iiState,jjSys,jjState,rate)
      end subroutine define_track_bot_rate
!----------------------------------------------------------------------------
      subroutine define_track_bot_rate_interface(BoxNumber,BoxNumberXY, &
                                  D3STATE,N3,D2STATE,N2,Depth,iiSys,iiState,jjSys,jjState,rate)
      use mem, ONLY: iiPel,iiBen
      implicit none

      integer,intent(IN)            ::BoxNumber
      integer,intent(IN)            ::BoxNumberXY
      integer,intent(IN)           :: N3
      integer,intent(IN)           :: N2
      real(RLEN),intent(IN)        :: D3STATE(1:N3)
      real(RLEN),intent(IN)        :: D2STATE(1:N2)
      real(RLEN),intent(IN)        :: Depth
      integer,intent(IN)            ::iiSys
      integer,intent(IN)            ::iiState
      integer,intent(IN)            ::jjSys
      integer,intent(IN)            ::jjState
      real(RLEN),intent(IN)         ::rate
        
      integer                      :: ml
      integer                      :: mls
      integer                      :: m3=0
      integer                      :: ml3
      integer                      :: m2=0
      integer                      :: ml2
      integer                      :: m3s
      real(RLEN)                   :: trackrate

      if ( iiTrack<=2) return;
      if ( iiSys==iiPel .and. rate <= ZERO .and. &
                      iiState<=N3-ii3dTrack ) then
           ml=nr_3d_track(iiState);if ( ml<=0 ) return;
           m2=fix_3d_track_bot(iiState)
           if ( jjState > 0 ) then
             ml2=nr_2d_track(jjState);
             if (jjSYS==iiPel) then 
                write(LOGUNIT,'(''Error in Benthic-Pelagic flux for '',A  )') &
                                trim(var_names(stPelStateS-1+iiState))
                message="iiSys=jjSys=jjPel rate < 0.0"
                call BFM_ERROR("bencoup_track",message)
             endif
           endif
           if ( m2<=0 .or. m2== jjState) then
             flag_3d_track_bot(iiState)=1
             check_3d_track_bot(iiState)=check_3d_track_bot(iiState)+rate
           endif
           trackrate=rate* min(1.0,(p_small+D3STATE(ml))/(p_small+ D3STATE(iiState)))
           call flux(BoxNumber,iiPel,ml,ml,trackrate/Depth)
           if ( jjState.gt.0) then
             call flux(BoxNumberXY,iiBen,ml2,ml2,-trackrate)
           endif
      elseif ( iiSys==iiBen .and. rate > ZERO .and. &
                      iiState<=N2-ii2dTrack ) then
           ml=nr_2d_track(iiState); if (ml==0) return;
           m2=fix_3d_track_bot(jjState)
           if (jjSYS==iiBen) then 
              stop 'iiSys=jjSys=jjBen rate > 0.0'
           elseif ( m2<=0 .and.ml.gt.0 )  then
               flag_3d_track_bot(jjState)=1
               check_3d_track_bot(jjState)=check_3d_track_bot(jjState)+rate
           endif
           mls=nr_3d_track(jjState);
           trackrate=rate*  min(1.0,(p_small+D2STATE(ml))/(p_small+ D2STATE(iiState)))
           call flux(BoxNumberXY,iiBen,ml,ml,-trackrate)
           call flux(BoxNumber,iiPel,mls,mls, trackrate/Depth)
      elseif ( rate .ne.ZERO ) then
            !error : wrong combination of rates
      endif
      end subroutine define_track_bot_rate_interface
!----------------------------------------------------------------------------

      subroutine set_exchange_states_benpel(iiPelState,iiBenState,ben_check)

         implicit none

         integer,intent(IN)              :: iiPelState
         integer,intent(IN)              :: iiBenState
         integer,intent(IN),optional     :: ben_check

         !  This fix_*d_track_bot are only filled (one time of more) before
         !  the first time the labeled rates are called.
         if ( iiPelState <=0 ) then
            message="iiPelState<=0"
            call BFM_ERROR("set_exchange_states_benpel",message)
         endif
         if ( iiBenState ==0 ) then
            message="iiBenState==0"
            call BFM_ERROR("set_exchange_states_benpel",message)
         endif
                              fix_3d_track_bot(iiPelState)=iiBenState
         if ( iiBenState> 0 ) then
                 fix_2d_track_bot(iiBenState)=iiPelState
                 check_ben_to_pel(iiBenState)=1
                 if (present(ben_check)) check_ben_to_pel(iiBenState)=ben_check
         endif

      end subroutine set_exchange_states_benpel
!----------------------------------------------------------------------------

      subroutine definitive_loss_track(iiSys,BoxNumber,STATE,SINK,N,N3,PELBOTTOM,Depth)
      use mem, ONLY: iiPel,iiBen

      implicit none
      integer,intent(IN)           :: iiSys
      integer,intent(IN)           :: BoxNumber
      integer,intent(IN)           :: N
      integer,intent(IN)           :: N3
      real(RLEN),intent(IN)        :: STATE(1:N)
      real(RLEN),intent(IN)        :: SINK(1:N,1:N)
      real(RLEN),intent(INOUT)     :: PELBOTTOM(1:N3)
      real(RLEN),intent(IN)        :: Depth

      integer                      :: mt,ml,m3
      real(RLEN)                   :: dev_loss

      select case (iiSys)
        case (iiPel)
          do mt=1,N-ii3dTrack
             ml=nr_3d_track(mt);if (ml.gt.0) then
               if ( PELBOTTOM(mt) < ZERO ) then
                 dev_loss=(SINK(mt,mt)*SEC_PER_DAY +PELBOTTOM(mt)/Depth) *min(1.0,(1.0D-80+STATE(ml))/(1.0D-80+STATE(mt)))
                  if (dev_loss .gt.1.0D-04 ) then
                     call flux(BoxNumber,iiPel,ml,ml,  -dev_loss)
                  elseif (dev_loss .lt.-1.0D-04 ) then
                    write(LOGUNIT,'(''Error in Pelagic definitive loss flux for '',A  )') &
                              trim(var_names(stPelStateS-1+mt))
                    write(LOGUNIT,'(''definitive loss  rate var (/m2/day):'',G13.5)') dev_loss*Depth
                    write(LOGUNIT,'(''SINK of this state var (/m2/day):'',G13.5)') SINK(mt,mt)*Depth*SEC_PER_DAY
                    write(LOGUNIT,'(''flux to bottom  of this state var (/m2/day):'',G13.5)') PELBOTTOM(mt)
                    write(LOGUNIT,'(''value of this state var (/m3):'',G13.5)') STATE(mt)
                    write(LOGUNIT,'(''value of '',A,''  (/m3):'',G13.5)') & 
                              trim(var_names(stPelStateS-1+ml)), STATE(ml) 
                    call BFM_ERROR("definitive_loss_track","SINK is smaller then bottom-water flux")
                  endif
               endif
             endif
          enddo
        case (iiBen)
          do mt=1,N-ii2dTrack
            ml=nr_2d_track(mt);if (ml.gt.0) then
              m3=fix_2d_track_bot(mt)
              if ( m3 .gt.0)  then
                if ( PELBOTTOM(m3) > ZERO ) then
                  dev_loss=(SINK(mt,mt)*SEC_PER_DAY -PELBOTTOM(m3)) *min(1.0,(1.0D-80+STATE(ml))/(1.D0-80+STATE(mt)))
                  if ( check_ben_to_pel(mt).eq.0) then
                      call flux(BoxNumber,iiBen,ml,ml,  -dev_loss)
                  elseif (dev_loss .gt.1.0D-04 ) then
                      call flux(BoxNumber,iiBen,ml,ml,  -dev_loss)
                  elseif (dev_loss .lt.-1.0D-04 ) then
                     write(LOGUNIT,'(''Error in Benthic definitive loss flux for '',A  )') &
                               trim(var_names(stBenStateS-1+mt))
                    write(LOGUNIT,'(''definitive loss  rate var (/m2/day):'',G13.5)') dev_loss
                    write(LOGUNIT,'(''SINK of this state var (/m2/day):'',G13.5)') SINK(mt,mt)*SEC_PER_DAY
                    write(LOGUNIT,'(''flux to bottom  of this state var (/m2/day):'',G13.5)') PELBOTTOM(m3)
                    write(LOGUNIT,'(''value of this state var (/m2):'',G13.5)') STATE(mt)
                    write(LOGUNIT,'(''value of '',A,''  (/m2):'',G13.5)')  &
                              trim(var_names(stBenStateS-1+ml)),STATE(ml)
                     call BFM_ERROR("definitive_loss_track","SINK is smaller then bottom-water flux")
                  endif
                endif
              endif
            endif
          enddo
      end select
      return
      end subroutine definitive_loss_track

      subroutine benpelcoup_track(BoxNumber,BoxNumberXY,D3STATE,N3,D2STATE,N2, PELBOTTOM,Depth)
      use mem, ONLY: iiPel,iiBen
      implicit none

      integer,intent(IN)           :: BoxNumber
      integer,intent(IN)           :: BoxNumberXY
      integer,intent(IN)           :: N3
      integer,intent(IN)           :: N2
      real(RLEN),intent(IN)        :: D3STATE(1:N3)
      real(RLEN),intent(IN)        :: D2STATE(1:N2)
      real(RLEN),intent(INOUT)     :: PELBOTTOM(1:N3)
      real(RLEN),intent(IN)        :: Depth

      integer                      ::mt,ml,mf,m3,m2,ml3,ml2
      real(RLEN)                   :: r

      do mt=1,N3-ii3dTrack
        ml=nr_3d_track(mt);if (ml >0) then
          mf=flag_3d_track_bot(mt)
          m2=fix_3d_track_bot(mt)
          ! mf==1: instead of automatic assigning rates for labeled states
          ! it is controlled by the routine define_track_bot_rate
          if ( mf.ne.1 .and. m2== -1) then
              write(LOGUNIT,'(''Error in Benthic-Pelagic flux for '',A  )') &
                                trim(var_names(stPelStateS-1+mt))
              call BFM_ERROR("bencoup_track","No tracked fluxes are defined")
          elseif ( m2 == -1) then
            r= abs(check_3d_track_bot(mt)- PELBOTTOM(mt))/(p_small + abs(PELBOTTOM(mt)))
            if ( r .gt.1.0D-04 ) then
              write(LOGUNIT,'(''Error in Benthic-Pelagic flux for '',A  )') &
                                trim(var_names(stPelStateS-1+mt))
              write(LOGUNIT,*) "PELBOTTOM=",PELBOTTOM(mt)
              write(LOGUNIT,*) "check==",check_3d_track_bot(mt)
              call BFM_ERROR("bencoup_track","Not all tracked fluxes are defined")
            endif
          elseif ( mf==1 ) then
            ! distribution amoung sources of these fluxes are already explicity 
            ! defined with define_track_bot_rate
            r= abs(check_3d_track_bot(mt)- PELBOTTOM(mt))/(p_small + abs(PELBOTTOM(mt)))
            if ( r .gt.1.0D-04 ) then
              write(LOGUNIT,'(''Error in Benthic-Pelagic flux: '',A, '' with '',A )') &
                                trim(var_names(stPelStateS-1+mt)),trim(var_names(stBenStateS-1+m2))
              write(LOGUNIT,*) "PELBOTTOM=",PELBOTTOM(mt)
              write(LOGUNIT,*) "check==",check_3d_track_bot(mt)
              write(message,*) 'Wrong use of define_track_bot_rate:not all sources are known!'
              call BFM_ERROR("bencoup_track",message)
            endif
          endif
          if (PELBOTTOM(mt) .ne. ZERO .and.  m2==0)  then
              ! flux found for which no destination is found!
              write(LOGUNIT,'(''Error in Benthic-Pelagic flux: '',A)') &
                              trim(var_names(stPelStateS-1+mt))
              write(message,*) 'Missing rate pelagic/rate definition'
              call BFM_ERROR("bencoup_track",message)
          elseif (PELBOTTOM(mt) < ZERO .and. mf==0 )  then
            ! automatic calculation of the labeled concentration
            PELBOTTOM(ml)=PELBOTTOM(mt)*min(1.0,(p_small+D3STATE(ml))/(p_small+ D3STATE(mt)))
            call flux(BoxNumber,  iiPel,ml,ml,   PELBOTTOM(ml)/Depth)
            ml2=nr_2d_track(m2);
            call flux(BoxNumberXY,iiBen,ml2,ml2,-PELBOTTOM(ml))
          endif
        endif
      enddo
      do mt=1,N2-ii2dTrack
        ml=nr_2d_track(mt);if (ml.ne.0) then
          m3=fix_2d_track_bot(mt)
          if ( m3 .gt.0)  then
            mf=flag_3d_track_bot(m3)
            ml3=nr_3d_track(m3)
            !only automatic  (mf==0)
            if ( PELBOTTOM(m3) > ZERO .and.  mf ==0 ) then
                PELBOTTOM(ml3)=PELBOTTOM(m3)*min(1.0,(p_small+D2STATE(ml))/(p_small+ D2STATE(mt)))
               call flux(BoxNumber,  iiPel,ml3,ml3, PELBOTTOM(ml3)/Depth)
               call flux(BoxNumberXY,iiBen,ml,ml,  -PELBOTTOM(ml3))
            endif
          endif
        endif
      enddo

      check_3d_track_bot=ZERO;
      end subroutine benpelcoup_track

      subroutine check_track_states(mode,where)
      use mem, ONLY: iiPel,iiBen
      use mem,only:D3STATE,NO_D3_BOX_STATES,NO_BOXES,nr_3d_track,ii3dTrack,&
              D2STATE,NO_D2_BOX_STATES,NO_BOXES_XY,  nr_2d_track,ii2dTrack

        implicit none
        integer,intent(IN)           :: mode
        character(len=*),intent(IN)  :: where

        integer                      ::m
        integer                      ::k
        logical                      ::l=.false.
 
        if (iiTrack==0) return;
        if ( p_check_track==0 .and. mode.gt.0 ) return
        if ( mode > 0 ) then
          k=p_check_track; m= k/mode;l=(m==1)
          if (m>1) then    
             k=k-100 ;m= k/mode; l=(m==1)
          endif
          if (m>1) then    
             k=k-10; l=(k==1)
          endif
          if ( .not.l )  return
        endif

        !write(LOGUNIT,*) 'D3_check_track_rates, mode =',mode,' ',where
        if (mode.eq.-1 .or. l) &
        call check_d_track_states(iiPel,where,D3STATE,NO_D3_BOX_STATES,NO_BOXES,  &
                                              nr_3d_track,ii3dTrack)
        !write(LOGUNIT,*) 'D2_track_rates'
        if (mode.eq.-2 .or. l) &
        call check_d_track_states(iiBen,where,D2STATE,NO_D2_BOX_STATES,NO_BOXES_XY,  &
                                              nr_2d_track,ii2dTrack)
      end subroutine check_track_states

      subroutine check_d_track_states(iiSys,where,STATE, n_states,n_boxes,nr_track,n)
        use mem, ONLY: iiPel,iiBen

        implicit none
        integer,intent(IN)           :: iiSys
        integer,intent(IN)           :: n_states
        integer,intent(IN)           :: n_boxes
        integer,intent(IN)           :: n
        integer,intent(IN)           :: nr_track(1:n)
        character(len=*),intent(IN)  :: where
!JM        real(RLEN),intent(IN)        :: STATE(1:n_states,1:n_boxes)
        real(RLEN),intent(IN)        :: STATE(1:n_boxes,1:n_states)

        integer                      ::i
        integer                      ::mt
        integer                      ::ml
        integer                      ::n_rstates
        integer                      ::start
        real(RLEN),dimension(n_boxes) :: t,st

        n_rstates=n_states-n
        do mt=1,n_rstates
          ml=nr_track(mt); 
          if (ml.gt.0) then
!JM            t(1:n_boxes)= (STATE(mt,1:n_boxes) - STATE(ml,1:n_boxes)) 
!JM            st(1:n_boxes)=t(1:n_boxes) /(1.D-80 +STATE(mt,1:n_boxes))
            t(1:n_boxes)= (STATE(mt,1:n_boxes) - STATE(1:n_boxes,ml)) 
            st(1:n_boxes)=t(1:n_boxes) /(1.D-80 +STATE(1:n_boxes,ml))
            if (minval(st(1:n_boxes)) .lt. -1.00.and.&
                          minval(t(1:n_boxes)).lt.-2.0D-04) then     
              do i=1,n_boxes
                if ( st(i) .lt.-1.00.and.t(i).lt.-2.0E-04 ) then
                   start=0;;if (iiSys.eq.iiBen) start=stBenStateS-1
                   write(LOGUNIT,'('' layer='',i4)') i
                   write(LOGUNIT,*) 'The call of this routine is ', where
                   write(LOGUNIT,'('' Error: conc. of '',A,'' larger than '',A)') &
                           trim(var_names(start+ml)), trim(var_names(start+mt))
                   write(LOGUNIT,*) "values :", STATE(ml,i),STATE(mt,i) 
                   write(message,*) 'Wrong by definition'
                   call BFM_ERROR("check_track_states",message)
                endif
              enddo
            endif
          endif
        enddo
      end subroutine check_d_track_states

      subroutine set_track_rates(SOURCE,SINK,STATE, n_states,n_boxes,nr_track,n)
!       use mem,only:track_error,ppN3n
        implicit none

        integer,intent(IN)           :: n_states
        integer,intent(IN)           :: n_boxes
        integer,intent(IN)           :: n
        integer,intent(IN)           :: nr_track(1:n_states)
!JM        real(RLEN),intent(INOUT)     :: SOURCE(1:n_states,1:n_states,1:n_boxes)
!JM        real(RLEN),intent(INOUT)     :: SINK(1:n_states,1:n_states,1:n_boxes)
!JM        real(RLEN),intent(IN)        :: STATE(1:n_states,1:n_boxes)
        real(RLEN),intent(INOUT)     :: SOURCE(1:n_boxes,1:n_states,1:n_states)
        real(RLEN),intent(INOUT)     :: SINK(1:n_boxes,1:n_states,1:n_states)
        real(RLEN),intent(IN)        :: STATE(1:n_boxes,1:n_states)

        integer                      ::jt
        integer                      ::jl
        integer                      ::mt
        integer                      ::ml
        integer                      ::n_rstates
        real(RLEN),dimension(n_boxes) :: p
        real(RLEN),dimension(n_boxes) :: t

        real(RLEN),dimension(n_boxes) :: r,s
        real(RLEN),dimension(n_boxes) :: t2
        real(RLEN)                    :: old
      
        r=ZERO;s=ZERO;old=ZERO;
        ! Calculate number of "real" states : total number of states -number of labeldstates
        n_rstates=n_states-n
        do mt=1,n_rstates
           ml=nr_track(mt); 
           if (ml.gt.0) then
             do jt=1,n_rstates
                jl= nr_track(jt)
                if (jl.gt.0.and.ml.ne.jl)  then
!JM                  p(1:n_boxes)=min(1.0,(p_small+STATE(ml,1:n_boxes))/(p_small+STATE(mt,1:n_boxes)))
!JM                  t=SOURCE(jt,mt,1:n_boxes); t=t*p; 
!JM                  SOURCE(jl,ml,1:n_boxes)=t(1:n_boxes) ;
!JM                  t=  SINK(mt,jt,1:n_boxes) ; t=t*p; 
!JM                  SINK(ml,jl,1:n_boxes)  = t(1:n_boxes);
                  p(1:n_boxes)=min(1.0,(p_small+STATE(1:n_boxes,ml))/(p_small+STATE(1:n_boxes,mt)))
                  t=SOURCE(1:n_boxes,jt,mt); t=t*p; 
                  SOURCE(1:n_boxes,jl,ml)=t(1:n_boxes) ;
                  t=  SINK(1:n_boxes,mt,jt) ; t=t*p; 
                  SINK(1:n_boxes,ml,jl)  = t(1:n_boxes);
                endif
             enddo
           endif
        enddo
      end subroutine set_track_rates


      subroutine calculate_labeled_fraction(DIAGNOS, STATE, n_states,n_diags,n_boxes,nr_track,n,start)

        implicit none
        integer,intent(IN)           :: n_states
        integer,intent(IN)           :: n_diags
        integer,intent(IN)           :: n_boxes
        integer,intent(IN)           :: n
        integer,intent(IN)           :: nr_track(1:n_states)
        integer,intent(IN)           :: start
!JM        real(RLEN),intent(INOUT)     :: DIAGNOS(1:n_diags,1:n_boxes)
!JM        real(RLEN),intent(IN)        :: STATE(1:n_states,1:n_boxes)
        real(RLEN),intent(INOUT)     :: DIAGNOS(1:n_boxes,1:n_diags)
        real(RLEN),intent(IN)        :: STATE(1:n_boxes,1:n_states)

        integer                      ::j
        integer                      ::mt
        integer                      ::ml
        integer                      ::n_rstates

        ! Calculate number of "real" states : total number of states -number of labeldstates
        n_rstates=n_states-n;j=start
        do mt=1,n_rstates
           ml=nr_track(mt); 
           if (ml > 0) then
!JM              DIAGNOS(j,:)=STATE(ml,:)/(p_small+STATE(mt,:))
              DIAGNOS(:,j)=STATE(:,ml)/(p_small+STATE(:,mt))
              j=j+1
           endif
        enddo
      end subroutine calculate_labeled_fraction

      subroutine calculate_average_labeled_fraction(BoxNumberXY,DIAGNOS,STATE, &
                                          n_states,n_diags,n_boxes,Depth,nr_track,n,start)

        implicit none
        integer,intent(IN)           :: BoxNumberXY
        integer,intent(IN)           :: n_states
        integer,intent(IN)           :: n_diags
        integer,intent(IN)           :: n_boxes
!JM        real(RLEN),intent(INOUT)     :: DIAGNOS(1:n_diags,1:n_boxes)
!JM        real(RLEN),intent(IN)        :: STATE(1:n_states,1:n_boxes)
        real(RLEN),intent(INOUT)     :: DIAGNOS(1:n_boxes,1:n_diags)
        real(RLEN),intent(IN)        :: STATE(1:n_boxes,1:n_states)
        real(RLEN),intent(IN)        :: Depth(1:n_boxes)
        integer,intent(IN)           :: n
        integer,intent(IN)           :: nr_track(1:n_states)
        integer,intent(IN)           :: start

        integer                      ::j
        integer                      ::mt
        integer                      ::ml
        integer                      ::n_rstates

        ! Calculate number of "real" states : total number of states -number of labeldstates
        n_rstates=n_states-n;j=start
        do mt=1,n_rstates
           ml=nr_track(mt); 
           if (ml > 0) then
!JM              DIAGNOS(j,BoxNumberXY)=sum(DEPTH(:)*STATE(ml,:))/(p_small+sum(Depth(:)*STATE(mt,:)))
              DIAGNOS(BoxNumberXY,j)=sum(DEPTH(:)*STATE(:,ml))/(p_small+sum(Depth(:)*STATE(:,mt)))
              j=j+1
           endif
        enddo
      end subroutine calculate_average_labeled_fraction


      subroutine inittransportstatetypes_track(mode,STATE,STATETYPE, n_states,n_boxes,nr_track,n)
        implicit none
        integer,intent(IN)           :: mode
!JM        real(RLEN),intent(INOUT)     :: STATE(1:n_states,1:n_boxes)
        real(RLEN),intent(INOUT)     :: STATE(1:n_boxes,1:n_states)
        integer,intent(INOUT)        :: STATETYPE(1:n_states)
        integer,intent(IN)           :: n_states
        integer,intent(IN)           :: n_boxes
        integer,intent(IN)           :: n
        integer,intent(IN)           :: nr_track(1:n_states)

        integer                      ::i
        integer                      ::mt
        integer                      ::ml
        integer                      ::n_rstates

        if ( mode.eq.2) then
           ! allocation of an array with which certain checks can be set off
           allocate(check_ben_to_pel(1:n_states),stat=i)
        endif
        ! Calculate number of "real" states : total number of states -number of labeldstates
        ! only initialization where nr_track(;)<>0:
        !  where nr_track>0 : tracked state variable , wehre nr_track<0: tracer.
        n_rstates=n_states-n+1;
        do mt=n_rstates,n_states
           ml=nr_track(mt); 
           if (ml.gt.0) then
!JM                STATE(mt,1:n_boxes)=min(p_small,STATE(ml,1:n_boxes))
                STATE(1:n_boxes,mt)=min(p_small,STATE(1:n_boxes,ml))
                STATETYPE(mt)=STATETYPE(ml)
           elseif ( ml.eq.-1) then
!JM                STATE(mt,1:n_boxes)=p_small
                STATE(1:n_boxes,mt)=p_small
                STATETYPE(mt)=ALLTRANSPORT
           endif
        enddo
      end subroutine inittransportstatetypes_track

      subroutine fill_sfl_if_tracking( parallel,nr,numc,sfl)
        logical                    :: parallel
        integer,intent(IN)         :: nr
        real(RLEN),intent(INOUT)   :: sfl(0:numc)
        if (iiTrack ==0 .or. parallel ) return
        sfl(nr_3d_track(nr))=sfl(nr)
        return
      end subroutine fill_sfl_if_tracking

      subroutine check_ben_states_on_neg_values_for_track(nr,nr_states,nr_boxes,ig,jg,ccb)

        implicit none
        integer,intent(IN)          :: nr
        integer,intent(IN)          :: nr_states
        integer,intent(IN)          :: nr_boxes
        integer,intent(IN)          :: ig
        integer,intent(IN)          :: jg
!JM        real(RLEN),intent(INOUT)    :: ccb(1:nr_states,0:nr_boxes)
        real(RLEN),intent(INOUT)    :: ccb(0:nr_boxes,1:nr_states)

        integer                      ::i
        integer                      ::mr
        real(RLEN)                   ::r

        if (iiTrack==0) return;
        mr=nr_2d_track(nr)
        ! nr porints to the tracking state variable
        ! mr points to the normal state variable
        if ( nr.gt.nr_states-ii2dTrack) then
          do i=1,nr_boxes
!JM           r=ccb(mr,i)
           r=ccb(i,mr)
           if ( r.lt.ZERO ) then   
              write(LOGUNIT,*) 'Warning:',trim(var_names(StBenStateS-1+mr)),' is negative in', ig,jg
              write(LOGUNIT,*) 'Warning: therefor is the is ',trim(var_names(StBenStateS-1+nr)),&
                                ' set on the same negative value'
!JM              ccb(nr,i)=r;
              ccb(i,nr)=r;
!JM           elseif (r .lt. ccb(nr,i)) then
           elseif (r .lt. ccb(i,nr)) then
              write(message,'(''tracked '',A,'' is larger then'',A ,'' at '',i3,'','',i3)') &
                  trim(var_names(StBenStateS-1+nr)),trim(var_names(StBenStateS-1+mr))
              write(LOGUNIT,*) 'Warning:',trim(message)
              write(LOGUNIT,*) 'Warning: therefor is the is ',trim(var_names(StBenStateS-1+nr)),&
                                ' set on the samevalue'
!JM              ccb(nr,i)=ccb(mr,i);
              ccb(i,nr)=ccb(i,mr);
           endif
          enddo
        endif
      end subroutine check_ben_states_on_neg_values_for_track

!----------------------------------------------------------------------------
      subroutine set_all_exchanges_for_tracking_for_your_model
      ! This routine is a model dependent routine 
      ! All other routines can also used for other model set-ups within the concept of BFM
      ! routine to make a coupling between pelagic and benthic vars to be used
      ! in the tracking routine
      ! Is is assumed that only for a nutrient tracking is done!

        use mem,ONLY:  ppR6c,ppQ6c, ppR6n,ppQ6n, ppR6p,ppQ6p, &
                       ppR6s,ppQ6s, ppR1c,ppQ1c, ppR1n,ppQ1n, ppR1p,ppQ1p, &
                       ppN1p,ppK1p, ppN3n,ppK3n, ppN4n,ppK4n, ppN5s,ppK5s, &
                       ppN6r,ppK6r, ppO3c,ppG3c, ppO3h,ppG3h, &
                       ppY3c,ppY3n, ppY3p, &
                       ppPhytoPlankton,ppMesoZooPlankton,ppMicroZooPlankton, &
                       iiPhytoPlankton,iiMesoZooPlankton,iiMicroZooPlankton, &
                       iiC,iiN,iiP,iiS
        implicit none
        integer     :: i,j

        ! List state vars in case of 
        ! a pelagic source the content is distributed among several benthic state variables
        ! a pelagic sink the content is supplied from several benthic state variables

        call set_exchange_states_benpel(ppN1p,-1)
        call set_exchange_states_benpel(ppN4n,ppK4n)

        do i=1,iiPhytoPlankton
          call set_exchange_states_benpel(ppPhytoPlankton(i,iiC),-1)
          call set_exchange_states_benpel(ppPhytoPlankton(i,iiN),-1)
          call set_exchange_states_benpel(ppPhytoPlankton(i,iiP),-1)
          j=ppPhytoPlankton(i,iiS);if ( j>0) call set_exchange_states_benpel(j,-1)
        enddo

        do i=1,iiMicroZooPlankton
          j=ppMicroZooPlankton(i,iiN)
          if (j>0)call set_exchange_states_benpel(j,-1)
          j=ppMicroZooPlankton(i,iiP)
          if (j>0)call set_exchange_states_benpel(j,-1)
        enddo

        call set_exchange_states_benpel(ppR6c,ppQ6c)
        call set_exchange_states_benpel(ppR6n,ppQ6n)
        call set_exchange_states_benpel(ppR6p,ppQ6p)
        call set_exchange_states_benpel(ppR6s,ppQ6s)

        call set_exchange_states_benpel(ppR1c,-1)
        call set_exchange_states_benpel(ppR1n,-1)
        call set_exchange_states_benpel(ppR1p,ppQ1p)

        call set_exchange_states_benpel(ppN3n,ppK3n,ben_check=0)
        call set_exchange_states_benpel(ppN5s,ppK5s)
        call set_exchange_states_benpel(ppN6r,ppK6r)

#ifdef INCLUDE_PELCO2
#ifdef INCLUDE_BENCO2
        call set_exchange_states_benpel(ppO3c,ppG3c)
        call set_exchange_states_benpel(ppO3h,ppG3h)
#endif
#endif

      end subroutine set_all_exchanges_for_tracking_for_your_model
      subroutine calc_track_rates_for_your_model
#ifdef INCLUDE_TRACK
      use mem,only : jK3G4n,jK34K24n,jK23K13n,K3n,K24n
      use mem,only : jtrK3G4n,jtrK34K24n,jtrK13K3n,trK3n,trK24n

      implicit none

      jtrK3G4n=jK3G4n*(trK3n)/(p_small+K3n)
      jtrK13K3n=min(ZERO,jK23K13n*(trK3n)/(p_small+K3n))
      jtrK34K24n=min(ZERO,jK34K24n*(trK24n)/(p_small+K24n))
#endif
      end subroutine calc_track_rates_for_your_model

      end module track
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
