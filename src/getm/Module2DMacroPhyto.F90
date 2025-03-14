#include "cppdefs.h"
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Modeule2DMacroPhyto
!
! DESCRIPTION
!   Parameter values for the phytoplankton groups
!
!

!   This file is generated directly from OpenSesame model code, using a code
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  module mem_2DMacroPhyto
!
! !USES:

  use global_mem,only:DONE,RLEN,ZERO,NMLUNIT,NML_OPEN,NML_READ, &
                      LOGUNIT,error_msg_prn
  use mem,  ONLY: NO_BOXES_XY,ppMacroStructure,iiMacroStructure, &
     ppMacroContent,iiC,iiP,iiN,iiL,D2STATE
  use mem_Param,ONLY:CalcMacroPhyto
  use domain,only: imin,jmin,imax,jmax,lonc,latc,az,ioff,joff,arcd1
  use time,only:CalDat,JulDay
  use exceptions, only:getm_error
  use mem_MacroPhyto,only:farm_name_local=>farm_name, &
    i_fb_local=>i_fb,iMsc_local=>iMsc,init_flag_local=>init_flag,  &
    save_year_local=>save_year,save_julian_local=>save_julian, &
    save_status_local=>save_status,save_day0_local=>save_day0, &
    surface_local=>surface,run_1d, &
!do these later    p_start_local=>p_start, p_fin_local=>p_fin,&                !JM
    p_depth_local=>p_depth, p_lines_local=>p_lines, &           !JM
    p_length_m_local=>p_length_m, p_vdist_m_local=>p_vdist_m, & !JM
    p_plants_a_local=>p_plants_a                                !JM
  use mem_MacroPhyto,only:test_MacroPhytogroup_status,p_start,p_fin
  use bio_var,only:ccb
!
!
! !AUTHORS
!   Piet Ruardij
!
!
!
! !REVISION_HISTORY
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
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Default all is public
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  public

  integer,dimension(:,:),allocatable    :: D2CalcMacroPhyto
  integer,dimension(:),allocatable      :: save_year,save_julian,save_day0, &
                                           save_status
  real(RLEN),dimension(:),allocatable   :: surface
  integer                               :: n,m,n_seq,itype,iyear
  integer,private,parameter             ::nmax=20
  real(RLEN),dimension(:),allocatable   :: i_fb_i,iMsc_i
  integer,dimension(nmax)               :: ntype,active
  character(len=10),dimension(nmax)     :: farm_name
  real(RLEN),dimension(nmax)            :: position_x,position_y,i_fb,iMsc
!JM added:
  real(RLEN),dimension(:),allocatable :: p_depth_i,p_lines_i,p_length_m_i, &
                                         p_vdist_m_i,p_plants_a_i
!not yet    integer,dimension(:),allocatable ::                                        p_start_i,p_fin_i
  real(RLEN),dimension(nmax)            :: p_depth,p_lines,p_length_m, &
                                           p_vdist_m,p_plants_a
!not yet  integer,dimension(nmax)            :: p_start,p_fin



 contains

  subroutine Init2dMacroPhyto
  use global_mem,only:bfmnmlopen
    implicit none
    integer             ::rc,i,j,k,l,ipos_x,ipos_y
    character(len=20)   ::file='D2MacroPhyto.nml'
    logical             ::yes=.false.
    logical             ::test_on_macrophyt
    character(len=80)   ::msg=""
    character(len=20)   ::hlp

    namelist/InitialMacroPhyto/n,iyear,ntype,farm_name, &
      position_x,position_y,i_fb,iMsc, &
      p_depth,p_lines,p_length_m,p_vdist_m,p_plants_a !, &   !JM added
!not yet                             p_start,p_fin                                           !JM added

    test_on_macrophyt=.false.
    INQUIRE(file=file,exist=yes)
    if  ( CalcMacroPhyto .and.yes) then
      allocate(D2CalcMacroPhyto(I2DFIELD),stat=rc)
      if (rc /= 0) call getm_error('Init2DMacroPhyto', &
                        'Error allocating memory (D2CalcMacroPhyto)')
      D2CalcMacroPhyto=0
      open(NMLUNIT,file=file,status='old',action='read',err=100)
      read(NMLUNIT,nml=InitialMacroPhyto,err=101)
!     k=2
!     do while(k>0)
!       j=bfmnmlopen(NMLUNIT,k,rc,file=file)
!       if( j>=0) then
!         if (rc.gt.0) goto 100
!         read(NMLUNIT,nml=InitialMacroPhyto,iostat=rc)
!         if (rc>0.and.j.ne.1) goto 101
!       endif
!       close(NML)
!      enddo

      if (n>nmax) then
        write(msg,'(A,I2)')"Numbers of area with Macrophyto is limited to",nmax
        call getm_error("Init2dMacroPhyto", msg)
      endif
      allocate(i_fb_i(n),stat=rc)
      if (rc/= 0) call getm_error('Init2DMacroPhyto','allocating i_fb_i')
      allocate(iMsc_i(n),stat=rc)
      if (rc/= 0) call getm_error('Init2DMacroPhyto','allocating iMsc_i')
      allocate(surface(n),stat=rc)
      if (rc/= 0) call getm_error('Init2DMacroPhyto','allocating surface')
!JM added:
      allocate(p_depth_i(n),stat=rc)
      if (rc/= 0) call getm_error('Init2DMacroPhyto','allocating p_depth_i')
      allocate(p_lines_i(n),stat=rc)
      if (rc/= 0) call getm_error('Init2DMacroPhyto','allocating p_lines_i')
      allocate(p_length_m_i(n),stat=rc)
      if (rc/= 0) call getm_error('Init2DMacroPhyto','allocating p_length_m_i')
      allocate(p_vdist_m_i(n),stat=rc)
      if (rc/= 0) call getm_error('Init2DMacroPhyto','allocating p_vdist_m_i')
      allocate(p_plants_a_i(n),stat=rc)
      if (rc/= 0) call getm_error('Init2DMacroPhyto','allocating p_plants_a_i')
!not yet      allocate(p_start_i(n),stat=rc)
!not yet      if (rc/= 0) call getm_error('Init2DMacroPhyto','allocating p_start_i')
!not yet      allocate(p_fin_i(n),stat=rc)
!not yet      if (rc/= 0) call getm_error('Init2DMacroPhyto','allocating p_fin_i')

      ipos_x=0;ipos_y=0;n_seq=0;active=0
      do i=1,n
!JM needs to be done here
        ipos_x=0; ipos_y=0
        k=(jmin+jmax)/2
        do j=imin+1,imax
!JM          if (az(j-1,k)==1 .or.az(j,k)==1) then
          if (lonc(j-1,k)<position_x(i).and.position_x(i)< lonc(j,k))ipos_x=j
!JM          endif
        enddo
        k=(imin+imax)/2
        do j=jmin+1,jmax
!JM          if (az(k,j-1)==1 .or.az(k,j)==1) then
          if (latc(k,j-1)<position_y(i).and.position_y(i)< latc(k,j))ipos_y=j
!JM          endif
        enddo
        if (ipos_x>0 .and.ipos_y>0) then
          k=ipos_y
          do j=imin+1,imax
            if  (az(j-1,k)==1 .or.az(j,k)==1) then
!!JM          if ( lonc(j-1,k)<position_x(i).or.position_x(i)< lonc(j,k)) then
              if ( lonc(j-1,k)<position_x(i).and.position_x(i)< lonc(j,k)) then
                ipos_x=j
                if (abs(lonc(j-1,k)-position_x(i)) &
                  <abs(position_x(i)-lonc(j,k)).or.az(j,k).ne.1) ipos_x=ipos_x-1
               endif
             endif
          enddo
          k=ipos_x
          do j=jmin+1,jmax
            if (az(k,j-1)==1 .or.az(k,j)==1) then
              if ( latc(k,j-1)<position_y(i).and.position_y(i)< latc(k,j)) then
                ipos_y=j
                if (abs(latc(k,j-1)-position_y(i))< &
                  abs(position_y(i)-latc(k,j)).or.az(k,j).ne.1) ipos_y=ipos_y-1
              endif
            endif
          enddo
          ! active is set negative for all farms at initialization of model..`
          if ( D2CalcMacroPhyto(ipos_x,ipos_y).ne.0) then
            active(i)=-D2CalcMacroPhyto(ipos_x,ipos_y)
          else
            n_seq=n_seq +1;D2CalcMacroPhyto(ipos_x,ipos_y)=n_seq
            active(i)=-n_seq
          endif
          STDERR "============================================================="
#ifdef GETM_PARALLEL
          write(msg,'(''i,j,='',I4,''('',I2,'') '',I4,''('',I2,'')'')')&
               ipos_x+ioff,ipos_x,ipos_y+joff,ipos_y
#else
          write(msg,'(''i,j,='',I4,'' '',I4)')ipos_x,ipos_y
#endif
          hlp="will be initiated"
          if (j<0) hlp="are continued"
          STDERR "at gridpoint ",msg(1:len_trim(msg)), &
                       " Modelling of MacroPhyto culture is prepared "
          STDERR "Parameters:"
          write(stderr,nml=InitialMacroPhyto)
          STDERR "============================================================="
          test_on_macrophyt=.true.
        endif
        i_fb_i(i)=i_fb(i)
        iMsc_i(i)=iMsc(i)
        surface(i)=DONE/arcd1(ipos_x,ipos_y)
!JM added:
        p_depth_i(i)=p_depth(i)
        p_lines_i(i)=p_lines(i)
        p_length_m_i(i)=p_length_m(i)
        p_vdist_m_i(i)=p_vdist_m(i)
        p_plants_a_i(i)=p_plants_a(i)
!not yet        p_start_i(i)=p_start(i)
!not yet        p_fin_i(i)=p_fin(i)

      enddo
    endif
    ! CalcMacrophyt on .false. if no active MacroPhyt s present
    CalcMacroPhyto=test_on_macrophyt
    if (.NOT.test_on_macrophyt) then
      if (allocated(D2CalcMacroPhyto)) then
         deallocate(D2CalcMacroPhyto);deallocate(surface)
         deallocate(i_fb_i);deallocate(iMsc_i)
!JM added:
         deallocate(p_depth_i)
         deallocate(p_lines_i)
         deallocate(p_length_m_i)
         deallocate(p_vdist_m_i)
         deallocate(p_plants_a_i)
!not yet         deallocate(p_start_i)
!not yet         deallocate(p_fin_i)
      endif
    else
      !important feature!!!!!!!
      run_1d=.false
      allocate(save_year(n),stat=rc)
      if (rc/= 0) call getm_error('Init2DMacroPhyto','allocating save_year')
      allocate(save_julian(n),stat=rc)
      if (rc/= 0) call getm_error('Init2DMacroPhyto','allocating save_julian')
      allocate(save_day0(n),stat=rc)
      if (rc/= 0) call getm_error('Init2DMacroPhyto','allocating save_day0')
      allocate(save_status(n),stat=rc)
      if (rc/= 0) call getm_error('Init2DMacroPhyto','allocating save_status')
      save_year=0;save_julian=0;save_day0=0;save_status=0

!     call CalDat(julianday,iye,i,j)
      do i=1,imax
        do j=1,jmax
          n_seq=D2CalcMacroPhyto(i,j)
          if (n_seq>0) then
            do m=1,n
              l=active(m)
!             itype=ntype(m)
              if (abs(l)==n_seq) then
                save_status(m)=-3
              endif
            enddo
          endif
        enddo
      enddo
    endif
    !set parameters for initial biomass and numbers of floating bodies on zero
    i_fb_local=ZERO
    iMsc_local=ZERO
!JM added:
    p_depth_local=ZERO
    p_lines_local=ZERO
    p_length_m_local=ZERO
    p_vdist_m_local=ZERO
    p_plants_a_local=ZERO
!not yet    p_start_local=ZERO
!not yet    p_fin_local=ZERO
    return
   9 STDERR "No line found with content: n=<integer>"
 100 call error_msg_prn(NML_OPEN,"Init2dMacroPhyto.f90","D2MacroPhyto.nml")
 101 call error_msg_prn(NML_READ,"Init2dMacroPhyto.f90","InitialMacroPhyto")
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end subroutine Init2dMacroPhyto

  subroutine Do2dMacroPhyto(julianday,i,j,boxnumber)
    implicit none
    integer,intent(IN)      ::julianday
    integer,intent(IN)      ::i,j
    integer,intent(IN)      ::boxnumber
    integer                 ::l,k,pp_any,iye,day0

    ! subdomains with grid points where macrophytes may grow.
    if (allocated(D2CalcMacroPhyto)) then
      ! active grid point with macrophyts
      save_status_local(:,boxnumber)=0
      n_seq=D2CalcMacroPhyto(i,j)
      if (n_seq==0) return
!      write(LOGUNIT,*)'Farm present: i,j,julianday',i,j,julianday
      do m=1,n
        l=active(m)
        itype=ntype(m)
!write(LOGUNIT,*)'m,l,n_seq,save_status(m)',m,l,n_seq,save_status(m)
!write(LOGUNIT,*)'init_flag_local',init_flag_local(itype)
        if (abs(l)==n_seq) then
          ! save_status=0: no growing macrophytes
          ! save_status=1: growing macrophytes culture
          ! save_status=2: initialization of macrophyte cluture
          ! save_status=-2 wait on signal from d-model if farm is started.
          ! save_status=-3: potential farm planned/continued.
          ! ReSet parameters
          call Set1DMacroPhytParameters(itype)
          surface_local(boxnumber)=surface(m)
!          itype=ntype(m)
!          farm_name_local(itype)=""
!          i_fb_local(itype)=ZERO
!          iMsc_local(itype)=ZERO
!          p_depth_local(itype)=ZERO
!          p_lines_local(itype)=ZERO
!          p_length_m_local(itype)=ZERO
!          p_vdist_m_local(itype)=ZERO
!          p_plants_a_local(itype)=ZERO
!not yet          p_start_local(itype)=ZERO
!not yet          p_fin_local(itype)=ZERO
          if (l.ne.0) then

            if (save_status(m)==-3) then
!write(LOGUNIT,*)"Init2dmacrophyto in if 1"
              !Check on new culature-----------------------------------
              call JulDay(iyear,1,1,day0)
              ! if k<p_start the model run is arrived in the next year
              k=p_fin(itype); if (k<p_start(itype)) k=k+365
!not yet              k=p_fin_i(m); if (k<p_start_i(m)) k=k+365
              !Continuing with existing culture
              !Start with new culture
              save_status(m)=-2
              write(LOGUNIT,*) 'Before check: julianday,day0,p_start,k', &
                      julianday,day0,p_start(itype),k
              if ((julianday< day0+k).and.  &
                   (day0+p_start(itype)< julianday + 2))save_status(m)=1
              write(LOGUNIT,*) 'julianday,day0,p_start', &
                                           julianday,day0,p_start(itype)
!not yet        (day0+p_start_i(m)< julianday))save_status(m)=1
!not yet        write(LOGUNIT,*) julianday,day0,p_start_i(m)
              write(LOGUNIT,*) ' &
                 Do2dMacroPhyto save_status, itype, m, l=', &
                                                save_status(m),itype,m,l
              if (save_status(m)==1) then
                STDERR LINE
                STDERR "Preparation is made for modelling macrophyt-farm ",m
                STDERR LINE

                call Set1DMacroPhytParameters(itype,m)   !JM do here!
write(LOGUNIT,*)"First call test_MacroPhytogroup_status"
write(LOGUNIT,*)"save_status(m),save_status_local(itype,boxnumber)",save_status(m),save_status_local(itype,boxnumber)
                call test_MacroPhytogroup_status(itype,boxnumber,iout=k)  !JM do here!
                surface_local(boxnumber)=surface(m)
              endif
            endif !if save_status(m)==-3

            save_status_local(itype,boxnumber)=save_status(m)
            if (l<0.and.save_status(m)==1) then
!write(LOGUNIT,*)"Init2dmacrophyto in if 2"
              !Situation at a start of a run where a active macroPhytfarm is
              ! active(k==1). Here will only be tested if indeed
              ! imacrophyt-biomasss present.
              l=-l
              !Reset saved parameters of test_MacroPhyto_status............
              active(m)=l
              pp_any=ppMacroStructure(itype,iiC)
              if (D2STATE(pp_any,boxnumber)>ZERO)  then
                STDERR LINE
                STDERR "Continuation with modelling MacroPhyto cultures"
                STDERR LINE
              endif

                call Set1DMacroPhytParameters(itype,m)   !JM do also here!
write(LOGUNIT,*)"First andahalf call test_MacroPhytogroup_status"
write(LOGUNIT,*)"save_status(m),save_status_local(itype,boxnumber)",save_status(m),save_status_local(itype,boxnumber)
                call test_MacroPhytogroup_status(itype,boxnumber,ifirstdayofmonth=1)  !JM do here!
              surface_local(boxnumber)=surface(m)
              save_status_local(itype,boxnumber)=save_status(m)
              save_julian_local(itype,boxnumber)=save_julian(m)
            elseif(l>0) then
!write(LOGUNIT,*)"Init2dmacrophyto in if 2 first else"
              pp_any=ppMacroStructure(itype,iiC)
              write(LOGUNIT,*) 'Msc;',D2STATE(pp_any,:)
              !a farm is present
              save_status_local(itype,boxnumber)=save_status(m)
              save_year_local(itype,boxnumber)=save_year(m)
              save_julian_local(itype,boxnumber)=save_julian(m)
              save_day0_local(itype,boxnumber)=save_day0(m)
            elseif (l<0.and.save_status(m)==-2) then
!write(LOGUNIT,*)"Init2dmacrophyto in if 2 second else"
              save_status(m)=0
              save_status_local(itype,boxnumber)=save_status(m)
              ! Error with initial values
              ! Initial values of macrophyts present in restart file
              ! of farm which is not in use or not restarted.
              ! Initial values are reset on 0
              pp_any=ppMacroStructure(itype,iiC)
              if (D2STATE(pp_any,boxnumber)>ZERO ) then
                STDERR LINE
                STDERR "Resetting  biomas MacroPhyto cultures"
                STDERR LINE
                D2STATE(pp_any,boxnumber) = ZERO
                pp_any=ppMacroContent(itype,iiC);D2STATE(pp_any,boxnumber)=ZERO
                pp_any=ppMacroContent(itype,iiN);D2STATE(pp_any,boxnumber)=ZERO
                pp_any=ppMacroContent(itype,iiP);D2STATE(pp_any,boxnumber)=ZERO
                pp_any=ppMacroContent(itype,iiL);D2STATE(pp_any,boxnumber)=ZERO
              endif
            endif
            ! test on status of new of existing culture
            !Be sure that in case of a new culture with the call in the BFM
            !the last test is completed.
!write(LOGUNIT,*)"Second call test_MacroPhytogroup_status"
!write(LOGUNIT,*)"save_status(m),save_status_local(itype,boxnumber)",save_status(m),save_status_local(itype,boxnumber)
            call test_MacroPhytogroup_status(itype,boxnumber,iout=k)
!write(LOGUNIT,*)'after test_MacroPhytogroup_status: k',k
!            if (k.eq.2)save_julian_local(itype,boxnumber)=save_julian(m)
            if (k.gt.0) then
              ! Set parameters
              call Set1DMacroPhytParameters(itype,m)
!              farm_name_local(itype)=farm_name(m)
!              i_fb_local(itype)=i_fb_i(m)
!              iMsc_local(itype)=iMsc_i(m)
!!JM added:
!              p_depth_local(itype)=p_depth_i(m)
!              p_lines_local(itype)=p_lines_i(m)
!              p_length_m_local(itype)=p_length_m_i(m)
!              p_vdist_m_local(itype)=p_vdist_m_i(m)
!              p_plants_a_local(itype)=p_plants_a_i(m)
!not yet              p_start_local(itype)=p_start_i(m)
!not yet              p_fin_local(itype)=p_fin_i(m)

              surface_local(boxnumber)=surface(m)
              if (k.eq.2) then
                 save_julian_local(itype,boxnumber)=save_julian(m)
write(LOGUNIT,*)"Third call test_MacroPhytogroup_status"
write(LOGUNIT,*)"save_status(m),save_status_local(itype,boxnumber)",save_status(m)
                 call test_MacroPhytogroup_status(itype,boxnumber)
              endif
              active(m)=abs(l)
            endif
          endif
        endif
      enddo
    endif
  end subroutine Do2dMacroPhyto

  subroutine Fin2dMacroPhyto(i,j,boxnumber)
    implicit none
    integer,intent(IN)      ::i
    integer,intent(IN)      ::j
    integer,intent(IN)      ::boxnumber
    integer                 ::l,k

    if (allocated(D2CalcMacroPhyto)) then
      n_seq=D2CalcMacroPhyto(i,j)
      if (n_seq==0) return
      do m=1,n
        l=active(m)
        itype=ntype(m)
        if (abs(l)==n_seq) then
          if (l.gt.0) then
             save_year(m)=save_year_local(itype,boxnumber)
             save_julian(m)=save_julian_local(itype,boxnumber)
             save_day0(m)=save_day0_local(itype,boxnumber)
             save_status(m)=save_status_local(itype,boxnumber)
          endif
        endif
      enddo
    endif
  end subroutine Fin2dMacroPhyto

  subroutine Set1DMacroPhytParameters(itype,seqnr)
  implicit none
  integer,intent(IN) :: itype
  integer,intent(IN),optional :: seqnr

  logical :: mode
  mode=present(seqnr)

  select case (mode)
     case (.false.)
       init_flag_local(itype)=.false.
       farm_name_local(itype)=""
       i_fb_local(itype)=ZERO
       iMsc_local(itype)=ZERO
       p_depth_local(itype)=ZERO
       p_lines_local(itype)=ZERO
       p_length_m_local(itype)=ZERO
       p_vdist_m_local(itype)=ZERO
       p_plants_a_local(itype)=ZERO
     case (.true.)
       farm_name_local(itype)=farm_name(seqnr)
       i_fb_local(itype)=i_fb_i(seqnr)
       iMsc_local(itype)=iMsc_i(seqnr)
!JM added:
       p_depth_local(itype)=p_depth_i(seqnr)
       p_lines_local(itype)=p_lines_i(seqnr)
       p_length_m_local(itype)=p_length_m_i(seqnr)
       p_vdist_m_local(itype)=p_vdist_m_i(seqnr)
       p_plants_a_local(itype)=p_plants_a_i(seqnr)
!not yet              p_start_local(itype)=p_start_i(m)
!not yet              p_fin_local(itype)=p_fin_i(m)
  end select
  end subroutine Set1DMacroPhytParameters

  end module mem_2DMacroPhyto
