!$Id$
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: coupling_getm_bfm_rivers
!
! !INTERFACE:
   module coupling_getm_bfm_rivers
!
! !BFM:
!    This is special routine for BFM  
!    Ths routine include code for:
!       4. river input calculations:
!            -set limits on input when no data are available
!            -set boundary condition for tracking-state variables.
!
! !USES:
!JM    use domain, only: imin,imax,jmin,jmax,kmax,ioff,joff
!JM    use domain, only: az,au,av
!JM    use domain, ONLY:H,imin,imax,jmin,jmax
    use exceptions, only:getm_error

!JM    character(len=80),private                 :: msg,sub
! !PUBLIC DATA MEMBERS:
   character(len=64),allocatable,dimension(:) :: tracked_river_name
   character(len=64)                          :: tracer_type
   REALTYPE                                   :: lon_min 
   REALTYPE                                   :: lat_min=-9999.0 
   REALTYPE                                   :: lon_max 
   REALTYPE                                   :: lat_max 
#ifdef INCLUDE_DIAGNOS_PRF
   integer                                    :: diag_end_sections(3,7)
#else
   integer                                    :: diag_end_sections(3,6)
#endif
   integer                                    :: start_tracking_in_jul
!JM   logical,public                             :: read_poro=.false.
!JM   integer,parameter                          :: DIAG_RESET=0,DIAG_ADD=1,&
!JM                                                  DIAG_AVERAGE=2,DIAG_INFO=3

!JM   public fill_diagn_bfm_vars, init_pel_co2, &
!          getm_bfm_bennut_calc_initial, set_2d_grid_parameters, &
!          init_2d_grid,unlabeled_var_index, &
!          assign_river_inputs_to_bio_states,check_reset_tracking, &
!          check_3d_track_rates,make_uv_flux_output,make_river_flux_output, &
!          output_2d_grid_parameters
   public assign_river_inputs_to_bio_states,make_river_flux_output

!EOP
!-----------------------------------------------------------------------

   contains
!-----------------------------------------------------------------------

!BOP

  integer function  assign_river_inputs_to_bio_states(ncid,river_name,river_nr,state_nr,iout)
#ifdef BFM_GOTM
     use netcdf
     use string_functions, ONLY: getseq_number
!JM     use bio_var, only: numc,var_names
     use bio_var, only: numc
     use bfm_output, only: var_names
     use mem, ONLY: ii3dTrack,nr_3d_track,iiTrack
     use trace_bdy,only: trace
     use time, only:JulDay
#endif


     implicit none
     integer,intent(IN)           ::ncid
     integer,intent(IN)           ::river_nr
     integer,intent(IN)           ::state_nr
     integer,intent(OUT)          ::iout
     character(len=*),intent(IN)  ::river_name

#ifdef BFM_GOTM
     integer                                    :: nr_real_states
     integer                                    :: unit=311
     integer                                    :: rc
     integer                                    :: n
     integer                                    :: m
     integer                                    :: err=0
     integer                                    :: nriver=0
     character(len=128)                         :: name=""
     character(len=20)                          :: river_part
!JM     REAL_4B                                    :: att_wrk
     REAL*4                                    :: att_wrk
     REALTYPE                                   :: lon_river
     REALTYPE                                   :: lat_river
     logical                                    :: use_tracer_type=.false.
     logical                                    :: exist=.false.
     integer                                    :: yy,mm,dd
     character(len=1)                           :: c1,c2
  
     namelist /trace_nml/ tracer_type

     nr_real_states=numc
     river_part="";if (len_trim(river_name)> 0 .and.river_nr.gt.0 ) river_part =trim(river_name)//'_'
     if ( iiTrack.ge.1) nr_real_states=numc-ii3dTrack;
     if ( ncid==0 ) then
       if ( index(var_names(1),'@riv_')>0 ) trace=.true.
       if (trace) then
         inquire(file='tracer_info.dat',exist=exist)
         tracer_type=''
         if (exist) then
              open(unit,file='tracer_info.dat',action='read',status='old',iostat=err)
              if ( err.ne.0) &
                 call getm_error("assign_river_inputs_to_bio_states","error when opening tracer_info.dat")
              read(unit,nml=trace_nml,iostat=err)
              if ( err.ne.0) &
                 call getm_error("assign_river_inputs_to_bio_states", &
                                        "error when reading namelist from tracer_info.dat")
              close(unit)
         endif
         use_tracer_type=len_trim(tracer_type).ne.0
         LEVEL3 'TRACER_RUN'
         LEVEL3 'tracer_type=',tracer_type
       elseif ( iiTrack.ge.1.and.ii3dTrack.gt.1 ) then
         LEVEL3 'TRACK_RUN'
         open(unit,file='track_info.dat',action='read',status='old',iostat=err)
         if ( err.ne.0) &
           call getm_error("assign_river_inputs_to_bio_states","error when opening track_info.dat")
         read(unit,*) name
         ! get a date from the first line in track_info.dat
         if( scan(name,'/-') .eq.0) &
          call getm_error("assign_river_inputs_to_bio_states","No valid date on first line of track_info.dat") 
          ! reclaculate the date in julian days  and put in integer_var n
          read(name,'(i4,a1,i2,a1,i2)') yy,c1,mm,c2,dd; call JulDay(yy,mm,dd,n)
         start_tracking_in_jul=n 
         read(unit,*) nriver
         if ( nriver >0 ) then
           allocate(tracked_river_name(nriver),stat=rc) ! NetCDF name of river which will be tracked
           if (rc /= 0) stop 'assign_river_inputs_to_bio_states:Error allocating memory (tracked_river_name)'
           do n=1,nriver
              read(unit,*) tracked_river_name(n)
              LEVEL4 "A nutrient or the tracer in river ",trim(tracked_river_name(n))," will be tracked"
           enddo
         else
           read(unit,*) lon_min,lat_min
           read(unit,*) lon_max,lat_max
           write(name,1000) lon_min,lat_min,lon_max,lat_max
 1000      format('All rivers in the area ',F5.1,',',F6.1,' to ',F5.1,',',F6.1,' will be tracked')
           LEVEL4  trim(name)
         endif
         close(unit)
       endif
     elseif ( iiTrack.eq.0.or.state_nr<=nr_real_states) then
        iout=-9999
        if (.not. trace ) then
          name=trim(river_part)//trim(var_names(state_nr))
        else
          n=state_nr-river_nr;if ( n.eq.0) iout=0;
          if (n .eq.0.or.(.not.use_tracer_type)) goto 90 ; !OK: no tracer  variable
          name=trim(river_part)//trim(tracer_type)
        endif
        err =  nf90_inq_varid(ncid,trim(name),iout)
        if (err .ne. NF90_NOERR) iout = -9999
        if ( iout.ne. -9999 ) LEVEL4 trim(river_name),': ',trim(var_names(state_nr))
     elseif ( ii3dTrack>1) then
        ! tracking only-----------------------------------------
        m=0;n=nr_3d_track(state_nr) ;iout=-2; 
        ! n> 0 : tracked state variable, n =-1: tracer.
        if ( n.gt. 0 .or.n.eq.-1 ) then
          if (allocated(tracked_river_name)) then
            ! check if river is track if this is the case  m will 
            !     get the sequence number of the river.
            m=size(tracked_river_name)
            m=getseq_number(river_name,tracked_river_name(1:m),m,.TRUE.)
            err=NF90_NOERR
          elseif (lat_min> -9998.0) then
            ! if river point is in the area  m will get the value 1
            err =  nf90_inq_varid(ncid,trim(river_name),m)
            if (err .eq. NF90_NOERR) then
              err=nf90_get_att(ncid,m,'longitude',att_wrk);
              if (err .eq. NF90_NOERR) then
                lon_river=att_wrk
                err=nf90_get_att(ncid,m,'latitude',att_wrk);
                if (err .eq. NF90_NOERR) then
                  lat_river=att_wrk;m=0
                  if (lat_river-lat_max< -1.0D-6 .and. lat_river-lat_min >1.0D-6 &
                  .and. lon_river-lon_max< -1.0D-6 .and. lon_river-lon_min >1.0D-6 ) m=1
                endif
              endif
            endif
          endif
        endif
        if ( err.eq.NF90_NOERR) then
          if (m>0) then
            if (n>0) then
              ! tracked state variable
              name=trim(river_part)//trim(var_names(n))
              err =  nf90_inq_varid(ncid,trim(name),iout)
              if ( err.eq.NF90_NOERR) then
                ! there is an timeseries found which will be used also for the tracked variable
                write(name,'(A,''-'',A,'' in '',A,'' will be tracked by adding this input to '',A)') &
                 trim(river_name), trim(var_names(n)),trim(river_name),trim(var_names(state_nr))
                 LEVEL4 trim(name) 
              else
               ! there is a tracked river found however no timeseries
               ! instead of a timeseries of the limits will be used which are used for the
               ! the full non_traced state variable.
               iout=-1
              endif
            elseif ( n.eq.-1) then
              ! sepecial river tracer state variable (tr_t)
              ! iout=0 controls river input such that the conv. of the river is set on zero.
              write(name,'(''Water from '',A,'' will be traced by setting the input for '',A,'' on 1'')') &
                   trim(river_name), trim(river_name)
              LEVEL4 trim(name) 
              iout=0
            endif
          endif
        endif
     endif

  90  assign_river_inputs_to_bio_states=err;
#else
     iout=0
     assign_river_inputs_to_bio_states=0;
#endif
  end function assign_river_inputs_to_bio_states

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: 
!
! !INTERFACE:
      subroutine  make_river_flux_output(height,dt,rSurface,rivcon,n,ig,jg)
    
!
! !USES:
#ifdef BFM_GOTM
     use variables_bio_3d, only: cut,cvt,ccb3d_out,counter_ave      !getm
     use domain, only: imin,imax,jmin,jmax,kmax,ioff,joff   !getm
     use domain, only: dxu,dxv,dyu,dyv                          !getm
     use domain, only: az                                       !getm
!JM     use bio_var,only: var_ids,var_ave,stBenFluxS,stBenFluxE    !bfm
     use bfm_output,only: var_ids,var_ave,stBenFluxS,stBenFluxE    !bfm
     use mem, only: flx_calc_nr,flx_CalcIn,flx_t,flx_states, &  !bfm
              flx_ostates,flx_SS,flx_cal_ben_start,flx_option   !bfm
     use time, only:secsprday
#endif
!
! !INPUT PARAMETERS:
     implicit none
     REALTYPE,intent(IN)         :: height
     REALTYPE,intent(IN)         :: dt
     REALTYPE,intent(IN)         :: rSurface !reciproke of surface
     REALTYPE,intent(IN)         :: rivcon(n)
     
     integer,intent(IN)          :: n
     integer,intent(IN)          :: ig
     integer,intent(IN)          :: jg
!
!
!
#ifdef BFM_GOTM
! LOCAL PARAMETERS:
      integer                                 :: i
      integer                                 :: j
      integer                                 :: l
      integer                                 :: k
      REALTYPE                                :: r
      REALTYPE                                :: s
!
!
!EOP
!-------------------------------------------------------------------------
!BOC
      ! diag_end_sections is initiliazed in the call to fill_diagn_bfm_vars(0,......
      ! this routine is called in init_getm_bio
      ! OUTPUT of is written in a "benthic" var. (D2FIELD) .. (one number per river)
      ! 
      k=diag_end_sections(2,3)
      j=flx_cal_ben_start;
      do i=stBenFluxS,stBenFluxE
       j=j+1
       if ( var_ids(i)==0) then
         !shortcut: no action
       else 
         k=k+1;
         r=_ZERO_;
         if (flx_option(j).eq.50.or.flx_option(j).eq.52) then
           if ( dt.eq._ZERO_ ) then
             s=_ZERO_
           elseif ( flx_option(j).eq.50) then
             ! Calculate vol per second 
             s=height/(dt*rSurface)
           else
             ! Calculate  input per m3 per day
             s=height/dt*secsprday
           endif
           ! filling of output vars
           ! filling hapens at every delta if ave==.true
           do l=flx_calc_nr(j-1)+1,flx_calc_nr(j)
             select case (flx_states(l))
               case (0)     ; r=s
               case default ; r=r+s* rivcon(flx_states(l))
             end select
           enddo
           select case (counter_ave(ig,jg))
             case (0)     ; ccb3d_out(k,ig,jg,1)=r
             case default ; ccb3d_out(k,ig,jg,1)=ccb3d_out(k,ig,jg,1)+r
           end select
         endif
       endif
      enddo
#endif
     end subroutine make_river_flux_output
!!EOC
!!-----------------------------------------------------------------------


end module coupling_getm_bfm_rivers

!-----------------------------------------------------------------------
!EOC
!-----------------------------------------------------------------------
! Copyright (C) 2008 - BFM                                              !
!-----------------------------------------------------------------------

