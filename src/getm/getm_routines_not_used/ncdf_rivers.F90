!$Id: ncdf_rivers.F90,v 1.7 2005-09-23 11:27:43 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: ncdf_river -
!
! !INTERFACE:
   module ncdf_river
!
! !DESCRIPTION:
! 
! !BFM
!  Extended to read info from netcdf files about BFM-state variables:
!        1. reading values  (timeseries) from netcdf file
!        2. when no boundary conditions are available make concentration
!           at boundary equal at adjacentent point in model domain and
!          if limits are set, limit this value between a minimum and a maximum
!        3. calculate states on basis of under states: H3c 
!           ( dissoled inorganic carbon)
!
! !USES:
   use netcdf
   use time, only: string_to_julsecs,time_diff,add_secs,in_interval
   use time, only: julianday,secondsofday,juln,secsn,timestep
   use time, only: write_time_string,timestr,timestep
   use rivers, only: nriver,river_data,river_name,river_flow,river_factor
   use rivers, only: ok,rriver,real_river_name,river_split
   use rivers, only: temp_missing,salt_missing
   use rivers, only: use_river_temp,use_river_salt,river_temp,river_salt
#ifdef GETM_BIO
   use bio, only: bio_calc
   use bio_var, only: numc
   use bfm_output, only: var_names
   use rivers, only: river_bio

#ifdef BFM_GOTM
   use rivers, only: river_bio,river_max,river_min,river_func
#endif
#endif
   IMPLICIT NONE
!
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_river_input_ncdf,get_river_data_ncdf
!
! !PRIVATE DATA MEMBERS:
   REALTYPE                            :: offset
   integer                             :: ncid,ndims,dims(2),unlimdimid,textr
   integer                             :: start(1),edges(1)
   integer                             :: timedim,time_id
   integer, allocatable                :: r_ids(:)
   integer, allocatable                :: salt_id(:)
   integer, allocatable                :: temp_id(:)
   integer, allocatable                :: r_salt(:)
   integer, allocatable                :: r_temp(:)
   REAL_4B, allocatable                :: river_times(:)
#ifdef GETM_BIO
   integer, allocatable                :: bio_id(:,:)
#ifdef BFM_GOTM
   integer, allocatable                :: missing_bio_id(:)
   integer, allocatable                :: r_bio(:,:)
   logical           ::trace,use_tracer_type
#endif
#endif
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_river_input_ncdf -
!
! !INTERFACE:
   subroutine init_river_input_ncdf(fn,nstart)
#ifdef BFM_GOTM
   use coupling_getm_bfm,ONLY: assign_river_inputs_to_bio_states,tracer_type,unlabeled_var_index
#endif
   IMPLICIT NONE
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: fn
   integer, intent(in)                 :: nstart
!
!
! !REVISION HISTORY:
!
!  See module for log.
!
! !LOCAL VARIABLES:
   integer                   :: i,j,m,n,k,ll,nn,mnl,ni
   integer                   :: err
   character(len=19)         :: tbuf
   integer                   :: j1,s1,j2,s2
   character(len=256)        :: time_units
   character(len=256)        :: bio_name
   REAL_4B,  dimension(2)   :: att_wrk
#ifdef BFM_GOTM
   character(len=256)        :: missing_river_name
   character(len=256)        :: use_function
#endif

!EOP
!-------------------------------------------------------------------------
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_river_input_ncdf() # ',Ncall
#endif

   LEVEL3 'init_river_input_ncdf'
#ifdef GETM_BIO
   LEVEL3 'Set limits and forcing fo river input'
#endif

   allocate(r_ids(rriver),stat=err)
   if (err /= 0) stop 'ncdf_river: Error allocating memory (r_ids)'

   allocate(r_salt(rriver),stat=err)
   if (err /= 0) stop 'ncdf_river: Error allocating memory (r_salt)'
   allocate(r_temp(rriver),stat=err)
   if (err /= 0) stop 'ncdf_river: Error allocating memory (r_temp)'

   allocate(salt_id(rriver),stat=err)
   if (err /= 0) stop 'ncdf_river: Error allocating memory (salt_id)'
   allocate(temp_id(rriver),stat=err)
   if (err /= 0) stop 'ncdf_river: Error allocating memory (temp_id)'

#ifdef BFM_GOTM
   allocate(r_bio(rriver,numc),stat=err)
   if (err /= 0) stop 'ncdf_river: Error allocating memory (r_bio)'
   allocate(bio_id(rriver,numc),stat=err)
   if (err /= 0) stop 'ncdf_river: Error allocating memory (bio_id)'
   bio_id = -1
   bio_id = -9999
   allocate(missing_bio_id(numc),stat=err)
   if (err /= 0) stop 'ncdf_river: Error allocating memory (missing_bio_id)'
   missing_bio_id = -9999
#endif

   err = nf90_open(fn,NF90_NOWRITE,ncid)
   if (err .ne. NF90_NOERR) go to 10

   err = nf90_inquire(ncid, unlimitedDimID = unlimdimid)
   if (err .NE. NF90_NOERR) go to 10

   err = nf90_inquire_dimension(ncid,unlimdimid,len = textr)
   if (err .ne. NF90_NOERR) go to 10

   err = nf90_inq_varid(ncid,"time",time_id)
   if (err .ne. NF90_NOERR) go to 10 

   ni=0;nn=0
   do n =1,nriver
     if ( ni.lt.n) then
        nn=nn+1;ni=ni+river_split(ni+1)
     endif
     if (ok(n) .gt. 0) then
        err = nf90_inq_varid(ncid,real_river_name(nn),r_ids(nn))
        if (err .ne. NF90_NOERR) go to 10
        r_salt(nn) = 0
        if ( use_river_salt ) then
           err =  nf90_inquire_attribute(ncid,r_ids(nn),'salt',xtype=ll)
           if (ll==NF90_SHORT) then
               err =  nf90_get_att(ncid,r_ids(nn),'salt',r_salt(nn))
           elseif ( ll==NF90_FLOAT) then
               r_salt(nn)=2;
               err =  nf90_get_att(ncid,r_ids(nn),'salt',att_wrk)
               LEVEL3 trim(real_river_name(nn)),': salinity is set on constant value of ',att_wrk(1)    
           endif
!        if (err .ne. NF90_NOERR) go to 10
         if (r_salt(nn) .eq. 1 ) then
            LEVEL3 'river salinity:    ',trim(real_river_name(nn))//trim('_salt')
            err =  nf90_inq_varid(ncid,trim(real_river_name(nn))//trim('_salt'),salt_id(nn))
            if (err .ne. NF90_NOERR) go to 10
           end if
        end if
        r_temp(nn) = 0
        if ( use_river_temp ) then
         err =  nf90_get_att(ncid,r_ids(nn),'temp',r_temp(nn))
!        if (err .ne. NF90_NOERR) go to 10
         if (r_temp(nn) .eq. 1 ) then
            LEVEL3 'river temperature: ',trim(real_river_name(nn))//trim('_temp')
            err =  nf90_inq_varid(ncid,trim(real_river_name(nn))//trim('_temp'),temp_id(nn))
            if (err .ne. NF90_NOERR) go to 10
         end if
      end if

#ifdef BFM_GOTM
      do m=1,numc
         err=assign_river_inputs_to_bio_states(ncid,real_river_name(nn),nn,m,ll)
         select case (ll)
           case (0)
              ! tracer
              river_bio(nn,m)=1.0;
            case (-1) 
              !track : river found which is tracked and for which no data are available
              ! assumed that concentration in river is equal to the concentration in entering 
              ! grid point corrected for min-max function
              river_bio(nn,m)=0.0
            case (-2) 
              ! track/trace variable found which is not tracked/traced
              ! there is no river input for this state variable
              river_bio(nn,m)=-9999.0
            case default
              bio_id(nn,m)=ll;      ! if bio_id>0 : forcing
              river_bio(nn,m)=0.0;  ! otherwise min_max functions will be used.
          end select
        end do
#endif
    end if
   end do
#ifdef BFM_GOTM
   do m=1,numc
     mnl=unlabeled_var_index(3,m)
     bio_name='limit_'//trim(var_names(mnl))
     err =  nf90_inq_varid(ncid,trim(bio_name),ll)
     if (err .eq. NF90_NOERR) then
        missing_river_name=''
        err=nf90_get_att(ncid,ll,'use',missing_river_name);
        if (err .eq. NF90_NOERR) then
           missing_river_name=missing_river_name(1:scan(missing_river_name,' '));
           bio_name=trim(missing_river_name)//'_'//trim(var_names(m))
           err =  nf90_inq_varid(ncid,trim(bio_name),missing_bio_id(m))
           if (err .ne. NF90_NOERR) missing_bio_id(m) = -2
           if ( missing_bio_id(m) .gt. -1 ) then
              write(bio_name,'(''If no data available for '',A,'' data of '',A,'' are used'')') & 
                                                   trim(var_names(m)),trim(missing_river_name) 
              LEVEL4 trim(bio_name)
           end if
        else
          use_function=''
          err=nf90_get_att(ncid,ll,'function',use_function);
          if (err .eq. NF90_NOERR) then
              river_func(m)=1;       
              write(time_units, &
                '(A,'': river inputs are derived from a function'')') &
                       trim(var_names(m))
              LEVEL4 trim(time_units)

          else
            err=nf90_get_att(ncid,ll,'valid_range',att_wrk);
            if (err .eq. NF90_NOERR) then
              river_min(m)=att_wrk(1);       
              river_max(m)=att_wrk(2);       
              write(time_units, '(A,'' limit min:'',F10.4,'' max.:'',F10.4)') &
                                                 trim(var_names(m)),river_min(m),river_max(m)
              LEVEL4 trim(time_units)
            endif
          endif
       endif
     else
        river_min(m)=0.0
        river_max(m)=0.0
        write(time_units, '(''Warning:no limit conditions found for '',A)') trim(var_names(m))
        LEVEL4 trim(time_units)
        write(time_units, '(A,'' limit min:'',F10.4,'' max.:'',F10.4)') &
                                                 trim(var_names(m)),river_min(m),river_max(m)
        LEVEL4 trim(time_units)
     endif
   enddo
#endif

   allocate(river_times(textr),stat=err)
   if (err /= 0) stop  &
      'init_river_input_ncdf: Error allocating memory (river_times)'

   err =  nf90_get_att(ncid,time_id,'units',time_units)
   if (err .ne. NF90_NOERR) go to 10
   call string_to_julsecs(time_units,j1,s1)
   err = nf90_get_var(ncid,time_id,river_times)
   if (err .ne. NF90_NOERR) go to 10

   offset = time_diff(julianday,secondsofday,j1,s1)
   if( offset .lt. river_times(1) ) then
      FATAL 'Model simulation starts before available river data'
      call write_time_string(julianday,secondsofday,tbuf)
      FATAL 'Simulation starts: ',tbuf
      call add_secs(j1,s1,nint(river_times(1)),j2,s2)
      call write_time_string(j2,s2,tbuf)
      FATAL 'River file starts: ',tbuf
      stop 'init_river_input_ncdf'
   else
      LEVEL3 'River offset time ',offset
   endif

   call add_secs(j1,s1,nint(river_times(textr)),j2,s2)
!   if( time_diff(j1,s1,j2,s2) .lt. _ZERO_ ) then
!      FATAL 'Not sufficient river data available'
!      stop 'init_river_input_ncdf'
!   endif

   call get_river_data_ncdf(nstart)

#ifdef DEBUG
   write(debug,*) 'Leaving init_river_input_ncdf()'
   write(debug,*)
#endif
   return
10 FATAL 'init_river_input_ncdf: ',nf90_strerror(err)
   stop
   end subroutine init_river_input_ncdf
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_river_data_ncdf - .
!
! !INTERFACE:
   subroutine get_river_data_ncdf(loop)
   IMPLICIT NONE
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: loop
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!
!  See module for log.
!
! !LOCAL VARIABLES:
   integer                   :: i,j,n,m,indx,err,k,nn,ni
   REALTYPE                  :: t
   REAL_4B                   :: x(1)
   logical, save             :: first=.true.
   integer, save             :: save_n=1,last_indx=-1
   REALTYPE, save            :: t_1,t_2,loop0
!EOP
!-------------------------------------------------------------------------
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'get_river_data_ncdf() # ',Ncall
#endif
   edges(1) = 1
#define NO_INTERPOL

#ifdef NO_INTERPOL
   if (first) then
      loop0=loop-1
      first = .false.
   endif
   t = (loop-loop0)*timestep
   do indx=save_n,textr
      if (river_times(indx) .ge. real(t + offset)) EXIT
   end do
   if (indx .gt. last_indx) then
      call write_time_string()
      LEVEL3 timestr, ': reading river data .... ',indx
      last_indx = indx
      start(1) = indx
      ni=0;nn=0
      do n =1,nriver
        if ( ni.lt.n) then
           nn=nn+1;ni=ni+river_split(ni+1)
        endif
        if (ok(n) .gt. 0) then
!         err = nf90_get_var(ncid,r_ids(nn),x,start=start,count=edges)
          err = nf90_get_var(ncid,r_ids(nn),x,start,edges)
          if (err .ne. NF90_NOERR) go to 10
          river_flow(nn) = river_factor*x(1)
          river_salt(nn) = salt_missing
          river_temp(nn) = temp_missing
          if ( r_salt(nn) .ge. 1 ) then
            if ( r_salt(nn) .eq. 1 ) then
               err = nf90_get_var(ncid,salt_id(nn),x,start=start,count=edges)
!              err = nf90_get_var(ncid,salt_id(nn),x,start,edges)
            elseif ( r_salt(nn) .eq. 2 ) then
               err =  nf90_get_att(ncid,r_ids(nn),'salt',x)
            endif
            if (err .ne. NF90_NOERR) go to 10
            river_salt(nn) = x(1)
          end if
          if ( r_temp(nn) .eq. 1 ) then
!            err = nf90_get_var(ncid,temp_id(nn),x,start,edges)
             err = nf90_get_var(ncid,temp_id(nn),x,start=start,count=edges)
             if (err .ne. NF90_NOERR) go to 10
             river_temp(nn) = x(1)
          end if
#ifdef BFM_GOTM
         if ( .not.trace  ) then
            do j=1,numc
               if (bio_id(nn,j) .gt. 0) then
                  err = nf90_get_var(ncid,bio_id(nn,j),x,start=start,count=edges)
                  if (err .ne. NF90_NOERR) go to 10
                  river_bio(nn,j) = x(1)
               elseif(missing_bio_id(j) .gt. 0) then
                  err = nf90_get_var(ncid,missing_bio_id(j),x,start=start,count=edges)
                  if (err .ne. NF90_NOERR) go to 10
                  river_bio(nn,j) = x(1)
               end if
            end do
         elseif ( use_tracer_type ) then
             if (bio_id(nn,n) .gt. 0) then
                err = nf90_get_var(ncid,bio_id(nn,n),x,start=start,count=edges)
                if (err .ne. NF90_NOERR) go to 10
                river_bio(nn,nn) = x(1)
             endif
         else
              river_bio(nn,nn) = 1.0 
         endif
#endif
       end if
     end do
   endif

#else
!AS this is not to be used !
   t = loop*timestep
   do indx=save_n,textr
      if (river_times(indx) .gt. real(t + offset)) EXIT
   end do
   ! First time through we have to initialize t_1
   if (first) then
      LEVEL3 'reading first river data - indx = ',indx
      first = .false.
      if (indx .gt. 1) then
         indx = indx-1
      end if
      save_n = indx
      start(1) = indx
      t_1 = river_times(indx) - offset
      t_2 = t_1

      do n =1,nriver
         if (ok(n) .gt. 0) then
            err = nf90_get_var(ncid,r_ids(n),x,start=start,count=edges)
            if (err .ne. NF90_NOERR) go to 10
            river_flow(n) = x(1)
         end if
      end do
   else
      if (indx .gt. save_n) then
         LEVEL3 'reading new river data - indx = ',indx
         save_n = indx
         t_1 = t_2
         t_2 = river_times(indx) - offset
      end if
   end if
#endif

#ifdef DEBUG
   write(debug,*) 'Leaving get_river_data_ncdf()'
   write(debug,*)
#endif
   return
10 FATAL 'get_river_data_ncdf: ',nf90_strerror(err)
   stop
   end subroutine get_river_data_ncdf
!EOC

!-----------------------------------------------------------------------

   end module ncdf_river

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
