!$Id: ncdf_3d_bdy.F90,v 1.10 2005-05-04 11:50:57 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  ncdf_3d_bdy - input in NetCDF format
!
! !INTERFACE:
   module ncdf_3d_bdy
!
! !DESCRIPTION:
! 
! !BFM
!    Extended to read info from netcdf files about BFM-state variables:
!        1. reading values  (timeseries) from netcdf file
!        2. when no boundary conditions are available make concentration
!           at boundary equal at adjacentent point in model domain and
!          if limits are set, limit this value between a minimum and a maximum
!        3. calculate states on basis of under states: H3c 
!           ( dissoled inorganic carbon)
!    JM Jan 2013: finished reading biological variables from non-climatology file
!
! !USES:
   use netcdf
   use domain, only: imin,imax,jmin,jmax,kmax,ioff,joff
   use domain, only: nsbv,NWB,NNB,NEB,NSB,bdy_index
   use domain, only: wi,wfj,wlj,nj,nfi,nli,ei,efj,elj,sj,sfi,sli
   use domain, only: H
   use m2d, only: dtm
   use variables_3d, only: hn
   use bdy_3d, only: T_bdy,S_bdy
   use time, only: string_to_julsecs,time_diff,add_secs
   use time, only: julianday,secondsofday,juln,secsn
   use time, only: write_time_string,timestr
#ifdef BFM_GOTM
   use domain, only: bdy_to,bdy_seq
   use bdy_3d, only: cc3d_bdy,cc3d_id,cc3d_bdy_max,cc3d_bdy_min,cc3d_bdy_func, &
                     cc3d_bdy_multi, &
                     no_gradient_ico_missing_values,limit_all_bio,area
   use bio, only: bio_calc
   use bio_var, only: numc
   use bfm_output, only: var_names
#endif
   IMPLICIT NONE
!
   private
!
   public                              :: init_3d_bdy_ncdf,do_3d_bdy_ncdf
!
! !PRIVATE DATA MEMBERS:
   integer                             :: ncid
   integer                             :: time_id,temp_id,salt_id
   integer                             :: start(4),edges(4)
   integer                             :: zax_dim,zax_len,zax_pos
   integer                             :: time_dim,time_len,time_pos
   logical                             :: climatology=.false.
   logical                             :: from_3d_fields
   integer                             :: n_bio
   REALTYPE                            :: offset
   REALTYPE                            :: test_no_gradient
   REAL_4B, allocatable,  dimension(:) :: bdy_times,wrk
   REAL_4B, allocatable,  dimension(:) :: zlev
   REALTYPE, allocatable, dimension(:)   :: ST_one,h_one
   REALTYPE, allocatable, dimension(:,:) :: T_old, T_new
   REAL_4B,  allocatable, dimension(:,:) :: T_wrk
   REALTYPE, allocatable, dimension(:,:) :: S_old, S_new
   REAL_4B,  allocatable, dimension(:,:) :: S_wrk
   REALTYPE, allocatable, dimension(:,:,:) :: T_bdy_clim,S_bdy_clim
#ifdef BFM_GOTM
   REAL_4B,  dimension(2)   :: att_wrk
   REAL_4B                  :: att_val
   REALTYPE, allocatable, dimension(:,:,:)   :: cc3d_bdy_old, cc3d_bdy_new
   REAL_4B,  allocatable, dimension(:,:,:)   :: cc3d_bdy_wrk
   REAL_4B,  allocatable, dimension(:,:)     :: cc3d_bdy_tmp
   REALTYPE, allocatable, dimension(:,:,:,:) :: cc3d_bdy_clim
   logical                  :: trace
#endif
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: ncdf_3d_bdy.F90,v $
!  Revision 1.10  2005-05-04 11:50:57  kbk
!  support for non-climatological 3D boundaries (S,T)
!
!  Revision 1.9  2004/04/06 16:32:29  kbk
!  TimeDiff --> time_diff
!
!  Revision 1.8  2003/12/16 16:50:41  kbk
!  added support for Intel/IFORT compiler - expanded TABS, same types in subroutine calls
!
!  Revision 1.7  2003/10/07 15:10:42  kbk
!  use zax_dim as argument to dim_len
!
!  Revision 1.6  2003/08/03 09:19:41  kbk
!  optimised reading of climatological boundary data
!
!  Revision 1.5  2003/05/05 15:44:20  kbk
!  reads boundary values from 3D fields as individual columns
!
!  Revision 1.4  2003/04/23 11:54:03  kbk
!  cleaned code + TABS to spaces
!
!  Revision 1.3  2003/04/07 16:19:52  kbk
!  parallel support
!
!  Revision 1.1.1.1  2002/05/02 14:01:49  gotm
!  recovering after CVS crash
!
!  Revision 1.1  2001/10/17 13:28:27  bbh
!  Initial import
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: init_3d_bdy_ncdf -
!
! !INTERFACE:
   subroutine init_3d_bdy_ncdf(fname)
!
! !DESCRIPTION:
!  kurt,kurt
!
! !USES:
#ifdef BFM_GOTM
    use string_functions, only:getseq_number,empty
    use coupling_getm_bfm,  only:unlabeled_var_index
    use variables_bio_3d, only: cc3d
#endif

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: fname
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  See log for module
!
! !LOCAL VARIABLES:
   character(len=256)        :: units
   character(len=19)         :: tbuf
   integer                   :: j1,s1,j2,s2,errint
   integer                   :: ndims, nvardims
   integer                   :: vardim_ids(4)
   integer, allocatable, dimension(:):: dim_ids,dim_len
   character(len=16), allocatable :: dim_name(:)
   integer                   :: rc,err
   integer                   :: i,j,k,l,m,n,id,ib,ll,ibnl
   integer                   :: bio_bs
   integer                   :: NaB(4)
   character(len=1)          :: gdir(4)
   character(len=20)         :: text
   character(len=256)        :: use_function
   logical                   :: limit_vars=.false.,ll_exist
   integer                   :: local_ncid
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   write(debug,*) 'ncdf_init_3d_bdy (NetCDF)'
   write(debug,*) 'Reading from: ',trim(fname)
#endif

   LEVEL3 'init_3d_bdy_ncdf'

   err = nf90_open(fname,NF90_NOWRITE,ncid)
   if (err .NE. NF90_NOERR) go to 10

   err = nf90_inquire(ncid, nDimensions = nDims)
   if (err .NE. NF90_NOERR) go to 10

   allocate(dim_ids(ndims),stat=rc)
   if (rc /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (dim_ids)'

   allocate(dim_len(ndims),stat=rc)
   if (rc /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (dim_len)'

   allocate(dim_name(ndims),stat=rc)
   if (rc /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (dim_name)'

   do n=1,ndims
      err = nf90_inquire_dimension(ncid,n,name=dim_name(n),len=dim_len(n))
      if (err .NE. NF90_NOERR) go to 10
      LEVEL4 n,dim_name(n), dim_len(n)
   end do

   if(ndims .eq. 4) then
!     We are reading boundary values from a full 3D field
!     We assume COARDS conventions
!     1 -> lon,x-axis
!     2 -> lat,y-axis
!     3 -> zax,levels
!     4 -> time
      LEVEL4 'boundary data from 3D fields'
      from_3d_fields=.true.
      zax_pos = 3
      time_pos = 4
   else
!     We are reading boundary values from a special boundary data file
!     The variables 'salt' and 'temp' must both exist and be spanned by
!     dimensions as:
!       1 -> zax,levels
!       2 -> bdy_points
!       3 -> time
      LEVEL4 'special boundary data file'
      from_3d_fields=.false.
      zax_pos = 1
      time_pos = 3
!     Note(BJB): This test may break backward compatibility,
!                so I leave it out for now:
      !if (ndims .NE. 3) stop 'init_3d_bdy_ncdf: Wrong number of dims 
   end if

!  We will use this information to actually find the dimension
!  index numbers in the data set.
!  Some of the tests will be repeated later (fixing is possible but not
!  high priority, BJB 2007-04-25).

   LEVEL4 ' ... checking variable "temp"'

   err = nf90_inq_varid(ncid,'temp',temp_id)
   if (err .NE. NF90_NOERR) go to 10

   err = nf90_inquire_variable(ncid,temp_id,ndims=nvardims)
   if (err .NE. NF90_NOERR) go to 10

   if (nvardims .NE. ndims) &
        stop 'init_3d_bdy_ncdf: Wrong number of dims in temp'

   err = nf90_inquire_variable(ncid,temp_id,dimids=vardim_ids)
   if (err .NE. NF90_NOERR) go to 10
   
   zax_dim  = vardim_ids(zax_pos)
   time_dim = vardim_ids(time_pos)

   ! The 'salt' part is only for error capture.
   LEVEL4 ' ... checking variable "salt"'

   err = nf90_inq_varid(ncid,'salt',salt_id)
   if (err .NE. NF90_NOERR) go to 10

   err = nf90_inquire_variable(ncid,salt_id,ndims=nvardims)
   if (err .NE. NF90_NOERR) go to 10

   if (nvardims .NE. ndims) &
        stop 'init_3d_bdy_ncdf: Wrong number of dims in salt'

   err = nf90_inquire_variable(ncid,salt_id,dimids=vardim_ids)
   if (err .NE. NF90_NOERR) go to 10

   if (zax_dim /= vardim_ids(zax_pos)) &
        stop 'init_3d_bdy_ncdf: Position of zax dimension of salt and temp differs'
   if (time_dim /= vardim_ids(time_pos)) &
        stop 'init_3d_bdy_ncdf: Position of time dimension of salt and temp differs'

   zax_len = dim_len(zax_dim)
   time_len = dim_len(time_dim)

   allocate(zlev(zax_len),stat=rc)
   if (rc /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (zlev)'

   err = nf90_inq_varid(ncid, dim_name(zax_dim), id)
   if (err .ne. NF90_NOERR) go to 10

   err = nf90_get_var(ncid,id,zlev)
   if (err .ne. NF90_NOERR) go to 10

!  a few sanity checks on the vertical axis for the 3D boundaries
   do n=1,zax_len
      if (zlev(n) .eq. NF90_FILL_REAL) then
         FATAL '3D boundary z-axis contains NF90_FILL_REAL values'
         FATAL 'proper interpolation cant be done'
         stop 'init_3d_bdy_ncdf'
      end if
   end do
!  not sure if this check is safe - kb
   if ( zlev(1) .ge. _ZERO_ .and. zlev(zax_len) .gt. _ZERO_ ) then
      LEVEL4 'converting positive z-axis (depth) values to negative'
      zlev = -_ONE_*zlev
   end if
!  check strict monotonicity
   do n=1,zax_len-1
      if ( .not. zlev(n) .gt. zlev(n+1) ) then
         FATAL '3D boundary z-axis not strict monotone: ',zlev(n),zlev(n+1)
         stop 'init_3d_bdy_ncdf'
      end if
   end do

   if( time_len .eq. 12) then
      climatology=.true.
      LEVEL4 'Assuming climatolgical 3D boundary conditions'
      LEVEL4 '# of times = ',time_len
   end if

   n_bio=0
#ifdef BFM_GOTM
   err = nf90_get_att(ncid,NF90_GLOBAL,'area',area)
   if (err .eq. NF90_NOERR) then
     LEVEL4 ' biological boundary data are meant for ',trim(area)
   else
      area=''
     LEVEL4 ' no area defined for biological data'
   endif
   LEVEL4 ' biological boundary data from 3D fields'
   n_bio=numc
   trace=.false.
   if (bio_calc) then 
     bio_bs=0
     do ib=1,n_bio
        ibnl=unlabeled_var_index(3,ib);ll=-1;
        err =  nf90_inq_varid(ncid,trim(var_names(ib)),cc3d_id(ib))
        if (err .eq. NF90_NOERR) then
           bio_bs=bio_bs+1
           LEVEL4 'full boundaries for '//trim(var_names(ib))
           ll=cc3d_id(ib)
           cc3d_bdy_func(ib)=0
        else
           err =  nf90_inq_varid(ncid,'limit_'//trim(var_names(ib)),ll)
           if (err .eq. NF90_NOERR) then
              cc3d_id(ib) = -1; limit_vars=.true.
           else
              cc3d_id(ib) = -99;ll=-1
           endif
        end if
        if ( ll >=0 ) then
            use_function='';text='function'
            err=nf90_get_att(ncid,ll,'function',use_function);
            if (err .eq. NF90_NOERR) then  
              cc3d_bdy_func(ib)=100000
              write(units, &
                '(''boundaries for '',A,'' are derived from a function'')') &
                       trim(var_names(ib))
              LEVEL4 trim(units)
            endif
            text='valid_range';
            err=nf90_get_att(ncid,ll,'valid_range',att_wrk);
            if (err .ne. NF90_NOERR) then
              if (cc3d_id(ib) <=0) goto 11
            else
              cc3d_bdy_min(ib)=att_wrk(1);       
              cc3d_bdy_max(ib)=att_wrk(2);       
            endif
        elseif (ibnl.ne.ib) then
            ! track state variable" always 0 on boundaries
            cc3d_bdy_min(ib)=0.0;
            cc3d_bdy_max(ib)=0.0;
            cc3d_id(ib) = -1
         endif
     enddo
     if ( limit_vars==.false. ) then
        inquire(file="bdy_limit.nc",exist=ll_exist)
        if (ll_exist) then
          LEVEL4 "bdy_limit.nc found!"
          err = nf90_open("bdy_limit.nc",NF90_NOWRITE,local_ncid)
          if (err .NE. NF90_NOERR) go to 10
          err = nf90_get_att(local_ncid,NF90_GLOBAL,'set_limits',text)
          if (text == 'all') limit_all_bio=.true.
          do ib=1,n_bio
            !check only variables for which no forcing data are present
            err =  nf90_inq_varid(local_ncid,trim(var_names(ib)),ll)
            if (err .eq. NF90_NOERR) then
              use_function='';units='';
              err=nf90_get_att(local_ncid,ll,'function',use_function);
              if (cc3d_id(ib)<0) then
                if (err .eq. NF90_NOERR) then
                  i=getseq_number(use_function,var_names,n_bio,.true.) 
                  if (i<=0) i=100000
                  if ( i>n_bio) then
                    cc3d_bdy_func(ib)=i
                    write(units, &
                   '(''boundaries for '',A,'' are derived from a function'')') &
                       trim(var_names(ib))
                  elseif ( cc3d_id(i)==-99.or.cc3d_id(i).eq.-1 ) then
                    cc3d_id(ib)=-99
                  else
                    cc3d_bdy_func(ib)=i
                    write(units, &
                   '(''boundaries for '',A,'' are derived from a function and '',A)') &
                       trim(var_names(ib)),trim(use_function)
                  endif
                  if (.NOT.EMPTY(units))LEVEL4 trim(units)
                endif
                !in case of definitions of valid-ranges in bdy_3d.nc and limit_bdy.nc
                ! the last one has the prioirty
                text='valid_range';
                err=nf90_get_att(local_ncid,ll,'valid_range',att_wrk);
                if (err .ne. NF90_NOERR) then
                   if (cc3d_id(ib) <=0) goto 11
                else
                  cc3d_bdy_min(ib)=att_wrk(1);       
                  cc3d_bdy_max(ib)=att_wrk(2);       
                  cc3d_id(ib)=-1
                endif
                err=nf90_get_att(local_ncid,ll,'multi',att_val);
                if (err .eq. NF90_NOERR) then
                  cc3d_bdy_multi(ib)=att_val;       
                endif
              elseif (err .eq. NF90_NOERR) then
                i=getseq_number(use_function,var_names,n_bio,.true.) 
                if (i.eq.ib) then
                  cc3d_bdy_func(ib)=i
                  err=nf90_get_att(local_ncid,ll,'multi',att_val);
                  if (err .eq. NF90_NOERR) cc3d_bdy_multi(ib)=att_val;       
                  if (cc3d_bdy_func(ib)==ib .and.cc3d_bdy_multi(ib) > 0.0) then 
                    write(units, &
                   '(''boundaries for '',A,'' are correct with '',G12.4)') &
                       trim(var_names(ib)), cc3d_bdy_multi(ib)
                   LEVEL4 trim(units)
                  endif
                endif
              endif
            endif
          enddo
        endif
        err=nf90_close(local_ncid)
     endif
     i=0
     do ib=1,n_bio
       if (cc3d_bdy_max(ib)>_ZERO_)  then
         if (.not.limit_all_bio.and.i==0) LEVEL4 &
          "next limits are used only for missing values in boundary timeseries &
           in case of (complete) missing boundaries"  
         write(units, &
            '(A,'' boundary limit min.value:'',F10.4,'' max.value:'',F10.4)')&
             trim(var_names(ib)),cc3d_bdy_min(ib),cc3d_bdy_max(ib)
             LEVEL4 trim(units)
         i=1
       endif
     enddo
     do ib=1,n_bio
       if (cc3d_bdy_func(ib)<-1)  then
         write(units, &
            '(''WARNING: NO boundaries for '',A)') trim(var_names(ib))
         LEVEL4 trim(units)
       endif
     enddo
     if ( bio_bs ==0 ) then
       i=getseq_number('@',var_names,n_bio,.false.) 
       if ( i>0) then
          trace=.true.
          LEVEL4 'trace_model: boundaries defined for every tracer'
          cc3d_id(1:n_bio)=-2
          bio_bs=n_bio
       endif
     endif

     allocate(cc3d_bdy(1:bio_bs,0:kmax,1:nsbv),stat=rc)
     if (rc /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (cc3d_clim)'

   endif
#endif



   if (climatology) then

      allocate(wrk(zax_len),stat=rc)
      if (rc /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (wrk)'

      allocate(T_bdy_clim(time_len,0:kmax,1:nsbv),stat=rc)
      if (rc /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (T_bdy_clim)'

      allocate(S_bdy_clim(time_len,0:kmax,1:nsbv),stat=rc)
      if (rc /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (S_bdy_clim)'

#ifdef BFM_GOTM
      if ( bio_calc) then 
        allocate(cc3d_bdy_clim(1:bio_bs,time_len,0:kmax,1:nsbv),stat=rc)
        if (rc /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (cc3d_bdy_clim)'
        if ( trace ) then
           cc3d_bdy_clim=0.0
           LEVEL4 'trace_model: all boundaries set on zero except: '
           ib=1
           m=0
           l=0
           k=bdy_index(1)
           NaB(1)= NWB; gdir(1)='W'
           NaB(2)= NNB; gdir(2)='N'
           NaB(3)= NEB; gdir(3)='E'
           NaB(4)= NSB; gdir(4)='S'
         ! find place of first state var which represent first open boundary..
           ib=getseq_number('@bndy_',var_names,n_bio,.false.)-1 
           do n=1,4
              do id=1,NaB(n)
                 l=l+1
                 i=ib+bdy_seq(l);
                 LEVEL4 var_names(i),' index=',i
                 k=bdy_index(l)
                 m=bdy_to(l)
                 cc3d_bdy_clim(i,1:time_len,0:kmax,k:m)=1.0;
              enddo
           enddo
        else
          cc3d_bdy_clim=-9999
        endif
      endif
#endif
!     we read each boundary column individually
!     here we can read from both a 3D field and from a
!     special boundary data file - only the arguments 'start' and 'edges'
!     varies in the calls to 'nf90_get_var()'
!     m counts the time
!     l counts the boundary number
!     k counts the number of the specific point
!     MUST cover the same area as in topo.nc
      if (from_3d_fields) then
        edges = 1
        edges(zax_pos) = zax_len
        start(zax_pos) = 1; edges(3) = dim_len(zax_dim);
         edges(4) = 1
      else
         start(1) = 1; edges(1) = dim_len(zax_dim);
         edges(2) = 1;
         edges(3) = 1
      end if

      do m=1,time_len
         start(time_pos) = m
         l = 0
         do n=1,NWB
            l = l+1
            k = bdy_index(l)
            i = wi(n)
            do j=wfj(n),wlj(n)
               if (from_3d_fields) then
                  start(1) = i+ioff ; start(2) = j+joff
               else
                  start(2) = k
               end if
               err = nf90_get_var(ncid,salt_id,wrk,start=start,count=edges)
               if (err .ne. NF90_NOERR) go to 10
               call interpol(0,zax_len,zlev,wrk,H(i,j),kmax,hn(i,j,:), &
                             S_bdy_clim(m,:,k),errint)
               err = nf90_get_var(ncid,temp_id,wrk,start=start,count=edges)
               if (err .ne. NF90_NOERR) go to 10
               call interpol(0,zax_len,zlev,wrk,H(i,j),kmax,hn(i,j,:), &
                             T_bdy_clim(m,:,k),errint)
#ifdef BFM_GOTM
               if ( bio_calc) then
                 bio_bs=0
                 do ib=1,n_bio
                   if ( cc3d_id(ib) > 0 ) then
                     bio_bs=bio_bs+1 
                     err = nf90_get_var(ncid,cc3d_id(ib),wrk,start=start,count=edges)
                     if (err .ne. NF90_NOERR) go to 10
                     call interpol(0,zax_len,zlev,wrk,H(i,j),kmax,hn(i,j,:), &
                             cc3d_bdy_clim(bio_bs,m,:,k),errint)
                   endif
                 enddo
               endif
#endif
               k = k+1
            end do
         end do

         do n = 1,NNB
            l = l+1
            k = bdy_index(l)
            j = nj(n)
            do i = nfi(n),nli(n)
               if (from_3d_fields) then
                  start(1) = i+ioff ; start(2) = j+joff
               else
                  start(2) = k
               end if
               err = nf90_get_var(ncid,salt_id,wrk,start=start,count=edges)
               if (err .ne. NF90_NOERR) go to 10
               call interpol(0,zax_len,zlev,wrk,H(i,j),kmax,hn(i,j,:), &
                             S_bdy_clim(m,:,k),errint)
               err = nf90_get_var(ncid,temp_id,wrk,start=start,count=edges)
               if (err .ne. NF90_NOERR) go to 10
               call interpol(0,zax_len,zlev,wrk,H(i,j),kmax,hn(i,j,:), &
                             T_bdy_clim(m,:,k),errint)
#ifdef BFM_GOTM
               if ( bio_calc) then
                 bio_bs=0
                 do ib=1,n_bio
                   if ( cc3d_id(ib) > 0 ) then
                     bio_bs=bio_bs+1
                     err = nf90_get_var(ncid,cc3d_id(ib),wrk,start=start,count=edges)
                     if (err .ne. NF90_NOERR) go to 10
                     call interpol(0,zax_len,zlev,wrk,H(i,j),kmax,hn(i,j,:), &
                             cc3d_bdy_clim(bio_bs,m,:,k),errint)
                   endif
                 enddo
               endif
#endif
               k = k+1
            end do
         end do

         do n=1,NEB
            l = l+1
            k = bdy_index(l)
            i = ei(n)
            do j=efj(n),elj(n)
               if (from_3d_fields) then
                  start(1) = i+ioff ; start(2) = j+joff
               else
                  start(2) = k
               end if
               err = nf90_get_var(ncid,salt_id,wrk,start=start,count=edges)
               if (err .ne. NF90_NOERR) go to 10
               call interpol(0,zax_len,zlev,wrk,H(i,j),kmax,hn(i,j,:), &
                             S_bdy_clim(m,:,k),errint)
               err = nf90_get_var(ncid,temp_id,wrk,start=start,count=edges)
               if (err .ne. NF90_NOERR) go to 10
               call interpol(0,zax_len,zlev,wrk,H(i,j),kmax,hn(i,j,:), &
                             T_bdy_clim(m,:,k),errint)
#ifdef BFM_GOTM
               if ( bio_calc) then
                 bio_bs=0
                 do ib=1,n_bio
                   if ( cc3d_id(ib) > 0 ) then
                     bio_bs=bio_bs+1
                     err = nf90_get_var(ncid,cc3d_id(ib),wrk,start=start,count=edges)
                     if (err .ne. NF90_NOERR) go to 10
                     call interpol(0,zax_len,zlev,wrk,H(i,j),kmax,hn(i,j,:), &
                             cc3d_bdy_clim(bio_bs,m,:,k),errint)
                   endif
                 enddo
               endif
#endif
               k = k+1
            end do
         end do

         do n = 1,NSB
            l = l+1
            k = bdy_index(l)
            j = sj(n)
            do i = sfi(n),sli(n)
               if (from_3d_fields) then
                  start(1) = i+ioff ; start(2) = j+joff
               else
                  start(2) = k
               end if
               err = nf90_get_var(ncid,salt_id,wrk,start=start,count=edges)
               if (err .ne. NF90_NOERR) go to 10
               call interpol(0,zax_len,zlev,wrk,H(i,j),kmax,hn(i,j,:), &
                             S_bdy_clim(m,:,k),errint)
               err = nf90_get_var(ncid,temp_id,wrk,start=start,count=edges)
               if (err .ne. NF90_NOERR) go to 10
               call interpol(0,zax_len,zlev,wrk,H(i,j),kmax,hn(i,j,:), &
                             T_bdy_clim(m,:,k),errint)
#ifdef BFM_GOTM
               if ( bio_calc) then
                 bio_bs=0
                 do ib=1,n_bio
                   if ( cc3d_id(ib) > 0 ) then
                     bio_bs=bio_bs+1
                     err = nf90_get_var(ncid,cc3d_id(ib),wrk,start=start,count=edges)
                     if (err .ne. NF90_NOERR) go to 10
                     call interpol(0,zax_len,zlev,wrk,H(i,j),kmax,hn(i,j,:), &
                             cc3d_bdy_clim(bio_bs,m,:,k),errint)
                   endif
                 enddo
               endif
#endif
               k = k+1
            end do
         end do
      end do
      err = nf90_close(ncid)

      test_no_gradient=minval(cc3d_bdy_clim(:,:,:,:))
      no_gradient_ico_missing_values= (test_no_gradient<_ZERO_)
      STDERR "test_no_gradient=",test_no_gradient
      STDERR "no_gradient_ico_missing_value=",no_gradient_ico_missing_values
      
   else

     if (from_3d_fields) then
         FATAL 'non-climatology bdy data only support special bdy data file'
         stop 'init_3d_bdy_ncdf'
      end if

      err = nf90_inq_varid(ncid,'time',time_id)
      if (err .NE. NF90_NOERR) go to 10

      err =  nf90_get_att(ncid,time_id,'units',units)
      if (err .NE. NF90_NOERR) go to 10

      allocate(bdy_times(time_len),stat=err)
      if (err /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (bdy_times)'

      err = nf90_get_var(ncid,time_id,bdy_times)
      if (err .NE. NF90_NOERR) go to 10

      call string_to_julsecs(units,j1,s1)
      offset = time_diff(julianday,secondsofday,j1,s1)
      if( offset .lt. bdy_times(1) ) then
         FATAL 'Model simulation starts before available boundary data'
         call write_time_string(julianday,secondsofday,tbuf)
         FATAL 'Simulation starts: ',tbuf
         call add_secs(j1,s1,nint(bdy_times(1)),j2,s2)
         call write_time_string(j2,s2,tbuf)
         FATAL 'Datafile starts:   ',tbuf
         stop 'init_3d_bdy_ncdf'
      else
         LEVEL3 'Boundary offset time ',offset
      end if

!     check if the bdy data file is long enough
      if( time_diff(juln,secsn,j1,s1) .gt. bdy_times(time_len) ) then
         FATAL 'Not enough 3D boundary data in file'
         call write_time_string(juln,secsn,tbuf)
         FATAL 'Simulation ends: ',tbuf
         call add_secs(j1,s1,nint(bdy_times(time_len)),j2,s2)
         call write_time_string(j2,s2,tbuf)
         FATAL 'Datafile ends:   ',tbuf
         stop 'init_3d_bdy_ncdf'
      end if

!      allocate(h_one(0:kmax),stat=err)
!    if (err /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (h_one)'
!    allocate(ST_one(0:kmax),stat=err)
!     if (err /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (ST_one)'
      allocate(T_old(0:kmax,1:nsbv),stat=err)
      if (err /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (T_old)'
      allocate(T_new(0:kmax,1:nsbv),stat=err)
      if (err /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (T_new)'
      allocate(T_wrk(zax_len,1:nsbv),stat=err)
      if (err /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (T_wrk)'

      allocate(S_old(0:kmax,1:nsbv),stat=err)
      if (err /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (S_old)'
      allocate(S_new(0:kmax,1:nsbv),stat=err)
      if (err /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (S_new)'
      allocate(S_wrk(zax_len,1:nsbv),stat=err)
      if (err /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (S_wrk)'

#ifdef BFM_GOTM
      if ( bio_calc) then 
        allocate(cc3d_bdy_old(1:bio_bs,0:kmax,1:nsbv),stat=rc)
        if (rc /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (cc3d_old)'
        allocate(cc3d_bdy_new(1:bio_bs,0:kmax,1:nsbv),stat=rc)
        if (rc /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (cc3d_bdy_new)'
        allocate(cc3d_bdy_wrk(1:bio_bs,zax_len,1:nsbv),stat=rc)
        if (rc /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (cc3d_bdy_wrk)'
        allocate(cc3d_bdy_tmp(zax_len,1:nsbv),stat=rc)
        if (rc /= 0) stop 'init_3d_bdy_ncdf: Error allocating memory (cc3d_bdy_tmp)'
      endif
#endif

!     Note(KK): We read in at once the data of all points
!               but only for the current time stag
      n = size(bdy_times)
      do i=1,n
         if(bdy_times(i) .ge. real(offset)) then
            EXIT
         end if
      end do

      if(i .gt. 1 .and. bdy_times(i) .gt. real(offset)) then
         i = i-1
      end if

      start(1) = 1; edges(1) = zax_len;
      start(2) = 1; edges(2) = nsbv;
      start(3) = i; edges(3) = 1

      err = nf90_get_var(ncid,temp_id,T_wrk,start=start,count=edges)
      if (err .ne. NF90_NOERR) go to 10

      err = nf90_get_var(ncid,salt_id,S_wrk,start=start,count=edges)
      if (err .ne. NF90_NOERR) go to 10

#ifdef BFM_GOTM
      bio_bs=0
      do ib=1,n_bio
        if ( cc3d_id(ib) > 0 ) then
          bio_bs=bio_bs+1
          err = nf90_get_var(ncid,cc3d_id(ib),cc3d_bdy_tmp, &
                                             start=start,count=edges)
          if (err .ne. NF90_NOERR) go to 11
          cc3d_bdy_wrk(bio_bs,:,:)=cc3d_bdy_tmp
        endif
      enddo
#endif
      l = 0
      do n=1,NWB
         l = l+1
         k = bdy_index(l)
         i = wi(n)
         do j=wfj(n),wlj(n)
            call interpol(0,zax_len,zlev,T_wrk(:,k),H(i,j),kmax,hn(i,j,:), &
                          T_new(:,k),errint)
            call interpol(0,zax_len,zlev,S_wrk(:,k),H(i,j),kmax,hn(i,j,:), &
                          S_new(:,k),errint)
!JM added:
#ifdef BFM_GOTM
            bio_bs=0
            do ib=1,n_bio
              if ( cc3d_id(ib) > 0 ) then
                bio_bs=bio_bs+1
                call interpol(1,zax_len,zlev,cc3d_bdy_wrk(bio_bs,:,k),H(i,j),&
                          kmax,hn(i,j,:),cc3d_bdy_new(bio_bs,:,k),errint)
                if (errint.gt.0) cc3d_bdy_new(bio_bs,:,k)=cc3d(i,j,k,bio_bs)
              endif
            enddo
#endif
            k = k+1
         end do
      end do

      do n = 1,NNB
         l = l+1
         k = bdy_index(l)
         j = nj(n)
         do i = nfi(n),nli(n)
            call interpol(0,zax_len,zlev,S_wrk(:,k),H(i,j),kmax,hn(i,j,:), &
                          S_new(:,k),errint)
            call interpol(0,zax_len,zlev,T_wrk(:,k),H(i,j),kmax,hn(i,j,:), &
                          T_new(:,k),errint)
!JM added:
#ifdef BFM_GOTM
            bio_bs=0
            do ib=1,n_bio
              if ( cc3d_id(ib) > 0 ) then
                bio_bs=bio_bs+1
                call interpol(1,zax_len,zlev,cc3d_bdy_wrk(bio_bs,:,k),H(i,j),&
                            kmax,hn(i,j,:),cc3d_bdy_new(bio_bs,:,k),errint)
                if (errint.gt.0) cc3d_bdy_new(bio_bs,:,k)=cc3d(i,j,k,bio_bs)
              endif
            enddo
#endif
            k = k+1
         end do
      end do

      do n=1,NEB
         l = l+1
         k = bdy_index(l)
         i = ei(n)
         do j=efj(n),elj(n)
            call interpol(0,zax_len,zlev,S_wrk(:,k),H(i,j),kmax,hn(i,j,:), &
                          S_new(:,k),errint)
            call interpol(0,zax_len,zlev,T_wrk(:,k),H(i,j),kmax,hn(i,j,:), &
                          T_new(:,k),errint)
!JM added:
#ifdef BFM_GOTM
            bio_bs=0
            do ib=1,n_bio
              if ( cc3d_id(ib) > 0 ) then
                bio_bs=bio_bs+1
                call interpol(1,zax_len,zlev,cc3d_bdy_wrk(bio_bs,:,k),H(i,j),&
                           kmax,hn(i,j,:),cc3d_bdy_new(bio_bs,:,k),errint)
                if (errint.gt.0) cc3d_bdy_new(bio_bs,:,k)=cc3d(i,j,k,bio_bs)
              endif
            enddo
#endif
            k = k+1
         end do
      end do

      do n = 1,NSB
         l = l+1
         k = bdy_index(l)
         j = sj(n)
         do i = sfi(n),sli(n)
            call interpol(0,zax_len,zlev,S_wrk(:,k),H(i,j),kmax,hn(i,j,:), &
                          S_new(:,k),errint)
            call interpol(0,zax_len,zlev,T_wrk(:,k),H(i,j),kmax,hn(i,j,:), &
                          T_new(:,k),errint)
!JM added:
#ifdef BFM_GOTM
            bio_bs=0
            do ib=1,n_bio
              if ( cc3d_id(ib) > 0 ) then
                bio_bs=bio_bs+1
                call interpol(1,zax_len,zlev,cc3d_bdy_wrk(bio_bs,:,k),H(i,j),&
                           kmax,hn(i,j,:),cc3d_bdy_new(bio_bs,:,k),errint)
                if (errint.gt.0) cc3d_bdy_new(bio_bs,:,k)=cc3d(i,j,k,bio_bs)
              endif
            enddo
#endif
            k = k+1
         end do
      end do
   end if



#ifdef DEBUG
   write(debug,*) 'Leaving init_3d_bdy_ncdf()'
   write(debug,*)
#endif
   return
10 FATAL 'init_3d_bdy_ncdf: ',nf90_strerror(err)
   stop
#ifdef BFM_GOTM
11 STDERR 'var_name=',var_names(ib), "text=",text
   FATAL 'init_3d_bdy_ncdf: ',nf90_strerror(err)
   stop
#endif
   end subroutine init_3d_bdy_ncdf
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: do_3d_bdy_ncdf -
!
! !INTERFACE:
   subroutine do_3d_bdy_ncdf(loop)
!
! !DESCRIPTION:
!  kurt,kurt
!
! !USES:
   use time, only: day,month,secondsofday,days_in_mon,leapyear,secsprday
#ifdef BFM_GOTM
   use modify_meteo_input,only:do_modify_temp_boundary
   use mem,only: ppZ4c
#ENDIF
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: loop
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   integer         :: err
   REALTYPE        :: rat
   integer         :: monthsecs,prev,this,next
   logical, save   :: first=.true.
   integer, save   :: loop0
   REALTYPE        :: t
   REALTYPE, save  :: t1=_ZERO_,t2=-_ONE_
   integer         :: i,j,k,l,n
   integer         :: ib,errint
   integer         :: bio_bs=0
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   write(debug,*) 'do_3d_bdy_ncdf (NetCDF) climatology=',climatology
#endif
   if ( climatology ) then
      if (time_len .eq. 12) then
         this = month
         monthsecs = secsprday*days_in_mon(leapyear,month)
         rat=((day-1)*secsprday+secondsofday)/float(monthsecs)
         next=this+1
         if (next .gt. time_len) next=1
         prev=this-1
         if (prev .eq. 0) prev=time_len
      else
         STDERR 'do_3d_bdy_ncdf: climatology time_len .ne. 12'
         stop
      end if

      S_bdy=(1.-rat)*0.5*(S_bdy_clim(prev,:,:)+S_bdy_clim(this,:,:))  &
         +     rat*0.5*(S_bdy_clim(next,:,:)+S_bdy_clim(this,:,:))
      T_bdy=(1.-rat)*0.5*(T_bdy_clim(prev,:,:)+T_bdy_clim(this,:,:))  &
         +     rat*0.5*(T_bdy_clim(next,:,:)+T_bdy_clim(this,:,:))
#ifdef BFM_GOTM
      call do_modify_temp_boundary(T_bdy,kmax,nsbv)
      if ( trace) then
         cc3d_bdy(1:n_bio,:,:)=(1.-rat)*0.5*(cc3d_bdy_clim(1:n_bio,prev,:,:)  &
                  +cc3d_bdy_clim(1:n_bio,this,:,:))   &
                  +rat*0.5*(cc3d_bdy_clim(1:n_bio,next,:,:) &
                  +cc3d_bdy_clim(1:n_bio,this,:,:))
      elseif ( bio_calc) then
        bio_bs=0
        do ib=1,n_bio
         if ( cc3d_id(ib) > 0 ) then
            bio_bs=bio_bs+1 
            cc3d_bdy(bio_bs,:,:)=(1.-rat)*0.5*(cc3d_bdy_clim(bio_bs,prev,:,:)  &
                                  +cc3d_bdy_clim(bio_bs,this,:,:))   &
                      +rat*0.5*(cc3d_bdy_clim(bio_bs,next,:,:) &
                               +cc3d_bdy_clim(bio_bs,this,:,:))
         endif
        enddo
      endif
#endif
   else

      if (first) then
         loop0=loop-1
      endif
      t = (loop-loop0)*dtm
!     call write_time_string()
!     STDERR timestr,':-reading-3D-boundary-data ...'
!     STDERR "first,dtm,t,t2,climatology",first,dtm,t,t2,loop0,loop,climatology

      if(t .gt. t2 .or. first) then

         if (first) then
            first = .false.
            t2=t
         else
            call write_time_string()
            LEVEL2 timestr,': reading 3D boundary data ...'
         end if

         n = size(bdy_times)
         do i=1,n
            if(bdy_times(i) .ge. real(t + offset)) then
               EXIT
            end if
         end do
         start(1) = 1; edges(1) = zax_len;
         start(2) = 1; edges(2) = nsbv;
         start(3) = i; edges(3) = 1

         t1=t2
         t2 = bdy_times(i) - offset

         T_old = T_new
         S_old = S_new

         err = nf90_get_var(ncid,temp_id,T_wrk,start=start,count=edges)
         if (err .ne. NF90_NOERR) go to 10

         err = nf90_get_var(ncid,salt_id,S_wrk,start=start,count=edges)
         if (err .ne. NF90_NOERR) go to 10

#ifdef BFM_GOTM
         cc3d_bdy_old=cc3d_bdy_new   !JM added

         bio_bs=0
         do ib=1,n_bio
           if ( cc3d_id(ib) > 0 ) then
             bio_bs=bio_bs+1
             err = nf90_get_var(ncid,cc3d_id(ib),cc3d_bdy_tmp,  &
                                           start=start,count=edges)
             if (err .ne. NF90_NOERR) go to 11
             cc3d_bdy_wrk(bio_bs,:,:)=cc3d_bdy_tmp
           endif
         enddo
#endif
         l = 0
         do n=1,NWB
           l = l+1
           k = bdy_index(l)
           i = wi(n)
           do j=wfj(n),wlj(n)
             call interpol(1,zax_len,zlev,T_wrk(:,k),H(i,j),kmax,hn(i,j,:), &
                           T_new(:,k),errint)
             call interpol(1,zax_len,zlev,S_wrk(:,k),H(i,j),kmax,hn(i,j,:), &
                             S_new(:,k),errint)
!JM added:
#ifdef BFM_GOTM
             bio_bs=0
             do ib=1,n_bio
               if ( cc3d_id(ib) > 0 ) then
                bio_bs=bio_bs+1
              call interpol(1,zax_len,zlev,cc3d_bdy_wrk(bio_bs,:,k),H(i,j),&
                            kmax,hn(i,j,:),cc3d_bdy_new(bio_bs,:,k),errint)
                 endif
               enddo
#endif
               k = k+1
            end do
         end do

         do n = 1,NNB
           l = l+1
           k = bdy_index(l)
           j = nj(n)
           do i = nfi(n),nli(n)
             call interpol(1,zax_len,zlev,S_wrk(:,k),H(i,j),kmax,hn(i,j,:), &
                            S_new(:,k),errint)
             call interpol(1,zax_len,zlev,T_wrk(:,k),H(i,j),kmax,hn(i,j,:), &
                             T_new(:,k),errint)
!JM added:
#ifdef BFM_GOTM
             bio_bs=0
             do ib=1,n_bio
               if ( cc3d_id(ib) > 0 ) then
                bio_bs=bio_bs+1
              call interpol(1,zax_len,zlev,cc3d_bdy_wrk(bio_bs,:,k),H(i,j),&
                          kmax,hn(i,j,:),cc3d_bdy_new(bio_bs,:,k),errint)
               endif
             enddo
#endif
               k = k+1
            end do
         end do

         do n=1,NEB
           l = l+1
           k = bdy_index(l)
           i = ei(n)
           do j=efj(n),elj(n)
             call interpol(1,zax_len,zlev,S_wrk(:,k),H(i,j),kmax,hn(i,j,:), &
                             S_new(:,k),errint)
             call interpol(1,zax_len,zlev,T_wrk(:,k),H(i,j),kmax,hn(i,j,:), &
                             T_new(:,k),errint)
!JM added:
#ifdef BFM_GOTM
             bio_bs=0
             do ib=1,n_bio
               if ( cc3d_id(ib) > 0 ) then
                bio_bs=bio_bs+1
                call interpol(1,zax_len,zlev,cc3d_bdy_wrk(bio_bs,:,k),H(i,j),&
                             kmax,hn(i,j,:),cc3d_bdy_new(bio_bs,:,k),errint)
                 endif
               enddo
#endif
               k = k+1
            end do
         end do

         do n = 1,NSB
           l = l+1
           k = bdy_index(l)
           j = sj(n)
           do i = sfi(n),sli(n)
             call interpol(1,zax_len,zlev,S_wrk(:,k),H(i,j),kmax,hn(i,j,:), &
                           S_new(:,k),errint)
             call interpol(1,zax_len,zlev,T_wrk(:,k),H(i,j),kmax,hn(i,j,:), &
                           T_new(:,k),errint)
!JM added:
#ifdef BFM_GOTM
             bio_bs=0
             do ib=1,n_bio
               if ( cc3d_id(ib) > 0 ) then
                 bio_bs=bio_bs+1
                 call interpol(1,zax_len,zlev,cc3d_bdy_wrk(bio_bs,:,k),H(i,j),&
                            kmax,hn(i,j,:),cc3d_bdy_new(bio_bs,:,k),errint)
                 endif
               enddo
#endif
               k = k+1
            end do
         end do
      endif

      T_bdy = T_old + (T_new - T_old)*(t-t1)/(t2-t1)
      S_bdy = S_old + (S_new - S_old)*(t-t1)/(t2-t1)
#ifdef BFM_GOTM
      forall(i=1:n_bio,j=0:kmax,k=1:nsbv,cc3d_bdy_new(i,j,k)<_ZERO_) 
        cc3d_bdy_new(i,j,k)=cc3d_bdy_old(i,j,k)
      end forall
      forall(i=1:n_bio,j=0:kmax,k=1:nsbv) 
        cc3d_bdy(i,j,k) = cc3d_bdy_old(i,j,k) + (cc3d_bdy_new(i,j,k) &
                                       - cc3d_bdy_old(i,j,k))*(t-t1)/(t2-t1)
       end forall
#endif

   end if

#ifdef DEBUG
   write(debug,*) 'Leaving do_3d_bdy_ncdf()'
   write(debug,*)
#endif
   return
10 FATAL 'do_3d_bdy_ncdf: ',nf90_strerror(err)
   stop
#ifdef BFM_GOTM
11 STDERR 'var_name=',var_names(ib)
   FATAL 'do_3d_bdy_ncdf: ',nf90_strerror(err)
   stop
#endif
   end subroutine do_3d_bdy_ncdf
!EOC

!-----------------------------------------------------------------------

! quick and dirty - should be merged with kbk_interpol.F90 and
! grid_interpol.F90

   subroutine interpol(mode,nlev,zlev,wrk,depth,kmax,zm,col,ierr)

! !INPUT PARAMETERS:
   integer, intent(in)       :: mode
   integer, intent(in)       :: nlev
   REAL_4B, intent(in)       :: zlev(nlev),wrk(nlev)
   REALTYPE, intent(in)      :: depth
   integer, intent(in)       :: kmax
   REALTYPE, intent(in)      :: zm(0:kmax)

! !OUTPUT PARAMETERS:
   REALTYPE, intent(INOUT)     :: col(0:kmax)
   integer,intent(OUT)         :: ierr

! !LOCAL VARIABLES:

   REALTYPE                 :: pcol(0:kmax)
   REALTYPE                 :: zmodel(1:kmax),rat
   REAL_4B                  :: missing_value=-9999.0
   REAL_4B                  :: test
   integer                  :: k,li,n,nn,flag

   test=missing_value +1.0;
   ierr=0
   zmodel(1) = -depth + 0.5*zm(1)
   do k=2,kmax
      zmodel(k) = zmodel(k-1) + 0.5*(zm(k-1)+zm(k))
   end do

   do k=kmax,1,-1
      if (zmodel(k) .ge. zlev(1)) then
           if ( wrk(1)> test ) then
              pcol(k) = wrk(1)
           else
              ierr=1;pcol(k)=-1.0D+80
           endif
      endif
   end do

   do k=1,kmax
      if (zmodel(k) .le. zlev(nlev)) then
           if ( wrk(nlev)> test )then
               pcol(k) = wrk(nlev)
           else
               pcol(k)=-1.0D+80;ierr=1
           endif
      endif
   end do

   nn=nlev
   do k=1,kmax
      if (zmodel(k) .gt. zlev(nlev) .and. zmodel(k) .lt. zlev(1)) then
         do while (zlev(nn).le.zmodel(k))
              nn=nn-1;
         enddo
         if ( wrk(nn+1)>test.and. wrk(nn)>test ) then
           rat = (zmodel(k)-zlev(nn+1))/(zlev(nn)-zlev(nn+1))
           pcol(k) = (_ONE_-rat)*wrk(nn+1)+rat*wrk(nn)
         else 
           pcol(k)=-1.0D+80;ierr=1
         endif
      end if
   end do
   pcol(0)=pcol(1)
   if (mode ==0 .or.(mode==1.and.ierr==0)) col=pcol
   end subroutine interpol


!-----------------------------------------------------------------------

   end module ncdf_3d_bdy

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
