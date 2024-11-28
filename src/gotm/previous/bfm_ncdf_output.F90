!$Id: $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !FUNCTION: special_dims_1d 
!
! !INTERFACE:
! !INTERFACE:
     integer function special_dims_1d &
         (mode,ncid,n,name,extname,units,lon_dim,lat_dim,time_dim,vars_id)
!
! !DESCRIPTION:
!
!
! !USES:
!  default: all is private.
!
! !PUBLIC MEMBER FUNCTIONS
! store arrays for average computations (pelagic and benthic)

!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:
!
!  $Log: $
!EOP
!-----------------------------------------------------------------------

! !DESCRIPTION:
! Here, the output of biogeochemical parameters either as ascii or as
! NetCDF files is managed.
!
! !USES:
     use bio_bfm, only: calc_sigma_depth
#ifdef NETCDF_FMT
     use bfm_output
     use ncdfout, only: set_attributes,store_data,check_err
     use netcdf
#endif
     IMPLICIT NONE
!
! !INPUT PARAMETERS:
     integer, intent(in)                 :: mode
     integer, intent(in)                 :: ncid
     integer, intent(in)                 :: n
     character(*), intent(in)            :: name
     character(*), intent(in)            :: extname
     character(*), intent(in)            :: units
     integer, intent(in)                 :: lon_dim
     integer, intent(in)                 :: lat_dim
     integer, intent(in)                 :: time_dim
     integer, intent(inout)              :: vars_id
!
! !REVISION HISTORY:
!  Original author(s): Piet Ruardij
!
! !LOCAL VARIABLES:
!    logical, save             :: first=.true.
!    integer, save             :: nn
!    integer                   :: iret
     REALTYPE,parameter        :: ddu=2.0
     REALTYPE                  :: zz !r,s
     integer                   :: dims_4(4)
     integer                   :: i,j,status,altZ_id,dim_altZ
     REALTYPE                  :: arr(0:n)
     character(len=30)         :: altZ,altZ_longname
     character(len=6)          :: dum,alt_unit
!    REAL_4B                   :: vals(2)
!EOP
     select case (mode)
     case (2)
       if ( index(extname,'__Z' )==1 ) then
          j=index(extname,':')-1
          read(extname(1:j),*) dum,altZ, zz,alt_unit, altZ_longname
          status = nf90_inq_dimid(ncid, altZ, dim_altZ)
          if (status.ne.NF90_NOERR) then
            status=nf90_def_dim(ncid,altZ,n,dim_altZ)
            if (status.eq.NF90_NOERR) then
               dims_4(1)=dim_altZ
               status = nf90_def_var(ncid,altZ,NF90_REAL,dims_4(1:1),altZ_id)
               if (status.eq.NF90_NOERR) then
                  i=len_trim(altZ_longname)
                  i=index(extname(1:j),altZ_longname(1:i))
                  status= set_attributes(ncid,altZ_id,long_name=extname(i:j), &
                                      units=alt_unit,missing_value=-9999.0D+00)
                  arr=_ZERO_
                  call calc_sigma_depth(n,ddu,zz,arr(1:n))
                  status = nf90_enddef(ncid=ncid)
                  call check_err(status)
                  status = store_data(ncid,altZ_id,Z_SHAPE,n,array=arr)
                  call check_err(status)
                  status = nf90_redef(ncid)
                  call check_err(status)
               endif
            endif
          endif
          dims_4(1)=lon_dim;dims_4(2)=lat_dim
          dims_4(3)=dim_altZ;dims_4(4)=time_dim
          status = nf90_def_var(ncid,name,NF90_REAL,dims_4(1:4),vars_id)
          status= set_attributes(ncid,vars_id,long_name=trim(extname(j+2:)))
          status= set_attributes(ncid,vars_id,units=units, &
                                                  missing_value=-9999.0D+00)
          special_dims_1d=1
       else
          special_dims_1d=0
       endif
      end select
     end function special_dims_1d

     integer function set_att_1d_averaged(ncid,vars_id)
#ifdef NETCDF_FMT
       use netcdf
#endif
     IMPLICIT NONE
     integer, intent(in)                 :: ncid
     integer, intent(inout)              :: vars_id
     integer                   :: i
     integer                   :: status
     REAL_4B                   :: vals(2)
       vals(1) = 1.0
       i = nf90_put_att(ncid,vars_id,'averaged',vals(1:1))
       status=1;if (i.eq.NF90_NOERR) status=0
       set_att_1d_averaged=status
     return
     end function set_att_1d_averaged
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: Storing the results
!
! !INTERFACE:
     subroutine bio_1d_save
!
! !DESCRIPTION:
! Here, the output of biogeochemical parameters either as ascii or as
! NetCDF files is managed.
!
! !USES:
     use bio_var,only: nlev,bio_model,bio_setup,h_l,dt,c1dimz,cc,ccb, &
             diag,diagb,numc ! ,numbc !, &
     use output, only: out_fmt,ts
     use bfm_output
#ifdef INCLUDE_DIAGNOS_PRF
     use bio_var,only:diagb_prf,numbc_prf
#endif
#ifdef NETCDF_FMT
     use netcdf
     use ncdfout, only: ncid
     use ncdfout, only: lon_dim,lat_dim,z_dim,time_dim,dims_4
     use ncdfout, only: define_mode,new_nc_variable,set_attributes,store_data
#endif
     IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
     logical, save             :: first=.true.
     integer, save             :: nn
     integer                   :: iret
     integer                   :: out_unit=67
     REALTYPE                  :: zz
     logical                   :: llcalc
     integer                   :: i,j,n

     integer,external          :: special_dims_1d
     integer,external          :: set_att_1d_averaged
!EOP
!-----------------------------------------------------------------------
!BOC
     select case (out_fmt)
       case (ASCII)
         if(first) then
            open(out_unit,file='bio.out',status='unknown')
            nn = ubound(cc(:,1),1)
            first = .false.
         end if
         write(out_unit,*)
         write(out_unit,*) trim(ts)
         zz = _ZERO_
         do i=nn,1,-1
           zz=zz+0.5*h_l(i)
           write(out_unit,'(F10.4,100(1x,E10.4E2))') zz,(cc(i,j),j=1,numc)
           zz=zz+0.5*h_l(i)
         end do
       case (NETCDF)
#ifdef NETCDF_FMT
#ifdef BFM_GOTM
        nn = ubound(cc(:,1),1)
        if (bio_model==6) then
          if (first) then
            first = .false.
            iret = define_mode(ncid,.true.)
            if (bio_setup/=2) then
              dims_4(1) = lon_dim; dims_4(2) = lat_dim
              dims_4(3) = z_dim  ; dims_4(4) = time_dim
              do n=stPelStateS,stPelFluxE
                if ( var_ids(n) /= 0 ) then
                  iret = new_nc_variable(ncid,var_names(n),NF90_REAL, &
                                            4,dims_4(1:4),var_ids(n))
                  iret = set_attributes(ncid,var_ids(n), &
                                units=var_units(n), &
                                long_name=var_long(n),missing_value=-9999.0D+00)
                  if (var_ave(n) )j= set_att_1d_averaged(ncid,var_ids(n))
                end if
              end do
            end if
            if (bio_setup>1) then ! define benthic variables
              dims_4(1) = lon_dim; dims_4(2) = lat_dim
              dims_4(3) = time_dim
              do n=stBenStateS,stBenFluxE
                if ( var_ids(n) /= 0 ) then
                  iret = new_nc_variable(ncid,var_names(n),NF90_REAL, &
                                       3,dims_4(1:3),var_ids(n))
                  iret = set_attributes(ncid,var_ids(n), &
                              units=var_units(n), &
                              long_name=var_long(n),missing_value=-9999.0D+00)
                endif
                 if (var_ave(n) )j= set_att_1d_averaged(ncid,var_ids(n))
              end do
            end if
#ifdef INCLUDE_DIAGNOS_PRF
            STDERR 'bio_1d_save numbc_prf=',numbc_prf
            do n=stPRFDiagS,stPRFDiagE
              j=0;if ( var_ids(n) /= 0 ) &
              j=special_dims_1d(2,ncid,numbc_prf,var_names(n),var_long(n), &
                       var_units(n),lon_dim,lat_dim,time_dim,var_ids(n))
!             write(LOGUNIT,*) 'special_dims_1d j=',j,var_ids(n),var_names(n)
              if (var_ave(n) )j= set_att_1d_averaged(ncid,var_ids(n))
            enddo
#endif
            iret = define_mode(ncid,.false.)
          end if

          do n=stPelStateS,stPelStateE
            if ( (var_ids(n)> 0) .and. (.not.var_ave(n) )) &
              iret = store_data(ncid,var_ids(n),XYZT_SHAPE,nlev, &
                                                  array=cc(1:nlev,n))
          end do
          i=0
          do n=stPelDiagS,stPelDiagE
            i=i+1
            if ( (var_ids(n)> 0).and. (.not.var_ave(n) ) ) &
              iret = store_data(ncid,var_ids(n),XYZT_SHAPE,nlev, &
                                                 array=diag(i,1:nlev))
          end do

          i=0
          do n=stPelFluxS,stPelFluxE
            i=i+1
            if ( (var_ids(n)> 0) .and. (.not.var_ave(n))) then
              call make_flux_output(1,i,0,nlev,dt,c1dimz,llcalc)
              iret = store_data(ncid,var_ids(n),XYZT_SHAPE,nlev, &
                                                        array=c1dimz)
            end if
          end do
          j=0
          do n=stPelStateS,stPelFluxE
            if ( (var_ids(n)> 0) .and.var_ave(n) ) then
              j=j+1
              iret = store_data(ncid,var_ids(n),XYZT_SHAPE,nlev, &
                                               array=cc_ave(1:nlev,j))
            endif
          end do

! storage of benthic variables
! stored as scalar: to be modified if benvar are arrays
          if (bio_setup>1) then
            i=0
            do n=stBenStateS,stBenStateE
              i=i+1
              if ( (var_ids(n)> 0) .and. (.not.var_ave(n))) &
                iret = store_data(ncid,var_ids(n),XYT_SHAPE,1, &
                                                      scalar=ccb(1,i))
            end do
            i=0
            do n=stBenDiagS,stBenDiagE
             i=i+1
             if ( (var_ids(n)> 0) .and. (.not.var_ave(n))) &
               iret = store_data(ncid,var_ids(n),XYT_SHAPE,1, &
                                                      scalar=diagb(i,1))
            end do
            i=0
            do n=stBenFluxS,stBenFluxE
              i=i+1
              if ( (var_ids(n)> 0) .and. (.not.var_ave(n))) then
                call make_flux_output(2,i,0,nlev,dt,c1dimz,llcalc)
                iret = store_data(ncid,var_ids(n),XYT_SHAPE,1, &
                                                      scalar=c1dimz(1))
              endif
            end do
            j=0
            do n=stBenStateS,stBenFluxE
              if ( (var_ids(n)> 0) .and. var_ave(n)) then
                j=j+1
                iret = store_data(ncid,var_ids(n),XYT_SHAPE,1,&
                                                scalar=ccb_ave(1,j))
              endif
            end do
#if INCLUDE_DIAGNOS_PRF
            i=0
            do n=stPRFDiagS,stPRFDiagE
              i=i+1
              if ( (var_ids(n)> 0).and. (.not.var_ave(n) ) ) &
               iret = store_data(ncid,var_ids(n),XYZT_SHAPE,numbc_prf,&
                                        array=diagb_prf(1:numbc_prf,i))
            end do
            do n=stPRFDiagS,stPRFDiagE
              i=i+1
              if ( (var_ids(n)> 0).and. var_ave(n) ) &
               iret = store_data(ncid,var_ids(n),XYZT_SHAPE,numbc_prf,&
                                       array=ccb_ave_prf(1:numbc_prf,i))
            end do
#endif
          end if
       end if !bio_model
#endif
#endif
     case default
         FATAL 'A non valid output format has been chosen'
         stop 'bio_1d_save'
     end select

     return
     end subroutine bio_1d_save
!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------


