!$Id$
#include "cppdefs.h"
!-----------------------------------------------------------------------
! !ROUTINE : save_bfm_2d
!
! !BFM
!   routine to save fields with the same size as as the bathymetry.     
!   all 2d-filed which are saved are special fileds: parameters fileds
!   bathemetry, p_poro and adsorption parameters for phosphate and silica

! !INTERFACE
      subroutine save_bfm_2d(varname,var,vr,units,long_name)

      use domain,only:imin,imax,jmin,jmax
#ifdef BFM_GOTM
      use exceptions,only:netcdf_error,getm_error

      use netcdf
      use ncdf_common
!AN JM obsolete in getm      use ncdf_2d,only:ws
      use ncdf_3d,only:ncid,hh_missing
!  use grid_ncdf
      use domain, only: az
#endif

      implicit none
      character(len=*),intent(IN)   :: varname
      REALTYPE,intent(IN)           :: var(E2DFIELD) 
      REALTYPE,intent(IN)           :: vr(2) 
      character(len=*),intent(IN)   :: units
      character(len=*),intent(IN)   :: long_name

#ifdef BFM_GOTM
!     REAL_4B, dimension(:), allocatable :: ws
      REALTYPE,dimension(E2DFIELD) :: ws !AN
      integer                            :: status
      integer                            :: var_id
      integer                            :: bat_id
      integer                           :: start(2) 
      integer                           :: edges(2) 
      integer                            :: f2_dims(2) 
      integer                            :: xlen,ylen
      REALTYPE                           :: fv,mv
      character(len=80)                  :: msg
#ifdef DEBUG
      STDERR 'save_bfm_2d varname=',varname
#endif
      msg="bathymetry"
      status = nf90_inq_varid(ncid,msg,bat_id)
      if (status .ne. NF90_NOERR) goto 10
!     STDERR 'status,ncid,bat_id',status,ncid,bat_id

      msg="dims bathymetry"
      status = nf90_inquire_variable(ncid,bat_id,dimids=f2_dims)
      if (status .ne. NF90_NOERR) goto 10
!     STDERR 'status,ncid,bat_id',status,ncid,f2_dims

      xlen = imax-imin+1
      ylen = jmax-jmin+1

      msg= "to definemode"
      status = nf90_redef(ncid)
      if (status .ne. NF90_NOERR) goto 10

!     STDERR 'status,ncid,bat_id',status,ncid,f2_dims
      msg="define " // varname
      status = nf90_def_var(ncid,varname,NF90_FLOAT,f2_dims,var_id)
      if (status .ne. NF90_NOERR) goto 10

!     STDERR 'status,ncid,var_id',status,ncid,var_id
      fv = hh_missing
      mv = hh_missing
      call set_attributes(ncid,var_id,                                 &
                          units=units,long_name=long_name, &
                          valid_range=vr,FillValue=fv,missing_value=mv)

!  Set netCDF slice information
      start(1) = 1
      start(2) = 1
      edges(1) = xlen
      edges(2) = ylen

!  save var

      status = nf90_enddef(ncid)
      if (status .ne. NF90_NOERR) goto 10

      call cnv_2d(imin,jmin,imax,jmax,az,var,hh_missing,               &
                  imin,jmin,imax,jmax,ws)

      msg="save "//varname
      status = nf90_put_var(ncid,var_id,ws(_2D_W_),start,edges)
      if (status .ne. NF90_NOERR) goto 10 

      return
  10  call netcdf_error(status,"save_bfm_2d",msg)
      return
  20  call getm_error("save_bfm_2d",msg)
#endif
      end subroutine save_bfm_2d

