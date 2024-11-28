!$Id: create_restart_ncdf.F90,v 1.2 2007-10-19 07:52:36 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Create a GETM NetCDFNetCDF  hotstart file
!
! !INTERFACE:
   subroutine create_restart_ncdf(fname,loop,runtype)
!
! !DESCRIPTION:
!  Creates a new NetCDF formatted file for storing variables necessary
!  to make a correct GETM hotstart. The created file contains dimensions
!  (xax, yax, zax) as well as the (empty) variables. Variables are named
!  corresponding to the names used in the Fortran files. Only the actual
!  domain is stored (i.e. not the halo-zones). This allows easy use of
!  'ncmerge' to stitch a number of hotstart files together to cover the
!  entire computational domain. See read\_restart\_ncdf() for use.
!
! !BFM
!    -code added to initialize the storing of BFM state variables.
!    - two integers (ioff_0 and joff_0) a are defined in which the values 
!      of ioff and joff for subdomain 000  
!      They are define as GLOBAL  attribute
!      This give the possibility to a complete (ncmerged) hoststart-file
!      which are composed of subdomain hotstart-files  from another
!      sub-domain composition.
!        
! !USES:
   use time,only:write_time_string
   use netcdf
   use ncdf_restart
   use domain, only: ioff,joff
   use domain, only: imin,imax,jmin,jmax,kmax
   use domain, only: vert_cord
#ifdef GETM_BIO
   use bio, only: bio_calc
   use bio_var, only: numc
#ifdef BFM_GOTM
   use bio_var, only: numbc
   use bfm_output, only: var_names, &
            stPelStateS,stBenStateS, stPelStateE,stBenStateE
#endif
#endif
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)        :: fname
   integer, intent(in)                 :: loop
   integer, intent(in)                 :: runtype
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
! !LOCAL VARIABLES:
   character(len=80)         :: history,tts
   character(len=80)         :: title
   character(len=80)         :: str 
   integer                   :: i,j,dum_id
   integer                   :: ioff_0,joff_0
!EOP
!-------------------------------------------------------------------------
!BOC
!  create netCDF file
   status = nf90_create(fname, NF90_CLOBBER, ncid)
   if (status .NE. NF90_NOERR) go to 10
   status = nf90_set_fill(ncid, NF90_NOFILL, i)
   if (status .NE. NF90_NOERR) go to 10

!  length of netCDF dimensions
#ifdef _WRITE_HOT_HALOS_
   xlen = (imax+HALO)-(imin-HALO)+1
   ylen = (jmax+HALO)-(jmin-HALO)+1
#else
   xlen = imax-imin+1
   ylen = jmax-jmin+1
#endif
   zlen = kmax+1

   status = nf90_def_dim(ncid, "xax", xlen, xdim_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_dim(ncid, "yax", ylen, ydim_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_dim(ncid, "zax", zlen, zdim_id)
   if (status .NE. NF90_NOERR) go to 10

#ifndef NO_3D
#ifdef GETM_BIO
#ifndef BFM_GOTM
   if (bio_calc) then
      status = nf90_def_dim(ncid, "biodim", numc, biodim_id)
      if (status .NE. NF90_NOERR) go to 10
   end if
#else
   call write_time_string(str=str)
   status= nf90_put_att(ncid,NF90_GLOBAL,"restart_date",str)
   if (status .NE. NF90_NOERR) go to 10

!  ioff_0=(ioff-imax*(ioff/imax));joff_0=(joff-jmax*(joff/jmax))
   status= nf90_put_att(ncid,NF90_GLOBAL,"ioff",ioff)
   if (status .NE. NF90_NOERR) go to 10

   status= nf90_put_att(ncid,NF90_GLOBAL,"joff",joff)
   if (status .NE. NF90_NOERR) go to 10

#endif
#endif
#endif

   status = nf90_def_var(ncid, "loop", nf90_int, loop_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_var(ncid, "julianday", nf90_int, julianday_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_var(ncid, "secondsofday", nf90_int, secondsofday_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_var(ncid, "timestep", nf90_double, timestep_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_var(ncid, "xax", nf90_double, (/ xdim_id /), xax_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_var(ncid, "yax", nf90_double, (/ ydim_id /), yax_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_var(ncid, "zax", nf90_double, (/ zdim_id /), zax_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_var(ncid, "z", nf90_double, &
                            (/ xdim_id, ydim_id /), z_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_var(ncid, "zo", nf90_double, &
                            (/ xdim_id, ydim_id /), zo_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_var(ncid, "U", nf90_double, &
                            (/ xdim_id, ydim_id /), U_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_var(ncid, "SlUx", nf90_double, &
                            (/ xdim_id, ydim_id /), SlUx_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_var(ncid, "Slru", nf90_double, &
                            (/ xdim_id, ydim_id /), Slru_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_var(ncid, "V", nf90_double, &
                            (/ xdim_id, ydim_id /), V_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_var(ncid, "SlVx", nf90_double, &
                            (/ xdim_id, ydim_id /), SlVx_id)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_def_var(ncid, "Slrv", nf90_double, &
                            (/ xdim_id, ydim_id /), Slrv_id)
   if (status .NE. NF90_NOERR) go to 10

#ifndef NO_3D
   if (runtype .ge. 2)  then
      status = nf90_def_var(ncid, "ssen", nf90_double, &
                               (/ xdim_id, ydim_id /), ssen_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "ssun", nf90_double, &
                               (/ xdim_id, ydim_id /), ssun_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "ssvn", nf90_double, &
                               (/ xdim_id, ydim_id /), ssvn_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "sseo", nf90_double, &
                               (/ xdim_id, ydim_id /), sseo_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "ssuo", nf90_double, &
                               (/ xdim_id, ydim_id /), ssuo_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "ssvo", nf90_double, &
                               (/ xdim_id, ydim_id /), ssvo_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "Uinto", nf90_double, &
                               (/ xdim_id, ydim_id /), Uinto_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "Vinto", nf90_double, &
                               (/ xdim_id, ydim_id /), Vinto_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "uu", nf90_double, &
                               (/ xdim_id, ydim_id, zdim_id /), uu_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "vv", nf90_double, &
                               (/ xdim_id, ydim_id, zdim_id /), vv_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "ww", nf90_double, &
                               (/ xdim_id, ydim_id, zdim_id /), ww_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "uuEx", nf90_double, &
                               (/ xdim_id, ydim_id, zdim_id /), uuEx_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "vvEx", nf90_double, &
                               (/ xdim_id, ydim_id, zdim_id /), vvEx_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "tke", nf90_double, &
                               (/ xdim_id, ydim_id, zdim_id /), tke_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "eps", nf90_double, &
                               (/ xdim_id, ydim_id, zdim_id /), eps_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "num", nf90_double, &
                               (/ xdim_id, ydim_id, zdim_id /), num_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "nuh", nf90_double, &
                               (/ xdim_id, ydim_id, zdim_id /), nuh_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "hn", nf90_double, &
                               (/ xdim_id, ydim_id, zdim_id /), hn_id)
      if (status .NE. NF90_NOERR) go to 10

#ifndef NO_BAROCLINIC
      status = nf90_def_var(ncid, "T", nf90_double, &
                               (/ xdim_id, ydim_id, zdim_id /), T_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "S", nf90_double, &
                               (/ xdim_id, ydim_id, zdim_id /), S_id)
      if (status .NE. NF90_NOERR) go to 10
#endif
#ifdef SPM
      status = nf90_def_var(ncid, "spm", nf90_double, &
                               (/ xdim_id, ydim_id, zdim_id /), spm_id)
      if (status .NE. NF90_NOERR) go to 10

      status = nf90_def_var(ncid, "spmpool", nf90_double, &
                               (/ xdim_id, ydim_id /), spmpool_id)
      if (status .NE. NF90_NOERR) go to 10
#endif
#ifdef GETM_BIO
      if (bio_calc) then
#ifdef BFM_GOTM
        status= nf90_put_att(ncid,NF90_GLOBAL,"FirstPelVar", &
                                            trim(var_names(stPelStates)))
        if (status .NE. NF90_NOERR) go to 10
        status= nf90_put_att(ncid,NF90_GLOBAL,"FirstBenVar", &
                                            trim(var_names(stBenStates)))
        if (status .NE. NF90_NOERR) go to 10
        status = nf90_def_var(ncid, "numc", nf90_int, numc_id)
        if (status .NE. NF90_NOERR) go to 10
        status = nf90_def_var(ncid, "numbc", nf90_int, numbc_id)
        if (status .NE. NF90_NOERR) go to 10
        status = nf90_def_var(ncid, "start_3d", nf90_int, id_var_start_3d)
        if (status .NE. NF90_NOERR) go to 10
        status = nf90_def_var(ncid, "start_2d", nf90_int, id_var_start_2d)
        if (status .NE. NF90_NOERR) go to 10
        j=stPelStateS
        do i=1,numc
          status = nf90_def_var(ncid, var_names(j), nf90_double, &
                             (/ xdim_id, ydim_id, zdim_id /), dum_id)
          if (status .NE. NF90_NOERR) go to 11
          status=nf90_put_att(ncid,dum_id,"_FillValue",-9999.0D+00)
          if (status .NE. NF90_NOERR) go to 10
          status=nf90_put_att(ncid,dum_id,"missing_value",-9999.0D+00)
          if (status .NE. NF90_NOERR) go to 10
          status=nf90_put_att(ncid,dum_id,"type","d3")
          if (status .NE. NF90_NOERR) go to 10
          j=j+1
        enddo
        j= stBenStateS
        do i=1,numbc
          status = nf90_def_var(ncid, var_names(j), nf90_double, &
                            (/ xdim_id, ydim_id /),dum_id)
          if (status .NE. NF90_NOERR) go to 11
          status=nf90_put_att(ncid,dum_id,"_FillValue",-9999.0D+00)
          if (status .NE. NF90_NOERR) go to 10
          status=nf90_put_att(ncid,dum_id,"missing_value",-9999.0D+00)
          if (status .NE. NF90_NOERR) go to 10
          status=nf90_put_att(ncid,dum_id,"type","d2")
          if (status .NE. NF90_NOERR) go to 10
          j=j+1
        enddo
#else
         status = nf90_def_var(ncid, "bio", nf90_double, &
                      (/ xdim_id, ydim_id, zdim_id, biodim_id /), bio_id)
         if (status .NE. NF90_NOERR) go to 10
#endif
      endif
#endif
    endif
#endif
!  globals
   title="GETM NetCDF hotstart file"
   status = nf90_put_att(ncid,NF90_GLOBAL,'title',trim(title))
   if (status .NE. NF90_NOERR) go to 10
#ifdef GETM_BIO
   if (bio_calc.and.runtype.ge.2) then
#ifdef BFM_GOTM
     status= nf90_put_att(ncid,NF90_GLOBAL,'FirstPelVar', &
                                            trim(var_names(stPelStates)))
     if (status .NE. NF90_NOERR) go to 10
     status= nf90_put_att(ncid,NF90_GLOBAL,'FirstBenVar', &
                                            trim(var_names(stBenStates)))
     if (status .NE. NF90_NOERR) go to 10
   endif
#endif
#endif

   history = 'Generated by GETM, ver. '//RELEASE
   status = nf90_put_att(ncid,NF90_GLOBAL,'history',trim(history))
   if (status .NE. NF90_NOERR) go to 10

   ! leave define mode
   status = nf90_enddef(ncid)
   if (status .NE. NF90_NOERR) go to 10

   status = nf90_sync(ncid)
   if (status .NE. NF90_NOERR) go to 10

   return

   10 FATAL 'create_restart_ncdf: ',nf90_strerror(status)
   stop 'create_restart_ncdf'
   return
   11 STDERR "var,=",trim(var_names(j)),j 
   FATAL 'create_restart_ncdf: ',nf90_strerror(status)
   stop 'create_restart_ncdf'
   return
   end subroutine create_restart_ncdf
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2007 - Karsten Bolding (BBH)                           !
!-----------------------------------------------------------------------
