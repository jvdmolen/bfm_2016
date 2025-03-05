!$Id: read_restart_ncdf.F90,v 1.6 2007-11-12 13:50:17 kb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Read variables from a GETM NetCDF hotstart file
!
! !INTERFACE:
   subroutine read_restart_bio_ncdf(mode,nr,iout,name)
!
! !BFM:
!   It is assumed that the id-number start at 1 and increases by one.
!   routine to read one pelagic or one benthic state variable:
!      1. name
!      2. values
!   id are calculated from a start point and the sequence number of the state variable
!   this routine is call from output_restart_bio
!
! !USES:
#ifdef BFM_GOTM
   use variables_bio_3d,only: ffp,ffb
   use domain,only:az,kmax
   use domain,only:imin,imax,jmin,jmax,kmax
   use halo_zones, only: update_2d_halo,update_3d_halo,wait_halo
   use halo_zones, only: D_TAG
   use netcdf
   use ncdf_restart
#endif
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)      :: mode
   integer, intent(in)      :: nr
!
! !OUTPUT PARAMETERS:
   integer, intent(out)                  :: iout
   character(*),intent(out)              :: name
!
! !DEFINED PARAMTERS:
!
! !REVISION HISTORY:
!  Original author(s): Piet Ruardij
!
! !LOCAL VARIABLES
#ifdef BFM_GOTM
   integer                 :: dum_id
   integer                 :: start_id
!EOP
!-----------------------------------------------------------------------
!BOC
   select case (mode)
     ! read number of pelagic state variables
     case(1)
       !pelagic/3d variables are located after the hydrodynamical variables
        !to the until the first benthic var.
        if (start_id_3d==0 ) then
           status = nf90_get_var(ncid,numc_id,iout)
        elseif ( start_id_2d >0)  then
          iout=start_id_2d-start_id_3d
        else
          status = nf90_inquire(ncid,nVariables=iout)
          iout=iout-start_id_3d+1
        endif
        STDERR "rest_bio_ncdf 3d:",start_id_3d
     ! read number of benthic state variables
     case(2)
        !benthic/2d variables are always located after the pelagic variables
        !to the end of the file
        if (start_id_2d==0 ) then
            status = nf90_get_var(ncid,numbc_id,iout)
        else
          status = nf90_inquire(ncid,nVariables=iout)
          iout=iout-start_id_2d+1
        endif
        STDERR "rest_bio_ncdf 2d:",start_id_2d
     case(3)
       !read name of a pelagic state variable
       !id is calulated from start_3d_id and the sequence number
       !of the state variable
       start_id=start_id_3d;if (start_id==0) &
           status = nf90_get_var(ncid,id_var_start_3d,start_id)
       iout=start_id-1+nr
       status=nf90_inquire_variable(ncid,iout,name=name)
     case(4)
       !read name of  a benthic state variable
       start_id=start_id_2d;if (start_id==0) &
         status = nf90_get_var(ncid,id_var_start_2d,start_id)
       iout=start_id-1+nr
       status=nf90_inquire_variable(ncid,iout,name=name)
     case(5)
       !read hotstart values of a pelagic state variable
       start_id=start_id_3d;if (start_id==0) &
         status = nf90_get_var(ncid,id_var_start_3d,start_id)
       iout=start_id-1+nr
       status=nf90_inquire_variable(ncid,iout,name=name)
!JM       status = &
!JM          nf90_get_var(ncid,iout,ffp(iloc:ilen,jloc:jlen,0:kmax),start,edges)
!LEVEL1 "read_restart_bio_ncdf read pelabic"
!LEVEL1 "start",start
!LEVEL1 "edges",edges
       status = &
          nf90_get_var(ncid,iout,ffp(iloc:ilen,jloc:jlen,0:kmax),start(1:3),edges(1:3))
       if (status .NE. NF90_NOERR) then
          LEVEL3 "read_restart_bio_ncdf(): setting",name,"=0"
          ffp=_ZERO_
       else
         call update_3d_halo(ffp,ffp,az, &
                          imin,jmin,imax,jmax,kmax,D_TAG)
         call wait_halo(D_TAG)
       endif
     case(6)
       !read hotstart values of a benthic state variable
       start_id=start_id_2d;if (start_id==0) &
         status = nf90_get_var(ncid,id_var_start_2d,start_id)
       iout=start_id-1+nr
       status=nf90_inquire_variable(ncid,iout,name=name)
!JM       status = &
!JM          nf90_get_var(ncid,iout,ffb(iloc:ilen,jloc:jlen),start,edges)
!LEVEL1 "read_restart_bio_ncdf read benthic"
!LEVEL1 "start",start
!LEVEL1 "edges",edges
       status = &
          nf90_get_var(ncid,iout,ffb(iloc:ilen,jloc:jlen),start(1:3),edges(1:3))
       if (status .NE. NF90_NOERR) then
          LEVEL3 "read_restart_bio_ncdf(): setting",name,"=0"
          ffb=_ZERO_
       endif
       return
     case(7)
      status = nf90_close(ncid)
   end select

   if (status .NE. NF90_NOERR) go to 10

   return

10 FATAL 'read_restart__bio_ncdf: ',nf90_strerror(status)

   stop
   return

#endif
   end subroutine read_restart_bio_ncdf
!EOC

!-----------------------------------------------------------------------
! Copyright (C) 2009 - Piet Ruardij (NIOZ)                           !
!-----------------------------------------------------------------------
