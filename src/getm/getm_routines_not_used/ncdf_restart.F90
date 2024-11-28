!$Id: ncdf_restart.F90,v 1.2 2007-10-19 07:52:36 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: Encapsulate netCDF restart quantities
!
! !INTERFACE:
   module ncdf_restart
!
! !DESCRIPTION:
!  This module and the related *\_restart\_ncdf() subroutines provide a
!  drop-in replacement for the binary file hotstart facility in GETM.
!  The main reason for using NetCDF formatted hotstart files instead of
!  binary format is the abillity to use standard tools (nco, ncmerge) is
!  a much easier way to to introduce a new subdomain decomposition for
!  an already running set-up - without having to start all over again.
!  See $read\_restart\_ncdf()$ for further explanation.\newline
!  This modules just contains variables shared accros the *\_restart\_ncdf()
!  routines.
!
!  !BFM
!    integers added to support better BFM 
!
! !USES:
   use output
   IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
   integer                             :: ncid=-1
   integer                             :: xdim_id=-1
   integer                             :: ydim_id=-1
   integer                             :: zdim_id=-1
   integer                             :: xax_id
   integer                             :: yax_id
   integer                             :: zax_id
   integer                             :: loop_id
   integer                             :: julianday_id
   integer                             :: secondsofday_id
   integer                             :: timestep_id
   integer                             :: z_id,zo_id
   integer                             :: U_id,zu_id
   integer                             :: SlUx_id,Slru_id
   integer                             :: V_id,zv_id
   integer                             :: SlVx_id,Slrv_id
#ifndef NO_3D
   integer                             :: ssen_id,ssun_id,ssvn_id
   integer                             :: sseo_id,ssuo_id,ssvo_id
   integer                             :: Uinto_id,Vinto_id
   integer                             :: uu_id,vv_id,ww_id
   integer                             :: uuEx_id,vvEx_id
   integer                             :: tke_id,eps_id
   integer                             :: num_id,nuh_id
   integer                             :: hn_id
#ifndef NO_BAROCLINIC
   integer                             :: T_id,S_id
#endif
#ifdef SPM
   integer                             :: spm_id,spmpool_id
#endif
#ifdef GETM_BIO
   integer                             :: il,ih,iloc,ilen
   integer                             :: jl,jh,jloc,jlen
#ifdef BFM_GOTM
   integer                             :: z2dim_id=-1
   integer                             :: numc_id
   integer                             :: numbc_id
   integer                             :: id_var_start_3d
   integer                             :: id_var_start_2d
   integer                             :: start_id_3d
   integer                             :: start_id_2d
   integer                             :: biodim_id

#else
   integer                             :: biodim_id
   integer                             :: bio_id
#endif

#endif
#endif

   integer                             :: xlen,ylen,zlen
   integer                             :: status
   integer                             :: start(5),edges(5)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------

   end module ncdf_restart

!-----------------------------------------------------------------------
! Copyright (C) 2007 - Karsten Bolding (BBH)                           !
!-----------------------------------------------------------------------
