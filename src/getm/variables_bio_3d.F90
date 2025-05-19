!$Id: variables_3d.F90,v 1.11 2006-03-17 17:19:54 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: variables_3d - global 3D related variables \label{sec-variables-3d}
!
! !INTERFACE:
   module variables_bio_3d
!

! !USES:
   use domain,     only: imin,imax,jmin,jmax,kmax
   IMPLICIT NONE

#ifdef BFM_GOTM

   REALTYPE                        :: ffp(I3DFIELD)
   REALTYPE                        :: ffb(I2DFIELD)
   REALTYPE                        :: play2d(I2DFIELD)
   REALTYPE                        :: play3d(I3DFIELD)

   logical                         :: counter_reset=.true.
   integer                         :: flag_out=0
   REALTYPE                        :: bio_missing=-9998.0
   integer                         :: n_cc3d_out,n_ccb3d_out
#ifdef INCLUDE_DIAGNOS_PRF
   integer                         :: n_ccb3d_prf
#endif
   REALTYPE                        :: p_poro_default
   REALTYPE                        :: p_pK1_ae_default,p_pK5_ae_default

   integer, allocatable            :: counter_ave(:,:)
   REALTYPE,allocatable            :: p_poro_2d(:,:)
   REALTYPE,allocatable            :: p_pK1_ae_2d(:,:),p_pK5_ae_2d(:,:)
   REALTYPE, allocatable,target    :: cc3d(:,:,:,:)
   REALTYPE, allocatable,target    :: ccb3d(:,:,:,:)
   REALTYPE, allocatable,target    :: cc3d_out(:,:,:,:)
   REALTYPE, allocatable,target    :: ccb3d_out(:,:,:,:)
#ifdef INCLUDE_DIAGNOS_PRF
   REALTYPE, allocatable,target    :: ccb3d_prf(:,:,:,:)
#endif
   REALTYPE,public,allocatable     :: adv3d_courant(:,:,:)
   REALTYPE,public,allocatable     :: adv3d_number(:,:,:)
   integer,public,allocatable      :: sw_CalcPhyto_2d(:,:,:)
   integer,public,allocatable      :: d3_pelvar_type(:)
!JM   REALTYPE, public, dimension(:,:,:), allocatable       :: cut,cvt    !JM added
   REALTYPE, public, dimension(I3DFIELD), target      :: cut,cvt    !JM added

#endif


end module variables_bio_3d

