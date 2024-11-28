!$Id: process_model.F90,v 1.9 2005-12-27 06:51:49 hb Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Initialise the bio module
!
! !INTERFACE:
   subroutine process_model(first,numc,nlev,cc,pp,dd)
!
! !DESCRIPTION:
! This routine is called from the ODE solvers implemented in
! {\tt ode\_solver.F90} and updates the sources and sinks for
! the individual biogeochemical models.
!
! !USES:
   use bio_var, only: bio_model ! ,h_l,t_l
#ifdef BIO_TEMPLATE
   use bio_template, only: do_bio_template
#endif
#ifdef BIO_NPZD
   use bio_npzd, only: do_bio_npzd
#endif
#ifdef BIO_IOW
   use bio_iow, only: do_bio_iow
#endif
#ifdef BIO_MAB
   use bio_mab, only: do_bio_mab
#endif
#ifdef BIO_SED
   use bio_sed, only: do_bio_sed
#endif
#ifdef BIO_FASHAM
   use bio_fasham, only: do_bio_fasham
#endif
#ifdef BFM_GOTM
   use bio_bfm, only: do_bio_bfm
#endif
   IMPLICIT NONE

!
! !INPUT PARAMETERS:
   logical, intent(in)                 :: first
   integer, intent(in)                 :: numc,nlev
   REALTYPE, intent(in)                :: cc(0:nlev,1:numc)
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE, intent(inout)             :: pp(0:nlev,1:numc,1:numc)
   REALTYPE, intent(inout)             :: dd(0:nlev,1:numc,1:numc)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
   select case (bio_model)
#ifdef BIO_TEMPLATE
      case (-1)
         call do_bio_template(numc,nlev)
#endif
#ifdef BIO_NPZD
      case (1)
         call do_bio_npzd(first,numc,nlev,cc,pp,dd)
#endif
#ifdef BIO_IOW
      case (2)
         call do_bio_iow(first,numc,nlev,cc,pp,dd,h_l,t_l)
#endif
#ifdef BIO_SED
      case (3)
         call do_bio_sed(nlev,pp,dd)
#endif
#ifdef BIO_FASHAM
      case (4)
         call do_bio_fasham(first,numc,nlev,cc,pp,dd)
#endif
#ifdef BIO_MAB
      case (5)
         call do_bio_mab(first,numc,nlev,cc,pp,dd,h_l,t_l)
#endif
#ifdef BFM_GOTM
      case (6)
         call do_bio_bfm(first,.false.)
#endif
      case default
         stop "bio: no valid biomodel specified in bio.inp !"
   end select
   return

   end subroutine process_model
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
