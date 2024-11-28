#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !SUBROUTINE: Initialise variable components
!
! !INTERFACE:
   subroutine init_cnps(c,n,p,s,nc,pc,sc)
!
! !DESCRIPTION:
!  This subroutine initialises the other internal components
!  of biogeochemical variables
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
    REALTYPE,dimension(:),intent(in)           :: c
    REALTYPE,intent(in),optional               :: nc,pc,sc
!
! !OUTPUT PARAMETERS:
    REALTYPE,dimension(:),intent(out),optional :: n
    REALTYPE,dimension(:),intent(out),optional :: p
    REALTYPE,dimension(:),intent(out),optional :: s
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!
!  Revision 1.1 Piet Ruardij
!   extended with 3 "if's"  to make possible to use values from namelist for
!   nc,pc,and sc which are not defined in the input file with values (.nml -file)
!LOCAL VARIABLES:
    REALTYPE                     :: nc_ratio,pc_ratio,sc_ratio
!
!EOP
!-----------------------------------------------------------------------
!BOC

    nc_ratio = 0.0126 ! Redfield
    if (present(nc)) then
      if (nc> 0.0) nc_ratio = nc
    end if

    pc_ratio = 0.7862e-3 ! Redfield
    if (present(pc)) then
      if ( pc> 0.0) pc_ratio = pc
    end if

    sc_ratio = 0.0145 ! Redfield
    if (present(sc)) then
      if ( sc > 0.0) sc_ratio = sc
    end if

    if (present(n)) n = nc_ratio*c
    if (present(p)) p = pc_ratio*c
    if (present(s)) s = sc_ratio*c

  end subroutine init_cnps
!EOC



!-----------------------------------------------------------------------

