!$Id: bio_var.F90,v 1.7 2005-12-02 20:57:27 hb Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bio_var --- declaration of biological variables
!
! !INTERFACE:
   module bio_var
!
! !DESCRIPTION:
!  Here all variables necessary for the biogeochemical models are
!  declared, mostly as allocatable variables.
!
! !USES:
!  default: all is public.
   public
!
! !PUBLIC DATA MEMBERS:
   integer                               :: bio_model
   integer                               :: numc,numcc
   REALTYPE                              :: I_0
   REALTYPE, dimension(:), allocatable   :: zlev
   REALTYPE, dimension(:), allocatable   :: par
   REALTYPE, dimension(:,:), allocatable,target :: cc,ws        !BFM

   REALTYPE, dimension(:,:), allocatable :: cc_before_transport !BFM
   logical, dimension(:),allocatable     :: llws
   integer                               :: surface_flux_method=0
   integer                               :: n_surface_fluxes=-1
   integer                               :: calc_init_bennut_states
   REALTYPE                              :: rel_max_sedi_rate=1.0
   REALTYPE                              :: sfl_N3n=0.0,sfl_N4n=0.0
   REALTYPE, dimension(:), allocatable   :: sfl_read
   REALTYPE, dimension(:), allocatable   :: sfl,bfl
   integer, dimension(:), allocatable    :: posconc
   logical, dimension(:), allocatable    :: mussels_inhale
   logical, dimension(:,:), allocatable  :: particle_active
   integer, dimension(:,:), allocatable  :: particle_indx
   REALTYPE, dimension(:,:), allocatable :: particle_pos

   REALTYPE, parameter                   :: secs_pr_day=86400.0

   integer                               :: nlev
   integer                               :: bio_setup =1        !BFM
   integer,public                        :: pelvar_save_all=0   !BFM
   REALTYPE, dimension(:), allocatable   :: c1dimz
   REALTYPE, dimension(:), pointer   :: nuh_l
   REALTYPE, dimension(:), pointer   :: h_l
   REALTYPE, dimension(:), allocatable   :: t_l
   REALTYPE, dimension(:), allocatable   :: c1dimnumc

#ifdef BFM_GOTM
   REALTYPE             :: bathy_dep
   REALTYPE             :: dt=-1.0
   integer(8)           :: julianday
   ! additional storage variables for benthic and diagnostics
   integer              :: numbc
   integer              :: numc_diag,numbc_diag
   integer              :: numc_flux,numbc_flux

#ifdef INCLUDE_DIAGNOS_PRF
   integer              :: nprf
   integer              :: numbc_prf
#endif

   ! parameter values for the attributes pelvar_type and benvar_type
   integer, parameter   :: SINKSOURCE=-1
   integer, parameter   :: NOTRANSPORT=0
   integer, parameter   :: HORTRANSPORT=10
   integer, parameter   :: ALLTRANSPORT=20
   integer, parameter   :: SILTTRANSPORT=25

   !additional BFM pelagic arrays
   REALTYPE, dimension(:),     allocatable         ::  SSt,RRa
   REALTYPE, dimension(:,:,:), allocatable, target ::  dd,pp

   !benthic BFM arrays
   REALTYPE, dimension(:,:),   allocatable, target :: ccb
   REALTYPE, dimension(:,:,:), allocatable, target :: ddb,ppb

   !diagnostic output arrays  (pelagic and benthic)
   REALTYPE, dimension(:,:),   allocatable, target,public :: diag
   REALTYPE, dimension(:,:),   allocatable, target,public :: diagb

   ! type and save attributes of pelagic and benthic variables
   integer, dimension(:),      allocatable, target :: pelvar_type
   integer, dimension(:),      allocatable, target :: pelvar_type_dyn
   integer, dimension(:),      allocatable, target :: benvar_type

   ! store array for adv_messages
   REALTYPE                              :: adv_courant=0.0
   REALTYPE                              :: adv_number=0.0
   REALTYPE,dimension(:), allocatable    :: adv1d_courant
   REALTYPE,dimension(:), allocatable    :: adv1d_number
#ifdef INCLUDE_DIAGNOS_PRF
   REALTYPE, dimension(:,:),   allocatable, target :: diagb_prf
#endif

#endif 
!BFM

!
!EOP
!-----------------------------------------------------------------------

   end module bio_var

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
