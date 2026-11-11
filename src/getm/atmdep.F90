#include "cppdefs.h"

#ifdef BFM_GOTM

module atmdep

   use domain, only: imin, imax, jmin, jmax, kmax, az
   use time, only: timestep
   use m3d, only: M
   use variables_3d, only: hn
   use variables_bio_3d, only: cc3d
   use bio, only: bio_calc
   use bfm_output, only: var_names, stPelStateS, stPelStateE
   use string_functions, only: getseq_number
   implicit none
   private

   integer, parameter, public :: max_atmdep_fields = 16

   logical, public :: use_atmdep = .false.
   integer, public :: n_atmdep_fields = 0

   character(len=256), public :: atmdep_file = 'atmdep_n.nc'
   character(len=64), public :: atmdep_lon_name = 'rlon'
   character(len=64), public :: atmdep_lat_name = 'rlat'
   character(len=64), public :: atmdep_flux_vars(max_atmdep_fields) = ''
   character(len=64), public :: atmdep_target_vars(max_atmdep_fields) = ''
   integer :: atmdep_target_ids(max_atmdep_fields) = 0

   REALTYPE, allocatable, public :: atmdep_flux(:,:,:)

   logical :: targets_resolved = .false.

   public :: do_atmdep, reset_atmdep_targets

contains

!-----------------------------------------------------------------------
! Clear cached BFM indices whenever the atmospheric-deposition setup changes.
!-----------------------------------------------------------------------
   subroutine reset_atmdep_targets()
      atmdep_target_ids = 0
      targets_resolved = .false.
   end subroutine reset_atmdep_targets

!-----------------------------------------------------------------------
! Translate readable namelist names into pelagic BFM state-array indices.
!-----------------------------------------------------------------------
   subroutine resolve_atmdep_targets()
      integer :: field
      integer :: registry_id

      if (targets_resolved) return

      if (.not. allocated(var_names)) then
         stop 'resolve_atmdep_targets: BFM variable registry is not initialised'
      end if
      if (stPelStateS < 1 .or. stPelStateE < stPelStateS .or. &
          stPelStateE > size(var_names)) then
         stop 'resolve_atmdep_targets: invalid BFM pelagic state range'
      end if

      do field=1,n_atmdep_fields
         registry_id = getseq_number(trim(atmdep_target_vars(field)), &
                                     var_names,stPelStateE,.true.,stPelStateS)
         if (registry_id == 0) then
            LEVEL3 'Unknown atmospheric deposition target: ', &
                   trim(atmdep_target_vars(field))
            stop 'resolve_atmdep_targets: target is not a BFM pelagic state'
         end if

         atmdep_target_ids(field) = registry_id - stPelStateS + 1
         LEVEL3 'Atmospheric deposition target: ', &
                trim(atmdep_target_vars(field)),' -> pelagic state ', &
                atmdep_target_ids(field)
      enddo

      targets_resolved = .true.
   end subroutine resolve_atmdep_targets

!-----------------------------------------------------------------------
! Add atmospheric-deposition fluxes after physical 3-D integration and
! before biological tracer transport in do_getm_bio.
!-----------------------------------------------------------------------
   subroutine do_atmdep(loop,do_3d)
      integer, intent(in) :: loop
      logical, intent(in) :: do_3d

      integer :: i,j,field,target_id
      REALTYPE :: atmdep_dt,surface_thickness

#ifdef DEBUG
      LEVEL4 'do_atmdep step ',loop
#endif
      if (.not. use_atmdep) return
      if (.not. do_3d) return
      if (.not. bio_calc) return

      call resolve_atmdep_targets()
      atmdep_dt = M*timestep

      do j=jmin,jmax
         do i=imin,imax
            if (az(i,j) == 1) then
               surface_thickness = hn(i,j,kmax)
               if (surface_thickness > _ZERO_) then
                  do field=1,n_atmdep_fields
                     target_id = atmdep_target_ids(field)
                     if (target_id < 1 .or. target_id > size(cc3d,4)) then
                        stop 'do_atmdep: target is outside the pelagic state array'
                     end if
                     cc3d(i,j,kmax,target_id) = cc3d(i,j,kmax,target_id) + &
                          atmdep_flux(i,j,field) * atmdep_dt / surface_thickness
                  enddo
               end if
            end if
         enddo
      enddo

   end subroutine do_atmdep

end module atmdep

#endif
