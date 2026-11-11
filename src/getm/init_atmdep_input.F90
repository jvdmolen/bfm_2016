#include "cppdefs.h"

#ifdef BFM_GOTM

!-----------------------------------------------------------------------
! Initialise atmospheric-deposition configuration and NetCDF input.
!-----------------------------------------------------------------------
subroutine init_atmdep_input(input_dir,n)

   use domain, only: imin, imax, jmin, jmax
   use atmdep, only: use_atmdep, atmdep_file, max_atmdep_fields, &
                    n_atmdep_fields, atmdep_flux_vars, &
                    atmdep_target_vars, atmdep_flux, &
                    atmdep_lon_name, atmdep_lat_name, &
                    reset_atmdep_targets
   use ncdf_atmdep, only: init_atmdep_input_ncdf
   implicit none

   character(len=*), intent(in) :: input_dir
   integer, intent(in) :: n

   integer :: rc,field
   integer, parameter :: atmdep_nml_unit = 97
   character(len=256) :: nml_file
   character(len=512) :: atmdep_path

   namelist /atmdep_nml/ use_atmdep, atmdep_file, &
                         n_atmdep_fields, atmdep_flux_vars, &
                         atmdep_target_vars, &
                         atmdep_lon_name, atmdep_lat_name

   LEVEL1 'init_atmdep_input'

   nml_file = trim(input_dir) // 'atmdep.nml'
   open(atmdep_nml_unit,status='old',action='read',file=trim(nml_file), &
        iostat=rc)
   if (rc /= 0) stop 'init_atmdep_input: could not open atmdep.nml'
   read(atmdep_nml_unit,NML=atmdep_nml)
   close(atmdep_nml_unit)

   call reset_atmdep_targets()

   if (.not. use_atmdep) then
      LEVEL2 'Atmospheric deposition disabled'
      return
   end if

   if (n_atmdep_fields < 1 .or. n_atmdep_fields > max_atmdep_fields) then
      stop 'init_atmdep_input: invalid n_atmdep_fields'
   end if

   do field=1,n_atmdep_fields
      if (len_trim(atmdep_flux_vars(field)) == 0) then
         stop 'init_atmdep_input: empty atmospheric deposition flux variable'
      end if
      if (len_trim(atmdep_target_vars(field)) == 0) then
         stop 'init_atmdep_input: empty atmospheric deposition target variable'
      end if
   enddo

   allocate(atmdep_flux(imin:imax,jmin:jmax,n_atmdep_fields),stat=rc)
   if (rc /= 0) stop 'init_atmdep_input: Error allocating atmdep_flux'
   atmdep_flux = _ZERO_

   LEVEL2 'Atmospheric deposition enabled'
   LEVEL3 'atmdep_file=          ',trim(atmdep_file)
   LEVEL3 'n_atmdep_fields=      ',n_atmdep_fields

   if (len_trim(atmdep_file) == 0) then
      atmdep_path = ''
   else if (atmdep_file(1:1) == '/') then
      atmdep_path = trim(atmdep_file)
   else
      atmdep_path = trim(input_dir) // trim(atmdep_file)
   end if

   call init_atmdep_input_ncdf(atmdep_path,n_atmdep_fields, &
                               atmdep_flux_vars,atmdep_flux, &
                               atmdep_lon_name,atmdep_lat_name,n)

end subroutine init_atmdep_input

!-----------------------------------------------------------------------
! Release atmospheric-deposition input resources.
!-----------------------------------------------------------------------
subroutine clean_atmdep_input()

   use atmdep, only: use_atmdep, atmdep_flux, reset_atmdep_targets
   use ncdf_atmdep, only: clean_atmdep_input_ncdf
   implicit none

   call clean_atmdep_input_ncdf()
   if (allocated(atmdep_flux)) deallocate(atmdep_flux)
   call reset_atmdep_targets()
   use_atmdep = .false.

end subroutine clean_atmdep_input

#endif
