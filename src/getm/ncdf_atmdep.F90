#include "cppdefs.h"

#ifdef BFM_GOTM

module ncdf_atmdep
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use netcdf
   use exceptions, only: getm_error, netcdf_error
   use domain, only: imin, imax, jmin, jmax
   use domain, only: ill, ihl, jll, jhl, ilg, ihg, jlg, jhg
   use domain, only: az, lonc, latc
   use grid_interpol, only: init_grid_interpol, do_grid_interpol
   use time, only: jul0, secs0, timestep, CalDat, JulDay, string_to_julsecs

   implicit none

   private
   public init_atmdep_input_ncdf, get_atmdep_data_ncdf
   public clean_atmdep_input_ncdf

   integer :: ncid = -1
   integer :: n_ncdf_fields = 0
   integer, allocatable :: atmdep_var_ids(:)

   integer :: nx_source = -1
   integer :: ny_source = -1
   integer :: grid_scan = 1
   integer :: read_start(3) = 1
   integer :: read_count(3) = 1
   integer :: n_time_records = 1
   integer :: time_origin_year = -1
   integer :: time_origin_month = 1
   integer :: time_origin_day = 1
   integer :: cached_records(2) = -1
   REALTYPE, allocatable :: cached_flux(:,:,:,:), work(:,:)

   ! File timestamps in seconds from the stored origin (not integer seconds).
   integer, parameter :: time_kind = selected_real_kind(15,307)
   real(time_kind), allocatable :: atmdep_times(:)
   real(time_kind), allocatable :: atmdep_time_bounds(:,:)
   integer :: time_origin_jul = 0, time_origin_secs = 0

   REALTYPE, allocatable :: atmdep_lon(:)
   REALTYPE, allocatable :: atmdep_lat(:)
   REALTYPE, allocatable :: beta(:,:)
   REALTYPE, allocatable :: ti(:,:)
   REALTYPE, allocatable :: ui(:,:)
   REALTYPE, allocatable :: interp_field(:,:)
   integer, allocatable :: gridmap(:,:,:)

   REALTYPE :: southpole(3) = (/0.0,-90.0,0.0/)

contains

   subroutine clean_atmdep_input_ncdf()
      integer :: err

      ! Resource ownership, not the enabled flag, determines cleanup.
      if (ncid >= 0) then
         err = nf90_close(ncid)
         call check_ncdf(err,'clean_atmdep_input_ncdf: close deposition file')
         ncid = -1
      endif

      if (allocated(atmdep_var_ids)) deallocate(atmdep_var_ids)
      if (allocated(atmdep_times)) deallocate(atmdep_times)
      if (allocated(atmdep_time_bounds)) deallocate(atmdep_time_bounds)
      if (allocated(cached_flux)) deallocate(cached_flux)
      if (allocated(work)) deallocate(work)
      if (allocated(atmdep_lon)) deallocate(atmdep_lon)
      if (allocated(atmdep_lat)) deallocate(atmdep_lat)
      if (allocated(beta)) deallocate(beta)
      if (allocated(ti)) deallocate(ti)
      if (allocated(ui)) deallocate(ui)
      if (allocated(interp_field)) deallocate(interp_field)
      if (allocated(gridmap)) deallocate(gridmap)

      n_ncdf_fields = 0
      cached_records = -1
      nx_source = -1
      ny_source = -1
      grid_scan = 1
      read_start = 1
      read_count = 1
      n_time_records = 1
      time_origin_year = -1
      time_origin_month = 1
      time_origin_day = 1
      time_origin_jul = 0
      time_origin_secs = 0
      southpole = (/0.0,-90.0,0.0/)
   end subroutine clean_atmdep_input_ncdf

   subroutine check_ncdf(err, where)
      integer, intent(in) :: err
      character(len=*), intent(in) :: where

      if (err .ne. NF90_NOERR) then
         call netcdf_error(err, trim(where), 'atmospheric deposition input')
      endif
   end subroutine check_ncdf

   subroutine init_atmdep_input_ncdf(atmdep_file, n_fields, atmdep_var_names, &
                                     atmdep_flux, lon_name, lat_name, nstart)
      character(len=*), intent(in) :: atmdep_file
      integer, intent(in) :: n_fields
      character(len=*), intent(in) :: atmdep_var_names(:)
      REALTYPE, intent(inout) :: atmdep_flux(:,:,:)
      character(len=*), intent(in) :: lon_name
      character(len=*), intent(in) :: lat_name
      integer, intent(in) :: nstart

      integer :: err
      integer :: field
      integer :: ndims
      integer :: dimids(NF90_MAX_VAR_DIMS), field_dimids(NF90_MAX_VAR_DIMS)

      LEVEL1 'init_atmdep_input_ncdf'
      LEVEL2 'Atmospheric deposition file: ',trim(atmdep_file)

      err = nf90_open(trim(atmdep_file), NF90_NOWRITE, ncid)
      call check_ncdf(err, 'init_atmdep_input_ncdf: open atmdep file')

      cached_records = -1
      southpole = (/0.0,-90.0,0.0/)

      n_ncdf_fields = n_fields
      if (allocated(atmdep_var_ids)) deallocate(atmdep_var_ids)
      allocate(atmdep_var_ids(n_ncdf_fields))
      atmdep_var_ids = -1

      do field = 1,n_ncdf_fields
         if (len_trim(atmdep_var_names(field)) == 0) then
            stop 'init_atmdep_input_ncdf: empty atmospheric deposition variable name'
         endif
         err = nf90_inq_varid(ncid, trim(atmdep_var_names(field)), &
                              atmdep_var_ids(field))
         call check_ncdf(err, 'init_atmdep_input_ncdf: find deposition variable')
      enddo

      err = nf90_inquire_variable(ncid, atmdep_var_ids(1),ndims=ndims,dimids=dimids)
      call check_ncdf(err, 'init_atmdep_input_ncdf: inquire deposition variable')
      if (ndims /= 3) then
         call getm_error('init_atmdep_input_ncdf', &
                         'deposition fields must have dimensions (x,y,time)')
      endif
      do field=2,n_ncdf_fields
         err = nf90_inquire_variable(ncid,atmdep_var_ids(field), &
                                     ndims=ndims,dimids=field_dimids)
         call check_ncdf(err,'init_atmdep_input_ncdf: compare deposition dimensions')
         if (ndims /= 3) then
            call getm_error('init_atmdep_input_ncdf', &
                            'all deposition fields must have three dimensions')
         endif
         if (any(field_dimids(:3) /= dimids(:3))) then
            call getm_error('init_atmdep_input_ncdf', &
                            'all deposition fields must share (x,y,time) dimensions')
         endif
      enddo
      err = nf90_inquire_dimension(ncid, dimids(1), len=nx_source)
      call check_ncdf(err, 'init_atmdep_input_ncdf: inquire x dimension')
      err = nf90_inquire_dimension(ncid, dimids(2), len=ny_source)
      call check_ncdf(err, 'init_atmdep_input_ncdf: inquire y dimension')
      call init_atmdep_time_axis(dimids(3))

      call init_source_grid_read(lon_name,lat_name,dimids(1),dimids(2))

      if (allocated(cached_flux)) deallocate(cached_flux)
      if (allocated(work)) deallocate(work)
      allocate(cached_flux(size(atmdep_flux,1),size(atmdep_flux,2),n_fields,2), &
               work(read_count(1),read_count(2)),stat=err)
      if (err /= 0) call getm_error('init_atmdep_input_ncdf','cannot allocate endpoint cache')

      ! Match the native SSTS reader: initialize at the preceding step.
      call get_atmdep_data_ncdf(nstart-1, atmdep_flux)

   end subroutine init_atmdep_input_ncdf

   subroutine init_atmdep_time_axis(time_dimid)
      integer, intent(in) :: time_dimid

      integer :: err
      integer :: time_var_id, bounds_var_id
      integer :: ndims, time_dims(NF90_MAX_VAR_DIMS), bounds_dims(NF90_MAX_VAR_DIMS)
      integer :: bounds_ndims, bounds_vertex_count, p, k, rc
      integer :: yy, mm, dd, hh, minute, ss, earliest_jul
      real(time_kind) :: seconds_per_unit, attribute_value, gap_tolerance
      character(len=NF90_MAX_NAME) :: time_dim_name, bounds_name
      character(len=256) :: time_units, bounds_units, origin_text
      character(len=32) :: calendar, bounds_calendar
      character(len=19) :: origin_date

      err = nf90_inquire_dimension(ncid, time_dimid, name=time_dim_name, &
                                   len=n_time_records)
      call check_ncdf(err, 'init_atmdep_input_ncdf: inquire time dimension')

      if (n_time_records < 1) then
         call getm_error('init_atmdep_input_ncdf', &
                         'atmospheric deposition time dimension is empty')
      endif

      err = nf90_inq_varid(ncid, trim(time_dim_name), time_var_id)
      if (err .ne. NF90_NOERR) then
         err = nf90_inq_varid(ncid, 'time', time_var_id)
      endif
      call check_ncdf(err, 'init_atmdep_input_ncdf: find time variable')

      err = nf90_inquire_variable(ncid,time_var_id,ndims=ndims,dimids=time_dims)
      call check_ncdf(err, 'init_atmdep_time_axis: time dimensions')
      if (ndims /= 1) then
         call getm_error('init_atmdep_time_axis','time must be one-dimensional')
      endif
      if (time_dims(1) /= time_dimid) then
         call getm_error('init_atmdep_time_axis','time must use the deposition time dimension')
      endif

      err = nf90_get_att(ncid, time_var_id, 'units', time_units)
      call check_ncdf(err, 'init_atmdep_input_ncdf: read time units')

      time_units = adjustl(lowercase(time_units))
      p = index(time_units,' since ')
      if (p <= 1) then
         call getm_error('init_atmdep_time_axis','expected time units: <unit> since <date>')
      endif
      select case (trim(time_units(:p-1)))
         case ('seconds','second','sec','s')
            seconds_per_unit = 1.0_time_kind
         case ('minutes','minute','min')
            seconds_per_unit = 60.0_time_kind
         case ('hours','hour','hr','h')
            seconds_per_unit = 3600.0_time_kind
         case ('days','day','d')
            seconds_per_unit = 86400.0_time_kind
         case default
            call getm_error('init_atmdep_time_axis','unsupported time unit: '//trim(time_units))
      end select

      calendar = 'standard'
      err = nf90_get_att(ncid,time_var_id,'calendar',calendar)
      if (err == NF90_ENOTATT) then
         calendar = 'standard'
      else
         call check_ncdf(err,'init_atmdep_time_axis: calendar')
      endif
      select case (trim(adjustl(lowercase(calendar))))
         case ('standard','gregorian','proleptic_gregorian')
         case default
            call getm_error('init_atmdep_time_axis','unsupported calendar: '//trim(calendar))
      end select
      calendar = adjustl(lowercase(calendar))

      ! These deposition records are interval means. The midpoint coordinate
      ! labels each record; the associated bounds define when it applies.
      bounds_name = ''
      err = nf90_get_att(ncid,time_var_id,'bounds',bounds_name)
      if (err == NF90_ENOTATT) then
         call getm_error('init_atmdep_time_axis', &
                         'time bounds are required for interval-mean deposition')
      endif
      call check_ncdf(err,'init_atmdep_time_axis: read time bounds name')
      err = nf90_inq_varid(ncid,trim(bounds_name),bounds_var_id)
      call check_ncdf(err,'init_atmdep_time_axis: find time bounds variable')
      err = nf90_inquire_variable(ncid,bounds_var_id,ndims=bounds_ndims, &
                                  dimids=bounds_dims)
      call check_ncdf(err,'init_atmdep_time_axis: time bounds dimensions')
      if (bounds_ndims /= 2) then
         call getm_error('init_atmdep_time_axis','time bounds must be two-dimensional')
      endif
      if (bounds_dims(2) /= time_dimid) then
         call getm_error('init_atmdep_time_axis','time bounds must use the deposition time dimension')
      endif
      err = nf90_inquire_dimension(ncid,bounds_dims(1),len=bounds_vertex_count)
      call check_ncdf(err,'init_atmdep_time_axis: time bounds vertex dimension')
      if (bounds_vertex_count /= 2) then
         call getm_error('init_atmdep_time_axis','time bounds must contain start and end')
      endif

      bounds_units = ''
      err = nf90_get_att(ncid,bounds_var_id,'units',bounds_units)
      if (err == NF90_NOERR) then
         if (trim(adjustl(lowercase(bounds_units))) /= trim(time_units)) &
            call getm_error('init_atmdep_time_axis','time and time-bounds units differ')
      else if (err /= NF90_ENOTATT) then
         call check_ncdf(err,'init_atmdep_time_axis: time bounds units')
      endif
      bounds_calendar = ''
      err = nf90_get_att(ncid,bounds_var_id,'calendar',bounds_calendar)
      if (err == NF90_NOERR) then
         if (trim(adjustl(lowercase(bounds_calendar))) /= trim(calendar)) &
            call getm_error('init_atmdep_time_axis','time and time-bounds calendars differ')
      else if (err /= NF90_ENOTATT) then
         call check_ncdf(err,'init_atmdep_time_axis: time bounds calendar')
      endif

      ! GETM accepts a whole-second date string. Normalize date-only input;
      ! reject unhandled fractional seconds/time zones rather than ignore them.
      origin_text = adjustl(time_units(p+7:))
      k = len_trim(origin_text)
      if (k > 0) then
         if (origin_text(k:k) == 'z') origin_text(k:k) = ' '
      endif
      if (len_trim(origin_text) == 10) origin_text = trim(origin_text)//' 00:00:00'
      if (len_trim(origin_text) /= 19) then
         call getm_error('init_atmdep_time_axis','origin must be YYYY-MM-DD [HH:MM:SS], optionally Z')
      endif
      if (origin_text(11:11) == 't') origin_text(11:11) = ' '
      if (origin_text(5:5) /= '-' .or. origin_text(8:8) /= '-' .or. &
          origin_text(11:11) /= ' ' .or. origin_text(14:14) /= ':' .or. &
          origin_text(17:17) /= ':') then
         call getm_error('init_atmdep_time_axis','invalid time origin format')
      endif
      origin_date = origin_text(:19)
      read(origin_date,'(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)',iostat=rc) yy,mm,dd,hh,minute,ss
      if (rc /= 0) call getm_error('init_atmdep_time_axis','invalid time origin numbers')
      if (yy < 1583 .or. mm < 1 .or. mm > 12 .or. dd < 1 .or. dd > 31 .or. &
          hh < 0 .or. hh > 23 .or. minute < 0 .or. minute > 59 .or. ss < 0 .or. ss > 59) then
         call getm_error('init_atmdep_time_axis','invalid origin or origin before supported year 1583')
      endif
      call string_to_julsecs(origin_date, time_origin_jul, time_origin_secs)
      call CalDat(time_origin_jul, time_origin_year, time_origin_month, &
                  time_origin_day)
      if (yy /= time_origin_year .or. mm /= time_origin_month .or. dd /= time_origin_day) then
         call getm_error('init_atmdep_time_axis','invalid calendar date in time origin')
      endif

      ! nf90_get_var does not automatically unpack scale_factor/add_offset.
      err = nf90_get_att(ncid,time_var_id,'scale_factor',attribute_value)
      if (err == NF90_NOERR) then
         if (attribute_value /= 1.0_time_kind) &
            call getm_error('init_atmdep_time_axis','packed time coordinates are unsupported')
      else if (err /= NF90_ENOTATT) then
         call check_ncdf(err,'init_atmdep_time_axis: scale_factor')
      endif
      err = nf90_get_att(ncid,time_var_id,'add_offset',attribute_value)
      if (err == NF90_NOERR) then
         if (attribute_value /= 0.0_time_kind) &
            call getm_error('init_atmdep_time_axis','packed time coordinates are unsupported')
      else if (err /= NF90_ENOTATT) then
         call check_ncdf(err,'init_atmdep_time_axis: add_offset')
      endif

      if (allocated(atmdep_times)) deallocate(atmdep_times)
      allocate(atmdep_times(n_time_records),stat=rc)
      if (rc /= 0) call getm_error('init_atmdep_time_axis','cannot allocate time coordinates')
      err = nf90_get_var(ncid,time_var_id,atmdep_times)
      call check_ncdf(err,'init_atmdep_time_axis: read time coordinates')
      if (any(.not. ieee_is_finite(atmdep_times))) &
         call getm_error('init_atmdep_time_axis','time coordinates must be finite')
      call reject_missing_times(time_var_id,'_FillValue')
      call reject_missing_times(time_var_id,'missing_value')
      if (any(abs(atmdep_times) >= real(NF90_FILL_DOUBLE,time_kind))) &
         call getm_error('init_atmdep_time_axis','invalid default fill value in time coordinates')

      atmdep_times = atmdep_times * seconds_per_unit
      if (any(.not. ieee_is_finite(atmdep_times))) &
         call getm_error('init_atmdep_time_axis','time conversion overflow')
      ! The supported calendars agree after the Gregorian transition.
      call JulDay(1583,1,1,earliest_jul)
      if (any(atmdep_times < -real(time_origin_jul-earliest_jul,time_kind)*86400.0_time_kind &
                            -real(time_origin_secs,time_kind))) &
         call getm_error('init_atmdep_time_axis','timestamps before 1583 are unsupported')
      do k=2,n_time_records
         if (atmdep_times(k) <= atmdep_times(k-1)) &
            call getm_error('init_atmdep_time_axis','time coordinates must be strictly increasing')
      enddo

      if (allocated(atmdep_time_bounds)) deallocate(atmdep_time_bounds)
      allocate(atmdep_time_bounds(2,n_time_records),stat=rc)
      if (rc /= 0) call getm_error('init_atmdep_time_axis','cannot allocate time bounds')
      err = nf90_get_var(ncid,bounds_var_id,atmdep_time_bounds)
      call check_ncdf(err,'init_atmdep_time_axis: read time bounds')
      if (any(.not. ieee_is_finite(atmdep_time_bounds))) &
         call getm_error('init_atmdep_time_axis','time bounds must be finite')
      call reject_missing_time_bounds(bounds_var_id,'_FillValue')
      call reject_missing_time_bounds(bounds_var_id,'missing_value')
      if (any(abs(atmdep_time_bounds) >= real(NF90_FILL_DOUBLE,time_kind))) &
         call getm_error('init_atmdep_time_axis','invalid default fill value in time bounds')
      atmdep_time_bounds = atmdep_time_bounds * seconds_per_unit
      if (any(.not. ieee_is_finite(atmdep_time_bounds))) &
         call getm_error('init_atmdep_time_axis','time-bounds conversion overflow')

      do k=1,n_time_records
         if (atmdep_time_bounds(1,k) >= atmdep_time_bounds(2,k)) &
            call getm_error('init_atmdep_time_axis','each time interval must have positive length')
         if (atmdep_times(k) < atmdep_time_bounds(1,k) .or. &
             atmdep_times(k) > atmdep_time_bounds(2,k)) &
            call getm_error('init_atmdep_time_axis','time coordinate lies outside its bounds')
         if (k > 1) then
            gap_tolerance = max(1.0e-6_time_kind, &
               64.0_time_kind*epsilon(1.0_time_kind)* &
               max(abs(atmdep_time_bounds(1,k)),abs(atmdep_time_bounds(2,k-1)), &
                   1.0_time_kind))
            if (abs(atmdep_time_bounds(1,k)-atmdep_time_bounds(2,k-1)) > gap_tolerance) &
               call getm_error('init_atmdep_time_axis','time intervals must be contiguous')
         endif
      enddo

      LEVEL3 'atmdep time records: ',n_time_records
      LEVEL3 'atmdep time origin year/month/day: ',time_origin_year, &
             time_origin_month,time_origin_day
      LEVEL3 'atmdep time seconds from origin: ',atmdep_times(1),atmdep_times(n_time_records)
      LEVEL3 'atmdep interval coverage: ',atmdep_time_bounds(1,1), &
             atmdep_time_bounds(2,n_time_records)

   end subroutine init_atmdep_time_axis

   function lowercase(text) result(lower)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: lower
      integer :: k, code
      lower = text
      do k=1,len(text)
         code = iachar(text(k:k))
         if (code >= iachar('A') .and. code <= iachar('Z')) lower(k:k) = achar(code+32)
      enddo
   end function lowercase

   subroutine reject_missing_times(varid,attribute)
      integer, intent(in) :: varid
      character(len=*), intent(in) :: attribute
      integer :: err, attribute_length
      real(time_kind) :: missing
      err = nf90_inquire_attribute(ncid,varid,attribute,len=attribute_length)
      if (err == NF90_ENOTATT) return
      call check_ncdf(err,'init_atmdep_time_axis: inquire '//attribute)
      if (attribute_length /= 1) &
         call getm_error('init_atmdep_time_axis','expected scalar '//attribute)
      err = nf90_get_att(ncid,varid,attribute,missing)
      if (err == NF90_ENOTATT) return
      call check_ncdf(err,'init_atmdep_time_axis: '//attribute)
      if (any(atmdep_times == missing)) &
         call getm_error('init_atmdep_time_axis','missing time coordinate: '//attribute)
   end subroutine reject_missing_times

   subroutine reject_missing_time_bounds(varid,attribute)
      integer, intent(in) :: varid
      character(len=*), intent(in) :: attribute
      integer :: err, attribute_length
      real(time_kind) :: missing
      err = nf90_inquire_attribute(ncid,varid,attribute,len=attribute_length)
      if (err == NF90_ENOTATT) return
      call check_ncdf(err,'init_atmdep_time_axis: inquire bounds '//attribute)
      if (attribute_length /= 1) &
         call getm_error('init_atmdep_time_axis','expected scalar bounds '//attribute)
      err = nf90_get_att(ncid,varid,attribute,missing)
      call check_ncdf(err,'init_atmdep_time_axis: bounds '//attribute)
      if (any(atmdep_time_bounds == missing)) &
         call getm_error('init_atmdep_time_axis','missing time bound: '//attribute)
   end subroutine reject_missing_time_bounds

   subroutine init_source_grid_read(lon_name,lat_name,x_dimid,y_dimid)
      character(len=*), intent(in) :: lon_name
      character(len=*), intent(in) :: lat_name
      integer, intent(in) :: x_dimid, y_dimid

      integer :: err
      integer :: id
      integer :: i, j, ndims
      integer :: axis_dimids(NF90_MAX_VAR_DIMS)
      logical :: ok
      REALTYPE :: tmp

      LEVEL2 'Initialising atmospheric deposition grid interpolation'

      if (nx_source < 2 .or. ny_source < 2) then
         call getm_error('init_source_grid_read', &
                         'source-grid longitude and latitude axes need at least two points')
      endif

      allocate(atmdep_lon(nx_source),stat=err)
      if (err /= 0) stop 'init_atmdep_input_ncdf: Error allocating atmdep_lon'
      err = nf90_inq_varid(ncid,trim(lon_name),id)
      call check_ncdf(err, 'init_atmdep_input_ncdf: find longitude variable')
      err = nf90_inquire_variable(ncid,id,ndims=ndims,dimids=axis_dimids)
      call check_ncdf(err,'init_source_grid_read: longitude dimensions')
      if (ndims /= 1 .or. axis_dimids(1) /= x_dimid) then
         call getm_error('init_source_grid_read', &
                         'longitude must be the one-dimensional deposition x axis')
      endif
      err = nf90_get_var(ncid,id,atmdep_lon)
      call check_ncdf(err, 'init_atmdep_input_ncdf: read longitude variable')

      allocate(atmdep_lat(ny_source),stat=err)
      if (err /= 0) stop 'init_atmdep_input_ncdf: Error allocating atmdep_lat'
      err = nf90_inq_varid(ncid,trim(lat_name),id)
      call check_ncdf(err, 'init_atmdep_input_ncdf: find latitude variable')
      err = nf90_inquire_variable(ncid,id,ndims=ndims,dimids=axis_dimids)
      call check_ncdf(err,'init_source_grid_read: latitude dimensions')
      if (ndims /= 1 .or. axis_dimids(1) /= y_dimid) then
         call getm_error('init_source_grid_read', &
                         'latitude must be the one-dimensional deposition y axis')
      endif
      err = nf90_get_var(ncid,id,atmdep_lat)
      call check_ncdf(err, 'init_atmdep_input_ncdf: read latitude variable')

      if (any(.not. ieee_is_finite(atmdep_lon)) .or. &
          any(.not. ieee_is_finite(atmdep_lat))) then
         call getm_error('init_source_grid_read','source-grid axes must be finite')
      endif
      do i=2,nx_source
         if (atmdep_lon(i) <= atmdep_lon(i-1)) then
            call getm_error('init_source_grid_read', &
                            'longitude axis must be strictly increasing')
         endif
      enddo
      if (.not. (all(atmdep_lat(2:) > atmdep_lat(:ny_source-1)) .or. &
                 all(atmdep_lat(2:) < atmdep_lat(:ny_source-1)))) then
         call getm_error('init_source_grid_read', &
                         'latitude axis must be strictly monotonic')
      endif

      call init_source_grid_mapping()

      if (atmdep_lat(1) .gt. atmdep_lat(2)) then
         LEVEL3 'Reverting atmospheric deposition lat-axis and setting grid_scan to 0'
         grid_scan = 0
         do j=1,ny_source/2
            tmp = atmdep_lat(j)
            atmdep_lat(j) = atmdep_lat(ny_source-j+1)
            atmdep_lat(ny_source-j+1) = tmp
         enddo
      endif

      allocate(beta(E2DFIELD),stat=err)
      if (err /= 0) stop 'init_atmdep_input_ncdf: Error allocating beta'
      beta = _ZERO_
      allocate(ti(E2DFIELD),stat=err)
      if (err /= 0) stop 'init_atmdep_input_ncdf: Error allocating ti'
      ti = -999.
      allocate(ui(E2DFIELD),stat=err)
      if (err /= 0) stop 'init_atmdep_input_ncdf: Error allocating ui'
      ui = -999.
      allocate(gridmap(E2DFIELD,1:2),stat=err)
      if (err /= 0) stop 'init_atmdep_input_ncdf: Error allocating gridmap'
      gridmap = -999
      allocate(interp_field(E2DFIELD),stat=err)
      if (err /= 0) stop 'init_atmdep_input_ncdf: Error allocating interp_field'
      interp_field = _ZERO_

      call init_grid_interpol(ill,ihl,jll,jhl,az,lonc,latc, &
                              atmdep_lon,atmdep_lat,southpole, &
                              gridmap,beta,ti,ui)

      ok = .true.
      do j=jmin,jmax
         do i=imin,imax
            if (az(i,j) .gt. 0 .and. &
                (ui(i,j) .lt. _ZERO_ .or. ti(i,j) .lt. _ZERO_)) then
               ok = .false.
               LEVEL3 'atmdep interpolation error at (i,j) ',i,j
            endif
         enddo
      enddo
      if (.not. ok) then
         call getm_error('init_atmdep_input_ncdf', &
                         'Some atmospheric deposition interpolation coefficients are not valid')
      endif

      read_start(1) = minval(gridmap(:,:,1),mask=(gridmap(:,:,1).gt.0))
      read_start(2) = minval(gridmap(:,:,2),mask=(gridmap(:,:,2).gt.0))
      read_count(1) = min(maxval(gridmap(:,:,1))+1,nx_source) - read_start(1) + 1
      read_count(2) = min(maxval(gridmap(:,:,2))+1,ny_source) - read_start(2) + 1
      read_start(3) = 1
      read_count(3) = 1

      where (gridmap(:,:,1) .gt. 0) gridmap(:,:,1) = gridmap(:,:,1) - read_start(1) + 1
      where (gridmap(:,:,2) .gt. 0) gridmap(:,:,2) = gridmap(:,:,2) - read_start(2) + 1

      if (grid_scan .eq. 0) then
         j = read_start(2)
         read_start(2) = ny_source - (read_start(2) + read_count(2) - 1) + 1
         read_count(2) = ny_source - j + 1 - read_start(2) + 1
      endif

      LEVEL3 'atmdep source x start/count: ',read_start(1),read_count(1)
      LEVEL3 'atmdep source y start/count: ',read_start(2),read_count(2)

   end subroutine init_source_grid_read

   subroutine init_source_grid_mapping()
      integer :: err, id, field
      character(len=NF90_MAX_NAME) :: mapping_name, field_mapping, mapping_type

      ! Follow the established ncdf_meteo/ncdf_ssts convention. If the
      ! deposition variable has no grid_mapping attribute, its one-dimensional
      ! axes are interpreted as ordinary geographic longitude and latitude.
      mapping_name = ''
      err = nf90_get_att(ncid,atmdep_var_ids(1),'grid_mapping',mapping_name)
      if (err == NF90_ENOTATT) then
         LEVEL3 'Atmospheric deposition grid: geographic lon/lat'
         return
      endif
      call check_ncdf(err,'init_source_grid_mapping: read grid_mapping attribute')

      do field=2,n_ncdf_fields
         field_mapping = ''
         err = nf90_get_att(ncid,atmdep_var_ids(field),'grid_mapping',field_mapping)
         if (err == NF90_ENOTATT .or. trim(field_mapping) /= trim(mapping_name)) then
            call getm_error('init_source_grid_mapping', &
                            'all deposition fields must use the same grid_mapping')
         endif
         call check_ncdf(err,'init_source_grid_mapping: compare grid_mapping attributes')
      enddo

      err = nf90_inq_varid(ncid,trim(mapping_name),id)
      call check_ncdf(err,'init_source_grid_mapping: find grid mapping variable')
      mapping_type = ''
      err = nf90_get_att(ncid,id,'grid_mapping_name',mapping_type)
      call check_ncdf(err,'init_source_grid_mapping: read grid_mapping_name')
      if (trim(lowercase(mapping_type)) /= 'rotated_latitude_longitude') then
         call getm_error('init_source_grid_mapping', &
                         'unsupported grid mapping: '//trim(mapping_type))
      endif

      ! CF stores the rotated north pole. GETM grid_interpol expects the
      ! equivalent south-pole description used by the native forcing readers.
      err = nf90_get_att(ncid,id,'grid_north_pole_latitude',southpole(1))
      call check_ncdf(err,'init_source_grid_mapping: grid_north_pole_latitude')
      err = nf90_get_att(ncid,id,'grid_north_pole_longitude',southpole(2))
      call check_ncdf(err,'init_source_grid_mapping: grid_north_pole_longitude')
      err = nf90_get_att(ncid,id,'north_pole_grid_longitude',southpole(3))
      if (err == NF90_ENOTATT) then
         southpole(3) = _ZERO_
      else
         call check_ncdf(err,'init_source_grid_mapping: north_pole_grid_longitude')
      endif

      if (any(.not. ieee_is_finite(southpole)) .or. &
          southpole(1) < -90.0 .or. southpole(1) > 90.0) then
         call getm_error('init_source_grid_mapping','invalid rotated-pole coordinates')
      endif

      if (southpole(2) >= _ZERO_) then
         southpole(2) = southpole(2) - 180.0
      else
         southpole(2) = southpole(2) + 180.0
      endif
      southpole(1) = -southpole(1)
      ! Retain the native reader convention; grid rotation about the pole is
      ! not used by this GETM interpolation implementation.
      southpole(3) = _ZERO_

      LEVEL3 'Atmospheric deposition grid: rotated latitude/longitude'
      LEVEL3 'atmdep south-pole lon/lat: ',southpole(2),southpole(1)
   end subroutine init_source_grid_mapping

   subroutine get_atmdep_data_ncdf(n, atmdep_flux)
      integer, intent(in) :: n
      REALTYPE, intent(inout) :: atmdep_flux(:,:,:)

      integer :: active_record, next_record, active_slot, next_slot, slot

      if (ncid < 0 .or. n_ncdf_fields < 1) then
         stop 'get_atmdep_data_ncdf: NetCDF reader not initialized'
      endif

      if (size(atmdep_flux,3) /= n_ncdf_fields) then
         stop 'get_atmdep_data_ncdf: field count changed after init'
      endif

      active_record = atmdep_record_for_time(n)
      next_record = min(active_record+1,n_time_records)
      active_slot = 0
      next_slot = 0
      do slot=1,2
         if (cached_records(slot) == active_record) active_slot = slot
         if (cached_records(slot) == next_record) next_slot = slot
      enddo

      ! Keep the active interval and prefetch the next interval. At a boundary,
      ! the formerly prefetched record becomes active without another read.
      if (active_slot == 0) then
         active_slot = 1
         if (next_slot == 1) active_slot = 2
         call read_atmdep_record(active_record,cached_flux(:,:,:,active_slot))
         cached_records(active_slot) = active_record
      endif
      if (next_record == active_record) then
         next_slot = active_slot
      else if (next_slot == 0) then
         next_slot = 3-active_slot
         call read_atmdep_record(next_record,cached_flux(:,:,:,next_slot))
         cached_records(next_slot) = next_record
      endif

      ! Each field is a mean over its time bounds. Hold that mean constant
      ! throughout the interval so its integrated deposition is preserved.
      atmdep_flux = cached_flux(:,:,:,active_slot)

   end subroutine get_atmdep_data_ncdf

   subroutine read_atmdep_record(record,atmdep_flux)
      integer, intent(in) :: record
      REALTYPE, intent(out) :: atmdep_flux(:,:,:)
      integer :: err, field

      read_start(3) = record
      atmdep_flux = _ZERO_

      do field = 1,n_ncdf_fields
         err = nf90_get_var(ncid, atmdep_var_ids(field), work, &
                            start=read_start, count=read_count)
         call check_ncdf(err, 'get_atmdep_data_ncdf: read deposition field')

         call flip_var(work)
         call do_grid_interpol(az,work,gridmap,ti,ui,interp_field)
         atmdep_flux(:,:,field) = interp_field(imin:imax,jmin:jmax)
      enddo

      LEVEL3 'atmdep time record: ',record
   end subroutine read_atmdep_record

   real(time_kind) function atmdep_seconds(n) result(seconds)
      integer, intent(in) :: n
      seconds = (real(jul0,time_kind)-real(time_origin_jul,time_kind))*86400.0_time_kind &
              + real(secs0,time_kind)-real(time_origin_secs,time_kind) &
              + real(n,time_kind)*real(timestep,time_kind)
   end function atmdep_seconds

   integer function atmdep_record_for_time(n)
      integer, intent(in) :: n
      integer :: left, right, middle
      real(time_kind) :: model_time
      character(len=256) :: msg

      ! Input runs before update_time(n). Use the native forcing convention
      ! n*timestep, relative to jul0/secs0, rather than the lagging calendar.
      ! Convert to real BEFORE multiplying: century-long seconds exceed int32.
      model_time = atmdep_seconds(n)
      if (.not. ieee_is_finite(model_time)) &
         call getm_error('atmdep_record_for_time','model time must be finite')
      if (model_time < atmdep_time_bounds(1,1) .or. &
          model_time >= atmdep_time_bounds(2,n_time_records)) then
         write(msg,'(a,es23.15,a,es23.15,a,es23.15)') &
            'time outside forcing coverage; seconds from file origin: model=',model_time, &
            ', first=',atmdep_time_bounds(1,1), &
            ', last=',atmdep_time_bounds(2,n_time_records)
         call getm_error('atmdep_record_for_time',trim(msg))
      endif

      if (n_time_records <= 1) then
         atmdep_record_for_time = 1
         return
      endif

      ! Find the latest interval start <= model_time. Intervals are [start,end),
      ! so an exact shared boundary selects the newly starting record.
      left = 1
      right = n_time_records
      do while (left < right)
         middle = left + (right-left+1)/2
         if (atmdep_time_bounds(1,middle) <= model_time) then
            left = middle
         else
            right = middle-1
         endif
      enddo
      atmdep_record_for_time = left

   end function atmdep_record_for_time

   subroutine flip_var(var)
      REALTYPE, intent(inout) :: var(read_count(1),read_count(2))

      select case (grid_scan)
         case (0)
            var = var(:,read_count(2):1:-1)
      end select

   end subroutine flip_var

end module ncdf_atmdep

#endif
