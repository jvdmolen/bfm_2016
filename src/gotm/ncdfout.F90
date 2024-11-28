!$Id: ncdfout.F90,v 1.19 2008-08-01 07:32:25 lars Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: ncdfout --- saving the results in NetCDF
!
! !INTERFACE:
   module ncdfout
   use netcdf

   IMPLICIT NONE

!
! !PUBLIC MEMBER FUNCTIONS:
   public define_mode, new_nc_variable, set_attributes, store_data
!
! !PUBLIC DATA MEMBERS:

!  netCDF file id
   integer, public          :: ncid=-1
   logical, public          :: first=.true.
   integer, public          :: set_no=0
!  dimension lengths
   integer, parameter        :: lon_len=1
   integer, parameter        :: lat_len=1
   integer                   :: depth_len
   integer                    :: time_len=NF90_UNLIMITED
   integer, public          :: start(4),edges(4)
!
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------


!
! !IROUTINE: Close files used for saving model results
!
! !INTERFACE:
   subroutine close_ncdf()
   IMPLICIT NONE
!
! !DESCRIPTION:
!  Closes the NetCDF file.
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  See ncdfout module
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: iret
!
!-------------------------------------------------------------------------
!BOC
   LEVEL1 'Output has been written in NetCDF'

   if (ncid .ne. -1) then
      iret = nf90_close(ncid)
      call check_err(iret)
   end if
   first=.true.
   set_no=0

   return
   end subroutine close_ncdf
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Begin or end define mode
!
! !INTERFACE:
   integer function define_mode(ncid,action)
!
! !DESCRIPTION:
!  Depending on the value of the argument {\tt action},
!  this routine put NetCDF in the `define' mode or not.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)       :: ncid
   logical, intent(in)       :: action
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  See ncdfout module
!
!EOP
!
! !LOCAL VARIABLES:
   integer         :: iret
!
!-----------------------------------------------------------------------
!BOC
   if(action) then
      iret = nf90_redef(ncid)
!kbk      call check_err(iret)
   else
      iret = nf90_enddef(ncid)
!kbk      call check_err(iret)
   end if
   define_mode = 0
   return
   end function define_mode
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Define a new NetCDF variable
!
! !INTERFACE:
   integer function new_nc_variable(ncid,name,data_type,n,dims,id)
!
! !DESCRIPTION:
!  This routine is used to define a new variable to store in a NetCDF file.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: ncid
   character(len=*), intent(in)        :: name
   integer, intent(in)                 :: data_type,n
   integer, intent(in)                 :: dims(1:n)
!
! !OUTPUT PARAMETERS:
   integer, intent(out)                :: id
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  See ncdfout module
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: iret
!
!-----------------------------------------------------------------------
!BOC
   iret = nf90_def_var(ncid,name,data_type,dims,id)
   call check_err(iret)
   new_nc_variable = iret
   return
   end function new_nc_variable
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set attributes for a NetCDF variable.
!
! !INTERFACE:
   integer function set_attributes(ncid,id,               &
                         units,long_name,                 &
                         valid_min,valid_max,valid_range, &
                         scale_factor,add_offset,         &
                         FillValue,missing_value,         &
                         C_format,FORTRAN_format)
!
! !DESCRIPTION:
!  This routine is used to set a number of attributes for
!  variables. The routine makes heavy use of the {\tt optional} keyword.
!  The list of recognized keywords is very easy to extend. We have
!  included a sub-set of the COARDS conventions.
!
! !USES:
!  IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)        :: ncid,id
   character(len=*), optional :: units,long_name
   REALTYPE, optional         :: valid_min,valid_max
   REALTYPE, optional         :: valid_range(2)
   REALTYPE, optional         :: scale_factor,add_offset
   REALTYPE, optional         :: FillValue,missing_value
   character(len=*), optional :: C_format,FORTRAN_format
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  See ncdfout module
!
! !LOCAL VARIABLES:
   integer                     :: len,iret
   REAL_4B                     :: vals_2(2),vals_1(1)
!
!EOP
!-----------------------------------------------------------------------
!BOC
   if(present(units)) then
      len = len_trim(units)
      iret = nf90_put_att(ncid,id,'units',units(1:len))
   end if

   if(present(long_name)) then
      len = len_trim(long_name)
      iret = nf90_put_att(ncid,id,'long_name',long_name(1:len))
   end if

   if(present(C_format)) then
      len = len_trim(C_format)
      iret = nf90_put_att(ncid,id,'C_format',C_format(1:len))
   end if

   if(present(FORTRAN_format)) then
      len = len_trim(FORTRAN_format)
      iret = nf90_put_att(ncid,id,'FORTRAN_format',FORTRAN_format(1:len))
   end if

   if(present(valid_min)) then
      vals_1(1) = valid_min
      iret = nf90_put_att(ncid,id,'valid_min',vals_1)
   end if

   if(present(valid_max)) then
      vals_1(1) = valid_max
      iret = nf90_put_att(ncid,id,'valid_max',vals_1)
   end if

   if(present(valid_range)) then
      vals_2(1) = valid_range(1)
      vals_2(2) = valid_range(2)
      iret = nf90_put_att(ncid,id,'valid_range',vals_2)
   end if

   if(present(scale_factor)) then
      vals_1(1) = scale_factor
      iret = nf90_put_att(ncid,id,'scale_factor',vals_1)
   end if

   if(present(add_offset)) then
      vals_1(1) = add_offset
      iret = nf90_put_att(ncid,id,'add_offset',vals_1)
   end if

   if(present(FillValue)) then
      vals_1(1) = FillValue
      iret = nf90_put_att(ncid,id,'_FillValue',vals_1)
   end if

   if(present(missing_value)) then
      vals_1(1) = missing_value
      iret = nf90_put_att(ncid,id,'missing_value',vals_1)
   end if

   set_attributes = 0
   return
   end function set_attributes
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Store values in a NetCDF file
!
! !INTERFACE:
   integer function store_data(ncid,id,var_shape,nlev, &
                               iscalar,iarray,scalar,array)
!
! !DESCRIPTION:
!  This routine is used to store a  variable in the NetCDF file.
!  The subroutine uses {\tt optional} parameters to find out which data
!  type to save.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: ncid,id,var_shape,nlev
   integer, optional                   :: iscalar
   integer, optional                   :: iarray(0:nlev)
   REALTYPE, optional                  :: scalar
   REALTYPE, optional                  :: array(1:nlev) !AN (1:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  See ncdfout module
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: iret,k,n=0
   integer                   :: idum(1:nlev)
   REAL_4B                   :: r4,dum(1:nlev)
!
!-----------------------------------------------------------------------
!BOC
   if (.not. present(iscalar) .and. .not. present(iarray) .and. &
       .not. present(scalar)  .and. .not. present(array) ) then
      FATAL 'At least one optional argument has to be passed to - store_data()'
      stop 'store_data'
   end if
   n = 0
   if(present(iscalar)) n = n+1
   if(present(iarray))  n = n+1
   if(present(scalar))  n = n+1
   if(present(array))   n = n+1
   if(n .ne. 1) then
      FATAL 'Only one optional argument must be passed to - store_data()'
      stop 'store_data'
   end if

   if (present(iscalar)) then
      select case (var_shape)
         case(POINT)
            iret = nf90_put_var(ncid,id,iscalar)
         case(T_SHAPE)
            start(1) = set_no; edges(1) = 1
            idum(1)=iscalar
            iret = nf90_put_var(ncid,id,start,edges,idum)
         case default
            FATAL 'A non valid - var_shape - has been passed in store_data()'
            stop 'store_data'
      end select
   else if (present(scalar)) then
      select case (var_shape)
         case(POINT)
            r4 = scalar
            iret = nf90_put_var(ncid,id,r4)
         case(T_SHAPE)
            start(1) = set_no; edges(1) = 1
            dum(1)=scalar
            iret = nf90_put_var(ncid,id,dum,start(1:1),edges(1:1))
         case(XYT_SHAPE)
            start(1) = 1;   edges(1) = lon_len
            start(2) = 1;   edges(2) = lat_len
            start(3) = set_no; edges(3) = 1
            dum(1)=scalar
            iret = nf90_put_var(ncid,id,dum,start(1:4),edges(1:4))
         case default
            FATAL 'A non valid - var_shape - has been passed in store_data()'
            stop 'store_data'
      end select
   else if (present(array)) then
      select case (var_shape)
         case(Z_SHAPE)
            k=1
            start(1) = 1;   edges(1) = nlev
         case(XYZT_SHAPE)
            k=4
            start(1) = 1;   edges(1) = lon_len
            start(2) = 1;   edges(2) = lat_len
            start(3) = 1;   edges(3) = nlev
            start(4) = set_no; edges(4) = 1
         case default
            FATAL 'A non valid - var_shape - has been passed in store_data()'
            stop 'store_data'
      end select
      dum=array!(1:nlev) !AN
      iret = nf90_put_var(ncid,id,dum,start(1:k),edges(1:k))
   else
   end if
   call check_err(iret)
   store_data = iret
   return
   end function store_data

   subroutine check_err(iret)
   integer,intent(IN)   :: iret
   if (iret .ne. NF90_NOERR) then
   print *, nf90_strerror(iret)
   stop
   endif
   end subroutine check_err
!EOC

!-----------------------------------------------------------------------

   end module ncdfout

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
