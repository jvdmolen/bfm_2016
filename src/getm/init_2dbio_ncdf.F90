!$Id: $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE:   init_2dbio_ncdf
! 
! !INTERFACE:
    subroutine init_2dbio_ncdf( mode,filename,name,ncben,status, var)
!
!
! !BFM
!   Routine meant to read 2d-fields of parameters used in the benthic model
!   of BFM ( porosities, phosphate and silica adsorption)
!
! !USES:
    use netcdf
    use exceptions, only: getm_error
    use domain,     only: iextr,jextr
    use domain,     only: imin,imax,jmin,jmax,kmax
    use domain,     only: ioff,joff
    use ncdf_topo,  only: ncdf_read_2d
!
! !INPUT PARAMETERS: 
    implicit none

    integer ,        intent(IN)          :: mode
    character(len=*),intent(IN)          :: filename
    character(len=*),intent(IN)          :: name
!
! !OUTPUT PARAMETERS:
    integer,         intent(INOUT)       :: ncben
    integer,         intent(INOUT)       :: status
    REALTYPE,intent(INOUT)               :: var(E2DFIELD)
! !LOCAL  VARIABLES:
    integer                              :: error=0
    integer                              :: i,k,l
    integer                              :: il,ih,jl,jh,ilocl,jlocl,iloch,jloch
    integer                              :: dimlen
    integer                              :: name_id
    integer,dimension(4)                 :: dimidsT
    character(len=80)                    :: text

! !DESCRIPTION:  read 2d-field with much less checks as done for bathymetry.
!          
!
!   01-07-2006 Piet Ruardij
! 
!EOP
!-------------------------------------------------------------------------
!BOC

    
    status=0
    if ( mode == 1 .or.mode == 11 ) then
       if (ncben /=0 ) then
         status=nf90_close(ncben)
         if (status .ne. NF90_NOERR) then
          text= "Error closing "//trim(filename)//"."
          goto 100
         endif
          k=len_trim(text)
          LEVEL2 text(1:k)
       endif
       status = nf90_open(filename,NF90_NOWRITE,ncben)
       if (status .ne. NF90_NOERR) then
          text= "Error opening "//trim(filename)//"."
          goto 100
       else 
          text= "Opening "//trim(filename)//"."
          k=len_trim(text)
          LEVEL2 text(1:k) 
       endif
    endif

    if ( mode  >= 0 ) then

!   Look for name
       status = nf90_inq_varid(ncben,name,name_id)
       if (status .eq. NF90_NOERR) then
          l=len_trim(name); k=len_trim(filename)
          write(text,'(''Found '',A,'' in '',A,''.'')') name(1:l),filename(1:k)
          k=len_trim(text)
          LEVEL2 text(1:k)
       else
          return
       endif
   
!   Is name a matrix?
       status = nf90_inquire_variable(ncben,name_id,ndims=dimlen)
       if (status .ne. NF90_NOERR) then
          write(text,'(''Could not get '',A,'' of '',A,''.'')')  &
                                              '''dimlen''',trim(filename)
          goto 100
       endif

       if (dimlen.lt.2) then
          text="name must have 2 dimensions."
          goto 100
       endif

!   Is the size of name consistent?
       status = nf90_inquire_variable(ncben,name_id,dimids=dimidsT)
       if (status .ne. NF90_NOERR) then
           write(text,'(''Could not get '',A,'' of '',A,''.'')')  &
                                              'dimensions',trim(filename)
            goto 100
       endif


       do i=1,dimlen
          status = nf90_inquire_dimension(ncben,dimidsT(i),len=il)
          if (status .ne. NF90_NOERR) then
            write(text,'(''Could not get dimlength'',i2,'' of '',A,''.'')') &
                                                          i,trim(filename)
          endif
    
          select case (i) 
            case(1)
                if (il.ne.iextr) error=i
                 !   Get i-dimension for dynamic allocation
            case(2)
                 if (il.ne.jextr) error=i
            case default 
               if (il.ne.1) error=i
         end select 
         if ( error .ne.0 ) then
           write(text, &
               '(''Length of dimension'',i2,'' in'',a,''inconsistent.'')') &
               error,name
           goto 100
         endif
       enddo

!  GLOBAL index range for variable to be read
      il    = max(imin+ioff,1);   ih    = min(imax+ioff,iextr)
      jl    = max(jmin+joff,1);   jh    = min(jmax+joff,jextr)
   
!  LOCAL index range for variable to be read
!  (different from GLOBAL range only for parallel runs)
      ilocl = max(imin-ioff,1);   jlocl = max(jmin-joff,1)
      iloch = ih-il+ilocl;        jloch = jh-jl+jlocl;

!  Read bathymetry
      call ncdf_read_2d(ncben,name_id,var(ilocl:iloch,jlocl:jloch),il,ih,jl,jh)
!     STDERR 'name',ilocl,iloch,jlocl,jloch,var(ilocl:iloch,jlocl:jloch)
   endif
   
   if ( mode .gt. 10.or. mode.lt.0 ) then
        status=nf90_close(ncben)
        if (status .ne. NF90_NOERR)  then
           text= "Error closing "//trim(filename)//"."
           goto 100
        endif
   endif

   if (status .eq. NF90_NOERR) status=0

   return

    100 call getm_error("ncdf_2dbio_ncdf()", text)   

    end
!EOC
!-----------------------------------------------------------------------
! Copyright by the GOTM-team and BFM-team under the GNU Public License 
! www.gnu.org
!-----------------------------------------------------------------------
