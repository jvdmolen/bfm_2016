!$Id: modify_meteo_input.F90,v 1.2 2003-04-23 12:04:08 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
!
! !MODULE:  modify_meteo_input - input specifications
!
! !INTERFACE:
   module modify_meteo_input
!
! !BFM:
!    routine to modify ECWMF-meteo input
!    Wind (only the velocities) and temperature can be modified 
!
! !USES:
   use meteo, only: t2,u10,v10,new_meteo
   use domain, only: imin,imax,jmin,jmax,az
   IMPLICIT NONE
!
! !PUBLIC DATA ME
   REALTYPE           :: t2_addition=0.0
   integer            :: uv10_way=1
   REALTYPE           :: uv10_factor=1.0
   INTEGER            :: modify_meteo=-1

! !LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------

  contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_modify_meteo_input - do modify meteo_input
!
! !INTERFACE:
   subroutine do_modify_meteo_input(file)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*),intent(IN)                    :: file
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  22Nov Author name Initial code
!
!
!EOP
!-------------------------------------------------------------------------
!BOC

     integer          ::file_id=283 
     integer          ::i 
     integer          ::j 
     REALTYPE         ::w,d,wm 
     namelist/modify/t2_addition,uv10_way,uv10_factor

     if (modify_meteo < 0) then
       modify_meteo=0
       uv10_way=1
       open(UNIT=file_id,FILE=file,ACTION='read',ERR=100)
       read(file_id,nml=modify,err=200)
       close(UNIT=file_id)
       LEVEL3 "-----meteo forcing modified ---------------------------------"
       LEVEL3 't2 increased with ',t2_addition
       select case  (uv10_way) 
         case (1) ;LEVEL3 'wind multiplied with  with ',uv10_factor
         case (2) ;LEVEL3 'exponential increase wind such that Beaufort7-->Beaufort8'
         case (3) ;LEVEL3 'power 2 increase wind such that Beaufort7-->Beaufort8'
         case (4) ;LEVEL3 'exponential increase wind such that Beaufort8-->Beaufort9'
       end select
       LEVEL3 "-------------------------------------------------------------"
       modify_meteo=1
     endif
     if (modify_meteo==1 .and. new_meteo ) then 
         do j=jmin,jmax
            do i=imin,imax
              if (az(i,j) .ge. 1) then
                 t2(i,j)=t2(i,j)+t2_addition
                 select case  (uv10_way) 
                    case (1)
                      u10(i,j)=u10(i,j) * uv10_factor
                      v10(i,j)=v10(i,j) * uv10_factor
                    case (2,3,4)
                      w=sqrt(u10(i,j)**2+v10(i,j)**2)
                      select case(uv10_way)
                       case (2) ;wm=exp (w*0.19273D+00)-1.0D+00
                       case (3) ;wm=2.0**(w*0.27406D+00)-1.0D+00
                       case (4) ;wm=w+exp(w*0.0792)-1.0D+00
                      end select
                      if ( w> 15.0) STDERR 'hw',i,j,w,wm
                      d=atan(v10(i,j)/(1.D-10+u10(i,j)))-1.0D+00
                      u10(i,j)=sign(abs(wm*cos(d)),u10(i,j))
                      v10(i,j)=sign(abs(wm*sin(d)),v10(i,j))
                end select
              end if
           end do
        end do
     endif
     return
 100 return
 200 stop 'reading data to modify meteo input' !AN added space
     
     end subroutine do_modify_meteo_input
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_modify_meteo_input - do modify meteo_input
!
! !INTERFACE:
   subroutine do_modify_temp_boundary(T,kmax,nsbv)
!
! !DESCRIPTION:
!
! !USES:

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   INTEGER,intent(IN)                        :: nsbv
   INTEGER,intent(IN)                        :: kmax
   REALTYPE,intent(INOUT)                    :: T(0:kmax,1:nsbv)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  22Nov Author name Initial code
!
!
!EOP
!-------------------------------------------------------------------------
     integer                            ::i,j
     REALTYPE                           ::r,s
     if ( t2_addition>0.0 ) then
       do j=1,nsbv
         r=max(0.0,T(kmax,j)-T(1,j))
         if ( r.gt.0.0 ) then
            do i=2,kmax
              s=T(i,j)-T(1,j)
              if (s.gt.0) T(i,j)=T(i,j) +min(1.0,s/r)*t2_addition ;
            enddo 
         endif
       enddo
     endif
     end subroutine do_modify_temp_boundary
!EOC

!-----------------------------------------------------------------------!BOC



!-----------------------------------------------------------------------!BOC

   end module modify_meteo_input

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Hans Burchard and Karsten Bolding               !
!-----------------------------------------------------------------------
