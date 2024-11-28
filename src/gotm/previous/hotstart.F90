!$Id: gotm.F90,v 1.42 2009-10-21 08:02:09 hb Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: hotstart --- the general framework \label{sec:hotstart}
!
! !INTERFACE:
   module hotstart
!
! !DESCRIPTION:
! This is 'where it all happens'. This module provides the internal
! routines {\tt hotstart_read()} to read and initialise hot model variables and
! {\tt hotstart_save()} to save hotstart variables
!
! !USES:
   use output,only:out_fn
   use meanflow
   use observations
   use time

   use airsea,      only: wind=>w,tx,ty,I_0,heat,precip,evap
#ifdef BFM_GOTM
   use airsea,      only: wx_wind=>u10,wy_wind=>v10,daylength             !BFM
#endif

   use airsea,      only: int_swr,int_heat,int_total
   use turbulence,  only: tke,eps,uu,vv,ww
   use turbulence,  only: num,nuh,nus
   use turbulence,  only: const_num,const_nuh
   use turbulence,  only: gamu,gamv,gamh,gams
   use turbulence,  only: kappa


#ifdef BIO
   use bio
!  use bio_fluxes
#ifdef BFM_GOTM
   use bfm_output,only: var_names
   use mem, only:D3STATE
   use mem, ONLY:D2STATE
#endif
#endif

   use bfm_output,only:stPelStates,stBenStates

   IMPLICIT NONE
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public hotstart_read, hotstart_save

   integer,public         ::sw_hotstart=-1,sw_hotstart_read=0,test_hotstart_time=1
   integer                ::yyyy_old=0
!
! !DEFINED PARAMETERS:
!
!
! !REVISION HISTORY:
!  Original author(s): Adrian Farcas, CEFAS, UK (2015-01)
!
!
!EOP
!
!!
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Delete hotstart variables \label{hotstart_read}
!
! !INTERFACE:
   subroutine hotstart_delete(mode)
!
! !DESCRIPTION:
!  This internal routine triggers the initialization of
!  model's hot variables

! !USES:
  IMPLICIT NONE
  integer,intent(IN)         ::mode
!
! !REVISION HISTORY:
!  Original author(s): Adrian Farcas, CEFAS, UK (2015-01)
!
!
!
!EOP
!
! !LOCAL VARIABLES:
   integer             ::     yyyy,mm,dd,il,j,k,yyyy_del
   logical             ::     file_exists
   Character(LEN=30)   ::     hotstartfile
!
!-----------------------------------------------------------------------
!BOC
  if (mode<=0) return
  call calendar_date(julianday,yyyy,mm,dd)
  yyyy_del=-yyyy

   il=0
   j=yyyy-mode
   STDERR "hotstart: test on delete",j,mode,il
   do while (il<10)
     k=dd+il
     write (hotstartfile,'(i4,''-'',i2.2,''-'',i2.2)') j,mm,k
     hotstartfile=trim(trim(out_fn)//"-"//trim(hotstartfile)//".hot")
     inquire( file=hotstartfile, exist=file_exists )
     STDERR hotstartfile,file_exists
     if (file_exists ) then
       open (unit=84,file=hotstartfile,status='unknown',form='unformatted')
       LEVEL2 'hotstartfile ',trim(hotstartfile) ,' deleted'
       close(unit=84,status='DELETE')
       yyyy_del=-yyyy_del
       return
     endif
     il=il+1
   end do
!  yyyy_del=-yyyy
!  if(yyyy_del>0) then
!      LEVEL2 'old hotstartfile ',trim(hotstartfile) ,' not deletedd'
!      LEVEL2 'according inquire statement no old hotstart file present!'
!  endif
   end subroutine hotstart_delete
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Read hotstart variables \label{hotstart_read}
!
! !INTERFACE:
   subroutine hotstart_read(mode)
!
! !DESCRIPTION:
!  This internal routine triggers the initialization of
!  model's hot variables

! !USES:
  IMPLICIT NONE
  integer,intent(IN)         ::mode
!
! !REVISION HISTORY:
!  Original author(s): Adrian Farcas, CEFAS, UK (2015-01)
!
!
!
!EOP
!
! !LOCAL VARIABLES:
   integer(8)          ::     hotjulianday
   integer             ::     hotsecondsofday,yyyy,mm,dd, &
                              minutes,hours
   logical             ::     file_exists
   Character(LEN=30)   ::     hotstartfile,date
!
!-----------------------------------------------------------------------
!BOC

   if (mode<0.or.mode>3) then
     STDERR 'sw_hotstart_read ouside range -1 to 10'      
     stop 'hotstart_read'
   endif        
   select case(mode)
   case(1)
     call calendar_date(julianday,yyyy,mm,dd)
     write (date,'(i4,''-'',i2.2,''-'',i2.2)') yyyy,mm,dd

      hotstartfile=trim(trim(out_fn)//"-"//trim(date)//".hot")
      inquire( file=hotstartfile, exist=file_exists )
      STDERR "Hotstart:"
      STDERR "  file_exists",file_exists,hotstartfile
      if (.not.file_exists ) then
        hotstartfile=trim("Keep-"//trim(date)//".hot")
        inquire( file=hotstartfile, exist=file_exists )
      endif
      STDERR "  file_exists",file_exists,hotstartfile
      STDERR "end hotstart:"
      yyyy_old=yyyy
   case(2)
     hotstartfile="Ecology.hot"
     inquire( file=hotstartfile, exist=file_exists )
      yyyy_old=yyyy
!  case(2)
!    hotstartfile="Pelagic.hot"
!    inquire( file=hotstartfile, exist=file_exists )
!  case(3)
!    hotstartfile="Benthic.hot"
!    inquire( file=hotstartfile, exist=file_exists )
   end select
   if (file_exists ) then
     open (unit=84,file=hotstartfile,status='unknown',form='unformatted')
     if (test_hotstart_time==1) then
       LEVEL1 'hotstart!',julianday
       LEVEL2 'reading_hot variables...'
       read (84) hotjulianday,hotsecondsofday
       LEVEL2  'hot read:',hotstartfile,hotjulianday,hotsecondsofday
       LEVEL2  'live:',julianday,secondsofday
       hours=secondsofday/3600
       minutes=(secondsofday-hours*3600)/60
       LEVEL2  'minutes;',minutes,'hours:',hours
       if ((hotjulianday.eq.julianday).and.&
         (hotsecondsofday.eq.secondsofday)) then
         read (84,err=10) D3STATE
         read (84,err=11) D2STATE
         read (84) t,s,u,v,h,tke,eps,num,nuh,nus,uu,vv,ww, &
                                    zeta,z0s,z0b,wind,wx_wind,wy_wind
         read (84) int_swr,int_heat,int_total
         close(84)
       else
         LEVEL1 'hotstart file variables mismatch:'
         LEVEL2  'julianday,secondsofday'
         LEVEL2  'hot read:',hotjulianday,hotsecondsofday
         hours=hotsecondsofday/3600
         minutes=(secondsofday-hours*3600)/60
         LEVEL2  'minutes;',minutes,'hours:',hours
         LEVEL2  'live:',julianday,secondsofday
         hours=secondsofday/3600
         minutes=(secondsofday-hours*3600)/60
         LEVEL2  'minutes;',minutes,'hours:',hours
         close(84)
         stop
       endif  
     elseif (test_hotstart_time==0) then  
       ! in own risk
       read (84) hotjulianday,hotsecondsofday
       read (84,err=10) D3STATE
       read (84,err=11) D2STATE
       close(84)
     endif
   endif 
     LEVEL2 'done.'
     return
 10  STDERR  &
       "error reading D3STATE:probably data generated from run without pelagic"
     stop
 11  STDERR  &
       "error reading D2STATE:probably data generated from run without benthic"
   end subroutine hotstart_read
!EOC



!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Write hotstart variables \label{hotstart_save}
!
! !INTERFACE:
   subroutine hotstart_save(mode)
!
! !DESCRIPTION:
!  This internal routine triggers the saving of
!  model's hot variables

! !USES:
  IMPLICIT NONE
  integer,intent(IN)     :: mode
!
! !REVISION HISTORY:
!  Original author(s): Adrian Farcas, CEFAS, UK (2015-01)
!
!
!
!EOP
!
! !LOCAL VARIABLES:
   integer            :: yyyy,mm,dd,hours,minutes
   character(LEN=30)  :: hotstartfile
!
!-----------------------------------------------------------------------
!BOC
   if ( mode<-1.or.mode>9) then
     STDERR 'sw_hotstart ouside range -1 to 10'      
     stop 'hotstart'
   endif
   call calendar_date(julianday,yyyy,mm,dd)
   if (mode>=0.and.yyyy==yyyy_old) return
   if (mode>=1) call hotstart_delete(mode)

   hours=secondsofday/3600
   minutes=(secondsofday-hours*3600)/60
   write (hotstartfile,'(i4,''-'',i2.2,''-'',i2.2)') yyyy,mm,dd
   LEVEL2 'hotstartfile:saving hot variables for ',hotstartfile
   LEVEL2 'mode=',mode,'minutes;',minutes,'hours:',hours
   hotstartfile=trim(trim(out_fn)//"-"//trim(hotstartfile)//".hot")
!  print *,'t var',t
!  print *,'zeta,z0s,z0b',zeta,z0s,z0b
   open (unit=84,file=hotstartfile,status='unknown',form='unformatted')
   write (84) julianday,secondsofday
   write (84) D3STATE
   write (84) D2STATE
   write (84) t,s,u,v,h,tke,eps,num,nuh,nus,uu,vv,ww, &
                                zeta,z0s,z0b,wind,wx_wind,wy_wind
   write (84) int_swr,int_heat,int_total
   close(84)
   yyyy_old=yyyy
   end subroutine hotstart_save
!EOC

!-----------------------------------------------------------------------



!EOC

!-----------------------------------------------------------------------

   end module hotstart

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
