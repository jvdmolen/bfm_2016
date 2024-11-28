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

   use airsea_driver,      only: wind=>w,tx,ty,I_0,heat,precip,evap
#ifdef BFM_GOTM
   use airsea_driver,      only: wx_wind=>u10,wy_wind=>v10,daylength             !BFM
#endif

   use airsea_driver,      only: int_swr,int_heat,int_total
   use turbulence,  only: tke,eps,uu,vv,ww
   use turbulence,  only: num,nuh,nus
   use turbulence,  only: const_num,const_nuh
   use turbulence,  only: gamu,gamv,gamh,gams
   use turbulence,  only: kappa


#ifdef BFM_GOTM
   use bfm_output,only: var_names
   use mem, only:D3STATE
   use mem, ONLY:D2STATE
!JM   use mem_Param, ONLY:p_CalcPelagicFlag,p_CalcBenthicFlag
   use mem_Param, ONLY:CalcPelagicFlag,CalcBenthicFlag
#endif

#ifdef BIO
   use bio
!  use bio_fluxes
#endif

   use bfm_output,only:stPelStates,stBenStates

   IMPLICIT NONE
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public hotstart_read, hotstart_save

   integer,public         ::sw_hotstart=1
   integer                ::yyyy_old=0,step=0
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
  call calendar_date(julianday,yyyy,mm,dd)
  if (step==0)step=max(step,dd)
!   step=-dd;return
! elseif(step<0)then; step=dd+step;endif
! ! This hotstart_read aims to check of  hotstart files of previous year
! ! can be deleted. if mode >0 the hotstart-file yyyy-mode will be deleted.
  yyyy_del=-yyyy

   il=step-1
   j=yyyy-mode
   STDERR "hotstart: test on delete",yyyy,mode,il,step
   do while (il>=-(step-1))
     k=dd+il
     write (hotstartfile,'(i4,''-'',i2.2,''-'',i2.2)') j,mm,k
     hotstartfile=trim(trim(out_fn)//"-"//trim(hotstartfile)//".hot")
     inquire( file=hotstartfile, exist=file_exists )
     if (file_exists ) then
       open (unit=84,file=hotstartfile,status='unknown',form='unformatted')
       LEVEL2 'hotstartfile ',trim(hotstartfile) ,' deleted'
       close(unit=84,status='DELETE')
       yyyy_del=-yyyy_del
       return
     endif
     il=il-1
   end do
   yyyy_del=-yyyy
   if(yyyy_del>0) then
       LEVEL2 'old hotstartfile ',trim(hotstartfile) ,' not deletedd'
       LEVEL2 'according inquire statement no old hotstart file present!'
   endif
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
   integer          ::     hotjulianday
   integer             ::     hotsecondsofday,yyyy,mm,dd, &
                              minutes,hours
   logical             ::     file_exists
   Character(LEN=30)   ::     hotstartfile,date
!
!-----------------------------------------------------------------------
!BOC

   if (mode==0) return
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
   if (file_exists ) then
     open(unit=84,file=hotstartfile,status='unknown',form='unformatted')
     LEVEL1 'hotstart!',julianday
     LEVEL2 'reading_hot variables...'
     read (84) hotjulianday,hotsecondsofday
         LEVEL2  'hot read:',hotstartfile,hotjulianday,hotsecondsofday
         LEVEL2  'live:',julianday,secondsofday
         hours=secondsofday/3600
         minutes=(secondsofday-hours*3600)/60
         LEVEL2  'minutes;',minutes,'hours:',hours
!    if (hotjulianday.eq.julianday)then
     if ((hotjulianday.eq.julianday).and.&
         (hotsecondsofday.eq.secondsofday)) then
         read (84) t,s,u,v,h,tke,eps,num,nuh,nus,uu,vv,ww, &
                                    zeta%value,z0s,z0b,wind,wx_wind%value,wy_wind%value
         read (84) int_swr,int_heat,int_total
         if (mode==3) then
            LEVEL2 'sw_hotstart=3 - reading ecology.hot'
            open(unit=85,file='ecology.hot',status='unknown',form='unformatted')
!JM            if (p_CalcPelagicFlag)    read (85,err=10) D3STATE
!JM            if (p_CalcBenthicFlag==3) read (85,err=11) D2STATE
            if (CalcPelagicFlag)    read (85,err=10) D3STATE
            if (CalcBenthicFlag==3) read (85,err=11) D2STATE
            close(85)
         else
!JM            if (p_CalcPelagicFlag)    read (84,err=10) D3STATE
!JM            if (p_CalcBenthicFlag==3) read (84,err=11) D2STATE
            if (CalcPelagicFlag)    read (84,err=10) D3STATE
            if (CalcBenthicFlag==3) read (84,err=11) D2STATE
         endif
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
!    print *,'t var',t
!    print *,'updated zeta,z0s,z0b,wind,wx_wind,wy_wind',zeta,z0s,z0b,wind,wx_wind,wy_wind
     LEVEL2 'done.'
     return
 10  STDERR  &
       "error reading D3STATE:probably data generated from run without pelagic"
     stop
 11  STDERR  &
       "error reading D2STATE:probably data generated from run without benthic"

   elseif (mode<=0) then
     LEVEL1 'no hotstart file, nothing changed!'
   endif


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
   LEVEL2 'julianday,secondsofday ',julianday,secondsofday
   write (84) t,s,u,v,h,tke,eps,num,nuh,nus,uu,vv,ww, &
                                zeta%value,z0s,z0b,wind,wx_wind%value,wy_wind%value
   write (84) int_swr,int_heat,int_total
!JM   if(p_CalcPelagicFlag)    write (84) D3STATE
   if(CalcPelagicFlag)    write (84) D3STATE
   open (unit=85,file='ecology.hot',status='unknown',form='unformatted')
   LEVEL2 'Saving Pelagic and Benthic state also in ecology.hot'
!JM   if(p_CalcPelagicFlag)    write (85) D3STATE
!JM   if(p_CalcBenthicFlag==3) write (85) D2STATE
   if(CalcPelagicFlag)    write (85) D3STATE
   if(CalcBenthicFlag==3) write (85) D2STATE
   close(85)
!JM   if(p_CalcBenthicFlag==3) write (84) D2STATE
   if(CalcBenthicFlag==3) write (84) D2STATE
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
