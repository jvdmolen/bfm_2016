!$Id: getm_error.F90,v 1.3 2004-04-06 16:54:33 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: error_functions ---
!
! !INTERFACE:
   MODULE gotm_error_msg
!
! !DESCRIPTION:
!   Aim: to get error-message  to the error message system of getm
! !USES:
   IMPLICIT NONE

!EOP
!
   public           :: set_d3_model_flag_for_gotm, gotm_error, &
                       output_gotm_error, set_warning_for_getm, &
                       get_warning_for_getm

   logical           :: d3_model_flag=.FALSE.
   logical          :: first 
   logical          :: warning=.FALSE.
   character(len=80):: hold_msg, hold_sub

!
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: gotm_error() - global error reporting routine
!
! !INTERFACE:
   subroutine gotm_error(sub,msg)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   character(len=*),intent(IN)         :: sub,msg
!
!
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'd3_model_flag=',d3_model_flag
   if ( d3_model_flag) then
     if ( .NOT.first ) then
         hold_msg=msg
         hold_sub=sub
         first=.TRUE.
     endif
   else
     FATAL "gotm_error:Called from: ",trim(sub)
     FATAL "gotm_error:Message:     ",trim(msg)
     FATAL "d3_model_flag",d3_model_flag
     if ( d3_model_flag) then
        stop "gotm_error()"
     else   
         first=.true.    
     endif
   endif

   return
   end subroutine gotm_error
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: set_d3_model_flag_for_gotm()  
!
! !INTERFACE:
   subroutine set_d3_model_flag_for_gotm(val)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   logical,intent(IN)                :: val
!
!
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'set_d3_model_flag_for_gotm val=',val
   if ( val) then
     LEVEL1 'Error messages to finish in GOTM will be handeld in GETM '
   endif
   d3_model_flag=val

  return
  end subroutine set_d3_model_flag_for_gotm
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: get_d3_model_flag_from_getm  
!
! !INTERFACE:
   subroutine get_d3_model_flag_from_getm(val)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   logical,intent(OUT)                :: val
!
!
!-----------------------------------------------------------------------
!BOC
   val=d3_model_flag

  return
  end subroutine get_d3_model_flag_from_getm
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: output_gotm_error() - global error reporting routine
!
! !INTERFACE:
   subroutine output_gotm_error(flag_msg,sub,msg)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   logical,intent(OUT)                          ::flag_msg
   character(len=*),intent(OUT)                 :: sub,msg
!
!
!-----------------------------------------------------------------------
!BOC
   if ( first ) then
     sub=hold_sub
     msg=hold_msg
     flag_msg=.true.
   else 
     flag_msg=.false.
   endif
   first=.false.
   return
   end subroutine output_gotm_error
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: set_warning_for_getm()
!
! !INTERFACE:
   subroutine set_warning_for_getm()
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
!
!-----------------------------------------------------------------------
!BOC
   warning=.true.

  return
  end subroutine set_warning_for_getm
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: get_warning_for_getm(flag_warning)
! !INTERFACE:
   subroutine get_warning_for_getm(flag_warning)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
! !INPUT PARAMETERS:
   logical,intent(OUT)                          ::flag_warning
!-----------------------------------------------------------------------
!BOC
   flag_warning=warning

  return
  end subroutine get_warning_for_getm
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: reset_warning_for_getm
! !INTERFACE:
   subroutine reset_warning_for_getm
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
! !INPUT PARAMETERS:
!-----------------------------------------------------------------------
!BOC
  warning=.false.
  return
  end subroutine reset_warning_for_getm
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: check_for_fatal_error()
!
! !INTERFACE:
   logical function check_for_fatal_error()
!
! !DESCRIPTION:
!       check_for_fatal_error cechk is a ftatl error is made (first=.True.)
!
! !USES:
   IMPLICIT NONE
   check_for_fatal_error=first
!
!
!-----------------------------------------------------------------------
!BOC

  return
   end function check_for_fatal_error



end module
!-----------------------------------------------------------------------
! Copyright (C) 2006 - Hans Burchard and Karsten Bolding (BBH)         !
!-----------------------------------------------------------------------
