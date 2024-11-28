!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model 
!
! FUNCTION
!   DefineSet
!
! FILE
!   DeFineset.f90
!
! DESCRIPTION
!   !
!  
!   This file is generated from f77 code, using a code generator which
!   transposes from the F77 code into F90
!   the code USES module file as used in the BFM model
!   F90 code generator written by P. Ruardij. 
!
! AUTHORS
! Piet Ruardij (rua@nioz.nl)  
!
! CHANGE_LOG
!   
!
! COPYING
!   
!   Copyright (C) 2004 P. Ruardij, the mfstep group, the ERSEM team 
!   (rua@nioz.nl, vichi@bo.ingv.it)
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation;
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      SUBROUTINE DefineSet(NUTR,mode,option,input,xinput,yinput)
        USE global_mem, ONLY:RLEN,LOGUNIT
        USE bennut_variables
        USE constants,only: LAYERS,DIFFUSION,POROSITY,ADSORPTION, &
          DEFINE,DOUBLE_DEFINE,PARAMETER_DEFINE,LABDA_1,LABDA_2, &
          SET_BOUNDARY
        USE bennut_constants
        USE bennut_interface,ONLY:input_para, transfer

        IMPLICIT  NONE
        integer,intent(IN) ::NUTR ! Specification
        integer,intent(IN) ::mode ! Specification
        integer,intent(IN) ::input ! Specification
        integer,intent(IN) ::option ! Specification
        REAL(RLEN),intent(IN) ::xinput ! Specification
        REAL(RLEN),intent(IN) ::yinput ! Specification
        integer ::i
        integer ::j
        integer ::j_1
        logical llError

!JM        integer,static  ::nr=0
        integer,save :: nr=0

        REAL(RLEN) ::r
        REAL(RLEN) ::s

        llError=.false.
        if (NUTR.ne.nutr_seq) stop 'error: wrong use of DefneSet'

        select case (mode)
        case(LAYERS)
          nr=0
          if (ns%status == LAYERS) then
            call input_para(0,option,input &
            ,xinput,yinput,ns%b(2),ns%equa-1,i)
            if (i == 1) ns%status= DIFFUSION
          else
            write(LOGUNIT,*) 'error DefineSet case: LAYERS'
            llError=.true.
          endif
        case(DIFFUSION)
          if (ns%status == DIFFUSION) then
            if (ns%b(ns%equa+1).lt.0.0D+00) ns%b(ns%equa+1)=1.0D30
            call input_para(1,option,input &
            ,xinput,yinput,ns%diff,ns%equa,i)
            if (i == 1) ns%status= POROSITY
          else
            write(LOGUNIT,*) 'error Defineset case: DIFFUSION'
            llError=.true.
          endif
        case(POROSITY)
          if (ns%status == POROSITY) then
            call input_para(1,option,input &
            ,xinput,yinput,ns%poro,ns%equa,i)
            if (i == 1) ns%status= ADSORPTION
          else
            write(LOGUNIT,*) 'error DefineSet case:POROSITY'
            llError=.true.
          endif
        case(ADSORPTION)
          if (ns%status == ADSORPTION) then
            call input_para(1,option,input &
            ,xinput,yinput,ns%ads,ns%equa,i)
            if (i == 1) ns%status= DEFINE
          else
            write(LOGUNIT,*) 'error DefineSet case:ADSORPTION'
            llError=.true.
          endif
        case(DEFINE,DOUBLE_DEFINE,PARAMETER_DEFINE)
          if (ns%status /= DEFINE.and.ns%equa /= 1) then 
            write(LOGUNIT,*) &
               'DefineSet:mode=DEFINE/DOUBLE_DEFINE:not allowed now'
            write(LOGUNIT,*) 'option (3d input parameter)=',option
            write(LOGUNIT,*) 'input (4d input parameter)=',input
            write(LOGUNIT,*) 'xinput (5d input parameter)=',xinput
            write(LOGUNIT,*) 'yinput (6d input parameter)=',yinput
            if ( ns%status==7) write(LOGUNIT,*) 'ns%ads=',ns%ads
            llError=.true.
          else
            i=option/10
            j=option-10*i
            if (i > ns%equa.or.(j < 0.or.j > 9)) then
              write(LOGUNIT,*) &
               'mode=DEFINE/DOUBLE_DEFINE/DOUBLE_DEINFE:not in interval'
               llError=.true.
            else
              nr=nr+1
              ns%coeffs(nr)%ia=input
              ns%coeffs(nr)%il=option
              ns%coeffs(nr)%labda=0.0D+00
              if ( input < 0 ) ns%coeffs(nr)%labda(1)=xinput
              if (mode /= DEFINE) ns%coeffs(nr)%labda(2)=yinput
              if (mode == PARAMETER_DEFINE) then
                r = transfer( LABDA_1,ns%coeffs(nr), xinput, ns%diff(i))
                s = transfer( LABDA_2,ns%coeffs(nr), yinput, ns%diff(i))
                ns%coeffs(nr)%labda(1) = r;
                ns%coeffs(nr)%labda(2) = s;
              endif
            endif
          endif
        case (SET_BOUNDARY)
          if (nr == ns%nn.and.(ns%status == DEFINE.or.ns%equa == 1)) then
            !wrap away all terms not necessary
            !the set of DE will consists nr coefficients (=terms)
            !and so and we have to define at most nn equations to solve the
            ! system.
            !lst()=0 : new definition:
            if (ns%lst(1) == 0) then
              ns%nn=nr
              do i=1,nr
                if ( any(ns%coeffs(i+1:nr)%il < ns%coeffs(i)%il)) then
                  write(LOGUNIT,*) &
                    'DefineSet: terms notin increasing sequence'
                  llError=.true.
                  continue
                else 
                  j=ns%coeffs(i)%il/10
                  if (ns%lst(j) == 0) then
                    ns%lst(j)=i
                    j_1=j-1
                    if (j_1.gt.0) then
                      ns%lfi(j_1)=i-1
                      if (ns%lst(j_1) == 0) then
                        write(LOGUNIT,*) &
                         'DefineSet:terms notin increasing sequence'
                        llError=.true.
                        continue
                      endif
                    endif
                  endif
                endif
              enddo
              ns%lfi(ns%equa)=ns%nn
            endif
            nn_boundaries=0
          endif
        end select
        if ( llError) then
          write(LOGUNIT,*) 'NUTR=',NUTR
          write(LOGUNIT,*) 'mode=',mode
          write(LOGUNIT,*) 'ns%status=',ns%status
          write(LOGUNIT,*) 'ns%equa=',ns%equa
          stop 'DefineSet'
        endif
        return
      end

