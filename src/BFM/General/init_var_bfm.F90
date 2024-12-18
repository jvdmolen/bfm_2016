#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Initialise BFM variables
!
! !INTERFACE:
   subroutine init_var_bfm(namlst,fname,unit,setup)
!
! !DESCRIPTION:
!  Allocate BFM variables and give initial values of
!  parameters and state variables
!
! !USES:
#ifndef NOT_STANDALONE
   use api_bfm
   use global_mem
#endif
#ifdef BFM_GOTM
   use global_mem,only:LOGUNIT,DONE
   use bio_var,only: var_ids,var_ave,var_names,nlev,numc,numbc, &
            calc_init_bennut_states,sw_BenSilt,NOTRANSPORT,ppb,ddb, &
            sfl_N3n,sfl_N4n,rel_max_sedi_rate,n_surface_fluxes, &
            surface_flux_method,stPelStateS,stPelStateE, &
            stBenStateS,stBenStateE,stPRFFluxE 
   use bio_bfm
#endif
#ifdef BFM_POM
   use api_pom
#endif
#ifdef BFM_OPA_OFFLINE
   use api_opa_offline
#endif
#ifdef BFM_OPA_PELAGOS
   use api_opa_pelagos
#endif
   use mem
   use mem_Phyto, ONLY: p_qnRc,p_qpRc,p_qsRc,p_qchlc,p_qnR2c
   use mem_BenBac, ONLY: p_qnHc=>p_qnc,p_qpHc=>p_qpc
   use mem_PelBac, ONLY: p_qnBc=>p_qnc,p_qpBc=>p_qpc
   use mem_MesoZoo, ONLY: p_qnMec=>p_qnc,p_qpMec=>p_qpc
   use mem_MicroZoo, ONLY: p_qnMic=>p_qnc,p_qpMic=>p_qpc
   use mem_BenPhyto, ONLY: p_useparams,CalculateBenPhyto
   use mem_PelBac, ONLY: p_version_PelBac=> p_version
   use mem_Silt,ONLY: p_SampleDepth
   use constants, ONLY: HOURS_PER_DAY,INTEGRAL, ZERO,PARAMETER,p_qnUc
   use mem_Param, ONLY: CalcPelagicFlag,CalcBenthicFlag,p_small, &
                        CalcPhytoPlankton,CalcMicroZooPlankton,          &
                        CalcMesoZooPlankton,CalcBenOrganisms,            &
                        CalcBenBacteria,CalcBacteria,CalcPelChemistry,   &
                        CalcBenPhyto,p_small,p_poro,p_d_tot,p_d_tot_2,p_1d_poro
   use string_functions, ONLY: getseq_number,empty,set_nco_style
   use mem_globalfun,   ONLY: IntegralExpDist_vector
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,          intent(in)        :: namlst
   character(len=*), intent(in)        :: fname
   integer,          intent(in)        :: unit
   integer,          intent(in)        :: setup

!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!
! !LOCAL VARIABLES:
   integer              :: icontrol,i,j,k,iiLastElement,nrphyto
   integer,parameter    :: NSAVE=200  ! Maximum no variables which can be saved
   REALTYPE             :: rhelp
   logical              :: nco_style_vars=.false.
   character(len=64),dimension(NSAVE):: var_save
   character(len=64),dimension(NSAVE):: ave_save
   REALTYPE  :: N1p0,N3n0,N4n0,N5s0,N6r0,  &
                P1c0,P2c0,P3c0,P4c0,P5c0,&
                Z2c0,Z3c0,Z4c0,Z5c0,Z6c0,B1c0,R9x0,R1c0,  &
                R2c0,R6c0,RZc0,O2o0,O3c0,O3h0, &
                P1l0,P2l0,P3l0,P4l0,P5l0

   REALTYPE, dimension(:), pointer  ::lcl_state
   REALTYPE, dimension(NO_BOXES_XY) ::alpha,r_xy1,r_xy2

   REALTYPE  :: Y1c0, Y2c0, Y3c0, Y4c0, Y5c0, &
                Q1c0, Q11c0, Q21c0, Q2c0, Q12c0, Q6c0, Yy3c0, &
                K3n0, G4n0, K15s0, BP1c0,H1c0,  &
                H2c0, H3c0, HNc0,Hac0, K1p0, K11p0, K21p0,     &
                K4n0, K14n0, K24n0, K6r0,K16r0,k26r0,k5s0, &
                D1m0, D2m0, D6m0, D7m0, D8m0, D9m0, G2o0,  &
                G3h0,G13h0,G23h0,G3c0,G13c0,G23c0, &
                p_qpQIc,p_qnQIc,p_qsQIc,p_poro0
   namelist /bfm_init_nml/ surface_flux_method,       &
                           n_surface_fluxes,          &
                           rel_max_sedi_rate,         &
                           sfl_N3n,sfl_N4n,sw_BenSilt,           &
                           N1p0,N3n0,N4n0,N5s0,N6r0,  &
                           P1c0,P2c0,P3c0,P4c0,P5c0,&
                           Z2c0,Z3c0,Z4c0,Z5c0,Z6c0,B1c0,R9x0,R1c0,  &
                           R2c0,R6c0,RZc0,O2o0,  &
                           P1l0,P2l0,P3l0,P4l0,P5l0

   namelist /bfm_save_nml/ nco_style_vars,var_save, ave_save

   namelist /bfm_ben_init_nml/calc_init_bennut_states,p_poro0, &
                           p_qpQIc,p_qnQIc,p_qsQIc, &
                           Y1c0, Y2c0, Y3c0,  Y4c0, Y5c0, BP1c0, &
                           Q1c0, Q11c0, Q21c0, Q2c0, Q12c0, Q6c0,          &
                           K3n0, G4n0, H1c0,          &
                           H2c0, H3c0, HNc0,Hac0, K1p0, K11p0, K21p0,   &
                           K4n0, k14n0, K24n0, K6r0,K16r0,K26r0,K5s0,    &
                           K15s0,&
                           D1m0, D2m0, D6m0, D7m0, D8m0, D9m0, G2o0, &
                           G3h0,G13h0,G23h0,G3c0,G13c0,G23c0
   interface
      subroutine init_cnps(c,n,p,s,nc,pc,sc)
         REALTYPE,dimension(:),intent(in)           :: c
         REALTYPE,intent(in),optional               :: nc,pc,sc
         REALTYPE,dimension(:),intent(out),optional :: n
         REALTYPE,dimension(:),intent(out),optional :: p
         REALTYPE,dimension(:),intent(out),optional :: s
      end subroutine init_cnps
   end interface
! COPYING
!
!   Copyright (C) 2006 P. Ruardij and Marcello Vichi
!   (rua@nioz.nl, vichi@bo.ingv.it)
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation;
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!EOP
!-----------------------------------------------------------------------
!BOC

   LEVEL2 'init_var_bfm'
   !---------------------------------------------
   ! Give reasonable initial values
   ! Overwritten by namelist parameters
   !---------------------------------------------
   surface_flux_method = -1
   n_surface_fluxes = 1

   !---------------------------------------------
   ! Pelagic variables
   !---------------------------------------------
   N1p0 = _ONE_
   N3n0 = _ONE_
   N4n0 = _ONE_
   N5s0 = _ONE_
   N6r0 = _ONE_
   O2o0 = 300.0
   O3c0 =  2101.0*12.0;
   O3h0 =  2275.0
   P1c0 = _ONE_
   P2c0 = _ONE_
   P3c0 = _ONE_
   P4c0 = _ONE_
   P5c0 = _ONE_
   P1l0 = _ONE_
   P2l0 = _ONE_
   P3l0 = _ONE_
   P4l0 = _ONE_
   P5l0 = _ONE_
   Z2c0 = _ONE_
   Z3c0 = _ONE_
   Z4c0 = _ONE_
   Z5c0 = _ONE_
   Z6c0 = _ONE_
   B1c0 = _ONE_
   R9x0 = _ONE_
   R1c0 = _ONE_
   R2c0 = _ONE_
   R6c0 = _ONE_
   RZc0 = 1.0D-80

   !---------------------------------------------
   ! Benthic variables
   !---------------------------------------------
   Y1c0  = _ONE_
   Y2c0  = _ONE_
   Y3c0  = _ONE_
   Yy3c0  = _ZERO_
   BP1c0  = _ONE_
   Y4c0  = _ONE_
   Y5c0  = _ONE_
   Q1c0  = _ONE_
   Q11c0 = _ONE_
   Q21c0 = _ONE_
   Q2c0  = _ONE_
   Q12c0 = _ONE_
   Q6c0  = _ONE_
   H1c0  = _ONE_
   H2c0  = _ONE_
   H3c0  = _ONE_
   HNc0  = _ONE_
   Hac0  = _ONE_
   K1p0  = _ONE_
   K11p0 = _ONE_
   K21p0 = _ONE_
   K3n0  = _ONE_
   G4n0  = _ONE_
   K4n0  = _ONE_
   K14n0 = _ONE_
   K24n0 = _ONE_
   K6r0  = _ONE_
   K5s0  = _ONE_
   K15s0  = _ONE_
   D1m0  = _ONE_
   D2m0  = _ONE_
   D6m0  = _ONE_
   D7m0  = _ONE_
   D8m0  = _ONE_
   D9m0  = _ONE_
   G2o0  = _ONE_
   G3c0  = _ZERO_
   G13c0  = _ZERO_
   G23c0  = _ZERO_
   G3h0  = _ZERO_
   G13h0  = _ZERO_
   G23h0  = _ZERO_

   !---------------------------------------------
   ! Open and read the namelist
   !---------------------------------------------
   icontrol=0
   open(namlst,file=fname,action='read',status='old',err=98)
   LEVEL2 "Read initial default initial values (3d model)"
   LEVEL2 "Read initial initial values (1d model)"
   read(namlst,nml=bfm_init_nml,err=99)
   write(LOGUNIT,nml=bfm_init_nml)
   p_qpQIc=-1.0;p_qnQIc=-1.0;p_qsQIc=-1.0
!  if (setup >=2 )  then
     read(namlst,nml=bfm_ben_init_nml,err=101)
     write(LOGUNIT,nml=bfm_ben_init_nml)
!  end if
   var_save=""
   ave_save=""
   var_ave=.false.
   LEVEL2 "Read variables names which will be saved in netcdf-file"
   read(namlst,nml=bfm_save_nml,err=100)
   close(namlst)
   icontrol=1
98 if ( icontrol == 0 ) then
     LEVEL3 'I could not open ',trim(fname)
     LEVEL3 'The initial values of the BFM variables are set to ONE'
     LEVEL3 'If thats not what you want you have to supply ',trim(fname)
   end if

   !---------------------------------------------
   ! Check variable to be saved and
   ! set the corresponding flag value in var_ids
   !---------------------------------------------
   if (nco_style_vars)  then
     do i=1,stPRFFluxE-1
      if (.not.( (i.ge.stPelStateS.and.i.le.stPelStateE).or. &
            (i.ge.stBenStateS.and.i.le.stBenStateE ) ) ) then
        var_names(i)=set_nco_style(var_names(i))
        continue
      endif
     enddo
   endif
   LEVEL2 "Variables selected for output"
   do i=1,NSAVE
      if (.NOT.empty(var_save(i))) then
            j=getseq_number(var_save(i),var_names,stPRFFluxE,.TRUE.)
            if ( j > 0 ) then  
              var_ids(j)=-1
              LEVEL3 trim(var_save(i))
            endif
      end if
      if ( .NOT.empty(var_save(i)) .AND. j==0 ) then
            LEVEL3 'Warning: variable ',trim(var_save(i)),' does not exist!'
      end if
   end do
   if (.NOT.empty(var_save(NSAVE))) then
      LEVEL3 "Warning:',NSAVE,'or more variables defined" 
      LEVEL3 "Warning:Change parameter NSAVE in init_var_bfm &
                   & or limit number of var_save variables"
      LEVEL3 "Warning:Potential memory problem!!"
   endif
   do i=1,NSAVE
      if (.NOT.empty(ave_save(i))) then
         j=getseq_number(ave_save(i),var_names,stPRFFluxE,.TRUE.)
         if ( .NOT.empty(ave_save(i)) .AND. j==0 ) then
            STDERR 'Warning: variable ',trim(ave_save(i)),' does not exist!'
         else if ( var_ids(j) <0 ) then
            STDERR 'Warning: Variable ',trim(ave_save(i)), &
               ' is already selected for output in var_save'
         else if ( j > 0 ) then
            var_ids(j)=-1
            var_ave(j)=.true.
            LEVEL3 trim(ave_save(i))
         end if
      end if
   end do
   if (.NOT.empty(ave_save(NSAVE))) then
      LEVEL3 "Warning:',NSAVE,'or more variables defined" 
      LEVEL3 "Warning:Change parameter NSAVE in init_var_bfm &
                            & or limit number of var_save variables"
      LEVEL3 "Warning:Potential memory problem!!"
   endif

   !---------------------------------------------
   ! Initialize BFM parameters
   !---------------------------------------------
   call Initialize

   !---------------------------------------------
   ! Initially set the number of sun hours
   ! equal to the number of hours in a day.
   !---------------------------------------------
   SUNQ = HOURS_PER_DAY

   !---------------------------------------------
   ! Initialise pelagic state variables
   ! also if using a benthic-only setup
   ! (for boundary conditions)
   !---------------------------------------------
      N1p = N1p0
      N3n = N3n0
      N4n = N4n0
      N5s = N5s0
      N6r = N6r0
      O2o = O2o0
      if ( ppO3c > 0 ) then
        D3STATE(:,ppO3c)=O3c0
        D3STATE(:,ppO3h)=O3h0 -N3n0
      endif
      P1c = P1c0
      if (p1l0 /= _ONE_) then
         P1l = P1l0
      else
         P1l = p1c0*p_qchlc(iiP1)
      end if
      P2c = P2c0
      if (p2l0 /= _ONE_) then
         P2l = P2l0
      else
         P2l = p2c0*p_qchlc(iiP2)
      end if
      P3c = P3c0
      if (p3l0 /= _ONE_) then
         P3l = P3l0
      else
         P3l = p3c*p_qchlc(iiP3)
      end if
      P4c = P4c0
      if (p4l0 /= _ONE_) then
         P4l = P4l0
      else
         P4l = p4c0*p_qchlc(iiP4)
      end if
      P5c = P5c0
      if (p5l0 /= _ONE_) then
         P5l = P5l0
      else
         P5l = p5c0*p_qchlc(iiP5)
      end if
      Pcc=1.0D-80
      P6c=1.0D-80
      P6l = P6c*p_qchlc(iiP6)
      Z2c = Z2c0
      Z3c = Z3c0
      Z4c = Z4c0
      Z2c = 1.0D-10
      Z5c = Z5c0
      Z6c = Z6c0
      B1c = B1c0
      R9x = R9x0
      R1c = R1c0
      R2c = R2c0 ; R2n=R2c*p_qnR2c
      R3c = 1.0D-10
      R6c = R6c0
      RZc = RZc0

      !---------------------------------------------
      ! Initialise other internal components
      ! with Redfield
      !---------------------------------------------
      do i=1,iiPhytoPLankton
        j= ppPhytoPlankton(i,iiN) 
        if ( j>0 ) then
          lcl_state=> PhytoPLankton(i,iiN)
          call init_cnps(c=PhytoPLankton(i,iiC),&
                           n=lcl_state,nc=p_qnRc(i))
        endif
        j= ppPhytoPlankton(i,iiP) 
        if ( j>0 ) then
           lcl_state=> PhytoPLankton(i,iiP)
           call init_cnps(c=PhytoPLankton(i,iiC), &
                             p=lcl_state,pc=p_qpRc(i))
        endif
        j= ppPhytoPlankton(i,iiS) 
        if ( j>0 ) then
           lcl_state=> PhytoPLankton(i,iiS)
           call init_cnps(c=PhytoPLankton(i,iiC), &
                        s=lcl_state,sc=p_qsRc(i))
        endif
      enddo
      do i=1,iiMesoZooPLankton
        j= ppMesoZooPlankton(i,iiN) 
        if ( j>  0 ) then
          lcl_state=> MesoZooPLankton(i,iiN)
          call init_cnps(c=MesoZooPLankton(i,iiC),&
                           n=lcl_state,nc=p_qnMec(i))
        endif
        j= ppMesoZooPlankton(i,iiP) 
        if ( j>  0 ) then
           lcl_state=> MesoZooPLankton(i,iiP)
           call init_cnps(c=MesoZooPLankton(i,iiC), &
                             p=lcl_state,pc=p_qpMec(i))
        endif
      enddo
      do i=1,iiMicroZooPLankton
        j= ppMicroZooPlankton(i,iiN) 
        if ( j>  0 ) then
           lcl_state=> MicroZooPlankton(i,iiN) 
           call init_cnps(c=MicroZooPLankton(i,iiC),&
                              n=lcl_state,nc=p_qnMic(i))
        endif
        j= ppMicroZooPlankton(i,iiP) 
        if ( j>  0 ) then
          lcl_state=> MicroZooPlankton(i,iiP) 
          call init_cnps(c=MicroZooPLankton(i,iiC),&
                               p=lcl_state,pc=p_qpMic(i))
        endif
      enddo
      call init_cnps(c=B1c,n=B1n,p=B1p,nc=p_qnBc,pc=p_qpBc)
      call init_cnps(c=R1c,n=R1n,p=R1p)
      call init_cnps(c=R6c,n=R6n,p=R6p,s=R6s)
   !---------------------------------------------
   ! Initialise benthic state variables
   !---------------------------------------------
   !MAV: need to always give init non-zero values
   ! because there are still part of the
   ! benthic system which are computed when setup=1
!   if (setup >=2) then
      p_poro=p_1d_poro
      write(LOGUNIT,*) 'p-poro',p_poro
      Y1c  = Y1c0
      Y2c  = Y2c0
      Y3c  = Y3c0
      Yy3c  = Yy3c0
      Ys3c  = ZERO
      Y4c  = Y4c0
      Y5c  = Y5c0
      BP1c  = BP1c0
      Q1c  = Q1c0
      Q11c = Q11c0
      Q21c = Q21c0
      Q2c  = Q2c0;Q2n=Q2c*p_qnR2c
      Q12c = Q12c0
      Q6c  = Q6c0
      H1c  = H1c0
      H2c  = H2c0
      H3c  = H3c0
      HNc  = HNc0
      Hac  = Hac0
      K1p  = K1p0
      Kp1p  = K1p
      K11p = K11p0
      K21p = K21p0
      K3n  = K3n0
      K13n  = K3n0
      K23n  = 0.1*K3n0
      Kp3n  = K3n
      G4n  = G4n0
      K4n  = K4n0
      Kp4n  = K4n
      K14n = K14n0
      K24n = K24n0
      K6r  = K6r0
      K16r  = K16r0
      K26r  = K26r0
      K5s  = K5s0
      K15s  = K15s0
      Dfm=0.1e-3
      Dcm=0.1e-3
!     irrenh_l=DONE
      write(LOGUNIT,*)'initial value DfmDcm=',Dfm,Dcm
      rhelp= CalculateBenPhyto(iiC,INTEGRAL,1,ZERO,Dfm(1),inverse=0.5D+00)
      Dlm=rhelp
      Kp4n=K4n*min(1.0,Dlm/D1m0)
      Kp3n=K3n
      write(LOGUNIT,*)'Dlm=',Dlm
      D1m  = D1m0
      D2m  = D2m0
      D6m  = D6m0
      D7m  = D7m0
      D8m  = D8m0
      D9m  = D9m0
      G2o  = G2o0
      DH2m=D6m
      DH3m=D6m
      DSm=-_ONE_

      if ( sw_BenSilt==1 .or.sw_BenSilt==2) then
        DSM=p_d_tot
        call CouplingBioSilt(-1,QSx)
        alpha=_ONE_/DSM
        r_xy1=p_SampleDepth
        r_xy2=p_d_tot
        QSx=r_xy1*QSx /IntegralExpDist_vector(-alpha, r_xy1)* &
              IntegralExpDist_vector( -alpha,r_xy2)

!       QSx=p_d_tot*QSx
      endif
      !---------------------------------------------
      ! Initialise organisms' internal components
      ! with Redfield
      !---------------------------------------------
      do i=1,iiBenPhyto
        nrphyto=p_useparams(i)
        j= ppBenPhyto(i,iiN) 
        if ( j>0 ) then
          lcl_state=> BenPhyto(i,iiN)
          call init_cnps(c=BenPhyto(i,iiC),&
                           n=lcl_state,nc=p_qnRc(nrphyto))
        endif
        write(LOGUNIT,*) 'bphyto ',BP1n(1),BP1c(1)
        j= ppBenPhyto(i,iiP) 
        if ( j>0 ) then
           lcl_state=> BenPhyto(i,iiP)
           call init_cnps(c=BenPhyto(i,iiC), &
                             p=lcl_state,pc=p_qpRc(nrphyto))
        endif
!       write(LOGUNIT,*) 'init qpc=',BP1p(1)/Bp1c(1)
        j= ppBenPhyto(i,iiS) 
        if ( j>0 ) then
           lcl_state=> BenPhyto(i,iiS)
           call init_cnps(c=BenPhyto(i,iiC), &
                        s=lcl_state,sc=p_qsRc(nrphyto))
        endif
!       write(LOGUNIT,*) 'init qsc=',BP1s(1)/Bp1c(1)
      enddo
      nrphyto=p_useparams(iiBP1)
      BP1l = BP1c0*p_qchlc(nrphyto)
      call init_cnps(c=Y1c,n=Y1n,p=Y1p)
      call init_cnps(c=Y2c,n=Y2n,p=Y2p)
      call init_cnps(c=Y3c,n=Y3n,p=Y3p)
      call init_cnps(c=Yy3c,n=Yy3n,p=Yy3p)
      call init_cnps(c=Y4c,n=Y4n,p=Y4p)
      call init_cnps(c=Y5c,n=Y5n,p=Y5p)
      call init_cnps(c=H1c,n=H1n,p=H1p,nc=p_qnHc, &
           pc=p_qpHc)
      call init_cnps(c=H2c,n=H2n,p=H2p,nc=p_qnHc, &
           pc=p_qpHc)
      ! H2 used for par.values H2 and H3 use same parameters!
      call init_cnps(c=H3c,n=H3n,p=H3p,nc=p_qnHc, &
           pc=p_qpHc)
      call init_cnps(c=HNc,n=HNn,p=HNp,nc=p_qnHc, &
           pc=p_qpHc)
!    end if

      !---------------------------------------------
      ! Initialise detritus' components  with Redfield
      !---------------------------------------------
      call init_cnps(c=Q1c,n=Q1n,p=Q1p,nc=p_qnQIc,pc=p_qpQIc)
      call init_cnps(c=Q11c,n=Q11n,p=Q11p,nc=p_qnQIc,pc=p_qpQIc)
      call init_cnps(c=Q21c,n=Q21n,p=Q21p,nc=p_qnQIc,pc=p_qpQIc)
      call init_cnps(c=Q6c,n=Q6n,p=Q6p,s=Q6s,nc=p_qnQIc,pc=p_qpQIc,sc=p_qsQIc)
      Qun=Q1n*0.05;Q1un=Q11n*0.05;Q2un=Q21n*0.05;Qpun=Qun
      Q1c=Q1c+Qun/p_qnUc;Q11c=Q11c+Q1un/p_qnUc;Q21c=Q21c+Q2un/p_qnUc;
#ifdef INCLUDE_MACROPHYT
       do j = 1,iiMacroStructure
        k=ppMacroStructure(j,iiC);D2STATE(:,k) = ZERO
       end do
       do j = 1,iiMacroContent
        k=ppMacroContent(j,iiC);D2STATE(:,k) = ZERO
        k=ppMacroContent(j,iiN);D2STATE(:,k) = ZERO
        k=ppMacroContent(j,iiP);D2STATE(:,k) = ZERO
        k=ppMacroContent(j,iiL);D2STATE(:,k) = ZERO
       end do
#endif

   !---------------------------------------------
   ! Check setup settings
   ! and finalize initialization
   !---------------------------------------------
   select case (setup)
      case (0)
      case (1) ! Pelagic only
         LEVEL2 "Pelagic-only setup (bio_setup=1), Switching off the benthic system"
         CalcBenthicFlag = 0
      case (2) ! Benthic only
         LEVEL2 "Benthic-only setup (bio_setup=2), Switching off the pelagic system"
         CalcPelagicFlag = .FALSE.
         CalcPhytoPlankton=.FALSE.
         CalcBacteria=.FALSE.
         CalcMesoZooPlankton=.FALSE.
         CalcMicroZooPlankton=.FALSE.
      case (3) ! Pelagic-Benthic coupling
         LEVEL2 "Pelagic-Benthic coupled setup (bio_setup=3)"
         if (CalcBenthicFlag == 0) &
            LEVEL3 'Warning, benthic system is switched off!'
         if (.NOT.CalcPelagicFlag) &
            LEVEL3 'Warning, pelagic system is switched off!'
   end select

   select case (CalcBenthicFlag)
     case (0)
        LEVEL3 "Benthic model is: not used"
     case (1)
        LEVEL3 "Benthic model is: simple nutrient return"
     case (2)
        LEVEL3 "Benthic model is: benthos + intermediate nutrient return"
     case (3)
        LEVEL3 "Benthic model is: benthos + Ruardij & Van Raaphorst"
! if element of H < 0.0: a call form getm!
        if ( ppG3c > 0 ) then 
           D2STATE(:,ppG3c)=G3c0
           D2STATE(:,ppG13c)=G13c0
           D2STATE(:,ppG23c)=G23c0
           D2STATE(:,ppG3h)=G3h0
           D2STATE(:,ppG13h)=G13h0
           D2STATE(:,ppG23h)=G23h0
           if ( G3c0== _ZERO_ ) D2STATE(:,ppG3c)=O3c0* p_poro(1)* D1m0;
           if ( G13c0== _ZERO_ )D2STATE(:,ppG13c)=O3c0* p_poro(1)*(D2m0-D1m0)
           if ( G23c0== _ZERO_ )D2STATE(:,ppG23c)=O3c0* p_poro(1)*(p_d_tot_2-D2m0)
           if ( G3h0== _ZERO_ )D2STATE(:,ppG3h)= O3h0* p_poro(1)* D1m0 
           if ( G13h0== _ZERO_ )D2STATE(:,ppG13h)=O3h0* p_poro(1)*(D2m0-D1m0)
           if ( G23h0== _ZERO_ )D2STATE(:,ppG23h)=O3h0* p_poro(1)*(p_d_tot_2-D2m0)
        endif
        STDERR 'Depth(1) calc_init_bennut_states=',Depth(1),calc_init_bennut_states
        if ( Depth(1) > 0.0 .and. calc_init_bennut_states >= 1) then
           LEVEL4 "Benthic nutrient State vars are calculated with an assumption of an"
           LEVEL4 "equilibrium between input (nutrient regeneratation) and output"
           LEVEL4 "(flux to water column and definitive loss processes)"
!          call InitBenthicNutrientDynamics(calc_init_bennut_states)
           call reset_diagonal(numbc,ddb)
           call reset_diagonal(numbc,ppb)
        endif
   end select

   !---------------------------------------------
   ! Initializing  silt
   !---------------------------------------------
LEVEL3 "init silt"
   InitializeModel=1
   call SiltDynamics
   call CalcCO2SatInField(nlev,numc,ERHO,ETW,ESW,D3STATE)

   InitializeModel=0
   !---------------------------------------------
   ! Zeroing of the switched off state variables
   !---------------------------------------------
   do j = 1,iiPhytoPlankton
      iiLastElement = iiL
      if (.NOT.CalcPhytoPlankton(j)) then
         do i = iiC,iiS
            k=ppPhytoPlankton(j,i) 
            if ( k.gt.0) then
              D3STATE(:,k) = p_small
              D3STATETYPE(k) = NOTRANSPORT
            endif
         end do
      end if
   end do
   do j = 1,iiMesoZooPlankton
      if (.NOT.CalcMesoZooPlankton(j)) then
         do i = iiC,iiP
            k=ppMesoZooPlankton(j,i) 
            if ( k>0) then
              D3STATE(:,k) = p_small
              D3STATETYPE(k) = NOTRANSPORT
            endif
         end do
      end if
   end do
   do j = 1,iiMicroZooPlankton
      if (.NOT.CalcMicroZooPlankton(j)) then
         do i = iiC,iiP
            k=ppMicroZooPlankton(j,i) 
            if ( k>0) then
              D3STATE(:,k) = p_small
              D3STATETYPE(k) = NOTRANSPORT
            endif
         end do
      end if
   end do
   if ( p_version_Pelbac .ne.4.and.p_version_PelBac.ne.5) then
      D3STATE(:,ppBac)=p_small
   else
      D3STATE(:,ppBac)=0.05*D3STATE(:,ppB1c)
   endif
   do j = 1,iiBenOrganisms
      iiLastElement = iiP
      if (.NOT.CalcBenOrganisms(j)) then
         do i = iiC,iiLastElement
            D2STATE(:,ppBenOrganisms(j,i)) = p_small
         end do
      end if
   end do
   do j = 1,iiBenBacteria
      iiLastElement = iiP
      if (.NOT.CalcBenBacteria(j)) then
         do i = iiC,iiLastElement
            D2STATE(:,ppBenBacteria(j,i)) = p_small
         end do
      end if
   end do
   do j = 1,iiBenPhyto
      iiLastElement = iiS
      if (.NOT.CalcBenPhyto(j)) then
         do i = iiC,iiLastElement
            D2STATE(:,ppBenPhyto(j,i)) = p_small
         end do
      end if
   end do
   if (.NOT.CalcBacteria) then
      B1c = p_small; B1n = p_small; B1p = p_small;
      D3STATETYPE(ppB1c) = NOTRANSPORT
      D3STATETYPE(ppB1n) = NOTRANSPORT
      D3STATETYPE(ppB1p) = NOTRANSPORT
   end if

   call InitTrack
   return

99  FATAL 'I could not read bfm_init_nml'
    stop 'init_var_bfm'
100 FATAL 'I could not read bfm_save_nml'
    stop 'init_var_bfm'
101 FATAL 'I could not read bfm_ben_init_nml'
    stop 'init_var_bfm'
102 FATAL 'I could not read bfm_ben_save_nml'
    stop 'init_var_bfm'

   end subroutine init_var_bfm
!EOC

