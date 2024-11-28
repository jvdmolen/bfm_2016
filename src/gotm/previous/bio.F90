#include"cppdefs.h"
!
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bio --- biological model \label{sec:bio}
!
! !INTERFACE:
   module bio
!
! !DESCRIPTION:
! This is the central module for all biogeochemical models.
! From here, after reading the namelist file {\tt bio.nml},
! the individual biogeochemical model is initialised, the memory
! is allocated, the advection and diffusion is called, the ODE solvers
! for the right hand sides are called, and simple Lagrangian particle
! calculations are managed.
!
! !USES:
   use bio_solver

#ifdef BFM_GOTM
   use bio_bfm, only : init_bio_bfm,var_info_bfm,assign_adv_rates
   use bio_bfm, only : reset_diagonal, allocate_memory_bfm, &
     test_on_all_negative_states, test_on_negative_states, &
     test_mass_conservation,do_nan_test_after_adv
!  use trace_bdy, only:init_trace_bdy,init_var_trace
   use controlled_messages,only:controlled_eval,controlled_output_point, &
                            controlled_output_status

#ifdef SPM
   use spm_bio, only : kd_spm !JMC added: kd by spm from interface to spm module
#endif

#endif

   use output, only: out_fmt

   use util

!
!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_bio, do_bio, end_bio,get_bio_updates,init_var_bio
   logical, public                     :: bio_calc=.false.
!
! !REVISION HISTORY:!
!  Original author(s): Hans Burchard & Karsten Bolding
!
!  $Log: bio.F90,v $
!  Revision 1.29  2005-12-27 11:23:04  hb
!  Weiss 1970 formula now used for surface oxygen saturation calculation
!  in bio_mab.F90
!
!  Revision 1.28  2005-12-27 06:51:49  hb
!  New biomodel bio_mab (bio_iow with additional sediment equation) added
!
!  Revision 1.27  2005-12-02 20:57:27  hb
!  Documentation updated and some bugs fixed
!
!  Revision 1.26  2005-11-18 10:59:35  kbk
!  removed unused variables - some left in parameter lists
!
!  Revision 1.25  2005/11/17 09:58:18  hb
!  explicit argument for positive definite variables in diff_center()
!
!  Revision 1.24  2005/10/11 08:43:44  lars
!  checked new transport routines
!
!  Revision 1.23  2005/09/19 21:07:00  hb
!  yevol replaced by adv_center and diff_center
!
!  Revision 1.22  2005/09/12 14:48:33  kbk
!  merged generic biological module support
!
!  Revision 1.21.2.1  2005/07/06 09:00:19  hb
!  moved init_bio() from do_bio() to time_loop-temporary no NPZD totn calculation
!
!  Revision 1.21  2004/08/18 11:34:14  hb
!  zlev now allocated from 0 to nlev
!
!  Revision 1.20  2004/08/02 11:44:12  kbk
!  bio module compiles and runs with GETM
!
!  Revision 1.19  2004/08/02 08:35:08  hb
!  no need to pass time information
!
!  Revision 1.18  2004/08/01 15:54:49  hb
!  call to light_fasham commented in again
!
!  Revision 1.17  2004/07/30 09:22:20  hb
!  use bio_var in specific bio models - simpliefied internal interface
!
!  Revision 1.16  2004/07/28 11:34:29  hb
!  Bioshade feedback may now be switched on or off, depending on bioshade_feedback set to .true. or .false. in bio.nml
!
!  Revision 1.15  2004/06/29 08:03:16  hb
!  Fasham et al. 1990 model implemented
!
!  Revision 1.14  2004/05/28 13:24:49  hb
!  Extention of bio_iow to fluff layer and surface nutrient fluxes
!
!  Revision 1.13  2004/04/13 09:18:54  kbk
!  size and temperature dependend filtration rate
!
!  Revision 1.12  2004/03/31 12:58:52  kbk
!  lagrangian solver uses - total_mussel_flux
!
!  Revision 1.11  2004/03/30 11:32:48  kbk
!  select between eulerian or lagrangian solver
!
!  Revision 1.10  2003/12/11 09:58:22  kbk
!  now compiles with FORTRAN_COMPILER=IFORT - removed TABS
!
!  Revision 1.9  2003/10/28 10:22:45  hb
!  added support for sedimentation only 1 compartment bio model
!
!  Revision 1.8  2003/10/16 15:42:16  kbk
!  simple mussesl model implemented - filter only
!
!  Revision 1.7  2003/10/14 08:00:09  hb
!  initialise sfl - no special treatment when cc(,) < 0
!
!  Revision 1.6  2003/09/16 12:11:24  hb
!  added new biological model - bio_iow
!
!  Revision 1.5  2003/07/23 12:27:31  hb
!  more generic support for different bio models
!
!  Revision 1.3  2003/04/05 07:01:41  kbk
!  moved bioshade variable to meanflow - to compile properly
!
!  Revision 1.2  2003/04/04 14:25:52  hb
!  First iteration of four-compartment geobiochemical model implemented
!
!  Revision 1.1  2003/04/01 17:01:00  hb
!  Added infrastructure for geobiochemical model
!
! !PRIVATE DATA MEMBERS:
!  from a namelist
   logical                   :: bio_eulerian=.true.
   REALTYPE                  :: cnpar=0.5
   integer                   :: w_adv_discr=P2_PDM
   integer                   :: ode_method=1
   integer                   :: split_factor=1
   logical                   :: bioshade_feedback=.true.
   logical                   :: bio_lagrange_mean=.true.
   integer                   :: bio_npar=10000

   integer                   :: mass_conservation_diff=1
   integer                   :: mass_conservation_adv=1
   logical                   :: nega_test_after_horadv=.false.
   logical                   :: nega_test_after_odv=.false.
   logical,public            :: mass_test_after_odv=.false.
   logical                   :: nan_test_after_adv=.false.
   logical                   :: nan_test_after_horadv=.false.
   logical,public            :: getm_stderr_control=.false.
   integer,public            :: warning_level=1
   logical,public            :: test_ben_on_zero_depth
   REALTYPE                  :: depth
   integer,public            :: ActualStep
   integer,public            :: ModelCallNr=1
   integer,public            :: BFMNeedCheck=0
   character(len=8),public   :: dir_bfmnml
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the bio module
!
! !INTERFACE:
   subroutine init_bio(namlst,fname,unit,nlev)
   use bio_var,only:bio_setup,bio_model,par,ALLTRANSPORT, &
     pelvar_type,pelvar_type_dyn
   use bfm_output,only:var_names,var_units, &
     var_long,stPelStateS,stBenStateS, stPelStateE,stBenStateE !,stPRFDiagS
!
! !DESCRIPTION:
! Here, the bio namelist {\tt bio.nml} is read and memory for the
! Lagrangian part of the model is allocated (note that the
! Lagrangian model up to now only works for the simple suspended matter model).
! If a Lagrangian particle method is chosen, particles are
! equidistantly distributed.
! The initial  Furthermore, information on the specific settings are
! written to standard output.
! Finally, the mussel module is called for initialisation.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: namlst
   character(len=*), intent(in)        :: fname
   integer, intent(in)                 :: unit
   integer, intent(in)                 :: nlev
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer                   :: rc=0
   integer                   ::n
   character(len=1)          :: text_tr
   namelist /bio_nml/ bio_calc,bio_model,bio_eulerian,  bio_setup,  &
          cnpar,w_adv_discr,ode_method,split_factor, &
          bioshade_feedback,bio_lagrange_mean,bio_npar
   namelist /bio_bfm_nml/ &
          mass_conservation_diff,mass_conservation_adv, &
          nega_test_after_horadv,nega_test_after_odv, &
          nan_test_after_adv, nan_test_after_horadv, &
          mass_test_after_odv,warning_level,test_ben_on_zero_depth, &
          dir_bfmnml
!EOP
!-----------------------------------------------------------------------
!BOC

   LEVEL1 'init_bio'

   bio_setup=1
   if ( bio_model == 7 ) then
      bio_setup =0
   endif

!  Open and read the namelist
   test_ben_on_zero_depth=.false.
   open(namlst,file=fname,action='read',status='old',err=98)
   read(namlst,nml=bio_nml,err=99)
   dir_bfmnml='BFM_NML'
   read(namlst,nml=bio_bfm_nml,err=199)
   close(namlst)

   if (bio_calc) then

      allocate(par(0:nlev),stat=rc)
      if (rc /= 0) STOP 'init_bio: Error allocating (par)'


      select case (bio_model)
#ifdef BFM_GOTM
      case (6)  ! The BFM (ERSEM) model
         call init_bio_bfm(nlev)
         call allocate_memory(nlev)
         call var_info_bfm()
         call allocate_memory_bfm(nlev)
!      case (7)  ! The trace model
!        call allocate_memory(nlev)                                 !BFM
!        call init_var_trace(bio_model)
#endif
      case default
         stop "bio: no valid biomodel specified in bio.nml !"
      end select


#ifdef BFM_GOTM
      if (bio_setup /= 2) then
        LEVEL2  'Transported? Pelagic variables:'
        pelvar_type_dyn=pelvar_type
        do n=stPelStateS,stPelStateE
        text_tr='';if (pelvar_type(n)>=ALLTRANSPORT) text_tr='y'
          LEVEL4 text_tr,'    ',trim(var_names(n)),' ',trim(var_units(n)), &
              ' ',trim(var_long(n))
        end do
      endif
      if (bio_setup >= 2) then
        LEVEL3 'Benthic variables:'
        do n=stBenStateS,stBenStateE
          LEVEL4 trim(var_names(n)),'  ',trim(var_units(n)), &
                                    '  ',trim(var_long(n))
        end do
      end if
#endif

      if ( bio_eulerian ) then
         LEVEL3 "Using Eulerian solver"
         select case (ode_method)
            case (1) ; LEVEL2 'Using euler_forward()'
            case (2) ; LEVEL2 'Using runge_kutta_2()'
            case (3) ; LEVEL2 'Using runge_kutta_4()'
            case (4) ; LEVEL2 'Using patankar()'
            case (5) ; LEVEL2 'Using patankar_runge_kutta_2()'
            case (6) ; LEVEL2 'Using patankar_runge_kutta_4()'
            case (7) ; LEVEL2 'Using modified_patankar()'
            case (8) ; LEVEL2 'Using modified_patankar_2()'
            case (9) ; LEVEL2 'Using modified_patankar_4()'
            case (10) ; LEVEL2 'Using emp_1()'
            case (11) ; LEVEL2 'Using emp_2()'
            case default
               stop "bio: no valid o`de_method specified in bio.nml!"
         end select
      end if

   STDERR "bio.F90 I"
   end if

   STDERR "bio.F90 II"
   return

98 LEVEL2 'I could not open bio.nml'
   LEVEL2 'If thats not what you want you have to supply bio.nml'
   LEVEL2 'See the bio example on www`.gotm.net for a working bio.nml'
   bio_calc = .false.
   return
99 FATAL 'I could not read bio_nml name list in bio.nml'
   stop 'stop init_bio'
199 FATAL 'I could not read bio_bfm_nml name list in bio.nml'
   stop 'stop init_bio'
   end subroutine init_bio
!EOC
!-----------------------------------------------------------------------
!BOP
! !IROUTINE: Initialise bio variables
!
! !INTERFACE:
   subroutine init_var_bio(namlst,unit,steps_in_outdelt)
!
! !DESCRIPTION:
! A call to this routine  initializes the variables of the bio module
! with meaningful values, depending on the water column properties set
! by a previous call to {\tt set\_env\_bio()}. The steps taken here
! include
! \begin{enumerate}
!  \item Allocating memory for the particle properties (if a particle
!  solver has been chosen),
!  \item Computing the position of the grid interfaces (needed for
!  averaging particle properties over grid cells)
!  \item Calling the initialization routines for the respective
!  bio modules
! \end{enumerate}
!
! If the bio module is used from a 3D code outside GOTM this routine
! should not be called. In this case, all initialization should be
! done inside the external calling program. Note in particular that
! the numbe of particles may change during the run, and variable memory
! needs to be allocated.
!
! !USES:
   use bio_var,only:bio_setup,bio_model
   use gotm_error_msg, only: get_d3_model_flag_from_getm
   IMPLICIT NONE
   integer, intent(in)                 :: namlst
   integer, intent(in)                 :: unit
   integer, intent(in)                 :: steps_in_outdelt
   logical                             ::d3_model_flag
!
!  Original author(s): Hans Burchard, Lars Umlauf, Karsten Bolding
!
!EOP

   select case (bio_model)

#ifdef BFM_GOTM
      case (6)  ! The BFM (ERSEM) model
         call init_var_bfm(namlst,'bio_bfm.nml', &
                                   unit,bio_setup,steps_in_outdelt)
         STDERR "after init_var_bfm: ontrolled_output_point"
         call get_d3_model_flag_from_getm(d3_model_flag)
         STDERR "d3_model_flag=",d3_model_flag
         if (.not.d3_model_flag) then
             STDERR "init_var_bio"
             !initialzation of the message system which produce info during the
             ! the run about progress, behaviour and dynamics of the model
             call controlled_output_point(-1,1,1)
         endif
!      case (7)  ! The trace model
!        call init_var_trace(bio_model)
#endif
      case default
         stop 'bio: no valid biomodel specified in bio.nml !'
      end select

   end subroutine init_var_bio
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Update the bio model \label{sec:do-bio}
!
! !INTERFACE:
#ifdef BFM_GOTM
   subroutine do_bio(ll_write_results,counter,d3_reset_count)
#endif
!
! !DESCRIPTION:
! This is the main loop for the biogeochemical model. Basically
! an operational split method is used, with first calculating the
! transport part, and than the reaction part.
! During the transport part, all sinks and sources are set to zero,
! and the surface fluxes are computed by calling the
! model specific surface flux subroutine. Then the mussel module
! is called.  For the Eulerian calculation, vertical advection
! (due to settling or rising or vertical migration) and vertical
! diffusion (due to mixing) and afterwards the light
! calculation (for the PAR) and the ODE solver for the right
! hand sides are called.
! It should be noted here that the PAR and the selfshading effect
! is calculated in a similar way for all biogeochemical models
! implemented in GOTM so far. In the temperature equation the
! absorption of solar radiation, $I(z)$, is the only source term,
! see equation (\ref{Iz}) section \ref{sec:temperature}.
! In (\ref{Iz}), a term $B(z)$ due to bioturbidity is used, which
! is calculated as a function of the biogeochemical particulate
! matter in the water column:
! \begin{equation}\label{B}
! B(z)=\exp\left(-k_c\int_z^0\left(\sum C_{turb}(\xi)\right)\,d\xi\right),
! \end{equation}
! where $k_c$ is the attenuation constant for self shading and
! $\sum C_{turb}$ is the sum of the biogeochemical particulate
! matter concentrations.
! The photosynthetically
! available radiation, $I_{PAR}$, follows from
! \begin{equation}
!   \label{light}
!   I_{PAR}(z)=I_0
! (1-a)\exp\left(\frac{z}{\tilde\eta_2}\right)
!   B(z).
! \end{equation}
!
!
! !USES:
   use bio_var,only:nlev,h_l,nuh_l,bio_model,bio_setup,cc,numc,numcc,dt, &
           pelvar_type, adv1d_courant,adv1d_number,ws,c1dimz,sfl,bfl, &
           posconc, adv_courant,adv_number, &
           ALLTRANSPORT,SILTTRANSPORT,llws ! ,cc_before_transport, diagb
   use bfm_output,only:reset_count,write_results
#ifdef BFM_GOTM
   use mem, only: ppR6s,ppN3n,ppN1p !track_error,
   use gotm_error_msg, only:gotm_error,set_warning_for_getm , &
                            get_d3_model_flag_from_getm
   use bio_bfm,only: output_warning_shiftings, do_nan_test_after_adv
#endif
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   logical,intent(in)                  :: ll_write_results
   integer,intent(inout),optional      :: counter
   logical,intent(in),optional         :: d3_reset_count
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   REALTYPE                  :: Qsour(0:nlev),Lsour(0:nlev)
   REALTYPE                  :: RelaxTau(0:nlev)
   REALTYPE                  :: dt_eff
   REALTYPE                  :: totsys_old(1:5,2,1)
   integer                   :: j
   integer                   :: split
   logical                   :: llsumh,llsumh2
#ifdef BFM_GOTM
   character(len=160)                 :: msg
   integer                   :: i
   integer                   :: kt=0,local_counter
   logical                   :: d3_model_flag
!  integer,save               ::follow=0
#endif
!EOP
!-----------------------------------------------------------------------
!BOC
   if (bio_calc) then
      write_results=ll_write_results
      Qsour    = _ZERO_
      Lsour    = _ZERO_
      RelaxTau = 1.e15

      call get_d3_model_flag_from_getm(d3_model_flag)
      ! Somtimes it happens that a concentration become negative
      ! after calculations of the 3d-transport
      llsumh2=(sum(h_l(1:nlev)) .gt. 0.21D+00)
      llsumh=(sum(h_l(1:nlev)) .gt. 0.21D+00)
!     if (ModelCallNr==2) STDERR 'I P6c(1)=', cc(1:3,ppP6c)
      if ((bio_model==6).and.d3_model_flag) then
        if (nega_test_after_horadv) & 
          call test_on_all_negative_states(d3_model_flag, llsumh2,bio_setup, &
                      h_l, numcc,nlev,warning_level, 'hor_adv', cc, kt )
        if ( kt.lt.0) then
          call gotm_error('do_bio', 'negative state value'); return
        elseif (kt.gt.0) then
!         call set_warning_for_getm
        endif
      endif

      if ( llsumh ) then
        do split=1,split_factor
          dt_eff=dt/float(split_factor)
          call test_mass_conservation(mass_test_after_odv,1,'',dt,totsys_old)
          if (.not.d3_model_flag)  ActualStep=ActualStep+1
          if ( bio_setup >0 ) then
            call ode_solver(ode_method,numc,nlev,dt_eff,cc)
          endif
          if (BFMNeedCheck>0) then
             STDERR "After ode_solver:"
             STDERR  cc(:,ppR6s)
          endif
          if (bio_model==6.and.nega_test_after_odv) then
            call test_on_all_negative_states (.true., llsumh2, &
              bio_setup, h_l, numcc,nlev,warning_level, 'ode_solver', cc, kt )
            if ( kt.lt.0) then
              call gotm_error('do_bio', 'negative state value'); return
            elseif (kt.gt.0) then
              call set_warning_for_getm
            endif
          endif
          call test_mass_conservation(mass_test_after_odv,3,'ode_solver', &
                                                   dt,totsys_old)
        end do !split_factor
      endif

      if (bio_eulerian) then
#ifdef BFM_GOTM
        local_counter=0;if (present(counter)) local_counter=counter
        if (present(d3_reset_count)) reset_count=d3_reset_count
        if (bio_model==6) then
          call assign_adv_rates
          do j=1,numcc
            if (bio_setup /= 2 ) then
              ! inclusive SILTTRANSPORT
              if (pelvar_type(j)>=ALLTRANSPORT) then
                if (warning_level >2) then
                  adv_courant=adv1d_courant(j)
                  adv_number=adv1d_number(j)
                endif
                c1dimz=cc(:,j)
!               if (modelCallNr==5.and.j.eq.ppR6s) &
!             STDERR "ModelCallNr,j,llsumh,llws(j)",ModelCallNr,j,llsumh,llws(j)
                if (llsumh ) then
                  if (llws(j)) then
                    call adv_center_bfm(j,nlev,dt,h_l,h_l,ws(:,j),flux,flux, &
                    _ZERO_,_ZERO_,w_adv_discr,mass_conservation_adv,c1dimz)
!-----------------JM added
                    msg='bio+advection'
                    call test_on_negative_states ( j,llsumh2,h_l,nlev,  &
                               trim(msg),c1dimz, kt,msg, i )
                    if (nan_test_after_adv) then
                      kt=do_nan_test_after_adv(j,nlev,cc(:,j), c1dimz,ws(:,j))
                      if (kt<0) then
                        STDERR 'oldcc=',cc(:,j)
                        STDERR 'newcc=',c1dimz
                        STDERR 'ws=',ws
                      endif
                    endif
                    msg='bio+advection+vdiff'
!-----------------JM end addition
                  else
                    msg='bio+vdiff'
                  endif
                 ! do diffusion step
!JM added for ppR9x: SILTTRANSPORT: Dirichlet bottom boundary condition
                  if (pelvar_type(j)>=SILTTRANSPORT.and.llws(j)) then
                    c1dimz(1)=cc(1,j)
                    call diff_center_bfm(nlev,dt,cnpar,posconc(j),h_l, &
                      Neumann,Dirichlet,sfl(j),cc(1,j),nuh_l,Lsour,Qsour, &
                      RelaxTau,c1dimz,mass_conservation_diff,c1dimz)
                  else
                    call diff_center_bfm(nlev,dt,cnpar,posconc(j),h_l, &
                      Neumann,Neumann,sfl(j),bfl(j),nuh_l,Lsour,Qsour, &
                      RelaxTau,c1dimz,mass_conservation_diff,c1dimz)
                  endif
                  call test_on_negative_states ( j,llsumh2,h_l,nlev,  &
                       trim(msg),c1dimz, kt,msg, local_counter )
                else
                  call diff_center_bfm(nlev,dt,cnpar,posconc(j),h_l, &
                   Neumann,Neumann,_ZERO_,_ZERO_,nuh_l,Lsour,Qsour, &
                   RelaxTau,c1dimz,mass_conservation_diff,c1dimz)
                  call test_on_negative_states ( j,llsumh2,h_l,nlev,  &
                       "bio+vert.diff" ,c1dimz,kt,msg,local_counter)
                endif
                if (present(counter))then
                  if (local_counter.gt.counter.and.kt.gt.0) kt=0
                  counter=local_counter
                endif
                if ((pelvar_type(j)>=ALLTRANSPORT.and.warning_level>1) &
                     .and.kt.ne.0) then
                   STDERR trim(msg)
                   call set_warning_for_getm
!               elseif (present(counter))then
!                 counter=counter+1
                endif
!               if (ModelCallNr==2.and.j==ppP6c) STDERR 'II P6c(1)=', c1dimz(1:3)
                cc(:,j)=c1dimz
                if (pelvar_type(j)>=ALLTRANSPORT.and.warning_level>2) then
                  adv1d_courant(j)=adv_courant
                  adv1d_number(j)=adv_number
                endif
              end if
            end if
          enddo
          if (.not.d3_model_flag) &
                call output_warning_shiftings(warning_level,local_counter)
          if ( kt.lt.0) then
            call gotm_error('do_bio', 'negative state value')
            return
          endif
          call CalcSiltResuspension()
          ! Output collected warnings during last timestep
          j=0;if (write_results)j=1
          !setup of the message system which produce info during the
          !run about progress, behaviour and dynamics of the model by counting
          ! the number of calls to  the routine controlled_messages in the code 
          call controlled_eval(j,'',control=0)
          if ((.not.d3_model_flag).and.(controlled_output_status()==-1)) then
            !first time initialization: allocate memory on basis of the 
            !  messages counted in controlled_output  
            call controlled_output_point(0,0,0)
            call controlled_output_point(1,1,1)
          endif
          ! Clear all collected warnings
          j=0;if(write_results)j=1; call controlled_eval(j,'',control=1)
          BFMNeedCheck=0
        elseif (bio_model==7) then
          do j=1,numcc
!           do diffusion step
            c1dimz=cc(:,j)
            call diff_center(nlev,dt,cnpar,posconc(j),h_l,Neumann,Neumann,&
                sfl(j),bfl(j),nuh_l,Lsour,Qsour,RelaxTau,c1dimz,c1dimz)
            call test_on_negative_states ( j,llsumh2,h_l,nlev,  &
                                      "bio+vert.diff" ,c1dimz,kt,msg )
            if (kt.gt.0) STDERR msg(1:len_trim(msg))
            cc(:,j)=c1dimz
          end do
        else ! other bio_model
#endif
          do j=1,numcc
!           do advection step
            if ( bio_model /= 7 ) call adv_center(nlev,dt,h_l,h_l,ws(:,j),&
                 flux, flux,_ZERO_,_ZERO_,w_adv_discr,cc(:,j))

!           do diffusion step
            call diff_center(nlev,dt,cnpar,posconc(j),h_l,Neumann,Neumann,&
                sfl(j),bfl(j),nuh_l,Lsour,Qsour,RelaxTau,cc(:,j),cc(:,j))
          end do
#ifdef BFM_GOTM
         end if !bio_model
#endif
      end if
   end if
   return
   end subroutine do_bio
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finish the bio calculations
!
! !INTERFACE:
   subroutine end_bio
!
! !DESCRIPTION:
!  Nothing done yet --- supplied for completeness.
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!
!EOP
!-----------------------------------------------------------------------
!BOC

   return
   end subroutine end_bio
!EOC
!!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Allocate memory for biological variables
!
! !INTERFACE:
   subroutine allocate_memory(nlev)
!
! !DESCRIPTION:
! Here, the memory for the global biogeochemical parameters
! such as concentrations, settling velocities, surface and bottom
! boundary fluxes, and various other parameters is allocated.
!
! !USES:
   use bio_var,only:numc, &
         cc,nuh_l,h_l,c1dimz,c1dimnumc,ws, &
         llws,bfl,sfl,posconc
#ifdef BFM_GOTM
   use bio_var,only:adv1d_courant,adv1d_number,pp,dd,pelvar_type, &
     pelvar_type_dyn
#endif
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,  intent(in)                :: nlev
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer                   :: rc=0
!EOP
!-----------------------------------------------------------------------
!BOC
   allocate(cc(0:nlev,1:numc),stat=rc)
   if (rc /= 0) STOP 'init_bio: Error allocating (cc)'
   cc=_ZERO_

   allocate(nuh_l(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_bio: Error allocating (nuh_l)'
   nuh_l=_ZERO_

   allocate(h_l(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_bio: Error allocating (h_l)'
   h_l=_ZERO_

   allocate(c1dimz(0:nlev),stat=rc)
   if (rc /= 0) STOP 'init_bio: Error allocating (c1dimz)'
   c1dimz=_ZERO_

   allocate(c1dimnumc(1:numc),stat=rc)
   if (rc /= 0) STOP 'init_bio: Error allocating (c1dimnumc)'
   c1dimnumc=_ZERO_

   allocate(ws(0:nlev,1:numc),stat=rc)

   if (rc /= 0) STOP 'init_bio: Error allocating (ws)'
   ws=_ZERO_

   allocate(llws(1:numc),stat=rc)
   if (rc /= 0) STOP 'init_bio: Error allocating (llws)'
   llws=.false.

   allocate(sfl(1:numc),stat=rc)
   if (rc /= 0) STOP 'init_bio: Error allocating (sfl)'
   sfl=_ZERO_

   allocate(bfl(1:numc),stat=rc)
   if (rc /= 0) STOP 'init_bio: Error allocating (bfl)'
   bfl=_ZERO_

   allocate(posconc(1:numc),stat=rc)
   if (rc /= 0) STOP 'init_bio: Error allocating (posconc)'
   posconc=1

#ifdef BFM_GOTM
   allocate(adv1d_courant(1:numc),stat=rc)
   if (rc /= 0) STOP 'init_bio: Error allocating (adv1d_courant)'
   adv1d_courant=_ZERO_

   allocate(adv1d_number(1:numc),stat=rc)
   if (rc /= 0) STOP 'init_bio: Error allocating (adv1d_number)'
   adv1d_number=_ZERO_

   allocate(pp(0:nlev,1:numc,1:numc),stat=rc)
   if (rc /= 0) STOP 'init_bio: Error allocating (pp)'
   pp=_ZERO_

   allocate(dd(0:nlev,1:numc,1:numc),stat=rc)
   if (rc /= 0) STOP 'init_bio: Error allocating (dd)'
   dd=_ZERO_

   ! allocate variable holding type and save attributes
   allocate(pelvar_type(1:numc),stat=rc)
   if (rc /= 0) STOP 'init_bio: Error allocating (pelvar_type)'
   pelvar_type = 0
   ! allocate variable holding type and save attributes
   allocate(pelvar_type_dyn(1:numc),stat=rc)
   if (rc /= 0) STOP 'init_bio: Error allocating (pelvar_type_dyn)'
   pelvar_type_dyn = 0

#endif
   return
   end subroutine allocate_memory

!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Return updated bio variable
!
! !INTERFACE:
   subroutine get_bio_updates(nlev,bioshade)
   use mem,   ONLY: xEPS, Depth
   use global_interface,   ONLY: CalcVerticalExtinction

   implicit none
   integer, intent(IN)                  :: nlev
   REALTYPE, intent(out),optional       :: bioshade(0:nlev)!

   integer                              :: i
   if (present(bioshade).and.bioshade_feedback) then
     ! 0= Special calculation of vertical extinction only controlled by
     !   biological constituents for use in gotm
     ! 1= full calulation used to caluclate bilolgical vertical extinction
     ! 2= Special calculation of vertical extinction only controlled by
     ! biological constituents and silt when silt is included in the BFM part.

     call  CalcVerticalExtinction(0)
     bioshade(1)=1.0
     do i=nlev,2,-1
       bioshade(i-1) = bioshade(i)*exp(-xEPS(i)*Depth(i))
     end do
     bioshade(1:nlev) =  bioshade(1:nlev)*exp(-xEPS(1:nlev)*Depth(1:nlev)*0.5)
   endif
   end subroutine get_bio_updates

!-----------------------------------------------------------------------
!BOP
!

!

!-----------------------------------------------------------------------


   end module bio

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!----------------------------------------------------------------------
