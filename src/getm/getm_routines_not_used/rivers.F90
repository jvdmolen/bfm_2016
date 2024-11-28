!$Id: rivers.F90,v 1.11 2006-03-01 15:54:08 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE:  rivers \label{sec-rivers}
!
! !INTERFACE:
   module rivers
!
! !DESCRIPTION:
!
!  This module includes support for river input. Rivers are treated the same
!  way as meteorology, i.e.\ as external module to the hydrodynamic model 
!  itself.
!  The module follows the same scheme as all other modules, i.e.\
!  {\tt init\_rivers}
!  sets up necessary information, and {\tt do\_rivers} updates
!  the relevant variables.
!  {\tt do\_river} is called in {\tt getm/integration.F90} 
!  between the {\tt 2d} and {\tt 3d} routines as it only
!  updates the sea surface elevation (in {\tt 2d}) and sea surface elevation,
!  and
!  optionally salinity and temperature (in {\tt 3d}). 
!  At present the momentum of the river water is not include, the model
!  however has a direct response to the river water because of the
!  pressure gradient introduced.
!
!  !BFM
!   code added to read boundary condition for BFM state variables
!
!
! !USES:
   use domain, only: imin,jmin,imax,jmax,ioff,joff
#if defined(SPHERICAL) || defined(CURVILINEAR)
   use domain, only: H,az,kmax,arcd1
#else
   use domain, only: H,az,kmax,ard1
#endif
   use m2d, only: dtm
   use variables_2d, only: z
#ifndef NO_BAROCLINIC
   use m3d, only: calc_salt,calc_temp
   use variables_3d, only: hn,ssen,T,S,rho
#endif
#ifdef GETM_BIO
   use bio, only: bio_calc
   use bio_var, only: numc
#ifdef BFM_GOTM
   use bio_var, only: c1dimnumc,cc
   use variables_bio_3d, only: cc3d
   use variables_bio_3d, only: bio_missing
#else
   use variables_3d, only: cc3d
#endif
#endif
   IMPLICIT NONE
!
   private
!
! !PUBLIC DATA MEMBERS:
   public init_rivers, do_rivers, clean_rivers
#ifdef GETM_BIO
   public init_rivers_bio
#endif
   integer, public                     :: river_method=0,nriver=0,rriver=0,nriver_sa=0
   logical,public                      :: use_river_temp = .false.
   logical,public                      :: use_river_salt = .false.
   character(len=64), public           :: river_data="rivers.nc"
   character(len=64), public, allocatable  :: river_name(:)
   character(len=64), public, allocatable  :: real_river_name(:)
   integer, public, allocatable        :: ok(:)
   REALTYPE, public, allocatable       :: river_flow(:)
   REALTYPE, public, allocatable       :: river_salt(:)
   REALTYPE, public, allocatable       :: river_temp(:)
   integer, public                     :: river_ramp= -1
   REALTYPE, public                    :: river_factor= _ONE_
   REALTYPE, public,parameter          :: temp_missing=-9999.0
   REALTYPE, public,parameter          :: salt_missing=-9999.0
   integer,  public, allocatable       :: river_split(:)
#ifdef GETM_BIO
   REALTYPE, public, allocatable       :: river_bio(:,:)
   REALTYPE, public, allocatable       :: river_min(:),river_max(:)
   integer,public,allocatable          :: river_func(:)                !concentration is calculated using a function!
   REALTYPE, public, allocatable       :: river_add(:)              !to do sens.analysis
   REALTYPE, public, allocatable       :: river_frc(:)              !one can modify concentrations
   integer, public, allocatable        :: seq_pointer_river(:)      !on the rivers 
#endif
!
! !PRIVATE DATA MEMBERS:
   integer                   :: river_format=2
   character(len=64)         :: river_info="riverinfo.dat"
   character(len=64)         :: river_sa="river-sa.dat"
   integer, allocatable      :: ir(:),jr(:)
   REALTYPE, allocatable     :: irr(:)
   REALTYPE, allocatable     :: macro_height(:)
   REALTYPE, allocatable     :: flow_fraction(:)
   integer                   :: macro_steps
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_rivers
!
! !INTERFACE:
   subroutine init_rivers
!
! !DESCRIPTION:
!
! First of all, the namelist {\tt rivers} is read from getm.F90 and
! a number of vectors with the length of {\tt nriver} (number of
! rivers) is allocated. Then, by looping over all rivers, the 
! ascii file {\tt river\_info} is read, and checked for consistency.
! The number of used rivers {\tt rriver} is calculated and it is checked
! whether they are on land (which gives a warning) or not. When a river name
! occurs more than once in {\tt river\_info}, it means that its runoff
! is split among several grid boxed (for wide river mouths).
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !LOCAL VARIABLES:
   integer                   :: i,j,n,nn,ni,m
   integer                   :: rc=0
   integer                   :: unit = 25 ! kbk
   logical                   :: outside
   REALTYPE                  :: area
   character(len=80)         :: msg
   NAMELIST /rivers/ &
            river_method,river_info,river_format,river_data,river_ramp, &
            river_factor,use_river_salt,use_river_temp
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_rivers() # ',Ncall
#endif

   LEVEL1 'init_rivers()'
   read(NAMLST,rivers)

   select case (river_method)
      case (0)
         LEVEL3 'River runoff not included.'
      case (1,2)
         LEVEL2 'river_method= ',river_method
         LEVEL2 'river_data=   ',trim(river_data)
         LEVEL2 'river_format= ',river_format
         LEVEL2 'river_ramp=   ',river_ramp
         LEVEL2 'river_factor= ',river_factor
         LEVEL2 'use_river_temp= ',use_river_temp
         LEVEL2 'use_river_salt= ',use_river_salt
         open(unit,file=river_info,action='read',status='old',err=90)
         read(unit,*) nriver
         allocate(ir(nriver),stat=rc) ! i index of rivers
         if (rc /= 0) stop 'rivers: Error allocating memory (ir)'
         allocate(jr(nriver),stat=rc) ! j index of rivers
         if (rc /= 0) stop 'rivers: Error allocating memory (jr)'
         allocate(ok(nriver),stat=rc) ! valid river spec.
         if (rc /= 0) stop 'rivers: Error allocating memory (ok)'
         allocate(river_name(nriver),stat=rc) ! NetCDF name of river.
         if (rc /= 0) stop 'rivers: Error allocating memory (river_name)'
         allocate(river_flow(nriver),stat=rc) ! river flux
         if (rc /= 0) stop 'rivers: Error allocating memory (river_flow)'
         allocate(macro_height(nriver),stat=rc) ! height over a macro tims-step
         if (rc /= 0) stop 'rivers: Error allocating memory (macro_height)'
         allocate(river_temp(nriver),stat=rc) ! temperature of river water
         if (rc /= 0) stop 'rivers: Error allocating memory (river_temp)'
         allocate(river_salt(nriver),stat=rc) ! salinity of river water
         if (rc /= 0) stop 'rivers: Error allocating memory (river_salt)'
         allocate(river_split(nriver),stat=rc) ! split factor for river water
         if (rc /= 0) stop 'rivers: Error allocating memory (river_split)'
         allocate(flow_fraction(nriver),stat=rc) ! areafactor of data for river 
         if (rc /= 0) stop 'rivers: Error allocating memory (flow_fraction)'
         allocate(irr(nriver),stat=rc) ! integrated river runoff
         if (rc /= 0) stop 'rivers: Error allocating memory (irr)'

         ok = 0
         rriver = 0 ! number of real existing rivers...
         flow_fraction = _ZERO_
         macro_steps=0
         do n=1,nriver
            read(unit,*) ir(n),jr(n),river_name(n)
            river_name(n) = trim(river_name(n))
            i = ir(n)-ioff
            j = jr(n)-joff
            river_temp(n) = temp_missing
            river_salt(n) = salt_missing
            river_flow(n) = _ZERO_
            irr(n) = _ZERO_
            macro_height(n) = _ZERO_
!           calculate the number of used rivers, they must be 
!           in sequence !
            rriver = rriver +1
            if ( n .gt. 1 ) then
               if (river_name(n) .eq. river_name(n-1))  rriver = rriver-1
            end if
            outside= &
                 i .lt. imin .or. i .gt. imax .or.  &
                 j .lt. jmin .or. j .gt. jmax
            if( .not. outside) then
               if(az(i,j) .eq. 0) then
                  LEVEL3 trim(river_name(n)),':',ir(n),jr(n),' n=',n,'Outside: river on land!'
                  ok(n) = 0
               else
                  ok(n) = 1
                  flow_fraction(n) = _ONE_/ARCD1
                  write(msg,'(A,'':'',I4,'' ('',I3,'') '',I4,''( '',I3,'') n='',I4)') &
                            trim(river_name(n)),ir(n),i,jr(n),j,n
                  LEVEL3 trim(msg)
               end if
            else
              LEVEL3 trim(river_name(n)),':',ir(n),jr(n),' n=',n,' Outside!'
            end if
         end do
         close(unit)

!        calculate the number of used gridboxes, they must be 
!        in sequence !
         LEVEL3 'Number of unique rivers: ',rriver
         allocate(real_river_name(rriver),stat=rc) ! NetCDF name of river.
         if (rc /= 0) stop 'rivers: Error allocating memory (rivers)'
         river_split = 1    ! normal case
         do n=2,nriver
           if (river_name(n) .eq. river_name(n-1))  river_split(n)=river_split(n-1)+1
         end do
         ni= nriver
         do n=1,nriver
            if (ni .ge. 1) then
               if ( river_split(ni) .gt. 1 ) then  
                  do m=1,river_split(ni)
                     river_split(ni-m+1) =  river_split(ni)
                  end do
               end if
               ni = ni - river_split(ni)
            end if
         end do
         LEVEL3 'split:',river_split
!        now river_split contains the number of gridboxes used 
!        for a single river
         nn = 1
         ni = 1
         do n=1,nriver
            if (ni .le. nriver) then
               real_river_name(nn) = river_name(ni)
               if ( river_split(ni) .gt. 1 ) then
                  area = _ZERO_
                  do m=1,river_split(ni) 
                     area = area +  flow_fraction(ni+m-1)
                  end do
                  do m=1,river_split(ni)
                     if ( area .gt. _ZERO_ ) then
                        flow_fraction(ni+m-1) = flow_fraction(ni+m-1)/area
                     else
                        flow_fraction(ni+m-1) = _ZERO_
                     end if
                  end do
               else
                  flow_fraction(ni) = _ONE_
               end if
               nn = nn + 1  
               ni = ni + river_split(ni)
            end if 
            if (ok(n) .eq. 0) then
               flow_fraction(n) = _ZERO_
            end if
         end do

      case default
         FATAL 'A non valid river_method has been selected'
         stop 'init_rivers'
   end select
   return

90 LEVEL2 'could not open ',trim(river_info),' for reading info on rivers'
   stop 'init_rivers()'

#ifdef DEBUG
   write(debug,*) 'Leaving init_rivers()'
   write(debug,*)
#endif
   return
   end subroutine init_rivers
!EOC

#ifdef GETM_BIO
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_rivers_bio
!
! !INTERFACE:
   subroutine init_rivers_bio()
!
! !DESCRIPTION:
! First, memory for storing the biological loads from rivers is 
! allocated.
! The variable - {\tt river\_bio} - is initialised to  - {\tt bio\_missing}.
!
! !USES:
   use coupling_getm_bfm,ONLY: assign_river_inputs_to_bio_states
#ifdef BFM_GOTM
   use string_functions, ONLY: getseq_number
   use bfm_output, only: var_names
#else
   use bio_var, only: var_names
#endif

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !LOCAL VARIABLES:
   integer                   :: rc=0
   integer                   :: i,j,k,l,m,n,isum,nok
   character(len=64)         :: l_river_name=''
   character(len=64)         :: l_var_name=''
   REALTYPE                  :: add=_ZERO_,frc=_ZERO_
   integer                   :: unit = 25 ! kbk
   logical                   :: ll_exist=.false.
!EOP
!-------------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'init_rivers_bio() # ',Ncall
#endif

#ifdef BFM_GOTM
   LEVEL1 'init_rivers_bio()'

   select case (river_method)
      case(1,2)
         allocate(river_bio(rriver,numc),stat=rc)
         if (rc /= 0) stop 'rivers: Error allocating memory (river_bio)'
         river_bio = bio_missing

         allocate(river_max(numc),stat=rc)
         if (rc /= 0) stop 'rivers: Error allocating memory (river_max)'
         river_max=bio_missing
         allocate(river_min(numc),stat=rc)
         if (rc /= 0) stop 'rivers: Error allocating memory (river_min)'
         river_min=bio_missing
         allocate(river_func(numc),stat=rc)
         if (rc /= 0) stop 'rivers: Error allocating memory (river_func)'
         river_func=0

         k=0
         rc= assign_river_inputs_to_bio_states(0,'',0,0,k)
         if ( rc.ne.0) call errsns(rc)

         inquire(file=river_sa,exist=ll_exist)
         if ( ll_exist) then
           open(unit,file=river_sa,action='read',status='old',err=90)
           read(unit,*) nriver_sa
           allocate(river_frc(nriver_sa),stat=rc)
           if (rc/= 0) stop 'rivers: Error allocating memory(river_frc)'
           allocate(river_add(nriver_sa),stat=rc)
           if (rc/= 0) stop 'rivers: Error allocating memory(river_add)'
           allocate(seq_pointer_river(nriver_sa),stat=rc)
           if (rc/= 0) stop 'rivers: Error allocating memory(seq_pointer_river)'
           n=0
           do i=1,nriver_sa
              read(unit,*)  l_river_name,l_var_name, frc,add
              m=getseq_number(l_var_name,var_names,numc,.TRUE.)
              k=getseq_number(l_river_name,river_name,nriver,.TRUE.)
              if ( ok(k).gt. 0 ) then
                 j=m+ 100 *k
                 n=n+1 ; l=1
                 if ( n>1 ) then
                   do k=1,n-1
                     if ( j > seq_pointer_river(k) ) l=k+1
                   enddo
                   do k=n,l+1,-1
                     seq_pointer_river(k)= seq_pointer_river(k-1)
                     river_frc(k)= river_frc(k-1)
                     river_add(k)= river_add(k-1)
                   enddo
                 endif
                 seq_pointer_river(l)=j
                 river_frc(l)=frc
                 river_add(l)=add
                 LEVEL2 'Sensitivity analysis for river.shortname:',l_river_name
                 LEVEL2 'Multiplication fraction :', river_frc(l)
                 LEVEL2 'Load addition:', river_add(l)
              endif
           enddo
           nriver_sa=n
           close(unit)
         endif
   end select
#endif


#ifdef DEBUG
   write(debug,*) 'Leaving init_rivers_bio()'
   write(debug,*)
#endif
   return
90 LEVEL2 'could not open ',trim(river_sa),' for reading info on rivers'
   stop 'init_rivers_bio()'
   end subroutine init_rivers_bio
!EOC
!-----------------------------------------------------------------------
#endif

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  do_rivers - updating river points \label{sec-do-rivers}
!
! !INTERFACE:
   subroutine do_rivers(do_3d)
!
! !DESCRIPTION:
! 
! Here, the temperature, salinity, sea surface elevation and layer heights
! are updated in the river inflow grid boxes. Temperature and salinity
! are mixed with riverine values proportional to the old volume and the
! river inflow volume at that time step, sea surface elevation is simply
! increased by the inflow volume divided by the grid box area, and
! the layer heights are increased proportionally. 
!
! !USES:
   use coupling_getm_bfm,only:make_river_flux_output
   use global_interface,only:CalcRiverConcentration
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   logical, intent(in)                 :: do_3d
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !LOCAL VARIABLES:
   integer                   :: i,j,k,m,n,ni,nk
   integer,save              :: nn=0
   REALTYPE                  :: ramp=_ONE_
   REALTYPE                  :: rvol,height
   REALTYPE                  :: svol,tvol,vol,r_extra
   REALTYPE                  :: dt_thermo,db_thermo,cor_thermo,cor_dep
   integer                   :: it_thermo,ib_thermo,kcor
   integer                   :: sa_n,sa_l 
   REALTYPE                  :: rivinput,r_av
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'do_rivers() # ',Ncall
#endif

!  river spin-up
   ramp=_ONE_
   if (river_ramp .gt. 0 .and. nn .lt. river_ramp) then
      ramp=min( _ONE_ , nn*_ONE_/river_ramp)
      nn=nn+1
   end if


   select case (river_method)
      case(0)
      case(1,2)
#ifdef BFM_GOTM
       !code to do sensitivity analysis with rivers
       sa_n=0;if ( nriver_sa > 0 ) sa_n=1
#endif
       macro_steps=macro_steps+1
       ni=0;nk=0
       do n=1,nriver
         if ( ni.lt.n) then
           nk=nk+1;ni=ni+river_split(ni+1)
         endif
         if(ok(n) .gt. 0) then
           i = ir(n)-ioff; j = jr(n)-joff
           rvol = dtm * river_flow(nk) * flow_fraction(n)
           irr(n) = irr(n) + rvol
           height = rvol * ARCD1
           z(i,j) = z(i,j) + height
#ifndef NO_BAROCLINIC
           macro_height(n)=macro_height(n)+height
!              on macrotime step adjust 3d fields
!          call CalculateThermovar(kmax,0.01D+00,hn(i,j,1:kmax),  &
!          rho(i,j,1:kmax),ib_thermo,db_thermo,it_thermo,dt_thermo,cor_thermo)
           kcor=kmax
!          if (dt_thermo.lt.-10.0) then
!            kcor=it_thermo-1
!            cor_thermo=sum(hn(i,j,it_thermo:kmax))
!            cor_dep=macro_height(n)/(H(i,j) &
!                               +ssen(i,j))/(H(i,j)+ssen(i,j)-cor_thermo) !1/m
!          endif
           if (do_3d) then
             if (calc_salt) then
               if ( abs(river_salt(nk) - salt_missing ) > 1.0E-3 ) then
                 rivinput=river_salt(nk) 
!                if (kcor<kmax) STDERR 'thermo:',it_thermo,rivinput
                 if (kcor<kmax) & 
                 rivinput=rivinput +cor_dep*sum(hn(i,j,it_thermo:kmax) &
                                       *(rivinput-S(i,j,it_thermo:kmax)))
                 S(i,j,1:kcor) = (S(i,j,1:kcor)*(H(i,j)+ssen(i,j))   &
                              + rivinput*macro_height(n))      &
                                / (H(i,j)+ssen(i,j)+macro_height(n))
!                if (kcor<kmax) STDERR 'thermo2:',it_thermo,rivinput
               else
                 S(i,j,1:kmax) = S(i,j,1:kmax)*(H(i,j)+ssen(i,j))   &
                              / (H(i,j)+ssen(i,j)+macro_height(n))
               end if
             end if
             if (calc_temp .and. abs(river_temp(nk) -temp_missing)>1.0E-3 ) then
               rivinput=river_temp(nk) 
               if (kcor<kmax) & 
               rivinput=rivinput +cor_dep*sum(hn(i,j,it_thermo:kmax) &
                                       *(rivinput-T(i,j,it_thermo:kmax)))
               T(i,j,1:kcor) = (T(i,j,1:kcor)*(H(i,j)+ssen(i,j))   &
                               + rivinput*macro_height(n))   &
                               / (H(i,j)+ssen(i,j)+macro_height(n))
            end if
#ifdef BFM_GOTM
            if (bio_calc) then
            cc=cc3d(i,j,:,:)
            c1dimnumc=-1.0
              do m=1,numc
                !if function is used to caluclate the conc, of a constituent
                ! (see further: river_func ) river_bio(nk,m)=0.0
                !                            river_max(m)=-99999
                ! resulting in the fact that the program will NOT pass any of 
                ! if/elseifstatements in this loop.
                rivinput=-1.
                !if there is a timeseries river_bio(nk,m) >0.0
                if ( river_bio(nk,m) .gt.1.0D-80 ) then
                  c1dimnumc(m)=river_bio(nk,m)
                  !code to do sensitivity analysis with rivers
                  if ( sa_n > 0 ) then
                    sa_l=100*n+m
                    if ( sa_l == seq_pointer_river(sa_n)) then
                      rivinput=rivinput * river_frc(sa_n)
                      rivinput=rivinput+ &
                                  river_add(sa_n)*dtm *flow_fraction(n) *ARCD1
                      sa_n=sa_n+1;if ( sa_n > nriver_sa) sa_n=0
                    endif
                  endif
                  !end code to do sensitivity analysis with rivers
                elseif (abs(river_bio(nk,m)-bio_missing).lt.1.0E-3) then
                  !if river_bio== bio-missing:there is no boundary condition....
                  !all tracked states will not modified by a not-tracked river.
                  ! HOWEVER THERE IS A DILUTION OF THE CONCENTRATION!
                  rivinput=1.0D-80                   
                  c1dimnumc(m)=rivinput
                elseif (abs(river_max(m)-bio_missing).gt.1.0E-3) then
                  ! code to limit riverinput if no data are given:
                  ! it assumed that concentration on the river is the same as 
                  ! in grid point where the river enters but limited!
                  r_av= &
                      sum(cc3d(i,j,1:kcor,m)*hn(i,j,1:kcor))/sum(hn(i,j,1:kcor))
                  rivinput=max(river_min(m),min(river_max(m),r_av))
                  c1dimnumc(m)=rivinput
                endif
              end do
              do m=1,numc
                r_extra=cc3d(i,j,kcor,m) 
                !code to calculate riverconcentration of a state var on basis 
                !of concentrations of other constituents in the same river.
                if ( river_func(m) == 1) then
!                 call CalcRiverConcentration(m,1,numc,T(i,j,kmax), &
!                  o   cc3d(i,j,1:kmax,1:numc),c1dimnumc)
                  call CalcRiverConcentration(m,T(i,j,1:kmax), &
                      cc(1:kmax,1:numc),cc3d(i,j,1:kmax,1:numc))
                endif
                !here the concentrations at the river points are modified::
                if ( c1dimnumc(m) >= 0.0) then
                  rivinput=c1dimnumc(m)
                  if (kcor<kmax) & 
                  rivinput=c1dimnumc(m) +cor_dep*sum(hn(i,j,it_thermo:kmax) &
                                   *(c1dimnumc(m)-cc3d(i,j,it_thermo:kmax,m)))
                  cc3d(i,j,1:kcor,m) = (cc3d(i,j,1:kcor,m)* &
                     (H(i,j)+ssen(i,j)) + rivinput*macro_height(n)) &
                                / (H(i,j)+ssen(i,j)+macro_height(n))
                   r_extra=(cc3d(i,j,kcor,m)-r_extra)/(1.0D-80+r_extra)
!                  if (r_extra> 5.0) &
!                      STDERR "WarningRivers: large change in:",m,r_extra 
                endif
              end do
              call make_river_flux_output(macro_height(n), &
                dtm*macro_steps,ARCD1, c1dimnumc,numc,i,j)
              end if
#endif
!             Changes of total and layer height due to river inflow:
              hn(i,j,1:kmax) = hn(i,j,1:kmax)/(H(i,j)+ssen(i,j)) &
                             *(H(i,j)+ssen(i,j)+macro_height(n))
              ssen(i,j) = ssen(i,j)+macro_height(n)
              macro_height(n) = _ZERO_
            end if
#endif
         end if
       end do
       if (do_3d) macro_steps=0
      case default
         FATAL 'Not valid rivers_method specified'
         stop 'init_rivers'
   end select

#ifdef DEBUG
   write(debug,*) 'Leaving do_rivers()'
   write(debug,*)
#endif
   return
   end subroutine do_rivers
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  clean_rivers
!
! !INTERFACE:
   subroutine clean_rivers
!
! !DESCRIPTION:
!
! This routine closes the river handling by writing the integrated
! river run-off for each river to standard output.
!
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !LOCAL VARIABLES:
   integer                   :: i,j,n
   REALTYPE                  :: tot=_ZERO_
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef DEBUG
   integer, save :: Ncall = 0
   Ncall = Ncall+1
   write(debug,*) 'clean_rivers() # ',Ncall
#endif

   select case (river_method)
      case(0)
      case(1,2)
         do n=1,nriver
            if(ok(n) .gt. 0) then
               i = ir(n); j = jr(n)
               LEVEL2 trim(river_name(n)),':  ' ,irr(n)/1.e6, '10^6 m3'
               tot = tot+irr(n)
            end if
         end do
      case default
         FATAL 'Not valid rivers_method specified'
         stop 'init_rivers'
   end select

#ifdef DEBUG
   write(debug,*) 'Leaving do_rivers()'
   write(debug,*)
#endif
   return
   end subroutine clean_rivers
!EOC

!-----------------------------------------------------------------------

   end module rivers

!-----------------------------------------------------------------------
! Copyright (C) 2001 - Karsten Bolding and Hans Burchard               !
!-----------------------------------------------------------------------

