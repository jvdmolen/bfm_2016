#include "DEBUG.h"
#include "INCLUDE.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !MODULE: CO2System
!
! DESCRIPTION
!   !

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  module CO2System

!#ifdef INCLUDE_PELCO2
!
! !USES:
   use global_mem, ONLY:RLEN,LOGUNIT,DONE,ZERO_KELVIN


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  real(RLEN),public  ::  K0Ac !  carbonate equilibrium I : [Co2]=k0Ac PCO2
  real(RLEN),public  ::  K1Ac !  carbonate equilibrium II : [CO2]*K1=[H+][HCO3-]
  real(RLEN),public  ::  K2Ac !  carbopnate equilibrium III  [HCO3-].K2=[H+].[CO3--]
  real(RLEN),public  ::  KwAc !
  real(RLEN),public  ::  KbAc ! constant for Boron equilibrium
  real(RLEN),public  ::  KsAc ! constant for sulphate equilibrium
  real(RLEN),public  ::  KfAc ! constants for fluor equilibirum
  real(RLEN),public  ::  KpAc(3) ! constants for phosphate equilibirum
  real(RLEN),public  ::  KsiAc ! constants for silica eqilbrium
  real(RLEN),public  ::  ldic ! 
  real(RLEN),public  ::  lpCO2 ! 
  real(RLEN),public  ::  sN1p !
  real(RLEN),public  ::  sN5s !
  real(RLEN),public  :: scl
  real(RLEN),public  :: ta
  integer,public     :: way

  integer,parameter   :: ISTATIC=1
  integer,parameter   :: DYNAMIC=2
  real(RLEN),parameter:: HplusBASIS=2225.0
  real(RLEN),parameter:: DICBasis=25212.0

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !  functions 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  contains

!EOP
!-------------------------------------------------------------------------!
!BOP


!DESCRIPTION
! code was based on the code by Dickson in Version 2 of ''Handbook &
! of Methods
! for the Analysis of the Various Parameters of the Carbon Dioxide &
! System
! in Seawater'', DOE, 1994 (SOP No. 3, p25-26).

! Derive simple terms used more than once

! !subroutine: CalcCO2System
!
! DESCRIPTION
!   !
! INTERFACE

  function CalcCO2System(mode,salt,temp,rho,n1p,n5s,Ac,CO2,HCO3,CO3,CAc,pH_in, &
                                            DIC_in,pCO2_in,DIC_out,pCO2_out,pH_out)
!
! USES:
  use mem_globalfun,   ONLY: rtsafe
! !INPUT PARAMETERS:
  IMPLICIT NONE
  integer                          :: mode
  real(RLEN),intent(IN)            :: salt
  real(RLEN),intent(IN)            :: temp
  real(RLEN),intent(IN)            :: rho 
  real(RLEN),intent(IN)            :: n1p
  real(RLEN),intent(IN)            :: n5s
  real(RLEN),intent(INOUT)         :: Ac
  real(RLEN),intent(IN),optional   :: pH_in
  real(RLEN),intent(IN),optional   :: DIC_in
  real(RLEN),intent(IN),optional   :: pCO2_in
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
  real(RLEN),intent(OUT)            :: CO2
  real(RLEN),intent(OUT)            :: HCO3
  real(RLEN),intent(OUT)            :: CO3
  real(RLEN),intent(OUT)            :: Cac
  real(RLEN),intent(OUT),optional   :: DIC_out
  real(RLEN),intent(OUT),optional   :: pCO2_out
  real(RLEN),intent(OUT),optional   :: pH_out
  integer                           :: CalcCO2System


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer     :: i
  real(RLEN)  :: hp
  real(RLEN)  :: hphp
  real(RLEN)  :: sqrts
  real(RLEN)  :: s15
  real(RLEN)  :: sis
  real(RLEN)  :: sqrtis
  real(RLEN)  :: tk
  real(RLEN)  :: dlogtk
  real(RLEN)  :: tmp1
  real(RLEN)  :: tmp2
  real(RLEN)  :: mkg_to_mmm3  ! mol/kg ---> mmol/m3 ( density to volume measure)
  real(RLEN)  :: d
  real(RLEN)  :: lpH
  logical     :: small_interval
  integer     :: error
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! The following global functions are used:rtsafe
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    tk  =   -ZERO_KELVIN+ temp
    dlogtk  =   log(tk)

    sis  =   19.924D+00* salt/( 1000.D+00- 1.005D+00* salt)
    sqrtis  =   sqrt(  sis)

    sqrts  =   sqrt(  salt)
    s15  =   salt* sqrts

    mkg_to_mmm3   = rho *1000.0;

    if (present(pH_in).and.present(DIC_in)) then
       way =3 
       ldic  = DIC_in/ mkg_to_mmm3
    elseif (present(DIC_in)) then
       !pH is estimated on basis of  Total Alkailinity and  DIC
       ldic  = DIC_in/ mkg_to_mmm3
       way =1 
    elseif (present(pCO2_in)) then
       !pH is estimated on basis of  Total Alkalinity and  DIC
       lpCO2  = pCO2_in
       way =2 
    endif


    !Ac: 1 mg-eq/l = 0,5 mmol/l CaCO. 3. = 50 mg/l CaCO. 3 
    ! mmol /m3 ---> mol/kg;
    ta    = Ac/mkg_to_mmm3    
    sN1p  = n1p/ mkg_to_mmm3
    sN5s  = n5s/ mkg_to_mmm3


    !---------------------equilibrium constants section--------------

    !   K0 from Weiss 1974
    !   unit: mol.kg-1.atm-1 okay!


    K0Ac = CalcK0Ac(salt,temp)

    !	K1Ac = [H][HCO3]/[H2CO3]
    !	K2Ac = [H][CO3]/[HCO3]

    !	Constants according to Mehrbach after Dickson and Millero 1987

    K1Ac = exp( - 2307.1266D+00/ tk+ 2.83655D+00- &
      1.5529413D+00* dlogtk+(- 4.0484D+00/ tk- 0.20760841D+00)* sqrts+ &
      0.08468345D+00* salt- 0.00654208D+00* s15+ log( DONE- 0.001005D+00* &
      salt))

    K2Ac = exp( - 3351.6106D+00/ tk- 9.226508D+00- &
      0.2005743D+00* dlogtk+(- 23.9722D+00/ tk- 0.10690177D+00)* sqrts+ &
      0.1130822D+00* salt- 0.00846934D+00* s15+ log( DONE- 0.001005D+00* &
      salt))

    !---------------------------OCMIP old----------------------------------
    ! Millero p.664 (1995) using Mehrbach et al. data on seawater scale
    !
    !    K1Ac=pow(10.0, (-1*(3670.7/ tk - 62.008 + 9.7944*dlogtk -
    !   &		salt * ( 0.0118 + 0.000116*salt))) );
    !
    !	K2Ac=pow(10.0, (-1.0*(1394.7/ tk + 4.777
    !                  - salt*( 0.0184 + 0.000118*salt))));
    !---------------------------OCMIP old &
    ! end----------------------------------



    error=0
    if ( mode ==  DYNAMIC ) then

        ! KbAc = [H][BO2]/[HBO2]
        ! Millero p.669 (1995) using data from Dickson (1990)

        KbAc = exp( (- 8966.90D+00+ sqrts*(- 2890.53D+00+ &
              sqrts*(- 77.942D+00+ sqrts*( 1.728D+00- 0.0996D+00* sqrts))))/ &
              tk+( 148.0248D+00+ 137.1942D+00* sqrts+ 1.62142D+00* salt)+(- &
              24.4344D+00- 25.085D+00* sqrts- 0.2474D+00* salt)* dlogtk+ &
              0.053105D+00* sqrts* tk)

        ! KpAc[1] = [H][H2PO4]/[H3PO4]

        ! DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)

        KpAc(1) = exp( - 4576.752D+00/ tk+ 115.525D+00- &
              18.453D+00* dlogtk+(- 106.736D+00/ tk+ 0.69171D+00)* sqrts+(- &
              0.65643D+00/ tk- 0.01844D+00)* salt)
            ! KpAc[2] = [H][HPO4]/[H2PO4]

            ! DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))

        KpAc(2) = exp( - 8814.715D+00/ tk+ 172.0883D+00- &
              27.927D+00* dlogtk+(- 160.340D+00/ tk+ 1.3566D+00)* sqrts+( &
              0.37335D+00/ tk- 0.05778D+00)* salt)

        ! KpAc[3] = [H][PO4]/[HPO4]

            ! DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)

        KpAc(3) = exp( - 3070.75D+00/ tk- 18.141D+00+( &
              17.27039D+00/ tk+ 2.81197D+00)* sqrts+(- 44.99486D+00/ tk- &
              0.09984D+00)* salt)

        !------------------------------------------------------------------------
        ! KsiAc = [H][SiO(OH)3]/[Si(OH)4]
        ! Millero p.671 (1995) using data from Yao and Millero (1995)

        KsiAc = exp( - 8904.2D+00/ tk+ 117.385D+00- &
           19.334D+00* dlogtk+ sqrtis*(- 458.79D+00/ tk+ 3.5913D+00+ sqrtis*( &
           188.74D+00/ tk- 1.5998D+00+ sis*(- 12.1652D+00/ tk+ 0.07871D+00)))+ &
              log( DONE- 0.001005D+00* salt))

        !----------------------------------------------------------------------
        ! KwAc = [H][OH]

        ! Millero p.670 (1995) using composite data
        !
        KwAc = exp( - 13847.26D+00/ tk+ 148.9652D+00- &
              23.6521D+00* dlogtk+( 118.67D+00/ tk- 5.977D+00+ 1.0495D+00* &
              dlogtk)* sqrts- 0.01615D+00* salt)

        !----------------------------------------------------------------------
        ! KsAc = [H][SO4]/[HSO4]

        ! Dickson (1990, J. chem. Thermodynamics 22, 113)

        KsAc = exp( - 4276.1D+00/ tk+ 141.328D+00- &
           23.093D+00* dlogtk+ sqrtis*(- 13856.0D+00/ tk+ 324.57D+00- &
           47.986D+00* dlogtk+ sqrtis*( 35474.0D+00/ tk- 771.54D+00+ &
           114.723D+00* dlogtk+ sqrtis*(- 2698.0D+00/ tk+  &
             sqrtis*( 1776.0D+00/ tk))))+ log( DONE- 0.001005D+00* salt))

        !-----------------------------------------------------------------------
        ! KfAc = [H][F]/[HF]

        ! Dickson and Riley (1979) -- change pH scale to total

        scl  =   salt/ 1.80655D+00
        KfAc = exp( 1590.2D+00/ tk- 12.641D+00+ &
              1.525D+00* sqrtis+ log( DONE- 0.001005D+00* salt)+ &
              log( DONE+( 0.1400D+00/ 96.062D+00)*( scl)/ KsAc))

        !-----------------equilibrium end--------------------------------------


        !-----Calculation of ta as function of pH------------------------------

        ! This routine expresses TA as a function of DIC, Ac and constants.
        ! It also calculates the derivative of this function with respect to
        ! Ac. It is used in the iterative solution for Ac. In the call
        ! ''x'' is the input value for Ac, ''fn'' is the calculated value &
        ! for TA
        ! and ''df'' is the value for dTA/dAc

        if ( way ==3 ) then
            ta=0.0;
            call CalcHplus(pH_in,Ac,lpH);
            Ac=Ac*mkg_to_mmm3
        else

          error = rtsafe( CalcHplus, 2.0D+00, 11.0D+00, 1.0D-6, lpH)

          if ( error == 0 ) then
            pH_out=lph
            hp = 10.0D+00**(-pH_out);
            hphp  =   hp* hp

            select case (way)
              case (1)
               CO2 = ldic* hphp/( hphp+ K1Ac* hp+ K1Ac* K2Ac)
               pCO2_out  =   CO2/ K0Ac
              case (2)
               CO2= lpCO2 *K0Ac
               ldic=CO2*(hphp+K1AC*hp +K1Ac*K2Ac)/hphp
               DIC_out = ldic  * mkg_to_mmm3
              end select

              HCO3  =   K1Ac* CO2/ hp
              CO3  =   K2Ac* HCO3/ hp

              ! transformation from mol/kg -----> mmmold/m3:
              CAc   =  mkg_to_mmm3* ( HCO3 +2*CO3 + 10**(-14)/hp)
              CO2  =   mkg_to_mmm3* CO2
              CO3  =   mkg_to_mmm3* CO3
              HCO3  =  mkg_to_mmm3* HCO3
          endif 
        endif
     endif


     if ( way.eq.3)  then
!       Calculate concentration of H+ and OH- ions.
        hp = 10**(-pH_in)
        hphp  =   hp* hp

!       Calculate fractions of carbonate ions that depends on pH:
        !pHCO3-
        tmp1 = (K1Ac * hp) / ((hphp) + (K1Ac * hp) + (K1Ac * K2Ac))
        !pCO3-
        tmp2= (K1Ac * K2Ac) / ((hphp) + (K1Ac * hp) + (K1Ac * K2Ac))
        CO2 = ldic* hphp/( hphp+ K1Ac* hp+ K1Ac* K2Ac)
        HCO3  =   K1Ac* CO2/ hp
        CO3  =   K2Ac* HCO3/ hp

!       Calculate Carbonate Alkalinity
        ta = (tmp1 + 2.0 * tmp2)*ldic + 10**(-14)/hp
        Cac    = ta*mkg_to_mmm3
     elseif ( mode == ISTATIC .or.( way==1.and.error >0) ) then
      if (Ac.gt.2.0*Dic_in) then 
        pH_out=15.0
        CO2=0.000001
        HCO3=0.000001
      else
        tmp1 = - 0.5D+00*( ta-ldic+ 4.0D+00* K2Ac/ K1Ac*( &
          2.0D+00* ta- ldic))/( DONE- 4.0D+00* K2Ac/ K1Ac)

        tmp2 = - K2Ac* (( 2.0D+00* ldic- ta))**(2.0D+00)/( &
          K1Ac*( DONE- 4.0D+00* K2Ac/ K1Ac))

        !pCO2:
        pCO2_out = ( tmp1+ sqrt( tmp1* tmp1- tmp2))/ K0Ac

        !CO2:
        CO2  =   K0Ac* pCO2_out

        !pH:
        pH_out = - log( K1Ac* CO2/(2.0D+00* ldic &
          - ta- 2.0D+00* CO2))/ log( 10.0D+00)

        !CO3:
        CO3  =   AC- DIC_in+ CO2

        !HCO3:
        HCO3  =   DIC_in- CO2- CO3

        !Cac
        CAc = (tmp1 + 2.0 * tmp2)*ldic + 10**(-14)/hp

        !transformation from mol/kg -----> mmmold/m3:
         CO2  =   mkg_to_mmm3* CO2
         CO3  =   mkg_to_mmm3* CO3
        HCO3  =   mkg_to_mmm3* HCO3
         CAc  =    mkg_to_mmm3* Cac
      ! of STATIC
      endif
    endif 

    if ( error > 0 ) then
        CalcCO2System=error;return;
    else
      CalcCO2System=0;
    endif

   end function CalcCO2System
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: CalcHplus
!
! DESCRIPTION
!   !

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  SUBROUTINE CalcHplus(lx,fn,df)
!
! !USES:
  ! The following Pelagic-states are used (NOT in fluxes): N1p, N5s, O3c
  ! The following global scalar vars are used: BoxNumber,BoxNumberXY
  ! The following Pelagic 1-d global boxvars are used: ERHO, Ac, ESW, K1Ac, &
  ! K2Ac, KsAc, KbAc, KwAc, KsiAc, KfAc
  ! The following Pelagic 2-d global boxvars  are used: KpAc
  ! The following global constants are used: RLEN

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Modules (use of ONLY is strongly encouraged!)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

! use global_mem, ONLY:RLEN,DONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

! !INPUT:
!
  real(RLEN),intent(IN)  :: lx
! !OUTPUT:
  real(RLEN),intent(OUT)  :: fn
  real(RLEN),intent(OUT)  :: df

!  
!
!EOP
!-------------------------------------------------------------------------!
!BOC
!
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  real(RLEN)  :: bt
  real(RLEN)  :: ft
  real(RLEN)  :: st
  real(RLEN)  :: k12
  real(RLEN)  :: k12p
  real(RLEN)  :: k123p
  real(RLEN)  :: c
  real(RLEN)  :: a
  real(RLEN)  :: a2
  real(RLEN)  :: da
  real(RLEN)  :: b
  real(RLEN)  :: b2
  real(RLEN)  :: db
  real(RLEN)  :: lfn
  real(RLEN)  :: ldf
  real(RLEN)  :: x
  real(RLEN)  :: x2
  real(RLEN)  :: x3
  real(RLEN)  :: xlog10

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! external functions
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  xlog10  =   2.30258509299D+00

  !--------------------- inputconcentrations-------------------------
  ! Calculate concentrations for borate, sulfate, and fluoride, phosphate
  ! and silicate
  ! Uppstrom (1974);
  bt  =   0.000232D+00* scl/ 10.811D+00
  ! 121 Morris & Riley (1966)
  st  =   0.14D+00* scl/ 96.062D+00
  ! Riley (1965)
  ft  =   0.000067D+00* scl/ 18.9984D+00
  !------------------------input conc. end------------------------
  ! Derive again simple terms used more than once


  x  =   (10.0D+00)**(- lx)
  x2  =   x* x
  x3  =   x2* x
  k12  =   K1Ac* K2Ac
  k12p  =   KpAc(1)* KpAc(2)
  k123p  =   k12p* KpAc(3)
  c  =   DONE+ st/ KsAc
  a  =   x3+ KpAc(1)* x2+ k12p* x+ k123p
  a2  =   a* a
  da  =   3.0D+00* x2+ 2.0D+00* KpAc(1)* x+ k12p
  b  =   x2+ K1Ac* x+ k12
  b2  =   b* b
  db  =   2.0D+00* x+ K1Ac


  ! This routine expresses TA as a function of DIC, htotal and constants.
  ! It also calculates the derivative of this function with respect to
  ! htotal. It is used in the iterative solution for htotal. In the call
  ! ''x'' is the input value for htotal, ''fn'' is the calculated value for TA
  ! and ''df'' is the value for dTA/dhtotal



  select case ( way ) 
    case (1,3)
      !
      !	fn = hco3+co3
      !	dfn = dhco3/dx+co3/dx
      !
       lfn = K1Ac* x* ldic/ b+ 2.0D+00* ldic* k12/ b+ bt/( &
          DONE+ x/ KbAc)
       ldf = (( K1Ac* ldic* b)- K1Ac* x* ldic* &
          db)/ b2- 2.0D+00* ldic* k12* db/ b2- bt/ KbAc/ &
          (DONE+ x/ KbAc)**(2.0D+00)
    case (2)
!          CO2= lpCO2 *K0Ac
!          ldic=CO2*(hphp+K1AC*hp +K1Ac*K2Ac)/hphp
       lfn=  k0Ac * ( lpCO2  + K1Ac * lpCO2/x  +    k12 * lpCO2/x2);
       ldf= -k0Ac * ( K1Ac * lpCO2/x2 + 2 *k12 * lpCO2/x3);
  end select
  !
  !	fn = fn+borate+oh+hpo4+2*po4+silicate+hfree+hso4+hf+h3po4-ta
  !	dfn = dfn/dx+dborate/dx+doh/dx+dhpo4/dx+2*dpo4/dx
  !                                +dsilicate/dx+dhfree/dx+dhso4/dx+dhf/dx+dh3po4/dx-ta
  !
  lfn= lfn + KwAc/ x+ sN1p* k12p* x/ a+ 2.0D+00* &
    sN1p* k123p/ a+ sN5s/( DONE+ x/ KsiAc)- x/ c- st/( DONE+ &
    KsAc/ x/ c)- ft/( DONE+ KfAc/ x)- sN1p* x3/ a- ta
  !
  !	df = dfn/dx
  !
  ldf=ldf - KwAc/ x2+( sN1p* k12p*( &
    a- x* da))/ a2- 2.0D+00* sN1p* k123p* da/ a2- sN5s/ &
    KsiAc/ (DONE+ x/ KsiAc)**(2.0D+00)- DONE/ c+ st* &
    (DONE+ KsAc/ x/ c)**(- 2.0D+00)*( KsAc/ c/ x2)+ ft* &
    (DONE+ KfAc/ x)**(- 2.0D+00)* KfAc/ x2- sN1p* x2*( &
    3.0D+00* a- x* da)/ a2

  ! give back amount mol/kg

  fn  =   lfn
  df  =   ldf*(- xlog10* x)

  ! ---pH results from iteration here required!---------------

  end subroutine CalcHplus
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: CalcK0Ac
!
! DESCRIPTION
!   !
    !   K0 from Weiss 1974
    !   unit: mol.kg-1.atm-1 okay!
!
! !INTERFACE
  function CalcK0Ac(salt,temp)
!
! !USES:
! !INPUT:
!
  real(RLEN),intent(IN)  :: salt
  real(RLEN),intent(IN)  :: temp
! !OUTPUT:
  real(RLEN) :: CalcK0Ac
!BOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!EOP

! !LOCAL
  real(RLEN) tk100

  tk100  =   (temp-ZERO_KELVIN)/ 100.0D+00

  CalcK0Ac = exp( 93.4517D+00/ tk100- 60.2409D+00+ 23.3585D+00* &
      log( tk100)+ salt*( 0.023517D+00+ tk100*(- 0.023656D+00+ 0.0047036D+00* tk100)))

  end function CalcK0Ac
!EOC
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!#endif

end module CO2System
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
