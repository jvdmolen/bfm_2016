!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL
!	   BFM - Biogeochemical Flux Model version 2.50-g
!
! MODULE
!   ModuleMem
!
! FILE
!   ModuleMem
!
! DESCRIPTION
!   Definition of Global Shared Memory
!  
!   This module contains all the structural definitions of the BFM
!   and sets up the memory layout.
!   It is automatically generated from the prototype file 
!   BFM/proto/ModuleMem.proto by including the information from 
!   BFM/General/GlobalDefsBFM.model
!   Do not directly edit this code because changes will be lost at
!   any new compilation.
!
! AUTHORS
!   Piet Ruardij & Marcello Vichi
!
! CHANGE_LOG
!   ---
!
! COPYING
!   
!   Copyright (C) 2006 P. Ruardij, the BFM team
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
!
!
!
!! IMPORTANT NOTE:
!! Do not change the lines starting with two comment characters "!" 
!! These lines are used by the parser to generate the final module file

!

#include"cppdefs.h"

#include "DEBUG.h"

      module mem
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Modules can optionally use (import) other modules
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        USE BFM_ERROR_MSG, ONLY: BFM_ERROR
        use global_mem, only:RLEN, ZERO
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        implicit none
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Default all is private
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        private
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! next variable can be used to track errors.
        ! By modifying the value in bio.F90 and introduction of output statements
        ! behind an if-statement on can more easy track error by llimiting
        ! the output only for certain cases........ Example:
        !  if ( track_error ===....) write( LOGUNIT,*)............ 
        integer,public          :: track_error=0

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! State variables Info
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! Definition of arrays which will hold all state variables and other
      ! global variables  used for exchange between submodels and/or output
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      real(RLEN),public,pointer,dimension(:,:) :: D3STATE
      real(RLEN),public,pointer,dimension(:,:,:) :: D3SOURCE
      real(RLEN),public,pointer,dimension(:,:,:) :: D3SINK
      integer,public,pointer,dimension(:) :: D3STATETYPE

      real(RLEN),public,pointer,dimension(:,:) :: D2STATE
      real(RLEN),public,pointer,dimension(:,:,:) :: D2SOURCE
      real(RLEN),public,pointer,dimension(:,:,:) :: D2SINK
      integer,public,pointer,dimension(:) :: D2STATETYPE

     real(RLEN),public,pointer,dimension(:,:) :: PELSURFACE

     real(RLEN),public,pointer,dimension(:,:) :: PELBOTTOM

     integer,public,pointer,dimension(:) :: iiPELSINKREF

      real(RLEN),public,pointer,dimension(:,:) :: D3DIAGNOS

      real(RLEN),public,pointer,dimension(:,:) :: D2DIAGNOS

      real(RLEN),public,pointer,dimension(:,:) :: D3DIAGNOS_PRF


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! GLOBAL system CONSTANTS
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      integer,parameter,public  ::NO_D3_BOX_STATES =56

      integer,parameter,public  ::NO_D2_BOX_STATES =101

      integer,parameter,public  ::NO_D3_BOX_DIAGNOSS =206

      integer,parameter,public  ::NO_D2_BOX_DIAGNOSS =469

      integer,parameter,public  ::NO_D3_BOX_DIAGNOSS_PRF=13

      integer,parameter,public  ::NO_D3_BOX_FLUX =10

      integer,parameter,public  ::NO_D2_BOX_FLUX =36

      integer,public  ::NO_BOXES
      integer,public  ::NO_BOXES_X
      integer,public  ::NO_BOXES_Y
      integer,public  ::NO_BOXES_Z
      integer,public  ::NO_STATES
      integer,public  ::NO_BOXES_XY
      integer,public  ::NO_BOXES_PRF=40

      integer,parameter,public  ::iiPel= 0
      integer,parameter,public  ::iiBen= 1000
      integer,parameter,public  ::iiReset= -1000
      integer,parameter,public  ::iiAdd= -1001
      integer,parameter,public  ::iiTotal= 0 
      integer,parameter,public  ::iiConsumption= -1 
      integer,parameter,public  ::iiProduction= 1 


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      !! GLOBAL definition of Pelagic (D3) state variables
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    integer,parameter,public :: ppR9x=1, ppO2o=2, ppN1p=3, ppN3n=4, ppN4n=5,&
     & ppN5s=6, ppN6r=7, ppB1c=8, ppB1n=9, ppB1p=10, ppBac=11,&
     & ppP1c=12, ppP1n=13, ppP1p=14, ppP1l=15, ppP1s=16, ppP2c=17,&
     & ppP2n=18, ppP2p=19, ppP2l=20, ppP3c=21, ppP3n=22, ppP3p=23,&
     & ppP3l=24, ppP4c=25, ppP4n=26, ppP4p=27, ppP4l=28, ppP5c=29,&
     & ppP5n=30, ppP5p=31, ppP5l=32, ppP5s=33, ppP6c=34, ppP6n=35,&
     & ppP6p=36, ppP6l=37, ppPcc=38, ppR1c=39, ppR1n=40, ppR1p=41,&
     & ppR2c=42, ppR2n=43, ppR3c=44, ppR6c=45, ppR6n=46, ppR6p=47,&
     & ppR6s=48, ppRZc=49, ppO3c=50, ppO3h=51, ppZ3c=52, ppZ4c=53, ppZ2c=54,&
     & ppZ5c=55, ppZ6c=56, ppP2s=0, ppP3s=0, ppP4s=0, ppP6s=0, ppR1s=0,&
     & ppR2p=0, ppR2s=0, ppR3n=0, ppR3p=0, ppR3s=0, ppZ3n=0, ppZ3p=0,&
     & ppZ4n=0, ppZ4p=0, ppZ2n=0, ppZ2p=0, ppZ5n=0, ppZ5p=0, ppZ6n=0, ppZ6p=0


    real(RLEN),public,dimension(:),pointer :: R9x, O2o, N1p, N3n, N4n, N5s,&
     & N6r, B1c, B1n, B1p, Bac, P1c, P1n, P1p, P1l, P1s, P2c, P2n, P2p,&
     & P2l, P3c, P3n, P3p, P3l, P4c, P4n, P4p, P4l, P5c, P5n, P5p, P5l,&
     & P5s, P6c, P6n, P6p, P6l, Pcc, R1c, R1n, R1p, R2c, R2n, R3c, R6c, R6n,&
     & R6p, R6s, RZc, O3c, O3h, Z3c, Z4c, Z2c, Z5c, Z6c


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      !! GLOBAL definition of Benthic (D2) state variables
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    integer,parameter,public :: ppBP1c=1, ppBP1n=2, ppBP1p=3, ppBP1l=4,&
     & ppBP1s=5, ppY1c=6, ppY1n=7, ppY1p=8, ppY2c=9, ppY2n=10,&
     & ppY2p=11, ppY4c=12, ppY4n=13, ppY4p=14, ppY5c=15,&
     & ppY5n=16, ppY5p=17, ppYy3c=18, ppYy3n=19, ppYy3p=20,&
     & ppY3c=21, ppY3n=22, ppY3p=23, ppYs3c=24, ppQ6c=25,&
     & ppQ6n=26, ppQ6p=27, ppQ6s=28, ppQ9x=29, ppQp9x=30,&
     & ppQSx=31, ppQun=32, ppQ1un=33, ppQ2un=34, ppQ2c=35,&
     & ppQ2n=36, ppQ12c=37, ppQ12n=38, ppQ1c=39, ppQ1n=40,&
     & ppQ1p=41, ppQ11c=42, ppQ11n=43, ppQ11p=44, ppQ21c=45,&
     & ppQ21n=46, ppQ21p=47, ppH1c=48, ppH1n=49, ppH1p=50,&
     & ppH2c=51, ppH2n=52, ppH2p=53, ppH3c=54, ppH3n=55,&
     & ppH3p=56, ppHNc=57, ppHNn=58, ppHNp=59, ppHac=60,&
     & ppKp1p=61, ppKp3n=62, ppKp4n=63, ppKn4n=64, ppKp5s=65,&
     & ppQpun=66, ppK1p=67, ppK11p=68, ppK21p=69, ppK4n=70, ppK14n=71,&
     & ppK24n=72, ppK5s=73, ppK15s=74, ppK6r=75, ppK16r=76, ppK26r=77,&
     & ppK3n=78, ppK13n=79, ppK23n=80, ppG2o=81, ppG4n=82, ppDfm=83,&
     & ppDcm=84, ppDlm=85, ppD1m=86, ppD2m=87, ppD6m=88, ppD7m=89,&
     & ppD8m=90, ppD9m=91, ppDSm=92, ppDH2m=93, ppDH3m=94, ppirri_bio=95,&
     & ppG3c=96, ppG3h=97, ppG13c=98, ppG13h=99, ppG23c=100, ppG23h=101,&
     & ppMcac=0, ppMcan=0, ppMcap=0, ppMcal=0, ppMsc=0


    real(RLEN),public,dimension(:),pointer :: BP1c, BP1n, BP1p, BP1l, BP1s,&
     & Y1c, Y1n, Y1p, Y2c, Y2n, Y2p, Y4c, Y4n, Y4p, Y5c, Y5n, Y5p,&
     & Yy3c, Yy3n, Yy3p, Y3c, Y3n, Y3p, Ys3c, Q6c, Q6n, Q6p, Q6s, Q9x,&
     & Qp9x, QSx, Qun, Q1un, Q2un, Q2c, Q2n, Q12c, Q12n, Q1c, Q1n, Q1p,&
     & Q11c, Q11n, Q11p, Q21c, Q21n, Q21p, H1c, H1n, H1p, H2c, H2n, H2p,&
     & H3c, H3n, H3p, HNc, HNn, HNp, Hac, Kp1p, Kp3n, Kp4n, Kn4n, Kp5s,&
     & Qpun, K1p, K11p, K21p, K4n, K14n, K24n, K5s, K15s, K6r, K16r, K26r,&
     & K3n, K13n, K23n, G2o, G4n, Dfm, Dcm, Dlm, D1m, D2m, D6m, D7m,&
     & D8m, D9m, DSm, DH2m, DH3m, irri_bio, G3c, G3h, G13c, G13h, G23c, G23h



    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ! Definition(s) of alternative Z-axis
    real(RLEN),public,dimension(:),pointer  :: seddepth
    real(RLEN),public,dimension(:),pointer  :: prmdepth
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! Constituent parameters:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    integer,parameter,public :: iiC=1, iiN=2, iiP=3, iiL=4, iiS=5, iiH=6,&
     & iiR=7, iiO=8, iiM=9


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! Group parameters:
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    integer,parameter,public :: iiPhytoPlankton=6, iiP1=1, iiP2=2, iiP3=3,&
     & iiP4=4, iiP5=5, iiP6=6
    integer,parameter,public :: iiPelDetritus=4, iiR1=1, iiR2=2, iiR3=3,&
     & iiR6=4
    integer,parameter,public     :: iiInorganic=1, iiO3=1
    integer,parameter,public :: iiMesoZooPlankton=3, iiZ3=1, iiZ4=2,&
     & iiZ2=3
    integer,parameter,public     :: iiMicroZooPlankton=2, iiZ5=1, iiZ6=2


    integer,parameter,public     :: iiBenPhyto=1, iiBP1=1
    integer,parameter,public :: iiBenOrganisms=4, iiY1=1, iiY2=2, iiY4=3,&
     & iiY5=4
    integer,parameter,public     :: iiSuspensionFeeders=2, iiYy3=1, iiY3=2
    integer,parameter,public     :: iiBenUrea=3, iiQu=1, iiQ1u=2, iiQ2u=3
    integer,parameter,public     :: iiBenCarboHydrates=2, iiQ2=1, iiQ12=2
    integer,parameter,public :: iiBenLabileDetritus=3, iiQ1=1, iiQ11=2,&
     & iiQ21=3
    integer,parameter,public :: iiBenBacteria=4, iiH1=1, iiH2=2, iiH3=3,&
     & iiHN=4
    integer,parameter,public :: iiBenthicPhosphate=3, iiK1=1, iiK11=2,&
     & iiK21=3
    integer,parameter,public :: iiBenthicAmmonium=3, iiK4=1, iiK14=2,&
     & iiK24=3
    integer,parameter,public     :: iiBenthicSilicate=2, iiK5=1, iiK15=2
    integer,parameter,public :: iiBenthicRedEq=3, iiK6=1, iiK16=2,&
     & iiK26=3
    integer,parameter,public :: iiBenthicNitrate=3, iiK3=1, iiK13=2,&
     & iiK23=3
    integer,parameter,public     :: iiBenthicCO2=3, iiG3=1, iiG13=2, iiG23=3


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      !  Global Variables
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      integer,public  :: BoxNumber
      integer,public  :: BoxNumberX
      integer,public  :: BoxNumberY
      integer,public  :: BoxNumberZ
      integer,public  :: BoxNumberXY

    real(RLEN),public                                    :: &
     LocalDelta,& !
     max_change_per_step,& !
     max_rate_per_step,& !
     Wind,& !Wind (m/s)
     SUNQ,& !Daylength in hours (h)
     ThereIsLight      !Forcing for day/night

    integer,public                                    :: &
     InitializeModel      ! Section -decr
     integer,public,dimension(:),allocatable           :: &
      CoupledtoBDc      !determine which Phytoplankton group is exchanged with Benthic Diatoms

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      !! GLOBAL definition of Pelagic (D3) variables which can be outputted in netcdf
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-

!    3d name                                                  description            unit
! ---------- ------------------------------------------------------------ ---------------
!        ETW                                                  Temperature               C
!        ESW                                                     Salinity             psu
!        EIR                                                   Irradiance         mE/m2/s
!       ERHO                                                      Density            g/m3
!        ESS                                           Suspended sediment           mg/m3
!      cmO2o                      Temperature-dependent oxygen saturation      mmol O2/m3
!       XO2o                                      Net O2 conc. = O2 - H2S      mmol O2/m3
!     eO2mO2                                   Relative Oxygen saturation               %
!       Chla                                                 Chlorophylla       mg Chl/m3
!     ws_out                                    settling velocity of Silt             m/s
!    flPTN6r                                     anaerobic mineralization   mmol S--/m3/d
!     sN4N3n                               first order nitrification rate              /d
!    flN4N3n                                           nitrification flux     mmol N/m3/d
!    flN3N4n                                   reverse-nitrification flux     mmol N/m3/d
!    flN3O4n                                  denitrification flux (sink)     mmol N/m3/d
!    flBaZTc                                        grazing on Nitrifiers     mmol C/m3/d
!      pNaNt                                 Archea part of PelNitrifiers               -
!    flR1O3c                             CO2 production by breakdown urea         mg C/m3
!    flR1N4n                             NH4 production by breakdown urea       mmol N/m3
!    flB1N4n                                      NH4 release by bacteria       mmol N/m3
!    flPIN4n                                 NH4 release by phytoplankton       mmol N/m3
!      rmB1c                                          bacterial mortality         mg C/m3
!    flR3R2c                                                   flux R3>R2       mg C/m3/d
! p_xfree_R2                                   TEP outside MacroAggregate               -
! fr_lim_BN_n           max.fract.of uptake limit for Nitrifiers (basis n)               -
! fr_lim_B1_n             max.fract.of uptake limit for Bacteria (basis n)               -
! fr_lim_B1_p             max.fract.of uptake limit for Bacteria (basis p)               -
! fl_xgrazing_B1c                                     grazing rate on Bacteria       mg C/m2/d
!        rml                                              Lysis mortality       mg C/m3/d
!      qpR6c                                             Quotum p/c in R6       mol P/g C
!      qnR6c                                             Quotum n/c in R6       mol N/g C
!      qsR6c                                            Quotum si/c in R6      mol Si/g C
!      qpB1c                                             Quotum p/c in B1       mol P/g C
!      qnB1c                                             Quotum n/c in B1       mol N/g C
!    pMIupZ4           part of food orginating from small food web for Z4               -
!    rnetPTc                            netto production of Phytoplankton       mg C/m3/d
!     flnDIn                             net uptake of dissolved nitrogen     mmol N/m3/d
!     flnDIp                             net uptake of dissolved nitrogen     mmol P/m3/d
!        Nun                                                         urea       mmol N/m3
!        Nup                                        rich labile Phosphate       mmol P/m3
!     sediR2                                       TEP sedimentation rate             m/d
!     sediR6                          Detritus sedimentation rate (C,N,P)             m/d
!    sediR6s                         Diatom frustiles  sedimentation rate             m/d
!     sediRZ                              FaecelPellet sedimentation rate             m/d
!     sediR9                                      Silt sedimentation rate             m/d
!     sediB1                                  Bacteria Sedimentation rate             m/d
!        POC                                          Particulate  Carbon          mgC/m3
!        PON                                        Particulate  Nitrogen        mmolN/m3
!        POP                                       Particulate  Phosphate        mmolP/m3
!        PTi                             index of most common Phyto group               -
!    limnuti                                    nutrient limitation index               -
!     O2o_vr                                           variance of Oxygen  (mmol O2/m3)^2
!     N1p_vr                                        variance of Phosphate   (mmol P/m3)^2
!     N3n_vr                                        variance of Phosphate   (mmol N/m3)^2
!    Chla_vr                                       variance of Chlorofyll     (mg C/m3)^2
!     R6c_vr                                         variance of Detritus     (mg C/m3)^2
!  xSizeMA_m                                          Size MacroAggregate               m
!  xSizeMA_d                               density-content MacroAggregate         afdw/m3
!      qR2P1                  proportion TEP to diatoms in MacroAggregate               -
!       xEPS                                 Total extinction coefficient             1/m
!     xEPS_0                       cdom+background extinction coefficient             1/m
!   xEPS_ESS                                  silt extinction coefficient             1/m
!   xEPS_Chl                                   Chl extinction coefficient             1/m
!       pCO2                                         partial CO2 pressure               -
!        CO2                                            CO2 concentration          mol/kg
!       HCO3                                          HCO3- concentration          mol/kg
!        CO3                                          CO3-- concentration          mol/kg
!         pH                                                           pH               -
!         Ac                                             Total alkalinity          mol/kg
!        CAc                                            Carbon alkalinity          mol/kg
!        DIC                                   Dissolved inorganic carbon          mol/kg

! iNPI(iiP1)                                     Nutrient depletion state               -
! iNPI(iiP2)                                     Nutrient depletion state               -
! iNPI(iiP3)                                     Nutrient depletion state               -
! iNPI(iiP4)                                     Nutrient depletion state               -
! iNPI(iiP5)                                     Nutrient depletion state               -
! iNPI(iiP6)                                     Nutrient depletion state               -
! sugPI(iiP1)                                 Specific Gross Prod. Diatoms     mg C/mg C/d
! sugPI(iiP2)                             Specific Gross Prod. Flagellates     mg C/mg C/d
! sugPI(iiP3)                       Specific Gross Prod. PicoPhytoPlankton     mg C/mg C/d
! sugPI(iiP4)                         Specific Gross Prod. Dinoflagellates     mg C/mg C/d
! sugPI(iiP5)             Specific Gross Prod. Resuspended Benthic Diatoms     mg C/mg C/d
! sugPI(iiP6)                    Specific Gross Prod. Phaeocystis colonies     mg C/mg C/d
! sunPI(iiP1)                                   Specific Net Prod. Diatoms     mg C/mg C/d
! sunPI(iiP2)                               Specific Net Prod. Flagellates     mg C/mg C/d
! sunPI(iiP3)                         Specific Net Prod. PicoPhytoPlankton     mg C/mg C/d
! sunPI(iiP4)                           Specific Net Prod. Dinoflagellates     mg C/mg C/d
! sunPI(iiP5)               Specific Net Prod. Resuspended Benthic Diatoms     mg C/mg C/d
! sunPI(iiP6)                      Specific Net Prod. Phaeocystis colonies     mg C/mg C/d
! sdoPI(iiP1)                                   Specific Mortality Diatoms     mg C/mg C/d
! sdoPI(iiP2)                               Specific Mortality Flagellates     mg C/mg C/d
! sdoPI(iiP3)                         Specific Mortality PicoPhytoPlankton     mg C/mg C/d
! sdoPI(iiP4)                           Specific Mortality Dinoflagellates     mg C/mg C/d
! sdoPI(iiP5)               Specific Mortality Resuspended Benthic Diatoms     mg C/mg C/d
! sdoPI(iiP6)                      Specific Mortality Phaeocystis colonies     mg C/mg C/d
! qpPc(iiP1)                                        Quotum P/C in Diatoms       mol P/g C
! qpPc(iiP2)                                    Quotum P/C in Flagellates       mol P/g C
! qpPc(iiP3)                              Quotum P/C in PicoPhytoPlankton       mol P/g C
! qpPc(iiP4)                                Quotum P/C in Dinoflagellates       mol P/g C
! qpPc(iiP5)                    Quotum P/C in Resuspended Benthic Diatoms       mol P/g C
! qpPc(iiP6)                           Quotum P/C in Phaeocystis colonies       mol P/g C
! qnPc(iiP1)                                        Quotum N/C in Diatoms       mol N/g C
! qnPc(iiP2)                                    Quotum N/C in Flagellates       mol N/g C
! qnPc(iiP3)                              Quotum N/C in PicoPhytoPlankton       mol N/g C
! qnPc(iiP4)                                Quotum N/C in Dinoflagellates       mol N/g C
! qnPc(iiP5)                    Quotum N/C in Resuspended Benthic Diatoms       mol N/g C
! qnPc(iiP6)                           Quotum N/C in Phaeocystis colonies       mol N/g C
! qsPc(iiP1)                                       Quotum Si/C in Diatoms      mol Si/g C
! qsPc(iiP2)                                   Quotum Si/C in Flagellates      mol Si/g C
! qsPc(iiP3)                             Quotum Si/C in PicoPhytoPlankton      mol Si/g C
! qsPc(iiP4)                               Quotum Si/C in Dinoflagellates      mol Si/g C
! qsPc(iiP5)                   Quotum Si/C in Resuspended Benthic Diatoms      mol Si/g C
! qsPc(iiP6)                          Quotum Si/C in Phaeocystis colonies      mol Si/g C
! qlPc(iiP1)                                      Quotum Chl/C in Diatoms       g Chl/g C
! qlPc(iiP2)                                  Quotum Chl/C in Flagellates       g Chl/g C
! qlPc(iiP3)                            Quotum Chl/C in PicoPhytoPlankton       g Chl/g C
! qlPc(iiP4)                              Quotum Chl/C in Dinoflagellates       g Chl/g C
! qlPc(iiP5)                  Quotum Chl/C in Resuspended Benthic Diatoms       g Chl/g C
! qlPc(iiP6)                         Quotum Chl/C in Phaeocystis colonies       g Chl/g C
! qpZc(iiZ3)                       Quotum P/C Carnivorous mesozooplankton       mol P/g C
! qpZc(iiZ4)                        Quotum P/C Omnivorous mesozooplankton       mol P/g C
! qpZc(iiZ2)                                Quotum P/C Filterfeederlarvae       mol P/g C
! qnZc(iiZ3)                       Quotum N/C Carnivorous mesozooplankton       mol N/g C
! qnZc(iiZ4)                        Quotum N/C Omnivorous mesozooplankton       mol N/g C
! qnZc(iiZ2)                                Quotum N/C Filterfeederlarvae       mol N/g C
! qp_mz(iiZ5)                                  Quotum P/C Microzooplankton       mol P/g C
! qp_mz(iiZ6)              Quotum P/C Heterotrophic nanoflagellates (HNAN)       mol P/g C
! qn_mz(iiZ5)                                  Quotum N/C Microzooplankton       mol N/g C
! qn_mz(iiZ6)              Quotum N/C Heterotrophic nanoflagellates (HNAN)       mol N/g C
! flPIR1n(iiP1)                                           Diatoms flux PI>R1     mmol N/m3/d
! flPIR1n(iiP2)                                       Flagellates flux PI>R1     mmol N/m3/d
! flPIR1n(iiP3)                                 PicoPhytoPlankton flux PI>R1     mmol N/m3/d
! flPIR1n(iiP4)                                   Dinoflagellates flux PI>R1     mmol N/m3/d
! flPIR1n(iiP5)                       Resuspended Benthic Diatoms flux PI>R1     mmol N/m3/d
! flPIR1n(iiP6)                              Phaeocystis colonies flux PI>R1     mmol N/m3/d
! flPIR1p(iiP1)                                           Diatoms flux PI>R1     mmol P/m3/d
! flPIR1p(iiP2)                                       Flagellates flux PI>R1     mmol P/m3/d
! flPIR1p(iiP3)                                 PicoPhytoPlankton flux PI>R1     mmol P/m3/d
! flPIR1p(iiP4)                                   Dinoflagellates flux PI>R1     mmol P/m3/d
! flPIR1p(iiP5)                       Resuspended Benthic Diatoms flux PI>R1     mmol P/m3/d
! flPIR1p(iiP6)                              Phaeocystis colonies flux PI>R1     mmol P/m3/d
! flPIR6n(iiP1)                                           Diatoms flux PI>R6     mmol N/m3/d
! flPIR6n(iiP2)                                       Flagellates flux PI>R6     mmol N/m3/d
! flPIR6n(iiP3)                                 PicoPhytoPlankton flux PI>R6     mmol N/m3/d
! flPIR6n(iiP4)                                   Dinoflagellates flux PI>R6     mmol N/m3/d
! flPIR6n(iiP5)                       Resuspended Benthic Diatoms flux PI>R6     mmol N/m3/d
! flPIR6n(iiP6)                              Phaeocystis colonies flux PI>R6     mmol N/m3/d
! flPIR6p(iiP1)                                           Diatoms flux PI>R6     mmol P/m3/d
! flPIR6p(iiP2)                                       Flagellates flux PI>R6     mmol P/m3/d
! flPIR6p(iiP3)                                 PicoPhytoPlankton flux PI>R6     mmol P/m3/d
! flPIR6p(iiP4)                                   Dinoflagellates flux PI>R6     mmol P/m3/d
! flPIR6p(iiP5)                       Resuspended Benthic Diatoms flux PI>R6     mmol P/m3/d
! flPIR6p(iiP6)                              Phaeocystis colonies flux PI>R6     mmol P/m3/d
! flPIR6s(iiP1)                                           Diatoms flux PI>R6    mmol Si/m3/d
! flPIR6s(iiP2)                                       Flagellates flux PI>R6    mmol Si/m3/d
! flPIR6s(iiP3)                                 PicoPhytoPlankton flux PI>R6    mmol Si/m3/d
! flPIR6s(iiP4)                                   Dinoflagellates flux PI>R6    mmol Si/m3/d
! flPIR6s(iiP5)                       Resuspended Benthic Diatoms flux PI>R6    mmol Si/m3/d
! flPIR6s(iiP6)                              Phaeocystis colonies flux PI>R6    mmol Si/m3/d
! fr_lim_PI_n(iiP1)               max.fract.of uptake limit for Diatoms (basisn)               -
! fr_lim_PI_n(iiP2)           max.fract.of uptake limit for Flagellates (basisn)               -
! fr_lim_PI_n(iiP3)     max.fract.of uptake limit for PicoPhytoPlankton (basisn)               -
! fr_lim_PI_n(iiP4)       max.fract.of uptake limit for Dinoflagellates (basisn)               -
! fr_lim_PI_n(iiP5) max.fract.of uptake limit for Resuspended Benthic Diatoms (basisn)               -
! fr_lim_PI_n(iiP6)  max.fract.of uptake limit for Phaeocystis colonies (basisn)               -
! fr_lim_PI_p(iiP1)               max.fract.of uptake limit for Diatoms (basisp)               -
! fr_lim_PI_p(iiP2)           max.fract.of uptake limit for Flagellates (basisp)               -
! fr_lim_PI_p(iiP3)     max.fract.of uptake limit for PicoPhytoPlankton (basisp)               -
! fr_lim_PI_p(iiP4)       max.fract.of uptake limit for Dinoflagellates (basisp)               -
! fr_lim_PI_p(iiP5) max.fract.of uptake limit for Resuspended Benthic Diatoms (basisp)               -
! fr_lim_PI_p(iiP6)  max.fract.of uptake limit for Phaeocystis colonies (basisp)               -
! fl_xgrazing_PIc(iiP1)                                      grazing rate on Diatoms       mg C/m2/d
! fl_xgrazing_PIc(iiP2)                                  grazing rate on Flagellates       mg C/m2/d
! fl_xgrazing_PIc(iiP3)                            grazing rate on PicoPhytoPlankton       mg C/m2/d
! fl_xgrazing_PIc(iiP4)                              grazing rate on Dinoflagellates       mg C/m2/d
! fl_xgrazing_PIc(iiP5)                  grazing rate on Resuspended Benthic Diatoms       mg C/m2/d
! fl_xgrazing_PIc(iiP6)                         grazing rate on Phaeocystis colonies       mg C/m2/d
! sediPI(iiP1)                                   Diatoms sedimentation rate             m/d
! sediPI(iiP2)                               Flagellates sedimentation rate             m/d
! sediPI(iiP3)                         PicoPhytoPlankton sedimentation rate             m/d
! sediPI(iiP4)                           Dinoflagellates sedimentation rate             m/d
! sediPI(iiP5)               Resuspended Benthic Diatoms sedimentation rate             m/d
! sediPI(iiP6)                      Phaeocystis colonies sedimentation rate             m/d
! sediMiZ(iiZ5)                          Microzooplankton sedimentation rate             m/d
! sediMiZ(iiZ6)      Heterotrophic nanoflagellates (HNAN) sedimentation rate             m/d
! sediMeZ(iiZ3)               Carnivorous mesozooplankton sedimentation rate             m/d
! sediMeZ(iiZ4)                Omnivorous mesozooplankton sedimentation rate             m/d
! sediMeZ(iiZ2)                        Filterfeederlarvae sedimentation rate             m/d
! PI_dw(iiP1)                                ash free dryweight of Diatoms         afdw/m3
! PI_dw(iiP2)                            ash free dryweight of Flagellates         afdw/m3
! PI_dw(iiP3)                      ash free dryweight of PicoPhytoPlankton         afdw/m3
! PI_dw(iiP4)                        ash free dryweight of Dinoflagellates         afdw/m3
! PI_dw(iiP5)            ash free dryweight of Resuspended Benthic Diatoms         afdw/m3
! PI_dw(iiP6)                   ash free dryweight of Phaeocystis colonies         afdw/m3
! eiPI(iiP1)                       Regulating factor for light in Diatoms               -
! eiPI(iiP2)                   Regulating factor for light in Flagellates               -
! eiPI(iiP3)             Regulating factor for light in PicoPhytoPlankton               -
! eiPI(iiP4)               Regulating factor for light in Dinoflagellates               -
! eiPI(iiP5)   Regulating factor for light in Resuspended Benthic Diatoms               -
! eiPI(iiP6)          Regulating factor for light in Phaeocystis colonies               -
! EPLi(iiP1)                                     Optimal light in Diatoms            W/m2
! EPLi(iiP2)                                 Optimal light in Flagellates            W/m2
! EPLi(iiP3)                           Optimal light in PicoPhytoPlankton            W/m2
! EPLi(iiP4)                             Optimal light in Dinoflagellates            W/m2
! EPLi(iiP5)                 Optimal light in Resuspended Benthic Diatoms            W/m2
! EPLi(iiP6)                        Optimal light in Phaeocystis colonies            W/m2

!  3_1d name                                                  description            unit
! ---------- ------------------------------------------------------------ ---------------
!     PrsM1p                                        Phosphate in sediment    mmol P/m3 pw
!      PrM1p                                        Phosphate in sediment    mmol P/m3 pw
!      PrM3n                                          Nitrate in sediment    mmol N/m3 pw
!      PrM4n                                         ammonium in sediment    mmol N/m3 pw
!     PreM5s                                         Silicate in sediment   mmol Si/m3 pw
!      PrM5s                                         Silicate in sediment   mmol Si/m3 pw
!      PrM6r                                        Red.Quiv. in sediment  mmol S--/m3 pw
!      PrQun                                                         Urea    mmol N/m3 pw
!      PrDIC                                              DIC in sediment    mmol C/m3 pw
!       PrAc                                               H+ in sediment   mmol H+/m3 pw
!       PrpH                                               pH in sediment               -
!      PrBDc                                              Benthic Diatoms         mg C/m3
!    PrBChlC                                      CChla ratio in sediment      g Chla/g C

     integer,parameter,public :: ppETW=1, ppESW=2, ppEIR=3,&
      & ppERHO=4, ppESS=5, ppcmO2o=6, ppXO2o=7,&
      & ppeO2mO2=8, ppChla=9, ppws_out=10, ppflPTN6r=11,&
      & ppsN4N3n=12, ppflN4N3n=13, ppflN3N4n=14, ppflN3O4n=15,&
      & ppflBaZTc=16, pppNaNt=17, ppflR1O3c=18, ppflR1N4n=19,&
      & ppflB1N4n=20, ppflPIN4n=21, pprmB1c=22, ppflR3R2c=23,&
      & ppp_xfree_R2=24, ppfr_lim_BN_n=25, ppfr_lim_B1_n=26,&
      & ppfr_lim_B1_p=27, ppfl_xgrazing_B1c=28, pprml=29, ppqpR6c=30,&
      & ppqnR6c=31, ppqsR6c=32, ppqpB1c=33, ppqnB1c=34,&
      & pppMIupZ4=35, pprnetPTc=36, ppflnDIn=37, ppflnDIp=38,&
      & ppNun=39, ppNup=40, ppsediR2=41, ppsediR6=42,&
      & ppsediR6s=43, ppsediRZ=44, ppsediR9=45, ppsediB1=46,&
      & ppPOC=47, ppPON=48, ppPOP=49, ppPTi=50, pplimnuti=51,&
      & ppO2o_vr=52, ppN1p_vr=53, ppN3n_vr=54, ppChla_vr=55, ppR6c_vr=56,&
      & ppxSizeMA_m=57, ppxSizeMA_d=58, ppqR2P1=59, ppxEPS=60, ppxEPS_0=61,&
      & ppxEPS_ESS=62, ppxEPS_Chl=63, pppCO2=64, ppCO2=65, ppHCO3=66,&
      & ppCO3=67, pppH=68, ppAc=69, ppCAc=70, ppDIC=71

     integer,public ::&
      & ppiNPI(iiPhytoPlankton),&
      & ppsugPI(iiPhytoPlankton),&
      & ppsunPI(iiPhytoPlankton),&
      & ppsdoPI(iiPhytoPlankton),&
      & ppqpPc(iiPhytoPlankton),&
      & ppqnPc(iiPhytoPlankton),&
      & ppqsPc(iiPhytoPlankton),&
      & ppqlPc(iiPhytoPlankton),&
      & ppqpZc(iiMesoZooPlankton),&
      & ppqnZc(iiMesoZooPlankton),&
      & ppqp_mz(iiMicroZooPlankton),&
      & ppqn_mz(iiMicroZooPlankton),&
      & ppflPIR1n(iiPhytoPlankton),&
      & ppflPIR1p(iiPhytoPlankton),&
      & ppflPIR6n(iiPhytoPlankton),&
      & ppflPIR6p(iiPhytoPlankton),&
      & ppflPIR6s(iiPhytoPlankton),&
      & ppfr_lim_PI_n(iiPhytoPlankton),&
      & ppfr_lim_PI_p(iiPhytoPlankton),&
      & ppfl_xgrazing_PIc(iiPhytoPlankton),&
      & ppsediPI(iiPhytoPlankton),&
      & ppsediMiZ(iiMicroZooPlankton), ppsediMeZ(iiMesoZooPlankton),&
      & ppPI_dw(iiPhytoPlankton), ppeiPI(iiPhytoPlankton),&
      & ppEPLi(iiPhytoPlankton)

     integer,parameter,public :: ppPrsM1p=1, ppPrM1p=2, ppPrM3n=3,&
      & ppPrM4n=4, ppPreM5s=5, ppPrM5s=6, ppPrM6r=7, ppPrQun=8,&
      & ppPrDIC=9, ppPrAc=10, ppPrpH=11, ppPrBDc=12, ppPrBChlC=13

     real(RLEN),public,dimension(:),pointer :: ETW, ESW, EIR, ERHO, ESS,&
      & cmO2o, XO2o, eO2mO2, Chla, ws_out, flPTN6r, sN4N3n, flN4N3n,&
      & flN3N4n, flN3O4n, flBaZTc, pNaNt, flR1O3c, flR1N4n, flB1N4n,&
      & flPIN4n, rmB1c, flR3R2c, p_xfree_R2, fr_lim_BN_n,&
      & fr_lim_B1_n, fr_lim_B1_p, fl_xgrazing_B1c, rml, qpR6c, qnR6c, qsR6c,&
      & qpB1c, qnB1c, pMIupZ4, rnetPTc, flnDIn, flnDIp, Nun, Nup, sediR2,&
      & sediR6, sediR6s, sediRZ, sediR9, sediB1, POC, PON, POP, PTi,&
      & limnuti, O2o_vr, N1p_vr, N3n_vr, Chla_vr, R6c_vr, xSizeMA_m,&
      & xSizeMA_d, qR2P1, xEPS, xEPS_0, xEPS_ESS, xEPS_Chl, pCO2, CO2, HCO3,&
      & CO3, pH, Ac, CAc, DIC

     real(RLEN),public,dimension(:,:),pointer :: iNPI, sugPI, sunPI, sdoPI,&
      & qpPc, qnPc, qsPc, qlPc, qpZc, qnZc, qp_mz,&
      & qn_mz, flPIR1n, flPIR1p, flPIR6n, flPIR6p, flPIR6s, fr_lim_PI_n,&
      & fr_lim_PI_p, fl_xgrazing_PIc, sediPI, sediMiZ, sediMeZ, PI_dw, eiPI,&
      & EPLi

     real(RLEN),public,dimension(:),pointer :: PrsM1p, PrM1p, PrM3n, PrM4n,&
      & PreM5s, PrM5s, PrM6r, PrQun, PrDIC, PrAc, PrpH, PrBDc, PrBChlC


   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   !!  GLOBAL definition of Benthic (D2) variables which can be outputted in netcdf
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!    2d name                                                  description            unit
! ---------- ------------------------------------------------------------ ---------------
!       EIRr                                             Light at surface               W
!     EUWIND                                                  EasternWind             m/s
!     EVWIND                                                 SouthernWind             m/s
!      dry_z                                  quotient RealDepth/MinDepth               -
!      ETAUB                                                 BottomStress             m/s
! EUCURR_LEVEL                                       WestEastSurfaceCurrent             m/s
! EVCURR_LEVEL                                     SouthNorthSurfaceCurrent             m/s
!      jrESS                                   resuspension rate  of silt         mg/m2/d
!   shiftDlm                                        Rate of change of Dpm             m/d
!   shiftD1m                                        Rate of change of D1m             m/d
!   shiftD2m                                        Rate of change of D2m             m/d
! G2_xavail_o                              flux corrected available oxygen   mmol O2/m2 pw
!     jO2Y2o              Uptake from water column of O2 by deposifeeders    mmol O2/m2/d
!    jcrrBTo                   Full rates- Corrected rates for low oxygen    mmol O2/m2/d
!      rrBTo                               Total benthic oxic respiration    mmol O2/m2/d
!      reK6o                             Anoxic respiration in oxic layer    mmol O2/m2/d
!      rrDTo                 Benthic anaerobic respiration in dentr.layer    mmol O2/m2/d
!      rrATo                          Total benthic anaerobic respiration    mmol O2/m2/d
!      reBTo                             Total benthic primary production    mmol O2/m2/d
!      ruBPc                               N uptake by benthic prim.prod.         mg C/m3
!      ruBPn                               N uptake by benthic prim.prod.     mmol N/m2/d
!      ruBPp                               N uptake by benthic prim.prod.     mmol P/m2/d
!      ruBPs                               N uptake by benthic prim.prod.    mmol Si/m2/d
!     ruBPn3                              NO3 uptake by benthic prim.prod     mmol N/m2/d
!      ruBTc                                         Total benthic uptake         mg C/m3
!      ruBTn                                         Total benthic uptake     mmol N/m2/d
!      ruBTp                                         Total benthic uptake     mmol P/m2/d
!      ruBTs                                         Total benthic uptake    mmol Si/m2/d
!      reBTc                         Total benthic aerobic mineralisation        mgC/m2/d
!      reBTn                         Total benthic aerobic mineralisation     mmol N/m2/d
!      reBTp                         Total benthic aerobic mineralisation     mmol P/m2/d
!      reBTs                     Total benthic aerobic silica dissolution    mmol Si/m2/d
!      reDTn             Benthic anaerobic mineralisation in denitr.layer     mmol N/m2/d
!      reDTp             Benthic anaerobic mineralisation in denitr.layer     mmol P/m2/d
!      reATc                       Total benthic anaerobic mineralisation        mgC/m2/d
!      reATn                       Total benthic anaerobic mineralisation     mmol N/m2/d
!      reATp                       Total benthic anaerobic mineralisation     mmol P/m2/d
!      reATs                   Total benthic anaerobic silica dissolution    mmol Si/m2/d
!  xEPS_Sedi  Total vertical extinction coefficient in top layer sediment             1/m
!     jQuBTn                                             urea consumption     mmol N/m2/d
!     jBTQun                                              urea production     mmol N/m2/d
!     jrrPTc                               rest respiration phytoplankton       mg C/m2/d
!     jrrMec                             rest respiration mesozooplankton       mg C/m2/d
!     jrrMic                            rest respiration microzooplankton       mg C/m2/d
!     irrenh                      Enhancement factor due to bioirrigation               -
!     turenh                       Enhancement factor due to bioturbation               -
!  pxturinD1        Part of bioturbation caused by organism in oxic layer               -
!     jG2K3o                          Oxygen consumption by nitrification    mmol O2/m2/d
!     jG2K7o                      ReOxidation of Red.Equiv. in oxic layer    mmol O2/m2/d
!        M1p                                      phosphate in oxic layer    mmol P/m3 pw
!       M11p                           phosphate in denitrification layer    mmol P/m3 pw
!       M21p                                    phosphate in anoxic layer    mmol P/m3 pw
!        M4n                                       ammonium in oxic layer    mmol N/m3 pw
!       M14n                            ammonium in denitrification layer    mmol N/m3 pw
!       M24n                                     ammonium in anoxic layer    mmol N/m3 pw
!        M3n                                         nitrate in porewater    mmol N/m3 pw
!        M5s                                       silicate in oxic layer   mmol Si/m3 pw
!       M15s                            silicate in denitrification layer   mmol Si/m3 pw
!        M6r                            reduction equivalent in porewater  mmol S--/m3 pw
!       Mp1p                                Phosphate in prim.prod. layer    mmol P/m3 pw
!       Mp3n                            Nitrate in prim.prod. layer layer    mmol N/m2 pw
!       Mp4n                                 Ammonium in prim.prod. layer    mmol N/m2 pw
!       Mp5s                                 Silicate in prim.prod. layer   mmol Si/m2 pw
! fr_lim_Ha_n    max.fract.of uptake limit for Nitrifying archea (basis n)               -
! fr_lim_Ha_o    max.fract.of uptake limit for Nitrifying archea (basis o)               -
!     cNIBTc            concentration of organisms which uptake nutrients            mg C
!      RI_Fc                               Detritus Food for Filterfeeder         mg C/m3
!      RI_Fn                               Detritus Food for Filterfeeder       mmol N/m3
!      RI_Fp                               Detritus Food for Filterfeeder       mmol N/m3
!      RI_Fs                               Detritus Food for Filterfeeder      mmol Si/m3
!     jPTY3c                  Phytoplankton filtered by SuspensionFeeders       mg C/m2/d
!     jPTY3n                  Phytoplankton filtered by SuspensionFeeders     mmol N/m2/d
!     jPTY3p                  Phytoplankton filtered by SuspensionFeeders     mmol P/m2/d
!     jPTZTn                          Phytoplankton grazed by Zooplankton     mmol N/m2/d
!     jPTZTp                         Phytoplankton grazed by Zoooplankton     mmol P/m2/d
!     jZIY3c               microzooplankton filtered by SuspensionFeeders       mg C/m2/d
!     jRIQIc detritus filtered by SuspensionFeeders sedimentated on sediment       mg C/m2/d
!     jRIQIs detritus filtered by SuspensionFeeders sedimentated on sediment    mmol Si/m2/d
!     jY3RIc                 total detritus excreted by SuspensionFeeders       mg C/m2/d
!     jY3RIs                 total detritus excreted by SuspensionFeeders    mmol Si/m2/d
!     jY3O3c      CO2 production by SuspensionFeeders released to pelagic       mg C/m2/d
!     jY3N4n             ammonium release by SuspensionFeeders to pelagic     mmol N/m2/d
!     jY3N1p            phosphate release by SuspensionFeeders to pelagic     mmol P/m2/d
!    jbotQ6c                        resuspension rate of benthic detritus       mg C/m2/d
!    jbotQ6n                        resuspension rate of benthic detritus     mmol N/m2/d
!    jbotQ6p                        resuspension rate of benthic detritus     mmol P/m2/d
!    jbotQ6s                        resuspension rate of benthic detritus    mmol Si/m2/d
!    jbotQ2c                             resuspension rate of benthic TEP       mg C/m2/d
!    jbotQ2n                             resuspension rate of benthic TEP     mmol N/m2/d
!  Depth_Ben                            depth of the layer above sediment               m
!    EIR_Ben                                                   irradiance         mE/m2/s
!    ETW_Ben                                                  temperature               C
!   ERHO_Ben                                                      density            g/m3
!    ESW_Ben                                                     salinity             psu
!    Ru_Benn                           urea concentration in water column       mmol N/M3
!    R3_Benc                                                 TEP in Pheao         mg C/m3
!    R3_Benn                                                 TEP in Pheao       mmol N/m3
!    R3_Benp                                                 TEP in Pheao       mmol P/m3
!    R2_Benc                                     TEP in diatom/aggregates         mg C/m3
!    R2_Benn                                     TEP in diatom/aggregates       mmol N/m3
!    O2o_Ben                                  oxygen conc. in the pelagic      mmol O2/m3
!  cmO2o_Ben                                           oxygen saturation.      mmol O2/m3
!    N1p_Ben                               phosphate conc. in the pelagic       mmol P/m3
!    N3n_Ben                                 nitrate conc. in the pelagic       mmol N/m3
!    N4n_Ben                                ammonium conc. in the pelagic       mmol N/m3
!    N5s_Ben                                silicate conc. in the pelagic      mmol Si/m3
!    N6r_Ben                             red. equiv. conc. in the pelagic     mmol S--/m3
! sediR6_Ben                                  Detritus sedimentation rate             m/d
! sediR2_Ben                                       EPS sedimentation rate             m/d
!   efilP6Y3  suspension-feeder-filtering limitation due to PhaeoColonies               -
!   efilPART  suspension-feeder-filtering limitation due to Siltparticles               -
!    ctfPm2c                  max. available phyto-food for Filterfeeders         mg C/m2
!   ctfZem2c           max. available mesozooo-food for SuspensionFeeders         mg C/m2
!   ctfZim2c           max. available microzoo-food for SuspensionFeeders         mg C/m2
!      sK4K3                               first order nitrification rate              /d
!     jK4K3n                                           nitrification flux     mmol N/m2/d
!     jKuK4n                                  urea flux for nitrification     mmol N/m2/d
!     jK3G4n                                         denitrification flux     mmol N/m2/d
!    jK23G4n                         denitrification flux in anoxic layer     mmol N/m2/d
!   jK31K21p                  flux of phosphate at lower benthic boundary     mmol P/m2/d
!   jK34K24n                   flux of ammonium at lower benthic boundary     mmol N/m2/d
!    jK13K3n                  flux of nitrate at oxygen penetration depth     mmol N/m2/d
!   jK23K13n                    flux of nitrate at lower benthic boundary     mmol N/m2/d
!   jK25K15s                   flux of silicate at lower benthic boundary    mmol Si/m2/d
!   jK36K26r             flux of red.equivalent at lower benthic boundary   mmol S--/m2/d
!    totPELc                                total mass present in pelagic         mg C/m2
!    totPELn                                total mass present in pelagic       mmol N/m2
!    totPELp                                total mass present in pelagic       mmol P/m2
!    totPELs                                total mass present in pelagic      mmol Si/m2
!    totBENc                                total mass present in benthos         mg C/m2
!    totBENn                                total mass present in benthos       mmol N/m2
!    totBENp                                total mass present in benthos       mmol P/m2
!    totBENs                                total mass present in benthos      mmol Si/m2
!    totSYSc                                 total mass present in system         mg C/m2
!    totSYSn                                 total mass present in system       mmol N/m2
!    totSYSp                                 total mass present in system       mmol P/m2
!    totSYSs                                 total mass present in system      mmol Si/m2
!     jtPelc                                            total flux in pel       mg C/m2/d
!     jtPeln                                            total flux in pel     mmol N/m2/d
!     jtPelp                                            total flux in pel     mmol P/m2/d
!     jtPels                                            total flux in pel    mmol Si/m2/d
!     jtBenc                                            total flux in ben       mg C/m2/d
!     jtBenn                                            total flux in ben     mmol N/m2/d
!     jtBenp                                            total flux in ben     mmol P/m2/d
!     jtBens                                            total flux in ben    mmol Si/m2/d
! jtotbenpelc                                    total flux from bento pel       mg C/m2/d
! jtotbenpeln                                    total flux from bento pel     mmol N/m2/d
! jtotbenpelp                                    total flux from bento pel     mmol P/m2/d
! jtotbenpels                                    total flux from bento pel    mmol Si/m2/d
!    jupPELc                            total CO2 uptake by Phytoplankton       mg C/m2/d
!   jminPELc                                 total pelagic mineralization       mg C/m2/d
!    jupBENc                  total CO2 uptake by benthic prim. producers       mg C/m2/d
!   jminBENc                                 total benthic mineralization       mg C/m2/d
!  totPELInc                                       total pelagic CO2 mass         mg C/m2
!  totBENInc                                       total benthic CO2 mass         mg C/m2
!    totPELh                                     total pelagic Alkalinity      mmol H+/m2
!    totBENh                                     total benthic Alkalinity      mmol H+/m2
!  jsdoMesoc            auto mortality loss (=grazing) of mesozooplankton       mg C/m2/d
!  jrsMicroc                            rest respiration microzooplankton       mg C/m2/d
!   jnetMeZc                             net pelagic second. mesoz. prod.       mg C/m2/d
!   jnetMiZc                             net pelagic second. microz prod.       mg C/m2/d
!    jnetB1c                              net pelagic production bacteria       mg C/m2/d
!    jnetY3c                                               net production filterfeeder mg C/m2/d
!    jspaY3c                         net spawning production filterfeeder       mg C/m2/d
!   jnetYy3c                           net larvae production filterfeeder       mg C/m2/d
!   jnetBPTc                               net benthic primary production       mg C/m2/d
!     jPLO3c                             total respiration Pelagic system       mg C/m2/d
!     jZIR6n                     detritus exretion meso+micro zooplankton     mmol N/m2/d
!     jZIR6p                     detritus exretion meso+micro zooplankton     mmol P/m2/d
!     jZIDIn                          loc exretion meso+micro zooplankton     mmol N/m2/d
!     jZIDIp                          loc exretion meso+micro zooplankton     mmol P/m2/d
!  jCaCO3Y3c                                         net shell production       mg C/m2/d
! Output2d_1                                                  Output_2d_1             any
! Output2d_2                                                  Output_2d_2             any
! Output2d_3                                                  Output_2d_3             any
! Output2d_4                                                  Output_2d_4             any
! jPelFishInput                   PotentialFoodAvailability for Pelagic Fish       mg C/m2/d
! jBenFishInput                   PotentialFoodAvailability for Benthic Fish       mg C/m2/d
!    SdTdzTh                                       dtdZ sum in thermoline               T
!       TMLd                                          Top Mix Layer Depth               m
!       BMLd                                   Bottom Mix Layer Top Depth               m
!     Hs_out                                      significant wave height               m
!     Tz_out                                    wave zero crossing period               s
!  u_orb_out                                        wave orbital velocity             m/s
!       TauW                                wave-induced bed-shear stress           m2/s2
!       TauC                             current-induced bed-shear stress           m2/s2
!     TauBed                    wave and current-induced bed-shear-stress            N/s2
!    eta_out                                                ripple height               m
!      DICae                              inorganic C in aerobic sediment          mol/kg
!      DICan                            inorganic C in anaerobic sediment          mol/kg
!    O3c_Ben                               inorg. C conc. in the pelagica         mg C/m3
!    O3h_Ben                               inorg. C conc. in the pelagica         mg C/m3
!    ctO3m2h                                  Total alkalinity in Pelagic          mol/m2
!      CAcae                        Carbon Alkalinity in aerobic sediment          mol/kg
!      CAcan                      Carbon Alkalinity in anaerobic sediment          mol/kg
!       Acae                         Total Alkalinity in aerobic sediment          mol/kg
!       Acan                       Total Alkalinity in anaerobic sediment          mol/kg
!       pHae                                          pH in aerobic layer               -
!       pHdn                                  pH in denitrification layer               -
!       pHan                                        pH in anaerobic layer               -
!     pCO2ae                                          pH in aerobic layer               -
!     pCO2an                                        pH in anaerobic layer               -
!      CO2ae                              CO2 concentration in oxic layer          mol/kg
!     HCO3ae                            HCO3- concentration in oxic layer          mol/kg
!      CO3ae                            CO3-- concentration in oxic layer          mol/kg
!    DIC_Ben                                   inorganic C in the palagic          mol/kg
!      CO2an                            CO2 concentration in anoxic layer          mol/kg
!     HCO3an                          HCO3- concentration in anoxic layer          mol/kg
!      CO3an                          CO3-- concentration in anoxic layer          mol/kg
!   jG33G23h                    alkalinity flux at lower benthic boundary    mmol H+/m2/d
!   jG33G23c                        flux of DIC at lower benthic boundary       mg C/m2/d

! jbotBPc(iiBP1)                         resuspension rate of Benthic Diatoms       mg C/m2/d
! sunBI(iiBP1)                            net prim. prod of Benthic Diatoms              /d
! sugBI(iiBP1)                          gross prim. prod of Benthic Diatoms              /d
! jBTQIc(iiQ1)                                           production of Urea       mg C/m2/d
! jBTQIc(iiQ11)                          production of Labile organic carbon       mg C/m2/d
! jBTQIc(iiQ21)                          production of Labile organic carbon       mg C/m2/d
! jBTQIn(iiQ1)                                           production of Urea     mmol N/m2/d
! jBTQIn(iiQ11)                          production of Labile organic carbon     mmol N/m2/d
! jBTQIn(iiQ21)                          production of Labile organic carbon     mmol N/m2/d
! jBTQIp(iiQ1)                                           production of Urea     mmol P/m2/d
! jBTQIp(iiQ11)                          production of Labile organic carbon     mmol P/m2/d
! jBTQIp(iiQ21)                          production of Labile organic carbon     mmol P/m2/d
! jQIBTc(iiQ1)                                          comsumption of Urea       mg C/m2/d
! jQIBTc(iiQ11)                         comsumption of Labile organic carbon       mg C/m2/d
! jQIBTc(iiQ21)                         comsumption of Labile organic carbon       mg C/m2/d
! jQIBTn(iiQ1)                                          comsumption of Urea     mmol N/m2/d
! jQIBTn(iiQ11)                         comsumption of Labile organic carbon     mmol N/m2/d
! jQIBTn(iiQ21)                         comsumption of Labile organic carbon     mmol N/m2/d
! jQIBTp(iiQ1)                                          comsumption of Urea     mmol P/m2/d
! jQIBTp(iiQ11)                         comsumption of Labile organic carbon     mmol P/m2/d
! jQIBTp(iiQ21)                         comsumption of Labile organic carbon     mmol P/m2/d
! jnetHIc(iiH1)                      net production Aerobic benthic bacteria       mg C/m2/d
! jnetHIc(iiH2)                    net production Anaerobic benthic bacteria       mg C/m2/d
! jnetHIc(iiH3)    net production Anaerobic benthic bacteria in anoxic layer       mg C/m2/d
! jnetHIc(iiHN)    net production Aerobic benthic nitrifying bacteria+archea       mg C/m2/d
! jugYIc(iiY1)                                  gross production Epibenthos       mg C/m2/d
! jugYIc(iiY2)                             gross production Deposit feeders       mg C/m2/d
! jugYIc(iiY4)                                 gross production Meiobenthos       mg C/m2/d
! jugYIc(iiY5)                           gross production Benthic predators       mg C/m2/d
! jnetYIc(iiY1)                                    net production Epibenthos       mg C/m2/d
! jnetYIc(iiY2)                               net production Deposit feeders       mg C/m2/d
! jnetYIc(iiY4)                                   net production Meiobenthos       mg C/m2/d
! jnetYIc(iiY5)                             net production Benthic predators       mg C/m2/d
! jmYIc(iiY1)                                         mortality Epibenthos       mg C/m2/d
! jmYIc(iiY2)                                    mortality Deposit feeders       mg C/m2/d
! jmYIc(iiY4)                                        mortality Meiobenthos       mg C/m2/d
! jmYIc(iiY5)                                  mortality Benthic predators       mg C/m2/d
! jmY3c(iiYy3)                           mortality Young Suspension feeders       mg C/m2/d
! jmY3c(iiY3)                         mortality (Adult) Suspension feeders       mg C/m2/d
! jrrY3c(iiYy3)                   total respiration Young Suspension feeders       mg C/m2/d
! jrrY3c(iiY3)                 total respiration (Adult) Suspension feeders       mg C/m2/d
! jrrYIc(iiY1)                                  rest respiration Epibenthos       mg C/m2/d
! jrrYIc(iiY2)                             rest respiration Deposit feeders       mg C/m2/d
! jrrYIc(iiY4)                                 rest respiration Meiobenthos       mg C/m2/d
! jrrYIc(iiY5)                           rest respiration Benthic predators       mg C/m2/d
! rugY3c(iiYy3)                   gross food uptake Young Suspension feeders       mg C/m2/d
! rugY3c(iiY3)                 gross food uptake (Adult) Suspension feeders       mg C/m2/d
! efsatY3(iiYy3)                     food saturation Young Suspension feeders               -
! efsatY3(iiY3)                   food saturation (Adult) Suspension feeders               -
! fr_lim_HI_n(iiH1) max.fract.of uptake limit for Aerobic benthic bacteria (basis n)               -
! fr_lim_HI_n(iiH2) max.fract.of uptake limit for Anaerobic benthic bacteria (basis n)               -
! fr_lim_HI_n(iiH3) max.fract.of uptake limit for Anaerobic benthic bacteria in anoxic layer (basis n)               -
! fr_lim_HI_n(iiHN) max.fract.of uptake limit for Aerobic benthic nitrifying bacteria+archea (basis n)               -
! fr_lim_HI_p(iiH1) max.fract.of uptake limit for Aerobic benthic bacteria (basis n)               -
! fr_lim_HI_p(iiH2) max.fract.of uptake limit for Anaerobic benthic bacteria (basis n)               -
! fr_lim_HI_p(iiH3) max.fract.of uptake limit for Anaerobic benthic bacteria in anoxic layer (basis n)               -
! fr_lim_HI_p(iiHN) max.fract.of uptake limit for Aerobic benthic nitrifying bacteria+archea (basis n)               -
! fr_lim_HI_o(iiH1)           max.fract.of uptake limit for Nitrifiers (basis o)               -
! fr_lim_HI_o(iiH2)           max.fract.of uptake limit for Nitrifiers (basis o)               -
! fr_lim_HI_o(iiH3)           max.fract.of uptake limit for Nitrifiers (basis o)               -
! fr_lim_HI_o(iiHN)           max.fract.of uptake limit for Nitrifiers (basis o)               -
! fr_lim_BPI_n(iiBP1)      max.fract.of uptake limit for Benthic Diatoms (basis n)               -
! fr_lim_BPI_p(iiBP1)      max.fract.of uptake limit for Benthic Diatoms (basis p)               -
! ZI_Fc(iiZ5)                                  Total Microzooplankton Food         mg C/m3
! ZI_Fc(iiZ6)              Total Heterotrophic nanoflagellates (HNAN) Food         mg C/m3
! ZI_Fn(iiZ5)                                  Total Microzooplankton Food       mmol N/m3
! ZI_Fn(iiZ6)              Total Heterotrophic nanoflagellates (HNAN) Food       mmol N/m3
! ZI_Fp(iiZ5)                                  Total Microzooplankton Food       mmol N/m3
! ZI_Fp(iiZ6)              Total Heterotrophic nanoflagellates (HNAN) Food       mmol N/m3
! jPIY3c(iiP1)                        Diatoms filtered by SuspensionFeeders       mg C/m2/d
! jPIY3c(iiP2)                    Flagellates filtered by SuspensionFeeders       mg C/m2/d
! jPIY3c(iiP3)              PicoPhytoPlankton filtered by SuspensionFeeders       mg C/m2/d
! jPIY3c(iiP4)                Dinoflagellates filtered by SuspensionFeeders       mg C/m2/d
! jPIY3c(iiP5)    Resuspended Benthic Diatoms filtered by SuspensionFeeders       mg C/m2/d
! jPIY3c(iiP6)           Phaeocystis colonies filtered by SuspensionFeeders       mg C/m2/d
! jP6Y3c(iiYy3)             Phaeocystis filtered by Young Suspension feeders       mg C/m2/d
! jP6Y3c(iiY3)           Phaeocystis filtered by (Adult) Suspension feeders       mg C/m2/d
! jZEY3c(iiZ3)    Carnivorous mesozooplankton filtered by SuspensionFeeders       mg C/m2/d
! jZEY3c(iiZ4)     Omnivorous mesozooplankton filtered by SuspensionFeeders       mg C/m2/d
! jZEY3c(iiZ2)             Filterfeederlarvae filtered by SuspensionFeeders       mg C/m2/d
! jPIQ6s(iiBP1)                  biogene silicate release by benthic diatoms    mmol Si/m2/d
! PI_Benc(iiP1)                           Diatoms Forcing for Benthic System         mg C/m3
! PI_Benc(iiP2)                       Flagellates Forcing for Benthic System         mg C/m3
! PI_Benc(iiP3)                 PicoPhytoPlankton Forcing for Benthic System         mg C/m3
! PI_Benc(iiP4)                   Dinoflagellates Forcing for Benthic System         mg C/m3
! PI_Benc(iiP5)       Resuspended Benthic Diatoms Forcing for Benthic System         mg C/m3
! PI_Benc(iiP6)              Phaeocystis colonies Forcing for Benthic System         mg C/m3
! PI_Benn(iiP1)                           Diatoms Forcing for Benthic System       mmol N/m3
! PI_Benn(iiP2)                       Flagellates Forcing for Benthic System       mmol N/m3
! PI_Benn(iiP3)                 PicoPhytoPlankton Forcing for Benthic System       mmol N/m3
! PI_Benn(iiP4)                   Dinoflagellates Forcing for Benthic System       mmol N/m3
! PI_Benn(iiP5)       Resuspended Benthic Diatoms Forcing for Benthic System       mmol N/m3
! PI_Benn(iiP6)              Phaeocystis colonies Forcing for Benthic System       mmol N/m3
! PI_Benp(iiP1)                           Diatoms Forcing for Benthic System       mmol N/m3
! PI_Benp(iiP2)                       Flagellates Forcing for Benthic System       mmol N/m3
! PI_Benp(iiP3)                 PicoPhytoPlankton Forcing for Benthic System       mmol N/m3
! PI_Benp(iiP4)                   Dinoflagellates Forcing for Benthic System       mmol N/m3
! PI_Benp(iiP5)       Resuspended Benthic Diatoms Forcing for Benthic System       mmol N/m3
! PI_Benp(iiP6)              Phaeocystis colonies Forcing for Benthic System       mmol N/m3
! PI_Benl(iiP1)                           Diatoms Forcing for Benthic System       mg Chl/m3
! PI_Benl(iiP2)                       Flagellates Forcing for Benthic System       mg Chl/m3
! PI_Benl(iiP3)                 PicoPhytoPlankton Forcing for Benthic System       mg Chl/m3
! PI_Benl(iiP4)                   Dinoflagellates Forcing for Benthic System       mg Chl/m3
! PI_Benl(iiP5)       Resuspended Benthic Diatoms Forcing for Benthic System       mg Chl/m3
! PI_Benl(iiP6)              Phaeocystis colonies Forcing for Benthic System       mg Chl/m3
! PI_Bens(iiP1)                           Diatoms Forcing for Benthic System      mmol Si/m3
! PI_Bens(iiP2)                       Flagellates Forcing for Benthic System      mmol Si/m3
! PI_Bens(iiP3)                 PicoPhytoPlankton Forcing for Benthic System      mmol Si/m3
! PI_Bens(iiP4)                   Dinoflagellates Forcing for Benthic System      mmol Si/m3
! PI_Bens(iiP5)       Resuspended Benthic Diatoms Forcing for Benthic System      mmol Si/m3
! PI_Bens(iiP6)              Phaeocystis colonies Forcing for Benthic System      mmol Si/m3
! ZE_Benc(iiZ3)       Carnivorous mesozooplankton Forcing for Benthic System         mg C/m3
! ZE_Benc(iiZ4)        Omnivorous mesozooplankton Forcing for Benthic System         mg C/m3
! ZE_Benc(iiZ2)                Filterfeederlarvae Forcing for Benthic System         mg C/m3
! ZE_Benn(iiZ3)       Carnivorous mesozooplankton Forcing for Benthic System       mmol N/m3
! ZE_Benn(iiZ4)        Omnivorous mesozooplankton Forcing for Benthic System       mmol N/m3
! ZE_Benn(iiZ2)                Filterfeederlarvae Forcing for Benthic System       mmol N/m3
! ZE_Benp(iiZ3)       Carnivorous mesozooplankton Forcing for Benthic System       mmol N/m3
! ZE_Benp(iiZ4)        Omnivorous mesozooplankton Forcing for Benthic System       mmol N/m3
! ZE_Benp(iiZ2)                Filterfeederlarvae Forcing for Benthic System       mmol N/m3
! puPIY3(iiP1)                   Diatoms availability for SuspensionFeeders               -
! puPIY3(iiP2)               Flagellates availability for SuspensionFeeders               -
! puPIY3(iiP3)         PicoPhytoPlankton availability for SuspensionFeeders               -
! puPIY3(iiP4)           Dinoflagellates availability for SuspensionFeeders               -
! puPIY3(iiP5) Resuspended Benthic Diatoms availability for SuspensionFeeders               -
! puPIY3(iiP6)      Phaeocystis colonies availability for SuspensionFeeders               -
! puZEY3(iiZ3) Carnivorous mesozooplankton availability for SuspensionFeeders               -
! puZEY3(iiZ4) Omnivorous mesozooplankton availability for SuspensionFeeders               -
! puZEY3(iiZ2)        Filterfeederlarvae availability for SuspensionFeeders               -
! puP6Y3(iiYy3) size related availability of PhaeoColonies  for Young Suspension feeders               -
! puP6Y3(iiY3) size related availability of PhaeoColonies  for (Adult) Suspension feeders               -
! sediPI_Ben(iiP1)                                   Diatoms sedimentation rate             m/d
! sediPI_Ben(iiP2)                               Flagellates sedimentation rate             m/d
! sediPI_Ben(iiP3)                         PicoPhytoPlankton sedimentation rate             m/d
! sediPI_Ben(iiP4)                           Dinoflagellates sedimentation rate             m/d
! sediPI_Ben(iiP5)               Resuspended Benthic Diatoms sedimentation rate             m/d
! sediPI_Ben(iiP6)                      Phaeocystis colonies sedimentation rate             m/d
! sediZE_Ben(iiZ3)               Carnivorous mesozooplankton sedimentation rate             m/d
! sediZE_Ben(iiZ4)                Omnivorous mesozooplankton sedimentation rate             m/d
! sediZE_Ben(iiZ2)                        Filterfeederlarvae sedimentation rate             m/d
! jnetPIc(iiP1)                            net primary production of Diatoms       mg C/m2/d
! jnetPIc(iiP2)                        net primary production of Flagellates       mg C/m2/d
! jnetPIc(iiP3)                  net primary production of PicoPhytoPlankton       mg C/m2/d
! jnetPIc(iiP4)                    net primary production of Dinoflagellates       mg C/m2/d
! jnetPIc(iiP5)        net primary production of Resuspended Benthic Diatoms       mg C/m2/d
! jnetPIc(iiP6)               net primary production of Phaeocystis colonies       mg C/m2/d

     integer,parameter,public :: ppEIRr=1, ppEUWIND=2,&
      & ppEVWIND=3, ppdry_z=4, ppETAUB=5, ppEUCURR_LEVEL=6,&
      & ppEVCURR_LEVEL=7, ppjrESS=8, ppshiftDlm=9,&
      & ppshiftD1m=10, ppshiftD2m=11, ppG2_xavail_o=12,&
      & ppjO2Y2o=13, ppjcrrBTo=14, pprrBTo=15, ppreK6o=16,&
      & pprrDTo=17, pprrATo=18, ppreBTo=19,&
      & ppruBPc=20, ppruBPn=21, ppruBPp=22,&
      & ppruBPs=23, ppruBPn3=24, ppruBTc=25,&
      & ppruBTn=26, ppruBTp=27, ppruBTs=28,&
      & ppreBTc=29, ppreBTn=30, ppreBTp=31,&
      & ppreBTs=32, ppreDTn=33, ppreDTp=34,&
      & ppreATc=35, ppreATn=36, ppreATp=37,&
      & ppreATs=38, ppxEPS_Sedi=39, ppjQuBTn=40,&
      & ppjBTQun=41, ppjrrPTc=42, ppjrrMec=43,&
      & ppjrrMic=44, ppirrenh=45, ppturenh=46,&
      & pppxturinD1=47, ppjG2K3o=48, ppjG2K7o=49,&
      & ppM1p=50, ppM11p=51, ppM21p=52,&
      & ppM4n=53, ppM14n=54, ppM24n=55,&
      & ppM3n=56, ppM5s=57, ppM15s=58,&
      & ppM6r=59, ppMp1p=60, ppMp3n=61,&
      & ppMp4n=62, ppMp5s=63, ppfr_lim_Ha_n=64,&
      & ppfr_lim_Ha_o=65, ppcNIBTc=66, ppRI_Fc=67,&
      & ppRI_Fn=68, ppRI_Fp=69, ppRI_Fs=70,&
      & ppjPTY3c=71, ppjPTY3n=72, ppjPTY3p=73,&
      & ppjPTZTn=74, ppjPTZTp=75, ppjZIY3c=76,&
      & ppjRIQIc=77, ppjRIQIs=78, ppjY3RIc=79,&
      & ppjY3RIs=80, ppjY3O3c=81, ppjY3N4n=82,&
      & ppjY3N1p=83, ppjbotQ6c=84, ppjbotQ6n=85,&
      & ppjbotQ6p=86, ppjbotQ6s=87, ppjbotQ2c=88,&
      & ppjbotQ2n=89, ppDepth_Ben=90, ppEIR_Ben=91,&
      & ppETW_Ben=92, ppERHO_Ben=93, ppESW_Ben=94,&
      & ppRu_Benn=95, ppR3_Benc=96, ppR3_Benn=97,&
      & ppR3_Benp=98, ppR2_Benc=99, ppR2_Benn=100,&
      & ppO2o_Ben=101, ppcmO2o_Ben=102, ppN1p_Ben=103,&
      & ppN3n_Ben=104, ppN4n_Ben=105, ppN5s_Ben=106,&
      & ppN6r_Ben=107, ppsediR6_Ben=108, ppsediR2_Ben=109,&
      & ppefilP6Y3=110, ppefilPART=111, ppctfPm2c=112,&
      & ppctfZem2c=113, ppctfZim2c=114, ppsK4K3=115,&
      & ppjK4K3n=116, ppjKuK4n=117, ppjK3G4n=118,&
      & ppjK23G4n=119, ppjK31K21p=120, ppjK34K24n=121,&
      & ppjK13K3n=122, ppjK23K13n=123, ppjK25K15s=124,&
      & ppjK36K26r=125, pptotPELc=126, pptotPELn=127,&
      & pptotPELp=128, pptotPELs=129, pptotBENc=130,&
      & pptotBENn=131, pptotBENp=132, pptotBENs=133,&
      & pptotSYSc=134, pptotSYSn=135, pptotSYSp=136,&
      & pptotSYSs=137, ppjtPelc=138, ppjtPeln=139,&
      & ppjtPelp=140, ppjtPels=141, ppjtBenc=142,&
      & ppjtBenn=143, ppjtBenp=144, ppjtBens=145,&
      & ppjtotbenpelc=146, ppjtotbenpeln=147, ppjtotbenpelp=148,&
      & ppjtotbenpels=149, ppjupPELc=150, ppjminPELc=151,&
      & ppjupBENc=152, ppjminBENc=153, pptotPELInc=154,&
      & pptotBENInc=155, pptotPELh=156, pptotBENh=157,&
      & ppjsdoMesoc=158, ppjrsMicroc=159, ppjnetMeZc=160,&
      & ppjnetMiZc=161, ppjnetB1c=162, ppjnetY3c=163,&
      & ppjspaY3c=164, ppjnetYy3c=165, ppjnetBPTc=166,&
      & ppjPLO3c=167, ppjZIR6n=168, ppjZIR6p=169,&
      & ppjZIDIn=170, ppjZIDIp=171, ppjCaCO3Y3c=172,&
      & ppOutput2d_1=173, ppOutput2d_2=174, ppOutput2d_3=175,&
      & ppOutput2d_4=176, ppjPelFishInput=177, ppjBenFishInput=178,&
      & ppSdTdzTh=179, ppTMLd=180, ppBMLd=181, ppHs_out=182,&
      & ppTz_out=183, ppu_orb_out=184, ppTauW=185, ppTauC=186,&
      & ppTauBed=187, ppeta_out=188, ppDICae=189, ppDICan=190,&
      & ppO3c_Ben=191, ppO3h_Ben=192, ppctO3m2h=193, ppCAcae=194,&
      & ppCAcan=195, ppAcae=196, ppAcan=197, pppHae=198,&
      & pppHdn=199, pppHan=200, pppCO2ae=201, pppCO2an=202,&
      & ppCO2ae=203, ppHCO3ae=204, ppCO3ae=205, ppDIC_Ben=206,&
      & ppCO2an=207, ppHCO3an=208, ppCO3an=209, ppjG33G23h=210,&
      & ppjG33G23c=211

     integer,public ::&
      & ppjbotBPc(iiBenPhyto),&
      & ppsunBI(iiBenPhyto),&
      & ppsugBI(iiBenPhyto),&
      & ppjBTQIc(iiBenLabileDetritus),&
      & ppjBTQIn(iiBenLabileDetritus),&
      & ppjBTQIp(iiBenLabileDetritus),&
      & ppjQIBTc(iiBenLabileDetritus),&
      & ppjQIBTn(iiBenLabileDetritus),&
      & ppjQIBTp(iiBenLabileDetritus),&
      & ppjnetHIc(iiBenBacteria),&
      & ppjugYIc(iiBenOrganisms),&
      & ppjnetYIc(iiBenOrganisms),&
      & ppjmYIc(iiBenOrganisms),&
      & ppjmY3c(iiSuspensionFeeders),&
      & ppjrrY3c(iiSuspensionFeeders),&
      & ppjrrYIc(iiBenOrganisms),&
      & pprugY3c(iiSuspensionFeeders),&
      & ppefsatY3(iiSuspensionFeeders),&
      & ppfr_lim_HI_n(iiBenBacteria),&
      & ppfr_lim_HI_p(iiBenBacteria),&
      & ppfr_lim_HI_o(iiBenBacteria),&
      & ppfr_lim_BPI_n(iiBenPhyto),&
      & ppfr_lim_BPI_p(iiBenPhyto),&
      & ppZI_Fc(iiMicroZooPlankton),&
      & ppZI_Fn(iiMicroZooPlankton),&
      & ppZI_Fp(iiMicroZooPlankton),&
      & ppjPIY3c(iiPhytoPlankton),&
      & ppjP6Y3c(iiSuspensionFeeders),&
      & ppjZEY3c(iiMesoZooPlankton),&
      & ppjPIQ6s(iiBenPhyto),&
      & ppPI_Benc(iiPhytoPlankton),&
      & ppPI_Benn(iiPhytoPlankton),&
      & ppPI_Benp(iiPhytoPlankton),&
      & ppPI_Benl(iiPhytoPlankton),&
      & ppPI_Bens(iiPhytoPlankton),&
      & ppZE_Benc(iiMesoZooPlankton),&
      & ppZE_Benn(iiMesoZooPlankton),&
      & ppZE_Benp(iiMesoZooPlankton),&
      & pppuPIY3(iiPhytoPlankton),&
      & pppuZEY3(iiMesoZooPlankton),&
      & pppuP6Y3(iiSuspensionFeeders),&
      & ppsediPI_Ben(iiPhytoPlankton), ppsediZE_Ben(iiMesoZooPlankton),&
      & ppjnetPIc(iiPhytoPlankton)

     real(RLEN),public,dimension(:),pointer :: EIRr, EUWIND, EVWIND, dry_z,&
      & ETAUB, EUCURR_LEVEL, EVCURR_LEVEL, jrESS, shiftDlm, shiftD1m,&
      & shiftD2m, G2_xavail_o, jO2Y2o, jcrrBTo, rrBTo, reK6o, rrDTo,&
      & rrATo, reBTo, ruBPc, ruBPn, ruBPp, ruBPs, ruBPn3,&
      & ruBTc, ruBTn, ruBTp, ruBTs, reBTc, reBTn, reBTp,&
      & reBTs, reDTn, reDTp, reATc, reATn, reATp, reATs,&
      & xEPS_Sedi, jQuBTn, jBTQun, jrrPTc, jrrMec, jrrMic, irrenh,&
      & turenh, pxturinD1, jG2K3o, jG2K7o, M1p, M11p, M21p,&
      & M4n, M14n, M24n, M3n, M5s, M15s, M6r,&
      & Mp1p, Mp3n, Mp4n, Mp5s, fr_lim_Ha_n, fr_lim_Ha_o,&
      & cNIBTc, RI_Fc, RI_Fn, RI_Fp, RI_Fs, jPTY3c,&
      & jPTY3n, jPTY3p, jPTZTn, jPTZTp, jZIY3c, jRIQIc,&
      & jRIQIs, jY3RIc, jY3RIs, jY3O3c, jY3N4n, jY3N1p,&
      & jbotQ6c, jbotQ6n, jbotQ6p, jbotQ6s, jbotQ2c, jbotQ2n,&
      & Depth_Ben, EIR_Ben, ETW_Ben, ERHO_Ben, ESW_Ben, Ru_Benn,&
      & R3_Benc, R3_Benn, R3_Benp, R2_Benc, R2_Benn, O2o_Ben,&
      & cmO2o_Ben, N1p_Ben, N3n_Ben, N4n_Ben, N5s_Ben, N6r_Ben,&
      & sediR6_Ben, sediR2_Ben, efilP6Y3, efilPART, ctfPm2c, ctfZem2c,&
      & ctfZim2c, sK4K3, jK4K3n, jKuK4n, jK3G4n, jK23G4n,&
      & jK31K21p, jK34K24n, jK13K3n, jK23K13n, jK25K15s, jK36K26r,&
      & totPELc, totPELn, totPELp, totPELs, totBENc, totBENn,&
      & totBENp, totBENs, totSYSc, totSYSn, totSYSp, totSYSs,&
      & jtPelc, jtPeln, jtPelp, jtPels, jtBenc, jtBenn,&
      & jtBenp, jtBens, jtotbenpelc, jtotbenpeln, jtotbenpelp, jtotbenpels,&
      & jupPELc, jminPELc, jupBENc, jminBENc, totPELInc, totBENInc, totPELh,&
      & totBENh, jsdoMesoc, jrsMicroc, jnetMeZc, jnetMiZc, jnetB1c, jnetY3c,&
      & jspaY3c, jnetYy3c, jnetBPTc, jPLO3c, jZIR6n, jZIR6p, jZIDIn,&
      & jZIDIp, jCaCO3Y3c, Output2d_1, Output2d_2, Output2d_3,&
      & Output2d_4, jPelFishInput, jBenFishInput, SdTdzTh, TMLd, BMLd,&
      & Hs_out, Tz_out, u_orb_out, TauW, TauC, TauBed, eta_out, DICae,&
      & DICan, O3c_Ben, O3h_Ben, ctO3m2h, CAcae, CAcan, Acae, Acan,&
      & pHae, pHdn, pHan, pCO2ae, pCO2an, CO2ae, HCO3ae, CO3ae,&
      & DIC_Ben, CO2an, HCO3an, CO3an, jG33G23h, jG33G23c

     real(RLEN),public,dimension(:,:),pointer :: jbotBPc, sunBI, sugBI,&
      & jBTQIc, jBTQIn, jBTQIp, jQIBTc, jQIBTn, jQIBTp,&
      & jnetHIc, jugYIc, jnetYIc, jmYIc, jmY3c, jrrY3c,&
      & jrrYIc, rugY3c, efsatY3, fr_lim_HI_n, fr_lim_HI_p, fr_lim_HI_o,&
      & fr_lim_BPI_n, fr_lim_BPI_p, ZI_Fc, ZI_Fn, ZI_Fp, jPIY3c, jP6Y3c,&
      & jZEY3c, jPIQ6s, PI_Benc, PI_Benn, PI_Benp, PI_Benl, PI_Bens,&
      & ZE_Benc, ZE_Benn, ZE_Benp, puPIY3, puZEY3, puP6Y3, sediPI_Ben,&
      & sediZE_Ben, jnetPIc


   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   !  boundary fluxes
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     real(RLEN),public,dimension(:),pointer :: jsurR9x, jsurO2o, jsurN1p,&
      & jsurN3n, jsurN4n, jsurN5s, jsurN6r, jsurB1c, jsurB1n,&
      & jsurB1p, jsurBac, jsurP1c, jsurP1n, jsurP1p, jsurP1l,&
      & jsurP1s, jsurP2c, jsurP2n, jsurP2p, jsurP2l, jsurP3c,&
      & jsurP3n, jsurP3p, jsurP3l, jsurP4c, jsurP4n, jsurP4p, jsurP4l,&
      & jsurP5c, jsurP5n, jsurP5p, jsurP5l, jsurP5s, jsurP6c, jsurP6n,&
      & jsurP6p, jsurP6l, jsurPcc, jsurR1c, jsurR1n, jsurR1p, jsurR2c,&
      & jsurR2n, jsurR3c, jsurR6c, jsurR6n, jsurR6p, jsurR6s, jsurRZc,&
      & jsurO3c, jsurO3h, jsurZ3c, jsurZ4c, jsurZ2c, jsurZ5c, jsurZ6c

     real(RLEN),public,dimension(:),pointer :: jbotR9x, jbotO2o, jbotN1p,&
      & jbotN3n, jbotN4n, jbotN5s, jbotN6r, jbotB1c, jbotB1n,&
      & jbotB1p, jbotBac, jbotP1c, jbotP1n, jbotP1p, jbotP1l,&
      & jbotP1s, jbotP2c, jbotP2n, jbotP2p, jbotP2l, jbotP3c,&
      & jbotP3n, jbotP3p, jbotP3l, jbotP4c, jbotP4n, jbotP4p, jbotP4l,&
      & jbotP5c, jbotP5n, jbotP5p, jbotP5l, jbotP5s, jbotP6c, jbotP6n,&
      & jbotP6p, jbotP6l, jbotPcc, jbotR1c, jbotR1n, jbotR1p, jbotR2c,&
      & jbotR2n, jbotR3c, jbotR6c, jbotR6n, jbotR6p, jbotR6s, jbotRZc,&
      & jbotO3c, jbotO3h, jbotZ3c, jbotZ4c, jbotZ2c, jbotZ5c, jbotZ6c


   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   !  Other 3d-Global Variables 
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     real(RLEN),public,dimension(:),allocatable           :: &
      OCDepth,& ! Cumulative/Column Depth measured from the surface
      Depth,& !depth of a layer
      ABIO_eps      ! the abiotic extinction coeff. calculated in silt models



   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   !  Other 2d-Global Variables 
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     real(RLEN),public,dimension(:),allocatable           :: &
      CO2_Ben,& !CO2 concentration in the pelagic mol/kg
      HCO3_Ben,& !HCO3- concentration in the peagic mol/kg
      CO3_Ben      !CO3-- concentration in the pelagic mol/kg

     integer,public,dimension(:),allocatable           :: &
      KPO4,& !dynamic profile
      KPO4sh,& !dynamic profile to estimate desorption/adsorption at sulphide horizon
      KNH4p,& !
      KNH4,& !
      KNO3,& !
      KNO3E,& !
      KRED,& !
      KSIO3eq,& !equilibirum profile
      KSIO3,& !dynamic profile to estimate desorption/adsorption at sulphide horizon
      KQ1c,& !
      KQun      ! Section -decr
      integer,public,dimension(:,:),allocatable         :: &
       sw_CalcPhyto      !var with which phytoplankton group can be set off in cas of low numbers Section -decr
       integer,public,dimension(:),allocatable           :: &
        PelBoxAbove,& !Pelagic boxes related to adjacent benthic box.
        KCO2,& !
        KHPLUS      !

   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   !  variables to generate flux_output 
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
           real(RLEN), public,dimension(:),allocatable ::flx_t
           integer,    public,dimension(:),allocatable ::flx_SS
           integer,    public,dimension(:),allocatable ::flx_states
           integer,    public,dimension(:),allocatable ::flx_ostates
           integer,    public,dimension(:),allocatable ::flx_calc_nr
           integer,    public,dimension(:),allocatable ::flx_CalcIn
           integer,    public,dimension(:),allocatable ::flx_option
           integer,    public                          ::flx_cal_ben_start
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   !! variables to  Track a nutrient throught the model
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
           ! iitrack>=1 : tracking procedures are active!
           integer,public     :: iiTrack=0
           integer,parameter,public     :: ii3dTrack=0
           integer,parameter,public     :: ii3dptTrack=0
           integer,parameter,public     :: ii3daptTrack=0
           integer,public,dimension(:),allocatable ::nr_3d_track
           integer,parameter,public     :: ii2dTrack=0
           integer,parameter,public     :: ii2dptTrack=0
           integer,public,dimension(:),allocatable ::nr_2d_track
           integer,public,dimension(:),allocatable&
            & ::flag_3d_track_bot
           integer,public,dimension(:),allocatable ::fix_3d_track_bot
           integer,public,dimension(:),allocatable ::fix_2d_track_bot
           real(RLEN),public,dimension(:),allocatable&
            & ::check_3d_track_bot

   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   !! SHARED GLOBAL FUNCTIONS (must be below contains)
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        public flux, flux_vector, Source, Source_D3_vector, Source_D2_vector , &  
               ResetSource, ResetSource_D3_vector, ResetSource_D2_vector , &  
               fixed_quota_flux_vector,sourcesink_flux_vector, &
               sourcesink_flux, set_for_state_fluxes_zero
       public ppPhytoPlankton, ppPelDetritus, ppInorganic,&
        & ppMesoZooPlankton, ppMicroZooPlankton, PhytoPlankton,&
        & PelDetritus, Inorganic, MesoZooPlankton, MicroZooPlankton

       public ppBenPhyto, ppBenOrganisms, ppSuspensionFeeders,&
        & ppBenUrea, ppBenCarboHydrates, ppBenLabileDetritus,&
        & ppBenBacteria, ppBenthicPhosphate, ppBenthicAmmonium,&
        & ppBenthicSilicate, ppBenthicRedEq, ppBenthicNitrate,&
        & ppBenthicCO2, BenPhyto, BenOrganisms,&
        & SuspensionFeeders, BenUrea, BenCarboHydrates,&
        & BenLabileDetritus, BenBacteria, BenthicPhosphate,&
        & BenthicAmmonium, BenthicSilicate, BenthicRedEq, BenthicNitrate,&
        & BenthicCO2


        contains

     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     !! Group Pelagic (D3) state functions
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

       function ppPhytoPlankton(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        integer ::ppPhytoPlankton
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        integer,dimension(6) ::&
         & referto=(/ppP1c,ppP2c,ppP3c,ppP4c,ppP5c,ppP6c/)
        integer,dimension(6) :: const_max=(/5,4,4,4,5,4/)
        integer,dimension(9) :: constituent_add=(/0,1,2,3,4,0,0,0,0/)

        if ( constituent <=const_max(n) ) then
         ppPhytoPlankton=referto(n)+ constituent_add(constituent)
        else
         ppPhytoPlankton=0
        endif

       END function

       function ppPelDetritus(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        integer ::ppPelDetritus
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        integer,dimension(4) :: referto=(/ppR1c,ppR2c,ppR3c,ppR6c/)
        integer,dimension(4) :: const_max=(/3,2,1,5/)
        integer,dimension(9) :: constituent_add=(/0,1,2,0,3,0,0,0,0/)

        if ( constituent <=const_max(n) ) then
         ppPelDetritus=referto(n)+ constituent_add(constituent)
        else
         ppPelDetritus=0
        endif

       END function

       function ppInorganic(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        integer ::ppInorganic
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        integer,dimension(1) :: referto=(/ppO3c/)
        integer,dimension(1) :: const_max=(/6/)
        integer,dimension(9) :: constituent_add=(/0,0,0,0,0,1,0,0,0/)

        if ( constituent <=const_max(n) ) then
         ppInorganic=referto(n)+ constituent_add(constituent)
        else
         ppInorganic=0
        endif

       END function

       function ppMesoZooPlankton(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        integer ::ppMesoZooPlankton
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        integer,dimension(3) :: referto=(/ppZ3c,ppZ4c,ppZ2c/)
        integer,dimension(3) :: const_max=(/1,1,1/)
        integer,dimension(9) :: constituent_add=(/0,1,2,0,0,0,0,0,0/)

        if ( constituent <=const_max(n) ) then
         ppMesoZooPlankton=referto(n)+ constituent_add(constituent)
        else
         ppMesoZooPlankton=0
        endif

       END function

       function ppMicroZooPlankton(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        integer ::ppMicroZooPlankton
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        integer,dimension(2) :: referto=(/ppZ5c,ppZ6c/)
        integer,dimension(2) :: const_max=(/1,1/)
        integer,dimension(9) :: constituent_add=(/0,1,2,0,0,0,0,0,0/)

        if ( constituent <=const_max(n) ) then
         ppMicroZooPlankton=referto(n)+ constituent_add(constituent)
        else
         ppMicroZooPlankton=0
        endif

       END function

       function PhytoPlankton(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        real(RLEN),dimension(:),pointer ::PhytoPlankton
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        PhytoPlankton => D3STATE(:,ppPhytoPlankton(n,constituent))

       END function

       function PelDetritus(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        real(RLEN),dimension(:),pointer ::PelDetritus
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        PelDetritus => D3STATE(:,ppPelDetritus(n,constituent))

       END function

       function Inorganic(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        real(RLEN),dimension(:),pointer ::Inorganic
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        Inorganic => D3STATE(:,ppInorganic(n,constituent))

       END function

       function MesoZooPlankton(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        real(RLEN),dimension(:),pointer ::MesoZooPlankton
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        MesoZooPlankton => D3STATE(:,ppMesoZooPlankton(n,constituent))

       END function

       function MicroZooPlankton(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        real(RLEN),dimension(:),pointer ::MicroZooPlankton
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        MicroZooPlankton => D3STATE(:,ppMicroZooPlankton(n,constituent))

       END function


     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     !! Group Benthic (D2) state functions
     !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

       function ppBenPhyto(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        integer ::ppBenPhyto
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        integer,dimension(1) :: referto=(/ppBP1c/)
        integer,dimension(1) :: const_max=(/5/)
        integer,dimension(9) :: constituent_add=(/0,1,2,3,4,0,0,0,0/)

        if ( constituent <=const_max(n) ) then
         ppBenPhyto=referto(n)+ constituent_add(constituent)
        else
         ppBenPhyto=0
        endif

       END function

       function ppBenOrganisms(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        integer ::ppBenOrganisms
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        integer,dimension(4) :: referto=(/ppY1c,ppY2c,ppY4c,ppY5c/)
        integer,dimension(4) :: const_max=(/3,3,3,3/)
        integer,dimension(9) :: constituent_add=(/0,1,2,0,0,0,0,0,0/)

        if ( constituent <=const_max(n) ) then
         ppBenOrganisms=referto(n)+ constituent_add(constituent)
        else
         ppBenOrganisms=0
        endif

       END function

       function ppSuspensionFeeders(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        integer ::ppSuspensionFeeders
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        integer,dimension(2) :: referto=(/ppYy3c,ppY3c/)
        integer,dimension(2) :: const_max=(/3,3/)
        integer,dimension(9) :: constituent_add=(/0,1,2,0,0,0,0,0,0/)

        if ( constituent <=const_max(n) ) then
         ppSuspensionFeeders=referto(n)+ constituent_add(constituent)
        else
         ppSuspensionFeeders=0
        endif

       END function

       function ppBenUrea(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        integer ::ppBenUrea
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        integer,dimension(3) :: referto=(/ppQun,ppQ1un,ppQ2un/)
        integer,dimension(3) :: const_max=(/2,2,2/)
        integer,dimension(9) :: constituent_add=(/0,0,0,0,0,0,0,0,0/)

        if ( constituent <=const_max(n) ) then
         ppBenUrea=referto(n)+ constituent_add(constituent)
        else
         ppBenUrea=0
        endif

       END function

       function ppBenCarboHydrates(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        integer ::ppBenCarboHydrates
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        integer,dimension(2) :: referto=(/ppQ2c,ppQ12c/)
        integer,dimension(2) :: const_max=(/2,2/)
        integer,dimension(9) :: constituent_add=(/0,1,0,0,0,0,0,0,0/)

        if ( constituent <=const_max(n) ) then
         ppBenCarboHydrates=referto(n)+ constituent_add(constituent)
        else
         ppBenCarboHydrates=0
        endif

       END function

       function ppBenLabileDetritus(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        integer ::ppBenLabileDetritus
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        integer,dimension(3) :: referto=(/ppQ1c,ppQ11c,ppQ21c/)
        integer,dimension(3) :: const_max=(/3,3,3/)
        integer,dimension(9) :: constituent_add=(/0,1,2,0,0,0,0,0,0/)

        if ( constituent <=const_max(n) ) then
         ppBenLabileDetritus=referto(n)+ constituent_add(constituent)
        else
         ppBenLabileDetritus=0
        endif

       END function

       function ppBenBacteria(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        integer ::ppBenBacteria
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        integer,dimension(4) :: referto=(/ppH1c,ppH2c,ppH3c,ppHNc/)
        integer,dimension(4) :: const_max=(/3,3,3,3/)
        integer,dimension(9) :: constituent_add=(/0,1,2,0,0,0,0,0,0/)

        if ( constituent <=const_max(n) ) then
         ppBenBacteria=referto(n)+ constituent_add(constituent)
        else
         ppBenBacteria=0
        endif

       END function

       function ppBenthicPhosphate(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        integer ::ppBenthicPhosphate
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        integer,dimension(3) :: referto=(/ppK1p,ppK11p,ppK21p/)
        integer,dimension(3) :: const_max=(/3,3,3/)
        integer,dimension(9) :: constituent_add=(/0,0,0,0,0,0,0,0,0/)

        if ( constituent <=const_max(n) ) then
         ppBenthicPhosphate=referto(n)+ constituent_add(constituent)
        else
         ppBenthicPhosphate=0
        endif

       END function

       function ppBenthicAmmonium(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        integer ::ppBenthicAmmonium
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        integer,dimension(3) :: referto=(/ppK4n,ppK14n,ppK24n/)
        integer,dimension(3) :: const_max=(/2,2,2/)
        integer,dimension(9) :: constituent_add=(/0,0,0,0,0,0,0,0,0/)

        if ( constituent <=const_max(n) ) then
         ppBenthicAmmonium=referto(n)+ constituent_add(constituent)
        else
         ppBenthicAmmonium=0
        endif

       END function

       function ppBenthicSilicate(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        integer ::ppBenthicSilicate
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        integer,dimension(2) :: referto=(/ppK5s,ppK15s/)
        integer,dimension(2) :: const_max=(/5,5/)
        integer,dimension(9) :: constituent_add=(/0,0,0,0,0,0,0,0,0/)

        if ( constituent <=const_max(n) ) then
         ppBenthicSilicate=referto(n)+ constituent_add(constituent)
        else
         ppBenthicSilicate=0
        endif

       END function

       function ppBenthicRedEq(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        integer ::ppBenthicRedEq
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        integer,dimension(3) :: referto=(/ppK6r,ppK16r,ppK26r/)
        integer,dimension(3) :: const_max=(/7,7,7/)
        integer,dimension(9) :: constituent_add=(/0,0,0,0,0,0,0,0,0/)

        if ( constituent <=const_max(n) ) then
         ppBenthicRedEq=referto(n)+ constituent_add(constituent)
        else
         ppBenthicRedEq=0
        endif

       END function

       function ppBenthicNitrate(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        integer ::ppBenthicNitrate
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        integer,dimension(3) :: referto=(/ppK3n,ppK13n,ppK23n/)
        integer,dimension(3) :: const_max=(/2,2,2/)
        integer,dimension(9) :: constituent_add=(/0,0,0,0,0,0,0,0,0/)

        if ( constituent <=const_max(n) ) then
         ppBenthicNitrate=referto(n)+ constituent_add(constituent)
        else
         ppBenthicNitrate=0
        endif

       END function

       function ppBenthicCO2(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        integer ::ppBenthicCO2
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        integer,dimension(3) :: referto=(/ppG3c,ppG13c,ppG23c/)
        integer,dimension(3) :: const_max=(/6,6,6/)
        integer,dimension(9) :: constituent_add=(/0,0,0,0,0,1,0,0,0/)

        if ( constituent <=const_max(n) ) then
         ppBenthicCO2=referto(n)+ constituent_add(constituent)
        else
         ppBenthicCO2=0
        endif

       END function

       function BenPhyto(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        real(RLEN),dimension(:),pointer ::BenPhyto
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        BenPhyto => D2STATE(:,ppBenPhyto(n,constituent))

       END function

       function BenOrganisms(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        real(RLEN),dimension(:),pointer ::BenOrganisms
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        BenOrganisms => D2STATE(:,ppBenOrganisms(n,constituent))

       END function

       function SuspensionFeeders(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        real(RLEN),dimension(:),pointer ::SuspensionFeeders
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        SuspensionFeeders =>&
         & D2STATE(:,ppSuspensionFeeders(n,constituent))

       END function

       function BenUrea(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        real(RLEN),dimension(:),pointer ::BenUrea
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        BenUrea => D2STATE(:,ppBenUrea(n,constituent))

       END function

       function BenCarboHydrates(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        real(RLEN),dimension(:),pointer ::BenCarboHydrates
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        BenCarboHydrates => D2STATE(:,ppBenCarboHydrates(n,constituent))

       END function

       function BenLabileDetritus(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        real(RLEN),dimension(:),pointer ::BenLabileDetritus
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        BenLabileDetritus =>&
         & D2STATE(:,ppBenLabileDetritus(n,constituent))

       END function

       function BenBacteria(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        real(RLEN),dimension(:),pointer ::BenBacteria
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        BenBacteria => D2STATE(:,ppBenBacteria(n,constituent))

       END function

       function BenthicPhosphate(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        real(RLEN),dimension(:),pointer ::BenthicPhosphate
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        BenthicPhosphate => D2STATE(:,ppBenthicPhosphate(n,constituent))

       END function

       function BenthicAmmonium(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        real(RLEN),dimension(:),pointer ::BenthicAmmonium
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        BenthicAmmonium => D2STATE(:,ppBenthicAmmonium(n,constituent))

       END function

       function BenthicSilicate(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        real(RLEN),dimension(:),pointer ::BenthicSilicate
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        BenthicSilicate => D2STATE(:,ppBenthicSilicate(n,constituent))

       END function

       function BenthicRedEq(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        real(RLEN),dimension(:),pointer ::BenthicRedEq
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        BenthicRedEq => D2STATE(:,ppBenthicRedEq(n,constituent))

       END function

       function BenthicNitrate(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        real(RLEN),dimension(:),pointer ::BenthicNitrate
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        BenthicNitrate => D2STATE(:,ppBenthicNitrate(n,constituent))

       END function

       function BenthicCO2(n,constituent)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        IMPLICIT NONE

        real(RLEN),dimension(:),pointer ::BenthicCO2
        integer, intent(IN) ::n
        integer, intent(IN) ::constituent

        BenthicCO2 => D2STATE(:,ppBenthicCO2(n,constituent))

       END function



          !With this routine all fluxes can be resetted on zero.
          ! This routine can be used when you decide during a run to exclude 
          ! a state variable in the calulationss
          subroutine set_for_state_fluxes_zero( iiSub,iistate)
            use constants, only: RLEN, ZERO,  SEC_PER_DAY
            use global_mem, only: LOGUNIT
            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            ! Implicit typing is never allowed
            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            implicit none

            integer,intent(IN) :: iiSub
            integer,intent(IN) :: iistate
          
            if ( iiSub== iiBen) then
              call ResetSource_D2_vector(iistate)
            else
              call ResetSource_D3_vector(iistate)
            endif
          end subroutine set_for_state_fluxes_zero

          subroutine flux_vector(iiSub,origin,destination,flux)

            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

            use constants, only: RLEN, ZERO,  SEC_PER_DAY
            use global_mem, only: LOGUNIT
            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            ! Implicit typing is never allowed
            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            implicit none

            integer,intent(IN) :: iiSub
            integer,intent(IN) :: origin
            integer,intent(IN) :: destination
            real(RLEN),intent(IN) :: flux(:)

            integer :: i
            character(len=8) :: D23

            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            !BEGIN compute
            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

            TESTNANVECTOR(flux,iiSub,origin,destination)
            CHECKFLUX(-1,iiSub,origin,destination)

            if ( origin /= destination )  then
              if ( minval(flux) < ZERO) then
                D23="Pelagic"
                if ( iiSub == iiBen) D23="Benthic"
                do i=1,size(flux)
                  if (flux(i)< 0.0D+00)  then
                    write(LOGUNIT,'(''Error flux_vector: negative flux at level:'',I4)') i
                    write(LOGUNIT,'(''In '',A,'':origin='',i4,'' destination='',i4)') &
                      D23, origin,destination
                    write(LOGUNIT,'(''In '',A,'':origin='',i4,'' &
                      destination='',i4)') D23, origin,destination
                    write(LOGUNIT,'(''flux='',(G16.8))') flux(i)
                    if ( iiSub== iiBen) then
                      write(LOGUNIT,*) "state value origin:",D2STATE(i,origin)
                      write(LOGUNIT,*) "state value destination:",D2STATE(i,destination)
                    else
                      write(LOGUNIT,*) "state value origin:",D3STATE(i,origin)
                      write(LOGUNIT,*) "state value destination:",D3STATE(i,destination)
                    endif
                  endif
                enddo
                call BFM_ERROR("flux_vector","negative flux")
              endif ! minval<0
              select case ( iiSub )
                case (iiPel)
                  D3SINK(:,origin,destination)  =  flux/SEC_PER_DAY
                  D3SOURCE(:,destination,origin)=  flux/SEC_PER_DAY
                case (iiBen)
                  D2SINK(:,origin,destination) =  flux/SEC_PER_DAY
                  D2SOURCE(:,destination,origin)   = flux/SEC_PER_DAY
              end select
            else
              select case ( iiSub )
                case (iiPel)
                  where (flux > 0.0D+00 )
                    D3SOURCE(:,origin,destination) =D3SOURCE(:,origin,destination) &
                      + flux/SEC_PER_DAY
                  elsewhere
                    D3SINK(:,destination,origin) =D3SINK(:,destination,origin) - &
                      flux/SEC_PER_DAY
                  end where
                case (iiBen)
                  where (flux > 0.0D+00 )
                    D2SOURCE(:,destination,origin) =D2SOURCE(:,destination,origin) &
                      + flux/SEC_PER_DAY
                  elsewhere
                    D2SINK(:,origin,destination) =D2SINK(:,origin,destination) - &
                      flux/SEC_PER_DAY
                  end where
              end select
            endif !origin <> destination

            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            !END compute
            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

            return
          end subroutine flux_vector

          subroutine testnan_vector(array,iiSub,origin,destination)
          use global_mem, only: LOGUNIT

            real(RLEN),intent(IN)    :: array(:)
            integer,intent(IN) :: iiSub
            integer,intent(IN) :: origin
            integer,intent(IN) :: destination
            integer:: i=0
            do i=1,size(array)
!JM              if (isnan(array(i))== .true. ) then
              if (isnan(array(i)).eqv. .true. ) then
                write(LOGUNIT,'(''at level:'',I4)') i
                write(LOGUNIT,'(''origin='',i4,'' destination='',i4)') &
                  origin,destination
                if ( iiSub== iiBen) then
                    write(LOGUNIT,*) "state value origin:",D2STATE(i,origin)
                    write(LOGUNIT,*) "state value destination:",D2STATE(i,destination)
                else
                    write(LOGUNIT,*) "state value origin:",D3STATE(i,origin)
                    write(LOGUNIT,*) "state value destination:",D3STATE(i,destination)
                endif
                STDERR 'Nan value in flux'
                stop 1002
              endif
            enddo
          end subroutine testnan_vector

          subroutine testnan(scalar,grid_nr,iiSub,origin,destination)
          use global_mem, only: LOGUNIT

            real(RLEN),intent(IN)    :: scalar
            integer,intent(IN) :: grid_nr
            integer,intent(IN) :: iiSub
            integer,intent(IN) :: origin
            integer,intent(IN) :: destination
!JM            if (isnan(scalar)== .true. ) then
            if (isnan(scalar).eqv. .true. ) then
            write(LOGUNIT,'(''origin='',i4,'' destination='',i4)') origin,destination
            if ( iiSub == iiBen)  then
                 write(LOGUNIT,*) "state value origin:",D2STATE(grid_nr,origin)
                 write(LOGUNIT,*) "state value destination:",D2STATE(grid_nr,destination)
              else 
                 write(LOGUNIT,*) "state value origin:",D3STATE(grid_nr,origin)
                 write(LOGUNIT,*) "state value destination:",D3STATE(grid_nr,destination)
            endif
            write(LOGUNIT,*) 'Nan value in scalar flux'
            stop 1003
          endif
        end subroutine testnan
        subroutine flux(grid_nr,iiSub,origin,destination,flow,error)

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

          use constants, only: RLEN, ZERO, SEC_PER_DAY
          use global_mem, only: LOGUNIT
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          ! Implicit typing is never allowed
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          implicit none

          integer,intent(IN)                 :: grid_nr
          integer,intent(IN)                 :: iiSub
          integer,intent(IN)                 :: origin
          integer,intent(IN)                 :: destination
          real(RLEN),intent(IN)              :: flow
          integer,intent(INOUT),optional     :: error

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          character(len=8)    :: D23
          !BEGIN compute
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

          TESTNAN(flow,grid_nr,iiSub,origin,destination)
          CHECKFLUX(grid_nr,iiSub,origin,destination)

          if ( origin /= destination ) then
            if ( flow < ZERO) then
              D23="Pelagic"
              if ( iiSub == iiBen) D23="Benthic"
              write(LOGUNIT,'(''In '',A,'':origin='',i4,'' destination='',i4)') &
                D23, origin,destination
              write(LOGUNIT,*) "Error in (scalar) vector  function: negative flux!"
              write(LOGUNIT,*) "origin,destination:", origin,destination
              write(LOGUNIT,*) flow
              if ( iiSub == iiBen)  then
                 write(LOGUNIT,*) "state value origin:",D2STATE(grid_nr,origin)
                 write(LOGUNIT,*) "state value destination:",D2STATE(grid_nr,destination)
              else 
                 write(LOGUNIT,*) "state value origin:",D3STATE(grid_nr,origin)
                 write(LOGUNIT,*) "state value destination:",D3STATE(grid_nr,destination)
              endif
              STDERR "Error in (scalar)flux function:negative flux !"
              call BFM_ERROR("flux","negative flux")
              if ( present(error)) error=1
            endif ! flow<0
            select case ( iiSub )
              case (iiPel)
                D3SINK(grid_nr,origin,destination)=  flow/SEC_PER_DAY
                D3SOURCE(grid_nr,destination,origin)= flow/SEC_PER_DAY
              case (iiBen)
                D2SINK(grid_nr,origin,destination)=  flow/SEC_PER_DAY
                D2SOURCE(grid_nr,destination,origin)= flow/SEC_PER_DAY
            end select
          else
            select case ( iiSub )
              case (iiPel)
                if (flow > 0.0 ) then
                  D3SOURCE(grid_nr,destination,origin)= D3SOURCE(grid_nr,destination,origin) &
                    +flow/SEC_PER_DAY
                else
                  D3SINK(grid_nr,origin,destination)= D3SINK(grid_nr,origin,destination) &
                    -flow/SEC_PER_DAY
                endif
              case (iiBen)
                if (flow > 0.0 ) then
                  D2SOURCE(grid_nr,destination,origin)= D2SOURCE(grid_nr,destination,origin) &
                    +flow/SEC_PER_DAY
                else
                  D2SINK(grid_nr,origin,destination)= D2SINK(grid_nr,origin,destination) &
                    -flow/SEC_PER_DAY
                endif
            end select
          endif
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          !END compute
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

          return
        end subroutine flux

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! vector function to get actual rate of change in the pelagic
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        function Source_D3_vector(iistate,mode)
          use constants, only: RLEN, ZERO, SEC_PER_DAY

          implicit none

          integer, intent(IN)             ::iistate
          integer, intent(IN)             ::mode
          real(RLEN) ::Source_D3_vector(size(D3SOURCE,DIM=1))

          ! Array in sum is by sum seen as 2D-array: DIM=1 and NOT 2
          select case (mode)
            case(iiConsumption)
              Source_D3_vector=sum(D3SINK(:,iistate,:),DIM=2)*SEC_PER_DAY
            case(iiTotal)
              Source_D3_vector=(sum(D3SOURCE(:,iistate,:),DIM=2)-&
                          sum(D3SINK(:,iistate,:),DIM=2))*SEC_PER_DAY
            case(iiProduction)
              Source_D3_vector=sum(D3SOURCE(:,iistate,:),DIM=2)*SEC_PER_DAY
          end select
        end function Source_D3_vector

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! vector function to get actual rate of change in the benthic
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        function Source_D2_vector(iistate,mode)
          use constants, only: RLEN, ZERO, SEC_PER_DAY

          implicit none

          integer, intent(IN) ::iistate
          integer, intent(IN) ::mode
          real(RLEN) ::Source_D2_vector(size(D2SOURCE,DIM=1))


          ! Array in sum is by sum seen as 2D-array: DIM=1 and NOT 2
          select case (mode)
            case(iiConsumption)
              Source_D2_vector=sum(D2SINK(:,iistate,:),DIM=2)*SEC_PER_DAY
            case(iiTotal)
              Source_D2_vector=(sum(D2SOURCE(:,iistate,:),DIM=2)-&
                          sum(D2SINK(:,iistate,:),DIM=2))*SEC_PER_DAY
            case(iiProduction)
              Source_D2_vector=sum(D2SOURCE(:,iistate,:),DIM=2)*SEC_PER_DAY
          end select
        end function Source_D2_vector
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! function to get actual rate of change
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        function source(iiSub,iibox,iistate)
          use constants, only: RLEN, ZERO, SEC_PER_DAY

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          ! Implicit typing is never allowed
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          implicit none

          real(RLEN)  ::Source
          integer, intent(IN)  ::iiSub
          integer, intent(IN)  ::iibox
          integer, intent(IN)  ::iistate
          if ( iiSub == iiPel )  then
            Source = (sum(D3SOURCE(iibox,iistate,:))- &
              sum(D3SINK(iibox,iistate,:)))*SEC_PER_DAY
          elseif ( iiSub == iiBen )  then
            Source = (sum(D2SOURCE(iibox,iistate,:))- &
              sum(D2SINK(iibox,iistate,:)))*SEC_PER_DAY
          endif
        end function source

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! vector to reset all rates of a state variable to zero
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        subroutine ResetSource_D3_vector(iistate)
          use constants, only: RLEN, ZERO, SEC_PER_DAY

          implicit none

          integer, intent(IN) ::iistate
          D3SOURCE(:,iistate,:)=ZERO
          D3SINK(:,iistate,:)=ZERO
        end subroutine ResetSource_D3_vector

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! vector to reset all rates of a state variable to zero
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        subroutine ResetSource_D2_vector(iistate)
          use constants, only: RLEN, ZERO, SEC_PER_DAY

          implicit none

          integer, intent(IN) ::iistate
          D2SOURCE(:,iistate,:)=ZERO
          D2SINK(:,iistate,:)=ZERO
        end subroutine ResetSource_D2_vector

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! vector to reset all rates of a state variable to zero
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        subroutine ResetSource(iiSub,iibox,iistate)
          use constants, only: RLEN, ZERO, SEC_PER_DAY

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          ! Implicit typing is never allowed
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          implicit none

          integer, intent(IN)  ::iiSub
          integer, intent(IN)  ::iibox
          integer, intent(IN)  ::iistate
          if ( iiSub == iiPel )  then
            D3SOURCE(iibox,iistate,:)=ZERO
            D3SINK(iibox,iistate,:)=ZERO
          elseif ( iiSub == iiBen )  then
            D2SOURCE(iibox,iistate,:)=ZERO
            D2SINK(iibox,iistate,:)=ZERO
          endif
        end subroutine ResetSource


        subroutine unicflux(grid_nr,iiSub,origin,destination)
        use constants, only: RLEN

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

          use global_mem, only: LOGUNIT
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          ! Implicit typing is never allowed
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          implicit none

          integer,intent(IN)    :: grid_nr
          integer,intent(IN)    :: origin
          integer,intent(IN)    :: iiSub
          integer,intent(IN)    :: destination

          real(RLEN) :: tot
          character(len=20):: type

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          !BEGIN compute
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


          select case ( iiSub )
            case (iiPel)
              type="D3"
              if ( grid_nr <=0  ) then
                tot=sum(D3SINK(:,origin,destination))
              else
                tot=D3SINK(grid_nr,origin,destination)
              endif
            case (iiBen)
              type="D2"
              if ( grid_nr <=0  ) then
                tot=sum(D2SINK(:,origin,destination))
              else
                tot=D2SINK(grid_nr,origin,destination)
              endif
            case (iiReset)
              D3SINK(:,:,:)=0.0D+00
              D2SINK(:,:,:)=0.0D+00
              return
          end select
          if ( tot > 0.0D+00  ) then
            write(LOGUNIT,'(''Double defintion '',A2,''-flux'')')type
            write(LOGUNIT,'(''origin:'',I3,'' destination:'',I3)') origin, destination
            if ( origin /= destination ) then
              STDERR 'double definition of fluxes'
              stop 1006
            endif
          endif
        !END compute
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        return
      end subroutine unicflux

      function fixed_quota_flux_vector(mode,iiSub,which,origin, &
                                          destination,flux,collect,name_routine)

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      use global_mem, only: LOGUNIT
      use constants, only: RLEN
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! Implicit typing is never allowed
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      implicit none

      integer            :: fixed_quota_flux_vector
      integer,intent(IN) :: mode
      integer,intent(IN) :: iiSub
      integer,intent(IN) :: which
      integer,intent(IN) :: origin
      integer,intent(IN) :: destination
      real(RLEN),intent(IN),dimension(:) :: flux
      real(RLEN),intent(INOUT),dimension(:) :: collect
      character(len=*),optional             :: name_routine

      character(len=32)                     :: h
      real(RLEN)                            :: r

      fixed_quota_flux_vector=0
      if ( origin> 0 .and.destination >0) then
         call flux_vector(iiSub,origin, destination,flux)
      else if ( origin > 0 ) then
         call flux_vector(iiSub,origin, origin,-flux)
      elseif ( destination > 0 ) then
         call flux_vector(iiSub,destination, destination,flux)
      elseif (iiSub < 0 ) then
         if ( mode==0)  return
         if ( sum(abs(flux)/(1.0D-80+abs(collect))-1.0D+00)> 1.0D-6) then
           !Check if we have to do with small numbers
           !If this is the case do not check!
           r=abs(sum(flux)) 
           if ( r> (1.0D-80 -r+abs(sum(collect))/1.0D-6 )) then
              h=''; if ( present(name_routine)) h=' in '//name_routine
              if ( iiSub==-iiN) then
                write(LOGUNIT,'(''Warning'',A,'': N:C quotumn not fixed'')'),trim(h)
                fixed_quota_flux_vector=1
              elseif (iiSub==-iiP) then
                write(LOGUNIT,'(''Warning'',A,'': P:C quotumn not fixed'')'),trim(h)
                fixed_quota_flux_vector=1
              endif
              return
           endif
         endif
      endif      
      if ( mode==1 ) then
        if ( (which == origin) .and.(origin.ne.destination)) then
           collect=collect-flux
        else
           collect=collect+flux
        endif
      endif
      end function fixed_quota_flux_vector

      subroutine sourcesink_flux(grid_nr,iiSub,origin,destination,flow)
            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            ! Implicit typing is never allowed
            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            implicit none

            integer,intent(IN) :: grid_nr
            integer,intent(IN) :: iiSub
            integer,intent(IN) :: origin
            integer,intent(IN) :: destination
            real(RLEN),intent(IN) :: flow

            if ( destination ==0 ) then
              call flux(grid_nr,iiSub,origin,origin,-flow)
            elseif ( origin ==0 ) then
              call flux(grid_nr,iiSub,destination,destination,flow)
            else
             call  flux(grid_nr,iiSub,origin,destination,flow)
            endif
      end subroutine sourcesink_flux

      subroutine sourcesink_flux_vector(iiSub,origin,destination,flux)
            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            ! Implicit typing is never allowed
            !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            implicit none

            integer,intent(IN) :: iiSub
            integer,intent(IN) :: origin
            integer,intent(IN) :: destination
            real(RLEN),intent(IN) :: flux(:)

            if ( destination ==0 ) then
              call flux_vector(iiSub,origin,origin,-flux)
            elseif ( origin ==0 ) then
              call flux_vector(iiSub,destination,destination,flux)
            else
             call  flux_vector(iiSub,origin,destination,flux)
            endif
      end subroutine sourcesink_flux_vector


      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! end of contain section
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

   end module mem

