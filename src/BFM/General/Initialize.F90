!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50-g
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: Initialize
!
! DESCRIPTION
!   Initialization of model
!   Allocation of memory for variables, reading of data files 

!   This file is generated directly from OpenSesame model code, using a code 
!   generator which transposes from the sesame meta language into F90.
!   F90 code generator written by P. Ruardij.
!   structure of the code based on ideas of M. Vichi.
!
! !INTERFACE
  SUBROUTINE Initialize
!
! USES:
  use global_mem, only: LOGUNIT
  use mem, only: InitializeModel,ppMicroZooplankton,ppMesoZooPlankton, &
    MicroZooplankton,MesoZooPlankton,iiMicroZooplankton,iiMesoZooPlankton, &
    iiBenPhyto,iiYy3,NO_BOXES,iiN,iiP,qpZc,qnZc,qp_mz,qn_mz,sw_CalcPhyto
!JM  use mem_Param
  use mem_Param, only: CalcBenthicFlag,CalcMesoZooPlankton,CalcSuspensionFeeders,CalcPhytoPlankton,CalcBenPhyto,InitParam
  use mem_Diffusion,ONLY:InitDiffusion
  use mem_WindOxReaeration_3,ONLY:InitWindOxReaeration_3
  use mem_PelGlobal,ONLY:InitPelGlobal
  use mem_PelChem,ONLY:InitPelChem
  use mem_PelBac,ONLY:InitPelBac
  use mem_MesoZoo,ONLY:InitMesoZoo
  use mem_MicroZoo,ONLY:InitMicroZoo
  use mem_Phyto,ONLY:InitPhyto
  use mem_Phaeo,ONLY:InitPhaeo
! use mem_PhotoAvailableRadiation,ONLY:InitPhotoAvailableRadiation
! use mem_LightAdaptation,ONLY:InitLightAdaptation
  use mem_BenOrganism,ONLY:InitBenOrganism
  use mem_FilterFeeder,ONLY:InitFilterFeeder
  use mem_BenBac,ONLY:InitBenBac
  use mem_BenNBac,ONLY:InitBenNBac
  use mem_Bioturbation,ONLY:InitBioturbation
  use mem_BenthicReturn1,ONLY:InitBenthicReturn1
  use mem_BenthicReturn2,ONLY:InitBenthicReturn2
  use mem_BenthicNutrient3,ONLY:InitBenthicNutrient3
  use mem_BenAmmonium,ONLY:InitBenAmmonium
  use mem_BenNitrate,ONLY:InitBenNitrate
  use mem_BenOxygen,ONLY:InitBenOxygen
  use mem_BenAnoxic,ONLY:InitBenAnoxic
  use mem_BenDenitriDepth,ONLY:InitBenDenitriDepth
  use mem_BenPhosphate,ONLY:InitBenPhosphate
  use mem_BenSilica,ONLY:InitBenSilica
  use mem_BenQ1Transport,ONLY:InitBenQ1Transport
  use mem_Settling,ONLY:InitSettling
  use mem_ControlBenPartNutrientBuffers,ONLY:InitControlBenPartNutrientBuffers
  use mem_CO2,ONLY:InitCO2
#ifdef INCLUDE_MACROPHYT
  use mem_MacroPhyto,ONLY:InitMacroPhyto,surface
#endif
#ifdef INCLUDE_DAAN
  use mem_Daan,ONLY:InitDaan
#endif
  use mem_BenCO2Transport,ONLY:InitBenCO2Transport
  use mem_BenAlkalinity,ONLY:InitBenAlkalinity
  use mem_Silt,ONLY:InitSilt
  use mem_BenPhyto,ONLY:InitBenPhyto,p_useparams
  use botflux,only:initbotflux

!  
!
! !AUTHORS
!   mfstep ERSEM team
!
! !REVISION_HISTORY
!   ---
!
! COPYING
!   
!   Copyright (C) 2006 P. Ruardij, the mfstep group, the ERSEM team 
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
!EOP
!-------------------------------------------------------------------------!
!BOC
!
!
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      implicit none
      integer             ::i
      integer             ::j


      InitializeModel=0
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! Allocate Memory for All global variables
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! Read all data files:(namelist files)
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      write(LOGUNIT,*) 'BEGIN of READING PARMETERS'
      call InitParam
      call InitDiffusion
      call InitWindOxReaeration_3
      call InitPelGlobal
      call InitPelChem
      call InitPelBac
      call InitMesoZoo
      call InitMicroZoo
      call InitPhyto
      call InitPhaeo
!     call InitPhotoAvailableRadiation
!     call InitLightAdaptation
      call InitBenOrganism
      call InitFilterFeeder
      call InitBenBac
      call InitBenNBac
      call InitBioturbation
      call InitBenthicReturn1
      call InitBenthicReturn2
      call InitBenthicNutrient3
      call InitBenAmmonium
      call InitBenNitrate
      call InitBenOxygen
      call InitBenAnoxic
      call InitBenDenitriDepth
      call InitBenPhosphate
      call InitBenSilica
      call InitBenQ1Transport
      call InitSettling
      call InitControlBenPartNutrientBuffers
#ifdef INCLUDE_MACROPHYT
      if (CalcMacroPhyto) then 
          call InitMacroPhyto
          surface=p_surface_1d
      endif
#endif
#ifdef INCLUDE_PELCO2
#ifdef INCLUDE_BENCO2
      call InitBenCO2Transport
      call InitBenAlkalinity
#endif
      call InitCO2
#endif
#ifdef INCLUDE_DAAN
      call InitDaan
#endif
      call InitSilt
      call initbotflux


!     if (iiBenPhyto==1.and.CalcBenPhyto(iiBenPhyto)) call InitBenPhyto
      call InitBenPhyto
      write(LOGUNIT,*) 'END of READING PARMETERS'
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      ! Read all other Init* files
      !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      call InitTransportStateTypes
      call CoupleInfoOnNSinkingToGotm
      
      !! Set/Reset Params which depend on each others:
      !! No filterfeeder larvae with young filterfeeders:
      CalcMesoZooPlankton(3)=CalcMesoZooPlankton(3) &
             .and.CalcBenthicFlag>0.and. CalcSuspensionFeeders(iiYy3)
      ! No resuspended diatoms with benthic diatoms:
      do i = 1 , iiBenPhyto
        j=p_useparams(i)
        CalcPhytoPlankton(j)=CalcPhytoPlankton(j) &
                      .and.CalcBenthicFlag>0.and. CalcBenPhyto(i)
      enddo

       sw_CalcPhyto=1
    END SUBROUTINE
!BOP
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model version 2.50
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
