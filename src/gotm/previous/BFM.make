#$Id: $
#
# Makefile to build the BFM model
#
#include $(GOTMDIR)/src/Rules.make
LIB     = $(LIBDIR)/libbio$(buildtype).a

#%.o: %.F90
#	${BFMSRC}/scripts/check_code -error -file $<
#	$(FC) $(F90FLAGS) $(EXTRA_FFLAGS) -c $< -o $@



# BFMDIR path is set in extras/bio 
# assuming that BFM is located at the same level of GOTM
#BFMSRC 	   = ../../../../bfm/src/BFM
#BFMGOTMSRC = ../../../../bfm/src/gotm
#BFMSHARE    = ../../../../bfm/src/share
BFMSRC 	   = $(BFMDIR)/src/BFM
BFMGOTMSRC = $(BFMDIR)/src/gotm
BFMSHARE    = $(BFMDIR)/src/share

DEFINES += -DBFM_GOTM
EMPTY_BFM=false
INCLUDE_NEW=true
INCLUDE_MACROPHYT=true
INCLUDE_TRACK=false
INCLUDE_PELCO2=true
INCLUDE_BENCO2=true
INCLUDE_DIAGNOS_PRF=true
INCLUDE_l1p0p=true
ifeq ($(INCLUDE_TRACK),true)
DEFINES += -DINCLUDE_TRACK
INCLUDE_PELCO2=false
INCLUDE_BENCO2=false
INCLUDE_DIAGNOS_PRF=false
endif
ifeq ($(INCLUDE_PELCO2),true)
DEFINES += -DINCLUDE_PELCO2
endif
ifeq ($(INCLUDE_BENCO2),true)
DEFINES += -DINCLUDE_BENCO2
endif
ifeq ($(INCLUDE_DIAGNOS_PRF),true)
DEFINES += -DINCLUDE_DIAGNOS_PRF
endif
ifeq ($(INCLUDE_MACROPHYT),true)
DEFINES += -DINCLUDE_MACROPHYT
endif
ifeq ($(INCLUDE_l1p0p),true)
DEFINES += -DINCLUDE_l1p0p
endif
ifeq ($(INCLUDE_NEW),true)
DEFINES += -DINCLUDE_NEW
endif


#DOCSRC	=  bio.F90 
#bio_var.F90 bio_template.F90 bio_npzd.F90 bio_iow.F90 \
#          bio_sed.F90 bio_fasham.F90 \
#          process_model.F90 ode_solvers.F90 bio_save.F90
#${LIB}(${BFMGOTMSRC}/adv_center_bfm.o)		\
#


OBJ   = \
${LIB}(bio_fluxes.o)		\
${LIB}(${BFMSRC}/General/init_cnps.o)			

ifeq ($(EMPTY_BFM),true)
	EMPTYBFMSRC=$(BFMDIR)/src/EMPTY_BFM
	OBJ   += \
	${LIB}(${BFMGOTMSRC}/gotm_error_msg.o)		\
	${LIB}(${EMPTYBFMSRC}/bio_bfm.o)		\
	${LIB}(${BFMGOTMSRC}/bfm_output.o)		\
	${LIB}(${BFMGOTMSRC}/bfm_ncdf_output.o)		\
	${LIB}(${EMPTYBFMSRC}/init_var_bfm.o)		\
	${LIB}(${EMPTYBFMSRC}/make_flux_output.o)	
else
	OBJ   += \
	${LIB}(${BFMSRC}/General/init_var_bfm.o)	
endif

OBJ   += \
${LIB}(${BFMGOTMSRC}/trace_bdy.o)		\
${LIB}(${BFMGOTMSRC}/process_model.o)		\
${LIB}(${BFMGOTMSRC}/bio_solver.o)		\
${LIB}(${BFMGOTMSRC}/adv_center_bfm.o)		\
${LIB}(${BFMGOTMSRC}/diff_center_bfm.o)		\
${LIB}(${BFMGOTMSRC}/bio.o)			\
${LIB}(${BFMGOTMSRC}/hotstart.o)

ifeq ($(EMPTY_BFM),true)

all: ${OBJ} 
	$(MOVE_MODULES_COMMAND)
else

BFM_MOD = \
	${LIB}(${BFMGOTMSRC}/bio_var.o)			\
	${LIB}(${BFMSRC}/Basis/string_functions.o)	\
	${LIB}(${BFMGOTMSRC}/gotm_error_msg.o)		\
	${LIB}(${BFMSRC}/Basis/ModuleGlobalMem.o)	\
	${LIB}(${BFMSRC}/Basis/ModuleConstants.o)	\
	${LIB}(${BFMSRC}/Basis/ModuleGlobFun.o)		\
	${LIB}(${BFMSRC}/Basis/bfm_error_msg.o)		\
	${LIB}(${BFMSRC}/Basis/ModuleMem.o)		\
	${LIB}(${BFMSRC}/Basis/ModuleParam.o)		\
	${LIB}(${BFMSRC}/Silt/wave.o)			\
	${LIB}(${BFMSRC}/Silt/spm_util.o)		\
	${LIB}(${BFMGOTMSRC}/bfm_output.o)		\
	${LIB}(${BFMGOTMSRC}/controlled_output.o)	\
	${LIB}(${BFMSRC}/Basis/BFMSinkingModule.o)	\
	${LIB}(${BFMSRC}/Basis/ModuleInterface.o)	\
	${LIB}(${BFMSRC}/Basis/PelBenInterActionModule.o)	\
	${LIB}(${BFMSRC}/Basis/CheckMassConservationModule.o)	\
	${LIB}(${BFMGOTMSRC}/bio_bfm.o)			\
	${LIB}(${BFMSRC}/General/LimitRates.o)		\
	${LIB}(${BFMSRC}/General/ModuleDiffusion.o)	\
	${LIB}(${BFMSRC}/General/RouseProfileModule.o)	\
	${LIB}(${BFMSRC}/General/NutrientControlModule.o)	\
	${LIB}(${BFMSRC}/General/CalcLim_pH.o)		\
	${LIB}(${BFMSRC}/Bennut/ModuleBenNutType.o)		\
	${LIB}(${BFMSRC}/Bennut/ModuleBenNutInterface.o)	\
	${LIB}(${BFMSRC}/Bennut/ModuleBenthicNutrient3.o)	\
	${LIB}(${BFMSRC}/Silt/ModuleSilt.o)		\
	${LIB}(${BFMSRC}/PelBen/botflux.o)		\
	${LIB}(${BFMSRC}/Ben/ModuleInfoForNextStep.o)		\
	${LIB}(${BFMSRC}/Ben/ModuleBenBac.o)		\
	${LIB}(${BFMSRC}/Ben/ModuleBenNBac.o)			\
	${LIB}(${BFMSRC}/Ben/ModuleBenOrganism.o)		\
	${LIB}(${BFMSRC}/Ben/ModuleFilterFeeder.o)		\
	${LIB}(${BFMSRC}/Ben/ModuleBioturbation.o)		\
	${LIB}(${BFMSRC}/Bennut/ModuleBenAmmonium.o)		\
	${LIB}(${BFMSRC}/Bennut/ModuleBenAnoxic.o)		\
	${LIB}(${BFMSRC}/Bennut/ModuleBenNitrate.o)		\
	${LIB}(${BFMSRC}/Bennut/ModuleBenNutConstants.o)	\
	${LIB}(${BFMSRC}/Bennut/ModuleBenNutVariables.o)	\
	${LIB}(${BFMSRC}/Oxygen/ModuleOxygen.o)			\
	${LIB}(${BFMSRC}/Bennut/ModuleBenPhosphate.o)		\
	${LIB}(${BFMSRC}/Bennut/ModuleBenSilica.o)		\
	${LIB}(${BFMSRC}/Bennut/ModuleBenthicReturn1.o)		\
	${LIB}(${BFMSRC}/Bennut/ModuleBenthicReturn2.o)		\
	${LIB}(${BFMSRC}/Bennut/BenGenericTransportModule.o)	\
	${LIB}(${BFMSRC}/Bennut/ModuleFindDepth.o)		\
	${LIB}(${BFMSRC}/Ben/TurbationDepthsModule.o)		\
	${LIB}(${BFMSRC}/CO2/CO2.o)				\
	${LIB}(${BFMSRC}/PelB/ModuleMesoZoo.o)			\
	${LIB}(${BFMSRC}/PelB/ModuleMicroZoo.o)			\
	${LIB}(${BFMSRC}/PelB/ModulePelBac.o)			\
	${LIB}(${BFMSRC}/PelB/ModulePelChem.o)			\
	${LIB}(${BFMSRC}/PelB/ModulePelGlobal.o)		\
	${LIB}(${BFMSRC}/PelB/ModulePhyto.o)			\
	${LIB}(${BFMSRC}/PelB/ModulePhaeo.o)			\
	${LIB}(${BFMSRC}/PelB/SizeRelatedProcessesModule.o)	\
	${LIB}(${BFMSRC}/Basis/ModuleBFMCodeForGetm.o)	\
	${LIB}(${BFMSRC}/PelBen/ModuleSettling.o)		\
	${LIB}(${BFMSRC}/PelBen/ModuleControlBenPartNutrientBuffers.o)	\
	${LIB}(${BFMSRC}/Ben/ModuleBenPhyto.o)			\
	${LIB}(${BFMSRC}/CO2/ModuleCO2Alk.o)			\
	${LIB}(${BFMSRC}/CO2/CO2System.o)			\

ifeq ($(INCLUDE_MACROPHYT),true)
BFM_MOD += \
	${LIB}(${BFMSRC}/PelB/ModuleMacroPhyto.o)				
endif

BFM_OBJ = \
	${LIB}(${BFMGOTMSRC}/bfm_ncdf_output.o)			\
	${LIB}(${BFMGOTMSRC}/GetDelta.o)			\
	${LIB}(${BFMGOTMSRC}/make_flux_output.o)		\
	${LIB}(${BFMSRC}/Basis/set_var_info_bfm.o)		\
	${LIB}(${BFMSRC}/Basis/check_if_in_output.o)		\
	${LIB}(${BFMSRC}/Basis/calc_sigma_depth.o)		\
	${LIB}(${BFMSRC}/Basis/Allocate3dMem.o)			\
	${LIB}(${BFMSRC}/Basis/Allocate2dMem.o)			\
	${LIB}(${BFMSRC}/Basis/AllocateOtMem.o)			\
	${LIB}(${BFMSRC}/Basis/InitBoxParams.o)			\
	${LIB}(${BFMSRC}/Basis/Ecology.o)			\
	${LIB}(${BFMSRC}/Basis/DefinitiveLossGain.o)		\
	${LIB}(${BFMSRC}/Basis/CompensateLoss.o)		\
	${LIB}(${BFMSRC}/General/Track.o)			\
	${LIB}(${BFMSRC}/General/InitTrack.o)			\
	${LIB}(${BFMSRC}/General/InitTransportStateTypes.o)	\
	${LIB}(${BFMSRC}/General/Initialize.o)			\
	${LIB}(${BFMSRC}/General/FindExtremeRates.o)		\
	${LIB}(${BFMSRC}/General/CalculateThermoVars.o)		\
	${LIB}(${BFMSRC}/General/CalcWaterProp.o)		\
	${LIB}(${BFMSRC}/General/TopPredLosses.o)		\
	${LIB}(${BFMSRC}/General/CalcAnoxicBacterialDynamics.o)	\
	${LIB}(${BFMSRC}/General/OutputAfterConservationTest.o)	\
	${LIB}(${BFMSRC}/Oxygen/WindOxReaeration_3.o)		\
	${LIB}(${BFMSRC}/Oxygen/CalcOxygenSaturation_3.o)	\
	${LIB}(${BFMSRC}/Oxygen/SetOxygenSaturation_3.o)	\
	${LIB}(${BFMSRC}/PelB/CalcChlorophylla.o)		\
	${LIB}(${BFMSRC}/PelB/CalcVerticalExtinction.o)		\
	${LIB}(${BFMSRC}/PelB/CalcLight.o)			\
	${LIB}(${BFMSRC}/PelB/MicroZoo.o)			\
	${LIB}(${BFMSRC}/PelB/MesoZoo.o)			\
	${LIB}(${BFMSRC}/PelB/MicroZoo.o)			\
	${LIB}(${BFMSRC}/PelB/PelBac.o)				\
	${LIB}(${BFMSRC}/PelB/PelChem.o)			\
	${LIB}(${BFMSRC}/PelB/PelGlobal.o)			\
	${LIB}(${BFMSRC}/PelB/Sedimentation.o)			\
	${LIB}(${BFMSRC}/PelB/PelagicSystem.o)			\
	${LIB}(${BFMSRC}/PelB/LimitNutrientUptake.o)		\
	${LIB}(${BFMSRC}/PelB/PhaeocystisCalc.o)		\
	${LIB}(${BFMSRC}/PelB/Phyto.o)				\
	${LIB}(${BFMSRC}/PelB/MixoDinoFlag.o)			\
	${LIB}(${BFMSRC}/PelBen/PelForcingForBen.o)		\
	${LIB}(${BFMSRC}/PelBen/Settling.o)			\
	${LIB}(${BFMSRC}/PelBen/WaterSedi_DetritusFlux.o)	\
	${LIB}(${BFMSRC}/PelBen/ControlBenPartNutrientBuffers.o)	\
	${LIB}(${BFMSRC}/PelBen/ResuspensionPartDetritus.o)	\
	${LIB}(${BFMSRC}/PelBen/Y3Z2Coup.o)			\
	${LIB}(${BFMSRC}/PelBen/ResuspensionBenPhyto.o)		\
	${LIB}(${BFMSRC}/PelBen/CouplingBioSilt.o)		\
	${LIB}(${BFMSRC}/PelBen/FilterLimPart.o)		\
	${LIB}(${BFMSRC}/Ben/CheckOnDryGridPoint.o)		\
	${LIB}(${BFMSRC}/Ben/BenBac.o)				\
	${LIB}(${BFMSRC}/Ben/BenNBac.o)				\
	${LIB}(${BFMSRC}/Ben/BenPhyto.o)			\
	${LIB}(${BFMSRC}/Ben/BenLimitNutrientUptake.o)		\
	${LIB}(${BFMSRC}/Ben/BenOrganism.o)			\
	${LIB}(${BFMSRC}/Ben/FastExcretion.o)			\
	${LIB}(${BFMSRC}/Ben/SizeRelatedMortY.o)		\
	${LIB}(${BFMSRC}/Ben/CompleteDepoProcesses.o)		\
	${LIB}(${BFMSRC}/Ben/DepositFeederDistribution.o)	\
	${LIB}(${BFMSRC}/Ben/GrowthStrategyY.o)			\
	${LIB}(${BFMSRC}/Ben/MaxInChannel.o)			\
	${LIB}(${BFMSRC}/Ben/PDiaInBenDynamics.o)		\
	${LIB}(${BFMSRC}/Ben/CalcSizeClassesYD.o)		\
	${LIB}(${BFMSRC}/Ben/BenthicSystem.o)			\
	${LIB}(${BFMSRC}/Ben/BenGlobal.o)			\
	${LIB}(${BFMSRC}/Ben/Bioturbation.o)			\
	${LIB}(${BFMSRC}/Ben/FilterFeeder.o)			\
	${LIB}(${BFMSRC}/Ben/CalcZeroOrderOxInLayer.o)		\
	${LIB}(${BFMSRC}/Ben/KeepAgeInfo.o)		\
	${LIB}(${BFMSRC}/Bennut/CalcBurialDetritus.o)		\
	${LIB}(${BFMSRC}/Bennut/RedistributeBenBacteria.o)	\
	${LIB}(${BFMSRC}/Bennut/BenGeneric1LayerProfile.o)	\
	${LIB}(${BFMSRC}/Bennut/BenAmmonium.o)			\
	${LIB}(${BFMSRC}/Bennut/BenAnoxic.o)			\
	${LIB}(${BFMSRC}/Bennut/BenDenitriDepth.o)		\
	${LIB}(${BFMSRC}/Bennut/BenNitrate.o)			\
	${LIB}(${BFMSRC}/Oxygen/BenOxygen.o)			\
	${LIB}(${BFMSRC}/Oxygen/CalcAvailOxInAerobic.o)		\
	${LIB}(${BFMSRC}/Bennut/BenPhosphate.o)			\
	${LIB}(${BFMSRC}/Bennut/BenDetritusTransport.o)		\
	${LIB}(${BFMSRC}/Bennut/BenSilica.o)			\
	${LIB}(${BFMSRC}/Bennut/BenthicNutrient2.o)		\
	${LIB}(${BFMSRC}/Bennut/BenthicNutrient3.o)		\
	${LIB}(${BFMSRC}/Bennut/ForcingBenthicNutrientModel.o)	\
	${LIB}(${BFMSRC}/Bennut/ForcingLocModel.o)		\
	${LIB}(${BFMSRC}/Bennut/BenthicReturn1.o)		\
	${LIB}(${BFMSRC}/Bennut/BenthicReturn2.o)		\
	${LIB}(${BFMSRC}/Bennut/BenProfiles.o)			\
	${LIB}(${BFMSRC}/Bennut/flux_at_deep_boundary.o)	\
	${LIB}(${BFMSRC}/Bennut/CalculateFromSet.o)		\
	${LIB}(${BFMSRC}/Bennut/CalculateMassFromSet_vector.o)	\
	${LIB}(${BFMSRC}/Bennut/CalculateSet.o)			\
	${LIB}(${BFMSRC}/Bennut/CalculateShift.o)		\
	${LIB}(${BFMSRC}/Bennut/FixProportionCoeff.o)		\
	${LIB}(${BFMSRC}/Bennut/CalculateTau.o)			\
	${LIB}(${BFMSRC}/Bennut/CompleteSet.o)			\
	${LIB}(${BFMSRC}/Bennut/InitializeSet.o)		\
	${LIB}(${BFMSRC}/Bennut/DefineSet.o)			\
	${LIB}(${BFMSRC}/Bennut/GetInfoFromSet.o)		\
	${LIB}(${BFMSRC}/Bennut/GetIntegerFromSet.o)		\
	${LIB}(${BFMSRC}/Bennut/GetInfoFromSet_vector.o)	\
	${LIB}(${BFMSRC}/Bennut/FindOptDepth.o)			\
	${LIB}(${BFMSRC}/Bennut/PrintSet.o)			\
	${LIB}(${BFMSRC}/Bennut/bess_exp.o)			\
	${LIB}(${BFMSRC}/Bennut/bessi0.o)			\
	${LIB}(${BFMSRC}/Bennut/bessi1.o)			\
	${LIB}(${BFMSRC}/Bennut/bessk0.o)			\
	${LIB}(${BFMSRC}/Bennut/bessk1.o)			\
	${LIB}(${BFMSRC}/Bennut/calculate_equation.o)		\
	${LIB}(${BFMSRC}/Bennut/FindLayerNr.o)			\
	${LIB}(${BFMSRC}/Bennut/funcalc.o)			\
	${LIB}(${BFMSRC}/Bennut/input_para.o)			\
	${LIB}(${BFMSRC}/Bennut/kfind.o)			\
	${LIB}(${BFMSRC}/Bennut/lubksb.o)			\
	${LIB}(${BFMSRC}/Bennut/ludcmp.o)			\
	${LIB}(${BFMSRC}/Bennut/CalculateFromLayer.o)		\
	${LIB}(${BFMSRC}/Bennut/qgaus_exp.o)			\
	${LIB}(${BFMSRC}/Bennut/re_store.o)			\
	${LIB}(${BFMSRC}/Bennut/set_max_sing.o)			\
	${LIB}(${BFMSRC}/Bennut/svbksb.o)			\
	${LIB}(${BFMSRC}/Bennut/svdcmp.o)			\
	${LIB}(${BFMSRC}/Bennut/transfer_term.o)		\
	${LIB}(${BFMSRC}/CO2/Alkalinity.o)			\
	${LIB}(${BFMSRC}/CO2/SurfaceCO2Processes.o)		\
	${LIB}(${BFMSRC}/CO2/PelCO2.o)				\
	${LIB}(${BFMSRC}/CO2/CalcCO2SatInField.o)		\
	${LIB}(${BFMSRC}/CO2/BenCO2Transport.o)			\
	${LIB}(${BFMSRC}/CO2/BenCO2Profiles.o)			\
	${LIB}(${BFMSRC}/CO2/BenpH.o)				\
	${LIB}(${BFMSRC}/CO2/BenAlkalinity.o)			\
	${LIB}(${BFMSRC}/Silt/Silt.o)				\
	${LIB}(${BFMSRC}/Silt/CalcSiltResuspension.o)

ifeq ($(INCLUDE_MACROPHYT),true)
BFM_OBJ += \
	${LIB}(${BFMSRC}/PelB/MacroPhyto.o)				
endif

ifeq ($(INCLUDE_l1p0p),true)
BFM_OBJ += \
	${LIB}(${BFMSRC}/General/Calculate1procentLight.o)				
endif

all: ${BFM_MOD} ${OBJ} ${BFM_OBJ}
	$(MOVE_MODULES_COMMAND)
	ls -l --time-style=+%ssec ${BFMSRC}/Basis/ModuleMem.F90 >\
 					 $(GOTMDIR)/src/mem_timestamp

$(BFM_MOD) : $(BFMSRC)/Basis/ModuleMem.F90

${BFMSRC}/Basis/ModuleMem.F90 : $(BFMSRC)/Basis/GlobalDefsBFM.model
	${BFMSRC}/scripts/GenerateGlobalBFMF90Code  $(DEFINES) \
	-read ${BFMSRC}/Basis/GlobalDefsBFM.model \
	-from ${BFMSRC}/proto -to ${BFMSRC}/Basis \
	-actions statemem alloc3dmem alloc2dmem allocOtmem netcdfmem \
	-to ${BFMSRC}/include -actions headermem  
endif

#------------------------------------------------------------------------------
# Copyright (C) 2006 - the GOTM-team and the BFM-team
#-----------------------------------------------------------------------
