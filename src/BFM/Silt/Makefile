## Makefile to build the BFM standalone
##
## uses makedepf90 to build dependencies "on the fly"
## http://www.helsinki.fi/~eedelman/makedepf90.html

include ./Rules.make

EXE	= bfm_standalone$(buildtype)
MAINOBJ	= ./standalone/standalone_main.o 
LIB	= $(LIBDIR)/libbfm$(buildtype).a
BFMSRC = ./BFM

BFM_MOD = \
	${LIB}(${BFMSRC}/General/ModuleConstants.o)			\
	${LIB}(${BFMSRC}/General/ModuleGlobFun.o)			\
	${LIB}(${BFMSRC}/General/bfm_error_msg.o)			\
	${LIB}(${BFMSRC}/General/ModuleMem.o)				\
	${LIB}(${BFMSRC}/General/ModuleParam.o)				\
	${LIB}(${BFMSRC}/General/ModuleInterface.o)			\
	${LIB}(${BFMSRC}/Light/ModuleLightAdaptation.o)			\
	${LIB}(${BFMSRC}/Light/ModulePhotoAvailableRadiation.o)		\
	${LIB}(${BFMSRC}/Oxygen/ModuleWindOxReaeration_3.o)		\
	${LIB}(${BFMSRC}/Bennut/ModuleBenthicNutrient3.o)		\
	${LIB}(${BFMSRC}/Bennut/ModuleBenAmmonium.o)			\
	${LIB}(${BFMSRC}/Bennut/ModuleBenAnoxic.o)			\
	${LIB}(${BFMSRC}/Bennut/ModuleBenDenitriDepth.o)		\
	${LIB}(${BFMSRC}/Bennut/ModuleBenNitrate.o)			\
	${LIB}(${BFMSRC}/Bennut/ModuleBenNutConstants.o)		\
	${LIB}(${BFMSRC}/Bennut/ModuleBenNutType.o)			\
	${LIB}(${BFMSRC}/Bennut/ModuleBenNutInterface.o)		\
	${LIB}(${BFMSRC}/Bennut/ModuleBenNutVariables.o)		\
	${LIB}(${BFMSRC}/Bennut/ModuleBenOxygen.o)			\
	${LIB}(${BFMSRC}/Bennut/ModuleBenPhosphate.o)			\
	${LIB}(${BFMSRC}/Bennut/ModuleBenQ1Transport.o)			\
	${LIB}(${BFMSRC}/Bennut/ModuleBenSilica.o)			\
	${LIB}(${BFMSRC}/Bennut/ModuleBenthicReturn1.o)			\
	${LIB}(${BFMSRC}/Bennut/ModuleBenthicReturn2.o)			\
	${LIB}(${BFMSRC}/Bennut/LimitRates.o)	 			\
	${LIB}(${BFMSRC}/CO2/CO2.o)			        \
	${LIB}(${BFMSRC}/PelB/ModuleMesoZoo.o)				\
	${LIB}(${BFMSRC}/PelB/ModuleMicroZoo.o)				\
	${LIB}(${BFMSRC}/PelB/ModulePelBac.o)				\
	${LIB}(${BFMSRC}/PelB/ModulePelChem.o)				\
	${LIB}(${BFMSRC}/PelB/ModulePelGlobal.o)			\
	${LIB}(${BFMSRC}/PelB/ModulePhyto.o)				\
	${LIB}(${BFMSRC}/PelB/ModulePhaeo.o)				\
	${LIB}(${BFMSRC}/PelBen/ModuleSettling.o)			\
	${LIB}(${BFMSRC}/PelBen/ModuleControlBenPartNutrientBuffers.o)	\
	${LIB}(${BFMSRC}/Ben/ModuleBenPhyto.o)				\
	${LIB}(${BFMSRC}/Ben/ModuleBenOrganism.o)			\
	${LIB}(${BFMSRC}/Ben/ModuleBenBac.o)				\
	${LIB}(${BFMSRC}/Ben/ModuleBenNBac.o)				\
	${LIB}(${BFMSRC}/Ben/ModuleFilterFeeder.o)			\
	${LIB}(${BFMSRC}/Ben/ModuleBioturbation.o)			\
	${LIB}(${BFMSRC}/CO2/CO2System.o)			        \
	${LIB}(${BFMSRC}/CO2/ModuleBenCO2Transport.o)		        \
	${LIB}(${BFMSRC}/CO2/ModuleBenAlkalinity.o)		        \
	${LIB}(${BFMSRC}/CO2/ModuleCO2.o)				\
 	${LIB}(${BFMSRC}/Silt/wave.o)         			        \
 	${LIB}(${BFMSRC}/Silt/spm_util.o)      			        \ 
	${LIB}(${BFMSRC}/Silt/ModuleSilt.o)        

	
BFM_OBJ = \
	${LIB}(${BFMSRC}/General/Track.o)				\
        ${LIB}(${BFMSRC}/General/InitTrack.o)				\
	${LIB}(${BFMGOTMSRC}/GetDelta.o)				\
	${LIB}(${BFMGOTMSRC}/D3toD1.o)					\
	${LIB}(${BFMGOTMSRC}/D2toD1.o)					\
	${LIB}(${BFMGOTMSRC}/make_flux_output.o)			\
	${LIB}(${BFMSRC}/General/set_var_info_bfm.o)			\
	${LIB}(${BFMSRC}/General/check_if_in_output.o)			\
	${LIB}(${BFMSRC}/General/AllocateMem.o)				\
	${LIB}(${BFMSRC}/General/InitBoxParams.o)			\
	${LIB}(${BFMSRC}/General/Initialize.o)				\
	${LIB}(${BFMSRC}/General/InitTransportStateTypes.o)		\
	${LIB}(${BFMSRC}/General/FindNaNInRates.o)			\
	${LIB}(${BFMSRC}/General/Ecology.o)				\
	${LIB}(${BFMSRC}/General/eTq.o)					\
	${LIB}(${BFMSRC}/General/CheckD2Rates.o)			\
	${LIB}(${BFMSRC}/General/CheckMassConservationNPS.o)		\
	${LIB}(${BFMSRC}/General/CheckMassConservationC.o) 	        \
	${LIB}(${BFMSRC}/General/CalcBndyConcentration.o)		\
	${LIB}(${BFMSRC}/General/CalcRiverConcentration.o)		\
	${LIB}(${BFMSRC}/General/ControlNutsBdy.o)			\
	${LIB}(${BFMSRC}/General/CalculateThermoVars.o)			\
	${LIB}(${BFMSRC}/Light/LightAdaptation.o)			\
	${LIB}(${BFMSRC}/Light/PhotoAvailableRadiation.o)		\
	${LIB}(${BFMSRC}/Oxygen/WindOxReaeration_3.o)			\
	${LIB}(${BFMSRC}/Oxygen/SurfaceO2Diffusion.o) 			\
	${LIB}(${BFMSRC}/Oxygen/CalcOxygenSaturation_3.o)		\
	${LIB}(${BFMSRC}/Oxygen/CalcSchmidtNumberOx.o)			\
	${LIB}(${BFMSRC}/PelB/CalcChlorophylla.o)			\
	${LIB}(${BFMSRC}/PelB/CalcVerticalExtinction.o)			\
	${LIB}(${BFMSRC}/PelB/MicroZoo.o)				\
	${LIB}(${BFMSRC}/PelB/MesoZoo.o)				\
	${LIB}(${BFMSRC}/PelB/MicroZoo.o)				\
	${LIB}(${BFMSRC}/PelB/PelBac.o)					\
	${LIB}(${BFMSRC}/PelB/PelChem.o)				\
	${LIB}(${BFMSRC}/PelB/PelGlobal.o)				\
	${LIB}(${BFMSRC}/PelB/PelagicSystem.o)				\
	${LIB}(${BFMSRC}/PelB/LimitNutrientUptake.o)			\
	${LIB}(${BFMSRC}/PelB/PhaeocystisCalc.o)			\
	${LIB}(${BFMSRC}/PelB/Phyto.o)					\
	${LIB}(${BFMSRC}/Ben/ResetTotMassVar.o)				\
	${LIB}(${BFMSRC}/PelBen/BentoPelCoup.o)				\
	${LIB}(${BFMSRC}/PelBen/PelForcingForBen.o)			\
	${LIB}(${BFMSRC}/PelBen/Sedimentation.o)			\
	${LIB}(${BFMSRC}/PelBen/Settling.o)				\
	${LIB}(${BFMSRC}/PelBen/ControlBenPartNutrientBuffers.o)	\
	${LIB}(${BFMSRC}/PelBen/RecalcPenetrationDepth.o)		\
	${LIB}(${BFMSRC}/PelBen/Y3Z2Coup.o)				\
	${LIB}(${BFMSRC}/PelBen/ResuspensionBenPhyto.o)			\
	${LIB}(${BFMSRC}/Ben/BenBac.o)					\
	${LIB}(${BFMSRC}/Ben/BenNBac.o)					\
	${LIB}(${BFMSRC}/Ben/BenPhyto.o)				\
	${LIB}(${BFMSRC}/Ben/BenLimitNutrientUptake.o)			\
	${LIB}(${BFMSRC}/Ben/BenOrganism.o)				\
	${LIB}(${BFMSRC}/Ben/BenthicSystem.o)				\
	${LIB}(${BFMSRC}/Ben/BenGlobal.o)				\
	${LIB}(${BFMSRC}/Ben/Bioturbation.o)				\
	${LIB}(${BFMSRC}/Ben/FilterFeeder.o)				\
	${LIB}(${BFMSRC}/Ben/ResetTotMassVar.o)				\
	${LIB}(${BFMSRC}/Ben/CorrectConcNearBed.o)			\
	${LIB}(${BFMSRC}/Bennut/BenAmmonium.o) 				\
	${LIB}(${BFMSRC}/Bennut/BenAnoxic.o) 				\
	${LIB}(${BFMSRC}/Bennut/BenDenitriDepth.o) 			\
	${LIB}(${BFMSRC}/Bennut/BenNitrate.o) 				\
	${LIB}(${BFMSRC}/Bennut/BenNitrogenShifting.o) 			\
	${LIB}(${BFMSRC}/Bennut/BenOxygen.o) 				\
	${LIB}(${BFMSRC}/Bennut/BenPhosphate.o) 			\
	${LIB}(${BFMSRC}/Bennut/BenQ1Transport.o) 			\
	${LIB}(${BFMSRC}/Bennut/BenSilica.o) 				\
	${LIB}(${BFMSRC}/Bennut/BenthicNutrient2.o) 			\
	${LIB}(${BFMSRC}/Bennut/BenthicNutrient3.o) 			\
	${LIB}(${BFMSRC}/Bennut/BenthicReturn1.o) 			\
	${LIB}(${BFMSRC}/Bennut/BenthicReturn2.o) 			\
	${LIB}(${BFMSRC}/Bennut/InitBenthicNutrient.o) 			\
	${LIB}(${BFMSRC}/Bennut/InitBenthicNutrient3.o) 		\
	${LIB}(${BFMSRC}/Bennut/BenProfiles.o) 		\
	${LIB}(${BFMSRC}/Bennut/CalculateFromSet.o) 			\
	${LIB}(${BFMSRC}/Bennut/CalculateSet.o) 			\
	${LIB}(${BFMSRC}/Bennut/CalculateShift.o) 			\
	${LIB}(${BFMSRC}/Bennut/FixProportionCoeff.o)			\
	${LIB}(${BFMSRC}/Bennut/CalculateTau.o) 			\
	${LIB}(${BFMSRC}/Bennut/CompleteSet.o) 				\
	${LIB}(${BFMSRC}/Bennut/InitializeSet.o) 			\
	${LIB}(${BFMSRC}/Bennut/DefineSet.o) 				\
	${LIB}(${BFMSRC}/Bennut/GetInfoFromSet.o) 			\
	${LIB}(${BFMSRC}/Bennut/PrintSet.o) 				\
	${LIB}(${BFMSRC}/Bennut/bess_exp.o) 				\
	${LIB}(${BFMSRC}/Bennut/bessi0.o) 				\
	${LIB}(${BFMSRC}/Bennut/bessi1.o) 				\
	${LIB}(${BFMSRC}/Bennut/bessk0.o) 				\
	${LIB}(${BFMSRC}/Bennut/bessk1.o) 				\
	${LIB}(${BFMSRC}/Bennut/calculate_equation.o) 			\
	${LIB}(${BFMSRC}/Bennut/calculatelayer.o) 			\
	${LIB}(${BFMSRC}/Bennut/funcalc.o) 				\
	${LIB}(${BFMSRC}/Bennut/input_para.o) 				\
	${LIB}(${BFMSRC}/Bennut/kfind.o) 				\
	${LIB}(${BFMSRC}/Bennut/lubksb.o) 				\
	${LIB}(${BFMSRC}/Bennut/ludcmp.o) 				\
	${LIB}(${BFMSRC}/Bennut/manage_coeff.o) 			\
	${LIB}(${BFMSRC}/Bennut/noutput.o) 				\
	${LIB}(${BFMSRC}/Bennut/qgaus_exp.o) 				\
	${LIB}(${BFMSRC}/Bennut/re_store.o) 				\
	${LIB}(${BFMSRC}/Bennut/set_max_sing.o) 			\
	${LIB}(${BFMSRC}/Bennut/svbksb.o) 				\
	${LIB}(${BFMSRC}/Bennut/svdcmp.o) 				\
	${LIB}(${BFMSRC}/Bennut/transfer.o) 				\
	${LIB}(${BFMSRC}/CO2/Alkalinity.o) 				\
	${LIB}(${BFMSRC}/CO2/PelCO2.o) 					\
	${LIB}(${BFMSRC}/CO2/CalcSchmidtNumberCO2.o) 			\
	${LIB}(${BFMSRC}/CO2/SurfaceCO2Processes.o) 			\
	${LIB}(${BFMSRC}/CO2/SurfaceCO2Diffusion.o) 			\
	${LIB}(${BFMSRC}/CO2/CalcCO2SatInField.o) 			\
	${LIB}(${BFMSRC}/CO2/BenCO2Transport.o) 			\
	${LIB}(${BFMSRC}/CO2/BenCO2Profiles.o)				\
	${LIB}(${BFMSRC}/CO2/BenpH.o) 					\
	${LIB}(${BFMSRC}/CO2/BenAlkalinity.o) 			        \
 	${LIB}(${BFMSRC}/Silt/Silt.o)        


APIBFM_MOD = \
	${LIB}(./share/api_bfm.o)			\
	${LIB}(./share/netcdf_bfm.o)			\
	${LIB}(./share/string_functions.o)		\
	${LIB}(./share/init_var_bfm.o)			\
	${LIB}(./share/init_cnps.o)			\
	${LIB}(./share/calcmean_bfm.o)			\
	${LIB}(${BFMSRC}/General/set_var_info_bfm.o)

STDALONE_MOD = \
	${LIB}(./standalone/time.o)	\
	${LIB}(./standalone/standalone.o)	\
	${LIB}(./standalone/integration-efw.o)	\
	${LIB}(./standalone/integration-lf.o)	\
	${LIB}(./standalone/integration-RK2.o)	\
	${LIB}(./standalone/D3toD1.o)		\
	${LIB}(./standalone/D2toD1.o)		\
	${LIB}(./standalone/ResetFluxes.o)	\
	${LIB}(./standalone/GetDelta.o)


all: ${BFM_MOD} ${APIBFM_MOD} ${STDALONE_MOD} ${BFM_OBJ} ${EXE}
	$(MOVE_MODULES_COMMAND)

$(EXE): $(MAINOBJ) 
	$(FC) $(FFLAGS) -o ../bin/$@ $(MAINOBJ) $(LIB) -L$(NETCDFLIBDIR) -lnetcdf

$(BFM_MOD) : $(BFMSRC)/General/ModuleMem.F90

${BFMSRC}/General/ModuleMem.F90 : $(BFMSRC)/General/GlobalDefsBFM.model
	${BFMSRC}/scripts/GenerateGlobalBFMF90Code  -read ${BFMSRC}/General/GlobalDefsBFM.model \
		-from ${BFMSRC}/proto -to ${BFMSRC}/General -actions statemem allocmem netcdfmem \
		-to ${BFMSRC}/include -actions headermem

distclean: clean

clean:
	$(RM) ../modules/$(FORTRAN_COMPILER)/*.{m,mod}
	$(RM) ../lib/$(FORTRAN_COMPILER)/lib*$(buildtype).a
	$(RM) ./share/*.o
	$(RM) ./standalone/*.o
	$(RM) ${BFMSRC}/*/*.o

depend:
	makedepf90 ${BFMSRC}/*/*.?90 ./share/*.?90 ./standalone/*.?90 > .depend

include .depend

