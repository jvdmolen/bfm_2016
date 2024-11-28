cmake_minimum_required(VERSION 3.0)

project(bfmgotm VERSION 20180626 LANGUAGES Fortran)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")

list(APPEND DEFINES -DNOT_STANDALONE -DBFM_GOTM )
include(IncludeOptions.txt OPTIONAL)
add_definitions( ${DEFINES} ) 

message(STATUS "BFM CMake: GETM_USE_BFM")
message(STATUS "BFM CMake: GETM_USE_BFM=${GETM_USE_BFM}")
if(NOT GETM_USE_BFM)
  option(GOTM_USE_NetCDF "Enable output in NetCDF format" ON)
  if(GOTM_USE_NetCDF)
    find_package(NetCDF REQUIRED)
    add_definitions(-DNETCDF_FMT -DREAL_4B=real\(4\))
message(STATUS "!!!!! NetCDF_INCLUDE_DIRS=${NetCDF_INCLUDE_DIRS}")
    include_directories("${NetCDF_INCLUDE_DIRS}")
    if (NetCDF_STATIC_MSVC_BUILD)
      # On Windows with a statically-compiled NetCDF library-compile all code against static runtime.
      # This MUST be done before any targets are added.
      add_compile_options("/libs:static")
    endif()
  endif(GOTM_USE_NetCDF)
endif()

set( BFMSRC "${CMAKE_CURRENT_SOURCE_DIR}/BFM" )
set( BFMGOTM "${CMAKE_CURRENT_SOURCE_DIR}/gotm" )
#JM set( GOTM "/home/rua/GOTM_SOURCES/gotm-4.1.0-wwz/src")

message(STATUS "GOTM_USE_BFM=${GOTM_USE_BFM}")
message(STATUS "DEFINES:${DEFINES}")
string (REPLACE ";" " " DEFINES "${DEFINES}")
message(STATUS "DEFINES:${DEFINES}")

#JM get_target_property(OUTVAR target COMPILE_DEFINITIONS)

# Global include directories
include_directories(BEFORE
                    "${PROJECT_SOURCE_DIR}/../include"
                    "${BFMSRC}/include"
                   )
#JM if (NOT GOTM_USE_BFM )
#JM add_library(gotm_basis
#JM              ${GOTM}/turbulence/turbulence.F90
#JM              ${GOTM}/util/util.F90
#JM             )
#JM endif()
message(STATUS "CMAKE_CURRENT_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}  ${CMAKE_CURRENT_BINARY_DIR}/../modules")

add_library(prebasis 
            ${BFMSRC}/General/ModuleGlobalMem.f90
            ${BFMSRC}/General/bfm_error_msg.F90
            ${BFMGOTM}/gotm_error_msg.F90
            ${BFMSRC}/General/ModuleConstants.f90
#JM            ${BFMSRC}/Basis/calc_sigma_depth.F90
          )
set(bfm_all ${bfm_all} prebasis)

add_custom_command(
    OUTPUT
       ${BFMSRC}/General/ModuleMem.F90
    DEPENDS
         ${BFMSRC}/General/GlobalDefsBFM.model
    COMMAND
#JM        ${BFMSRC}/scripts/GenerateGlobalBFMF90Code "${DEFINES}" -read ${BFMSRC}/General/GlobalDefsBFM.model -from ${BFMSRC}/proto -to ${BFMSRC}/General -actions statemem alloc3dmem alloc2dmem allocOtmem netcdfmem -to ${BFMSRC}/include -actions headermem
        ${BFMSRC}/scripts/GenerateGlobalBFMF90Code "${DEFINES}" -read ${BFMSRC}/General/GlobalDefsBFM.model -from ${BFMSRC}/proto -to ${BFMSRC}/General -actions statemem allocmem netcdfmem -to ${BFMSRC}/include -actions headermem
    COMMENT
    "Run GenerateGlobalBFMF90Code script"
    VERBATIM
    )

add_library(global_mod 
            ${BFMSRC}/General/ModuleMem.F90
            ${BFMSRC}/General/ModuleParam.F90
            ${BFMSRC}/General/ModuleGlobFun.f90
            ${BFMSRC}/General/ModuleInterface.F90
            ${BFMGOTM}/bio_var.F90
            ${BFMSRC}/share/string_functions.F90
            ${BFMSRC}/share/api_bfm.F90
#            ${BFMSRC}/share/init_var_bfm.F90
#            ${BFMSRC}/share/calcmean_bfm.F90
            ${BFMSRC}/share/init_cnps.F90
#            ${BFMSRC}/share/netcdf_bfm.F90
            ${BFMGOTM}/bfm_output.F90
            ${BFMGOTM}/controlled_output.F90
            ${BFMGOTM}/GetDelta.F90
            ${BFMGOTM}/diff_center_bfm.F90
            ${BFMGOTM}/calendar_date.F90
            ${BFMGOTM}/adv_center_bfm.F90
            ${BFMGOTM}/D2toD1.f90
            ${BFMGOTM}/D3toD1.f90
            ${BFMGOTM}/trace_bdy.F90    #JM added for coupling to getm
            ${BFMGOTM}/make_flux_output.F90
            ${BFMGOTM}/bio_bfm.F90
	    ${BFMSRC}/Ben/ModuleBenNBac.f90 #AN added
	    ${BFMSRC}/Ben/ModuleBenPhyto.F90 #JM added
	    ${BFMSRC}/Ben/ModuleFilterFeeder.f90 #JM added
	    ${BFMSRC}/Ben/ModuleBioturbation.f90 #AN added
	    ${BFMSRC}/General/ModuleDiffusion.f90 #AN added
	    ${BFMSRC}/Bennut/LimitRates.F90 #AN added
#JM	    ${BFMSRC}/Basis/BFMSinkingModule.F90 #AN added
#JM	    ${BFMSRC}/Basis/PelBenInterActionModule.F90 #AN added
#JM	    ${BFMSRC}/General/RouseProfileModule.F90 #AN added
	    ${BFMSRC}/Bennut/ModuleBenthicNutrient3.F90 #AN added
	    ${BFMSRC}/Ben/ModuleBenBac.f90 #AN added
	    ${BFMSRC}/Bennut/ModuleBenNitrate.f90 #AN added
	    ${BFMSRC}/Bennut/ModuleBenAmmonium.f90 #JM added
	    ${BFMSRC}/Bennut/ModuleBenPhosphate.f90 #JM added
	    ${BFMSRC}/Bennut/ModuleBenNutType.f90 #AN added
	    ${BFMSRC}/Bennut/ModuleBenNutInterface.f90 #AN added
	    ${BFMSRC}/PelBen/botflux.F90 #AN added
        ${BFMSRC}/Bennut/ModuleBenQ1Transport.F90 #AN added
            ${BFMSRC}/Bennut/ModuleBenNutVariables.f90
	    ${BFMSRC}/PelB/ModulePhyto.F90 #AN added
	    ${BFMSRC}/PelB/ModulePhaeo.f90 #AN added
	    ${BFMSRC}/PelB/ModulePelBac.f90 #AN added
	    ${BFMSRC}/PelB/ModulePelGlobal.f90 #AN added
	    ${BFMSRC}/PelB/ModuleMesoZoo.f90 #AN added
	    ${BFMSRC}/PelB/ModuleMicroZoo.f90 #AN added
	    ${BFMSRC}/PelB/ModulePelChem.f90 #JM added
#JM	    ${BFMSRC}/PelB/SizeRelatedProcessesModule.F90 #AN added
	    ${BFMSRC}/Bennut/ModuleBenSilica.f90 #AN added
#JM	    ${BFMSRC}/General/Calculate1procentLight.F90
	    ${BFMSRC}/Silt/ModuleSilt.f90
	    ${BFMSRC}/CO2/ModuleCO2.f90     # JM added: dependency CheckMassConservationC.F90
	    ${BFMSRC}/General/CheckMassConservationC.F90
	    ${BFMSRC}/General/CheckMassConservationNPS.F90
#JM        ${BFMSRC}/Ben/ModuleInfoForNextStep.F90


           )

#set_target_properties(global_mod PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(global_mod  prebasis turbulence )
set(bfm_all ${bfm_all} global_mod)

FUNCTION( dodir arg0 arg1 arg2 arg3  arg4)
  file(GLOB filelistname "${arg0}/${arg2}/.files")
  file(GLOB ${arg1} "${arg0}/${arg2}/*.[fF]90")

  foreach(item ${${arg1}} )
     set(search "[^/]+$" )
     string(REGEX MATCH "${search}" sname "${item}" )
     list(APPEND namesonly ${sname} )
  endforeach(item)

  if(filelistname STREQUAL "")
    file(WRITE "${arg0}/${arg2}/.files" "${namesonly}")
  else()
    file(READ "${filelistname}" filelist)

    if( NOT filelist STREQUAL namesonly)
      set( copy "${namesonly}")
      foreach(item ${namesonly} )
        list(FIND filelist "${item}" i)
        if ( i GREATER -1 )
          list(REMOVE_AT filelist ${i} )
          list(REMOVE_ITEM copy ${item} )
        endif()
      endforeach(item)
      message(STATUS "filelist<>${arg1}")
      message(STATUS "filelist${filelist}")
      if (NOT filelist STREQUAL "")
        message(STATUS "--------------------------------------------------------- " )
        message(STATUS "Warning: existing file(s) ${filelist} in ${arg0}/${arg2} removed?")
        message(STATUS "if yes, remove ${filelistname}" )
        message(STATUS "if none, restore files(s) ${copy} in ${arg0}/${arg2}" )
        message(STATUS "--------------------------------------------------------- " )
        message(FATAL_ERROR "")
      endif()
      if (NOT copy STREQUAL "")
        message(STATUS "--------------------------------------------------------- " )
        message(STATUS " Warning: newfile ${copy} in ${arg0}/${arg2} added?" )
        message(STATUS "if yes, remove ${filelistname}" )
        message(STATUS "if none, remove file(s) ${copy} in ${arg0}/${arg2}" )
        message(STATUS "--------------------------------------------------------- " )
        message(FATAL_ERROR "")
      endif()
    endif()
  endif()
  message(STATUS "files collected in ${arg0}/${arg2}")

   foreach(item ${${arg1}} )
     file(STRINGS "${item}" content)
     string (REPLACE ";" "\n" content "${content}")
     set(search "\n *(MODULE|module)" )
     string(REGEX MATCH "${search}" test "${content}" )
     if (test STREQUAL "" )
       list(APPEND subrfun ${item} )
     else()
       list(APPEND modul ${item} )
     endif()
    endforeach(item)
    set(${arg3} "${subrfun}" PARENT_SCOPE)
    set(${arg4} "${modul}" PARENT_SCOPE)
endfunction(dodir)

#dodir( ${BFMSRC} basis Basis basis_fun_list basis_mod_list)
#add_library(basis_fun ${basis_fun_list} )
#add_library(basis_mod ${basis_mod_list} )
#target_link_libraries(basis_mod  global_mod bio)
#target_link_libraries(basis_fun basis_mod prebasis)
#set(bfm_all ${bfm_all} basis_mod)
#set(bfm_all ${bfm_all} basis_fun)

dodir( ${BFMSRC} general General general_fun_list general_mod_list)
add_library(general_fun ${general_fun_list} )
add_library(general_mod ${general_mod_list} )
#target_link_libraries(general_fun bio pelben_fun global_mod basis_mod oxygen_mod pelagic_fun pelagic_mod benthic_mod pelben_mod co2dyn_mod siltdyn_mod co2dyn_fun )
target_link_libraries(general_fun bio pelben_fun global_mod oxygen_mod pelagic_fun pelagic_mod benthic_mod pelben_mod co2dyn_mod siltdyn_mod co2dyn_fun )
target_link_libraries(general_mod global_mod)
set(bfm_all ${bfm_all} general_fun)

dodir( ${BFMSRC} oxygen Oxygen oxygen_fun_list oxygen_mod_list)
add_library(oxygen_fun ${oxygen_fun_list} )
add_library(oxygen_mod ${oxygen_mod_list} )
target_link_libraries(oxygen_mod global_mod)
#target_link_libraries(oxygen_fun basis_mod oxygen_mod pelben_mod)
target_link_libraries(oxygen_fun oxygen_mod pelben_mod)
set(bfm_all ${bfm_all} oxygen_mod)
set(bfm_all ${bfm_all} oxygen_fun)

dodir( ${BFMSRC} pelagic PelB pelagic_fun_list pelagic_mod_list)
add_library(pelagic_fun ${pelagic_fun_list} )
add_library(pelagic_mod ${pelagic_mod_list} )
#target_link_libraries(pelagic_mod basis_mod prebasis gotm_fin)
target_link_libraries(pelagic_mod prebasis gotm_fin)
target_link_libraries(pelagic_fun general_fun general_mod pelagic_mod oxygen_fun gotm_fin  )
add_dependencies(pelagic_fun gotm_fin pelagic_mod)
set(bfm_all ${bfm_all} pelagic_mod)
set(bfm_all ${bfm_all} pelagic_fun)

dodir( ${BFMSRC} pelben PelBen pelben_fun_list pelben_mod_list)
add_library(pelben_fun ${pelben_fun_list} )
add_library(pelben_mod ${pelben_mod_list} )
target_link_libraries(pelben_mod global_mod prebasis)
target_link_libraries(pelben_fun benthic_mod benthic_fun siltdyn_mod pelagic_mod)
add_dependencies(pelben_fun pelagic_mod)
set(bfm_all ${bfm_all} pelben_mod)
set(bfm_all ${bfm_all} pelben_fun)

dodir( ${BFMSRC} benthic Ben benthic_fun_list benthic_mod_list)
add_library(benthic_fun ${benthic_fun_list} )
add_library(benthic_mod ${benthic_mod_list} )
target_link_libraries(benthic_mod global_mod pelagic_mod )
target_link_libraries(benthic_fun benthic_mod global_mod pelagic_mod pelben_fun bennut_mod bennut_fun pelagic_fun)
add_dependencies(benthic_fun benthic_mod)
set(bfm_all ${bfm_all} benthic_mod)
set(bfm_all ${bfm_all} benthic_fun)

if ( GOTM_USE_BFM )
#  target_link_libraries(benthic_fun benthic_mod bennut_mod co2dyn_mod pelben_mod turbulence )
  target_link_libraries(benthic_fun benthic_mod bennut_mod co2dyn_mod pelben_mod )
else()
  target_link_libraries(benthic_fun benthic_mod bennut_mod co2dyn_mod pelben_mod gotm_basis )
endif()

dodir( ${BFMSRC} bennut Bennut bennut_fun_list bennut_mod_list)
add_library(bennut_fun ${bennut_fun_list} )
add_library(bennut_mod ${bennut_mod_list} )
target_link_libraries(bennut_mod general_mod benthic_mod global_mod)
#target_link_libraries(bennut_fun co2dyn_mod pelben_mod bennut_mod basis_mod basis_fun pelben_fun )
target_link_libraries(bennut_fun co2dyn_mod pelben_mod bennut_mod pelben_fun )
add_dependencies(bennut_fun bennut_mod)
add_dependencies(bennut_mod benthic_mod)
set(bfm_all ${bfm_all} bennut_mod)
set(bfm_all ${bfm_all} bennut_fun)

dodir( ${BFMSRC} co2dyn CO2 co2dyn_fun_list co2dyn_mod_list)
add_library(co2dyn_fun ${co2dyn_fun_list} )
add_library(co2dyn_mod ${co2dyn_mod_list} )
target_link_libraries(co2dyn_mod global_mod)
#target_link_libraries(co2dyn_fun co2dyn_mod bennut_mod benthic_mod bennut_fun basis_fun)
target_link_libraries(co2dyn_fun co2dyn_mod bennut_mod benthic_mod bennut_fun)
set(bfm_all ${bfm_all} co2dyn_mod)
set(bfm_all ${bfm_all} co2dyn_fun)

dodir( ${BFMSRC} share share share_fun_list share_mod_list)
add_library(share_fun ${share_fun_list} )
add_library(share_mod ${share_mod_list} )
target_link_libraries(share_mod global_mod general_mod co2dyn_mod)
target_link_libraries(share_fun general_mod co2dyn_mod)
set(bfm_all ${bfm_all} share_mod)
set(bfm_all ${bfm_all} share_fun)

dodir( ${BFMSRC} siltdyn Silt siltdyn_fun_list siltdyn_mod_list)
add_library(siltdyn_fun ${siltdyn_fun_list} )
add_library(siltdyn_mod ${siltdyn_mod_list} )
target_link_libraries(siltdyn_mod global_mod)
#target_link_libraries(siltdyn_fun siltdyn_mod basis_mod )
target_link_libraries(siltdyn_fun siltdyn_mod )
set(bfm_all ${bfm_all} siltdyn_mod)
set(bfm_all ${bfm_all} siltdyn_fun)

set(bfm_up ${bfm_all} PARENT_SCOPE)

#make available for passing to GETM cmake
#set(global_mod global_mod PARENT_SCOPE)

# Install BFM library and include files, unless the host overruled this.
#if (NOT DEFINED BFM_NO_LIBRARY_INSTALL)
#   install(TARGETS ${bfm_all} 
#     EXPORT bfmConfig
#     DESTINATION lib)
#   install(DIRECTORY ${CMAKE_BINARY_DIR}/modules/\${BUILD_TYPE}/
#     DESTINATION include)
##  install(FILES ${PROJECT_SOURCE_DIR}/../include/fabm.h ${PROJECT_SOURCE_DIR}/drivers/${FABM_HOST}/fabm_driver.h
##    DESTINATION include)
#endif()
#
#export(EXPORT bfmConfig)
#install(EXPORT bfmConfig DESTINATION cmake)

