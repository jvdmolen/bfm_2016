
# Global include directories
INCLUDE_DIRECTORIES( "${NetCDF_INCLUDE_DIRS}")

include(${BFM_BASE}/src/IncludeOptions.txt OPTIONAL)
add_definitions( ${DEFINES} ) 
if (";${DEFINES};" MATCHES ";-DINCLUDE_MACROPHYT;" )
  option(MACROPHYT "Include macrophyte calculations" ON)
endif()

if(DEFINED ENV{GOTM_BUILD_DIR})
	set(GOTM_BUILD_DIR $ENV{GOTM_BUILD_DIR})
else()
	set(GOTM_BUILD_DIR "/export/lv6/user/anunes/home/build/_gotm")
endif()
message(STATUS GOTM_BUILD_DIR= ${GOTM_BUILD_DIR})


include_directories("${GOTM_BUILD_DIR}/modules" "${GOTM_BUILD_DIR}/bfm_src")

# coupling libraries

add_library(variables_bio_3d
            variables_bio_3d.F90
           )
#target_link_libraries(variables_bio_3d domian)

add_library(coupling_getm_bfm
            coupling_getm_bfm.F90
           )
target_link_libraries(coupling_getm_bfm global_mod variables_bio_3d bio gotm_fin)

add_library(coupling_getm_bfm_rivers
            coupling_getm_bfm_rivers.F90
           )
target_link_libraries(coupling_getm_bfm_rivers basis_fun global_mod variables_bio_3d)

add_library(restart_bfm
            output_restart_bio.F90
           )
target_link_libraries(restart_bfm global_mod coupling_getm_bfm bio bfm_ncdf)

add_library(read_restart_bfm
	    read_restart_ncdf.F90
            read_restart_bio_ncdf.F90
           )
target_link_libraries(read_restart_bfm global_mod)

add_library(output_2d_bfm
            init_2dbio_ncdf.F90
            save_bfm_2d.F90
           )
target_link_libraries(output_2d_bfm ${NetCDF_LIBRARIES})

add_library(modify_meteo_bfm 
            modify_meteo_input.F90
          )
target_link_libraries(modify_meteo_bfm global_mod)
target_link_libraries(modify_meteo_bfm ${GOTM_BUILD_DIR}/bfm_src/libglobal_mod.a)

if (MACROPHYT)
  add_library(macro_phyto_bfm
            Module2DMacroPhyto.F90
             )
  target_link_libraries(macro_phyto_bfm global_mod)
  target_link_libraries(macro_phyto_bfm ${GOTM_BUILD_DIR}/bfm_src/libglobal_mod.a)
  set(macro_phyto_bfm ${macro_phyto_bfm} PARENT_SCOPE)
endif()

set(variables_bio_3d ${variables_bio_3d} PARENT_SCOPE)
set(modify_meteo_bfm ${modify_meteo_bfm} PARENT_SCOPE)
set(coupling_getm_bfm ${coupling_getm_bfm} PARENT_SCOPE)
set(coupling_getm_bfm_rivers ${coupling_getm_bfm_rivers} PARENT_SCOPE)
set(restart_bfm ${restart_bfm} PARENT_SCOPE)
set(read_restart_bfm ${read_restart_bfm} PARENT_SCOPE)
set(output_2d_bfm ${output_2d_bfm} PARENT_SCOPE)

set(bfm_getmlibs variables_bio_3d modify_meteo_bfm coupling_getm_bfm coupling_getm_bfm_rivers restart_bfm read_restart_bfm output_2d_bfm)
set(bfm_getmlibs ${bfm_getmlibs} PARENT_SCOPE)



