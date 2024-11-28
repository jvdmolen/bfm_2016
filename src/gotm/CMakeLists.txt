
# Global include directories
INCLUDE_DIRECTORIES( "${NetCDF_INCLUDE_DIRS}")
include_directories(BEFORE
                    "${BFM_BASE}/include")

include(${BFM_BASE}/src/IncludeOptions.txt OPTIONAL)
add_definitions( ${DEFINES} ) 

if(DEFINED ENV{GOTM_BUILD_DIR})
	set(GOTM_BUILD_DIR $ENV{GOTM_BUILD_DIR})
	message(STATUS a)
else()
	set(GOTM_BUILD_DIR "/export/lv6/user/anunes/home/build/_gotm")
	message(STATUS b)
endif()
message(STATUS GOTM_BUILD_DIR= ${GOTM_BUILD_DIR})

include_directories("${GOTM_BUILD_DIR}/modules")

# coupling libraries

add_library(gotm_fin
              bio_bfm.F90
              process_model.F90
              bio_solver.F90
           )
#target_link_libraries(gotm_fin global_mod turbulence general_fun siltdyn_fun)
target_link_libraries(gotm_fin general_fun global_mod siltdyn_fun)
#target_include_directories(gotm_fin PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
#add_dependencies(gotm_fin bio)
set(bfm_cpl gotm_fin)

message( STATUS  "adding output_gotm library")
add_library(output_gotm
	    ncdfout.F90
	    output.F90
	    asciiout.F90
	    gotmout.F90
	    hotstart.F90
           )
#target_link_libraries( output_gotm global_mod meanflow turbulence observations airsea_driver)
target_link_libraries( output_gotm global_mod gotm_fin turbulence )
target_link_libraries( output_gotm  ${GOTM_BUILD_DIR}/meanflow/libmeanflow.a)
target_link_libraries( output_gotm  ${GOTM_BUILD_DIR}/airsea/libairsea_driver.a)
target_link_libraries( output_gotm  ${GOTM_BUILD_DIR}/observations/libobservations.a)
add_dependencies(gotm_fin output_gotm) 
#target_include_directories(output_gotm PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
#target_link_libraries( output_gotm global_mod)
set(bfm_cpl ${bfm_cpl} output_gotm)

message( STATUS  "adding bfm_ncdf library")
add_library(bfm_ncdf
              bfm_output_ncdf.F90
           )
#              gotm/make_flux_output.F90
#target_include_directories(bfm_ncdf PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries( bfm_ncdf global_mod output_gotm gotm_fin)
add_dependencies(bfm_ncdf gotm_fin)
set(bfm_cpl ${bfm_cpl} bfm_ncdf)

message( STATUS  "adding bio library")
add_library(bio
              bio.F90
           )
#target_include_directories(bio PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)
#target_include_directories(bio PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../modules>)
target_link_libraries( bio global_mod gotm_fin util )
#target_include_directories(bio PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
add_dependencies(bio bfm_ncdf )    #fake, but it makes cmake compile bio
set(bfm_cpl ${bfm_cpl} bio)

#set(gotm_fin ${variables_bio_3d} PARENT_SCOPE)
#set(modify_meteo_bfm ${modify_meteo_bfm} PARENT_SCOPE)
#set(coupling_getm_bfm ${coupling_getm_bfm} PARENT_SCOPE)
#set(coupling_getm_bfm_rivers ${coupling_getm_bfm_rivers} PARENT_SCOPE)
#set(restart_bfm ${restart_bfm} PARENT_SCOPE)
#set(read_restart_bfm ${read_restart_bfm} PARENT_SCOPE)
#set(output_2d_bfm ${output_2d_bfm} PARENT_SCOPE)
#set(bfm_getmlibs variables_bio_3d modify_meteo_bfm coupling_getm_bfm coupling_getm_bfm_rivers restart_bfm read_restart_bfm output_2d_bfm)
#set(bfm_up_getm ${bfm_all})
#set(bfm_up_getm ${bfm_up_getm} PARENT_SCOPE)
#set(bfm_getmlibs ${bfm_getmlibs} PARENT_SCOPE)
#set(bfm_up_files ${bfm_all_files})
#set(bfm_up_files ${bfm_up_files} PARENT_SCOPE)
#set(bfm_up ${bfm_up} PARENT_SCOPE)
set(bfm_cpl ${bfm_cpl} PARENT_SCOPE)


