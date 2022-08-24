#
# This is a modified/original file distributed in RESPACK code under GNU GPL ver.3.
# https://sites.google.com/view/kazuma7k6r
#

cmake_minimum_required(VERSION 2.8.10)
if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(FATAL_ERROR "cmake should be executed not for 'util/wan2respack' subdirectory, but for the top directory of REIWA90.")
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

set(SOURCES_UTIL util.f90)

add_custom_target(wan2respack ALL)
add_executable(gen_mk.x generate_sample_mk.f90 m_rd_dat_wfn.f90 ${SOURCES_UTIL})
add_executable(gen_wan.x generate_wan.f90 ${SOURCES_UTIL})
target_link_libraries(gen_mk.x ${LAPACK_LIBRARIES} m)
target_link_libraries(gen_wan.x ${LAPACK_LIBRARIES} m)
add_dependencies(wan2respack gen_mk.x gen_wan.x)

install(TARGETS gen_mk.x gen_wan.x RUNTIME DESTINATION bin)
install(FILES init.py wan2respack.py wan2respack_core.py wan2respack_pre.py qe2respack.py 
	DESTINATION bin
        PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE)
