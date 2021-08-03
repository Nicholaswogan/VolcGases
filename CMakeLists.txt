cmake_minimum_required(VERSION 3.4.0)
project(VolcGases LANGUAGES Fortran)

option(SKBUILD "Should be ON of being build by skbuild, 
and OFF of being build by regular cmake" OFF)

if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()
            
add_library(VolcGasesFort SHARED
            ${CMAKE_CURRENT_SOURCE_DIR}/src/VolcGasesFort.f90
            ${CMAKE_CURRENT_SOURCE_DIR}/src/minpack/dpmpar.f
            ${CMAKE_CURRENT_SOURCE_DIR}/src/minpack/enorm.f
            ${CMAKE_CURRENT_SOURCE_DIR}/src/minpack/fdjac2.f
            ${CMAKE_CURRENT_SOURCE_DIR}/src/minpack/lmder.f
            ${CMAKE_CURRENT_SOURCE_DIR}/src/minpack/lmder1.f
            ${CMAKE_CURRENT_SOURCE_DIR}/src/minpack/lmdif.f
            ${CMAKE_CURRENT_SOURCE_DIR}/src/minpack/lmpar.f
            ${CMAKE_CURRENT_SOURCE_DIR}/src/minpack/qrfac.f
            ${CMAKE_CURRENT_SOURCE_DIR}/src/minpack/qrsolv.f
            ${CMAKE_CURRENT_SOURCE_DIR}/src/VolcGasesFort_wrapper.f90)
            
set_target_properties(VolcGasesFort PROPERTIES PREFIX "")
set_target_properties(VolcGasesFort PROPERTIES SUFFIX ".so")

if (SKBUILD)
  install(TARGETS VolcGasesFort DESTINATION VolcGases)
else()
  install(TARGETS VolcGasesFort DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/VolcGases/)
endif()