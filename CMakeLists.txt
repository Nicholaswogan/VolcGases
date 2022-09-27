cmake_minimum_required(VERSION "3.14.0")
project(VOLCGASES LANGUAGES Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

option(SKBUILD "Should be ON of being build by skbuild, 
and OFF of being build by regular cmake" OFF)

if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

include(cmake/CPM.cmake)
            
add_subdirectory(src)
add_subdirectory(test)

if (SKBUILD)
  install(TARGETS volcgases DESTINATION VolcGases)
else()
  install(TARGETS volcgases DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/VolcGases/)
endif()