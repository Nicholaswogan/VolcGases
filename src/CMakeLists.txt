CPMAddPackage(
  NAME minpack
  VERSION 2.0.0
  GITHUB_REPOSITORY "fortran-lang/minpack"
  GIT_TAG "v2.0.0-rc.1"
  DOWNLOAD_ONLY ON
)
add_library(minpack
  ${minpack_SOURCE_DIR}/src/minpack.f90
)

add_library(volcgases SHARED
  volcgases_const.f90
  volcgases.f90
)
target_link_libraries(volcgases minpack)

# some extra optimizations
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  if ("${CMAKE_BUILD_TYPE}" MATCHES "Release")
    target_compile_options(minpack PRIVATE -funroll-loops -finline-functions)
    target_compile_options(volcgases PRIVATE -funroll-loops -finline-functions)
  endif()
endif()