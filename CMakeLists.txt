cmake_minimum_required(VERSION 3.0.2)

project(lightda-lenkf-rsm)

enable_language(Fortran)

# Prevent older CMake versions from setting include directories with -isystem,
# which prevents gfortran from finding Fortran module files
set(CMAKE_NO_SYSTEM_FROM_IMPORTED ON)

# Stop newer linkers on Linux using RUNPATH instead of RPATH, which causes
# problems loading indirect library dependencies.
# Note it may also be necessary to remove -Wl,--enable-new-dtags from linker
# flags elsewhere (for example this --enable-new-dtags may be set by recent
# versions of OpenMPI and must be removed afterward)

if($CMAKE_VERSION VERSION_GREATER 3.13.0)
  add_link_options($<$<CXX_COMPILER_ID:GNU>:-Wl,--disable-new-dtags>)
else()
  if(UNIX AND NOT APPLE)
    set (CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} -Wl,--disable-new-dtags)
  endif()
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/modules/")

if(DEFINED CMAKE_Fortran_COMPILER_VERSION)
  set(Fortran_MODULE_DIRECTORY include/fortran_modules/${CMAKE_Fortran_COMPILER_ID}/${CMAKE_Fortran_COMPILER_VERSION})
else()
  set(Fortran_MODULE_DIRECTORY include/fortran_modules/${CMAKE_Fortran_COMPILER_ID})
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${Fortran_MODULE_DIRECTORY})

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 9)
    # Compiler doesn't support overriding finalizers
  else()
    # Compiler supports overriding finalizers
    add_definitions(-DOVERRIDABLE_FINALIZERS)
  endif()
else()
  # Compiler supports overriding finalizers
  add_definitions(-DOVERRIDABLE_FINALIZERS)
endif()

find_package(lightda REQUIRED NO_POLICY_SCOPE)

find_package(LAPACK REQUIRED)

function(list_transform_prepend var prefix)
    set(temp "")
    foreach(f ${${var}})
        list(APPEND temp "${prefix}${f}")
    endforeach()
    set(${var} "${temp}" PARENT_SCOPE)
endfunction()

add_subdirectory(src)

include(CMakePackageConfigHelpers)

# Add all targets to the build-tree export set
export(TARGETS lightda-lenkf-rsm
  FILE "${PROJECT_BINARY_DIR}/lightda-lenkf-rsmTargets.cmake")

# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE lightda-lenkf-rsm)

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/lightda-lenkf-rsmConfig.cmake"
  INSTALL_DESTINATION lib/cmake/lightda-lenkf-rsm
  )

install(FILES
          "${CMAKE_CURRENT_BINARY_DIR}/lightda-lenkf-rsmConfig.cmake"
        DESTINATION lib/cmake/lightda-lenkf-rsm
)

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION ${Fortran_MODULE_DIRECTORY})

install(EXPORT lightda-lenkf-rsm
  FILE lightda-lenkf-rsmTargets.cmake
  DESTINATION lib/cmake/lightda-lenkf-rsm
)

# -------------------
# Build documentation
# -------------------

set ( SKIP_DOC_GEN FALSE CACHE BOOL
  "Disable building the API documentation with FORD" )
if ( NOT SKIP_DOC_GEN )
  find_program(FORD ford)

  find_package(PythonInterp)

  if(NOT PYTHON_VERSION_STRING VERSION_LESS 3.5)
    set(MARKDOWN_SUPPORTED_PYTHON true)
  endif()

  if(NOT MARKDOWN_SUPPORTED_PYTHON)
    message(WARNING "Python version ${PYTHON_VERSION_STRING} is not supported by Python-Markdown, not attempting to build documentation.")
  endif()

  if(FORD AND MARKDOWN_SUPPORTED_PYTHON)

    set(DOC_ROOT "${PROJECT_BINARY_DIR}/doc")
    set(DOC_DIR "${DOC_ROOT}/lightda-lenkf-rsm")
    set(PROJ_DIR "${PROJECT_SOURCE_DIR}")
    set(FORD_PROJECT_FILE "${CMAKE_CURRENT_BINARY_DIR}/lightda-lenkf-rsm.md")
    set(MACRO_FLAG "")
    string(REPLACE ";" "\ninclude:" FORD_PROJECT_INCLUDES "${PROJ_DIR}/src")
    set(FORD_PROJECT_EXCLUDE_DIRS "${MPI_Fortran_INCLUDE_DIRS};${PROJECT_BINARY_DIR};${PROJ_DIR}/build;${PROJ_DIR}/tests")
    string(REPLACE ";" "\nexclude_dir:" FORD_PROJECT_EXCLUDE_DIRS "${FORD_PROJECT_EXCLUDE_DIRS}")
    set(FORD_PROJECT_EXCLUDE_FILES "")
    string(REPLACE ";" "\nexclude:" FORD_PROJECT_EXCLUDE_FILES "${FORD_PROJECT_EXCLUDE_FILES}")

    configure_file(lightda-lenkf-rsm.md.in lightda-lenkf-rsm.md)

    # Pick the preprocessor to use based on the Fortran compiler
    if ( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel" )
      set ( FPP "fpp\n" )
    else ()
      set ( FPP "gfortran -E\n" ) # default to gfortran -E for gfortran and unsupported compilers
    endif ()
    file ( WRITE "${PROJECT_BINARY_DIR}/.PREPROCESSOR" "${FPP}" )

    # Dynamically generate the FORD outputs list
    message ( STATUS "Dynamically computing FORD output information..." )
    if ( NOT (DEFINED FORD_OUTPUTS_CACHED) )
      message ( STATUS "Running FORD to dynamically compute documentation outputs, this could take a while..." )
      execute_process ( COMMAND ${CMAKE_COMMAND} -E remove_directory ${DOC_DIR}
	COMMAND ${CMAKE_COMMAND} -E make_directory ${DOC_DIR}
	COMMAND "${FORD}" --debug -q ${MACRO_FLAG} -d "${PROJ_DIR}" -o "${DOC_DIR}" -p "${PAGE_DIR}" "${FORD_PROJECT_FILE}" OUTPUT_QUIET )
    else ()
      message ( STATUS "Re-using cached FORD outputs, rather than regenerating them" )
    endif()

    # Compile a list of FORD output files
    file ( GLOB_RECURSE FORD_OUTPUTS
      "${DOC_DIR}/*.*" )
    file ( GLOB_RECURSE FORD_CLEAN_OUTPUTS
      "${DOC_DIR}/*.*" )

    # Cache FORD outputs
    if ( (DEFINED FORD_OUTPUTS) AND ( NOT ( "${FORD_OUTPUTS}" STREQUAL "" ) ) )
      message ( STATUS "Caching FORD outputs" )
      set ( FORD_OUTPUTS_CACHED "${FORD_OUTPUTS}"
	CACHE STRING "variable containing FORD outputs to prevent rebuilding FORD docs" FORCE )
    endif ()
    message ( STATUS "Done dynamically computing FORD outputs." )

    foreach ( DOC_SRC_FILE ${PROJECT_BINARY_DIR}/lightda-lenkf-rsm.md ${SYSTEM_MPI_SOURCES} ${UTIL_SOURCES} ${ASSIM_SOURCES} ${ASSIM_PARALLEL_SOURCES} ${ADVECT1D_ASSIMILATE_INTERFACES_SOURCES} ${ADVECT1D_ASSIMILATE_SOURCES} ${ADVECT1D_SOURCES} ${PAGES_FILES} )
      list ( APPEND FORD_DEPENDS "${DOC_SRC_FILE}" )
    endforeach ()

    # Re-build FORD output if needed
    add_custom_command ( OUTPUT ${FORD_OUTPUTS_CACHED}
      COMMAND "${FORD}" --debug ${MACRO_FLAG} -d "${PROJ_DIR}" -o "${DOC_DIR}" "${FORD_PROJECT_FILE}"
      MAIN_DEPENDENCY "${FORD_PROJECT_FILE}"
      DEPENDS ${FORD_DEPENDS}
      COMMENT "Building HTML documentation for ${PROJECT_NAME} using FORD" )

    add_custom_target ( documentation ALL
      DEPENDS ${FORD_OUTPUTS_CACHED} )

    # Install documentation
    set ( INSTALL_API_DOCUMENTATION TRUE
      CACHE BOOL "Install FORD generated documentation?" )
    if ( INSTALL_API_DOCUMENTATION )
      install ( DIRECTORY "${DOC_ROOT}/" DESTINATION share/doc )
    else ()

      # FORD Not found
      message ( WARNING
	"FORD not found. Please set the CMake cache variable FORD to point to the installed FORD executable, and reconfigure or disable building the documentation. FORD can be installed from PYPI with `sudo pip install FORD` or from <https://github.com/cmacmackin/ford> If you do not wish to install FORD and build the JSON-Fortran documentation, then please set the CMake cache variable SKIP_DOC_GEN to TRUE." )

    endif ()
  endif ()

endif()

find_program(FPRETTIFY_EXECUTABLE NAMES "fprettify" DOC "Path to fprettify")

if(NOT FPRETTIFY_EXECUTABLE)
  message(STATUS "fprettify not found.")
else()
  set(DO_FPRETTIFY "${FPRETTIFY_EXECUTABLE}" "--diff" "-i" "2" "-l" "81")

  list_transform_prepend(ALL_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/)

  add_custom_target(
    fprettify
    COMMAND ${DO_FPRETTIFY}
    ${ALL_SOURCES}
    )

endif()

enable_testing()
add_subdirectory(tests)
