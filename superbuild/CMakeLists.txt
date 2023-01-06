cmake_minimum_required(VERSION 3.10)

include(ExternalProject)

project(lenkf-rsm)

set(
  system_mpi_GIT_URL
  "https://github.com/LightDA-assim/system_mpi.git" CACHE STRING
  "URL of system_mpi git repository")

set(
  fortran_exceptions_GIT_URL
  "https://github.com/LightDA-assim/fortran_exceptions.git" CACHE STRING
  "URL of fortran_exceptions git repository")

set(
  lightda_GIT_URL
  "https://github.com/LightDA-assim/lightda-core.git" CACHE STRING
  "URL of fortran_exceptions git repository")

ExternalProject_Add(
  system_mpi
  GIT_REPOSITORY ${system_mpi_GIT_URL}
  GIT_TAG main
  CMAKE_CACHE_ARGS
        -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/system_mpi)

ExternalProject_Add(
  fortran_exceptions
  GIT_REPOSITORY ${fortran_exceptions_GIT_URL}
  GIT_TAG main
  CMAKE_CACHE_ARGS
        -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/fortran_exceptions)

set(
  fhash_GIT_URL
  "https://github.com/jhaiduce/fhash.git" CACHE STRING
  "URL of fhash git repository")

ExternalProject_Add(
  fhash
  GIT_REPOSITORY ${fhash_GIT_URL}
  GIT_TAG add_cmake_build_files
  CMAKE_CACHE_ARGS
        -DPYTHON_EXECUTABLE:PATH=${PYTHON_EXECUTABLE}
        -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/fhash
	-DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE})

ExternalProject_Add(
  lightda
  GIT_REPOSITORY ${lightda_GIT_URL}
  GIT_TAG main
  CMAKE_CACHE_ARGS
    -Dsystem_mpi_DIR:STRING=${CMAKE_CURRENT_BINARY_DIR}/system_mpi/lib/cmake/system_mpi
    -Dfortran_exceptions_DIR:STRING=${CMAKE_CURRENT_BINARY_DIR}/fortran_exceptions/lib/cmake/fortran_exceptions
    -Dfhash_DIR:STRING=${CMAKE_CURRENT_BINARY_DIR}/fhash/lib/cmake/fhash
    -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/lightda
  DEPENDS system_mpi fortran_exceptions)

ExternalProject_Add(lightda-lenkf-rsm
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/.."
  CMAKE_CACHE_ARGS
    -Dsystem_mpi_DIR:STRING=${CMAKE_CURRENT_BINARY_DIR}/system_mpi/lib/cmake/system_mpi
    -Dfortran_exceptions_DIR:STRING=${CMAKE_CURRENT_BINARY_DIR}/fortran_exceptions/lib/cmake/fortran_exceptions
    -Dfhash_DIR:STRING=${CMAKE_CURRENT_BINARY_DIR}/fhash/lib/cmake/fhash
    -Dlightda_DIR:STRING=${CMAKE_CURRENT_BINARY_DIR}/lightda/lib/cmake/lightda
    -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/lightda-lenkf-rsm
  DEPENDS system_mpi fortran_exceptions lightda)
