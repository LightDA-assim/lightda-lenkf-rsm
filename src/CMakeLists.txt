set(LENKF_RSM_SOURCES lenkf_rsm.f90 lenkf_rsm_c.f90 lenkf_rsm_filter.F90)
add_library(lightda-lenkf-rsm SHARED ${LENKF_RSM_SOURCES})
target_link_libraries(lightda-lenkf-rsm PUBLIC lightda)
target_link_libraries(lightda-lenkf-rsm PRIVATE ${LAPACK_LIBRARIES})
target_include_directories(lightda-lenkf-rsm INTERFACE
  $<INSTALL_INTERFACE:${Fortran_MODULE_DIRECTORY}>)
