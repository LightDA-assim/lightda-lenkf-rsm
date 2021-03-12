set(LENKF_RSM_SOURCES lenkf_rsm.f90 lenkf_rsm_c.f90 lenkf_rsm_filter.F90)
add_library(lightda-lenkf-rsm SHARED ${LENKF_RSM_SOURCES})
target_link_libraries(lightda-lenkf-rsm PUBLIC lightda)
target_link_libraries(lightda-lenkf-rsm PRIVATE ${LAPACK_LIBRARIES})
target_include_directories(lightda-lenkf-rsm INTERFACE
  $<INSTALL_INTERFACE:${Fortran_MODULE_DIRECTORY}>)
install(TARGETS lightda-lenkf-rsm EXPORT lightda-lenkf-rsm
  LIBRARY DESTINATION lib ${CMAKE_INSTALL_LIBDIR})

set(ALL_SOURCES ${LENKF_RSM_SOURCES})

# Prepend path to sources and propagate to parent scope
list_transform_prepend(ALL_SOURCES "src/")
set(ALL_SOURCES ${ALL_SOURCES} PARENT_SCOPE)
