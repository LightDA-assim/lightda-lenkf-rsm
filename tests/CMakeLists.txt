
if(FPRETTIFY_EXECUTABLE)

  set(FPRETTIFY_CMD ${DO_FPRETTIFY})
  list(APPEND FPRETTIFY_CMD "${ALL_SOURCES}")
  string(REPLACE ";" "\;" FPRETTIFY_CMD "${FPRETTIFY_CMD}")
  add_test(NAME test_fprettify
    COMMAND ${CMAKE_COMMAND}
    -DCMD=${FPRETTIFY_CMD}
    -DFAIL_PATTERN=.
    -P ${CMAKE_CURRENT_SOURCE_DIR}/test_and_check_regex.cmake)
endif()
