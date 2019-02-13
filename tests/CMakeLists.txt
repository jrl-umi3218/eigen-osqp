set(BOOST_COMPONENTS unit_test_framework)
search_for_boost()

include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories(${Boost_INCLUDE_DIRS})

macro(addtest testname)
  add_executable(${testname} ${testname}.cpp)
  target_link_libraries(${testname} ${Boost_LIBRARIES} ${PROJECT_NAME})
  add_test(${testname}Unit ${testname})
  # Adding a project configuration file (for MSVC only)
  GENERATE_MSVC_DOT_USER_FILE(${testname})
endmacro()

addtest(CSCMatrixTest)
addtest(QPTest)