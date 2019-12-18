#
# Copyright 2018-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

find_package(Boost REQUIRED COMPONENTS unit_test_framework)

macro(addtest testname)
  add_executable(${testname} ${testname}.cpp)
  target_link_libraries(${testname} PUBLIC ${PROJECT_NAME} Boost::unit_test_framework)
  add_test(${testname}Unit ${testname})
  generate_msvc_dot_user_file(${testname})
endmacro()

addtest(CSCMatrixTest)
addtest(QPTest)
