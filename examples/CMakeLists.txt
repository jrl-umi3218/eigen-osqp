#
# Copyright 2018-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

macro(generate_example name)
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} PUBLIC ${PROJECT_NAME})
endmacro()

generate_example(dense_matrix_qp)
generate_example(sparse_matrix_qp)
