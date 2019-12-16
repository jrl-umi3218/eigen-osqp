macro(generate_example name)
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} ${PROJECT_NAME})
    target_include_directories(${name} PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../src>)
endmacro()

generate_example(dense_matrix_qp)
generate_example(sparse_matrix_qp)