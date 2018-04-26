set(SOURCES OSQP.cpp CSCMatrix.cpp)
set(HEADERS eigen-osqp/OSQP.h eigen-osqp/CSCMatrix.h)

add_definitions("-DEIGEN_OSQP_EXPORT")

add_library(${PROJECT_NAME} SHARED ${SOURCES} ${HEADERS})
target_link_libraries(${PROJECT_NAME} osqp)
set_target_properties(${PROJECT_NAME} PROPERTIES SOVERSION 0 VERSION 1.0.0)
if(${USE_LTO})
  target_compile_options(${PROJECT_NAME} PRIVATE "-flto")
  set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS "-flto")
endif()

set(INSTALL_PATH include/${PROJECT_NAME})
install(FILES ${HEADERS} DESTINATION ${INSTALL_PATH})
install(TARGETS ${PROJECT_NAME}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

