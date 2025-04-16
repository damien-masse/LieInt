# ==================================================================
#  codac / basics example - cmake configuration file
# ==================================================================

set(SRC_TESTS
    test_intExp1 
    test_lieSE2_0 
    test_lieSE2_1 
    test_lieSO3_0 
    test_contract
    integrationSE2/integrate_SE2_0
    integrationSE2/integrate_SE2_1
    integrationSO3/integrate_SO3_0
)

# Compilation

foreach(SRC_TEST ${SRC_TESTS})
  string(REPLACE "/" "_" EXECUTABLE ${SRC_TEST})
  add_executable(${EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${SRC_TEST}.cpp)
  target_include_directories(${EXECUTABLE} PUBLIC ${CMAKE_SOURCE_DIR} ${CODAC_INCLUDE_DIRS} ${IBEX_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIRS})
  target_link_libraries(${EXECUTABLE} ${CODAC_LIBRARIES} ${LIB})
endforeach(SRC_TEST ${SRC_TESTS})
