# ==================================================================
#  codac / basics example - cmake configuration file
# ==================================================================

  cmake_minimum_required(VERSION 3.0.2)
  project(LANGUAGES CXX)

# Adding IBEX

  # In case you installed IBEX in a local directory, you need 
  # to specify its path with the CMAKE_PREFIX_PATH option.
  # set(CMAKE_PREFIX_PATH "~/ibex-lib/build_install")
#
#  find_package(IBEX REQUIRED)
#  ibex_init_common() # IBEX should have installed this function
#  message(STATUS "Found IBEX version ${IBEX_VERSION}")
#  message(STATUS "Found IBEX inclusion ${IBEX_INCLUDE_DIRS}")
#
## Adding Eigen3
#
#  # In case you installed Eigen3 in a local directory, you need
#  # to specify its path with the CMAKE_PREFIX_PATH option, e.g.
#  # set(CMAKE_PREFIX_PATH "~/eigen/build_install")
#
  find_package(Eigen3 REQUIRED NO_MODULE)
  message(STATUS "Found Eigen3 version ${EIGEN3_VERSION}")
#
## Adding Codac
#
#  # In case you installed Codac in a local directory, you need 
#  # to specify its path with the CMAKE_PREFIX_PATH option.
#  # set(CMAKE_PREFIX_PATH "~/codac/build_install")
#
# set(CMAKE_PREFIX_PATH "/home/damien/Recherche/codac/codac/build_install")
  find_package(CODAC REQUIRED codac)
  message(STATUS "Found Codac version ${CODAC_VERSION}")
  message(STATUS "Found Codac inclusion ${CODAC_INCLUDE_DIRS}")
  message(STATUS "Found Codac flags ${CODAC_CFLAGS}")
  message(STATUS "Found Codac libraries ${CODAC_LIBRARIES}")

set(LIB intLie)
set(EXECS test_intExp1 test_lieSE2_0 test_lieSE2_1 test_lieSO3_0 test_contract)

# Compilation

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_compile_options(-Wall)

add_subdirectory(src)

add_subdirectory(tests)

