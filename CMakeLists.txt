# ==================================================================
#  codac / basics example - cmake configuration file
# ==================================================================

  cmake_minimum_required(VERSION 3.0.2)
  project(LANGUAGES CXX)

# Adding IBEX

  # In case you installed IBEX in a local directory, you need 
  # to specify its path with the CMAKE_PREFIX_PATH option.
  # set(CMAKE_PREFIX_PATH "~/ibex-lib/build_install")

  find_package(IBEX REQUIRED)
  ibex_init_common() # IBEX should have installed this function
  message(STATUS "Found IBEX version ${IBEX_VERSION}")
  message(STATUS "Found IBEX inclusion ${IBEX_INCLUDE_DIRS}")

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

#set(LIBS diff_inclusion)
set(EXECS intlie)

# Compilation

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_compile_options(-Wall)
#foreach(LIB1 ${LIBS})
#  add_library(${LIB1} STATIC diffincl.cpp diff_quadbox.cpp)
#  target_include_directories(${LIB1} PUBLIC ${CODAC_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
#  target_compile_options(${LIB1} PUBLIC ${CODAC_CFLAGS} ${PKG_CPAD_CFLAGS})
#  target_link_libraries(${LIB1} ${CODAC_LIBRARIES} Ibex::ibex ${PKG_CAPD_LDFLAGS} /usr/local/lib/libprim.a)
#endforeach(LIB1 ${LIBS})
foreach(EXECUTABLE ${EXECS})
  add_executable(${EXECUTABLE} ${EXECUTABLE}.cpp lieSE2.cpp)
  target_include_directories(${EXECUTABLE} PUBLIC ${CODAC_INCLUDE_DIRS} ${IBEX_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIRS})
  target_link_libraries(${EXECUTABLE} ${CODAC_LIBRARIES} Ibex::ibex ${LIBS} /usr/local/lib/libprim.a)
endforeach(EXECUTABLE ${EXECUTABLES})
