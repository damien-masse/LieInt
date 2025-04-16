# ==================================================================
#  LieInt | src configuation file (included from ..)
# ==================================================================

set(SRCS 
   intExp.cpp
   lieGroup.cpp
   lieSE2.cpp
   lieSO3.cpp
   intTFunc.cpp)

add_compile_options(-Wall)
add_library(${LIB} STATIC ${SRCS})
target_include_directories(${LIB} PUBLIC ${CODAC_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIRS} ${PKG_CAPD_INCLUDE_DIRS})
target_link_libraries(${LIB} ${CODAC_LIBRARIES} )
