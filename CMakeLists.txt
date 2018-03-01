
cmake_minimum_required(VERSION 2.6)
project(biot)
 
#Bring the headers, such as Student.h into the project
include_directories(${CMAKE_SOURCE_DIR})
include_directories(/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/qhull/2015.2/include )
include_directories(/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/getfem/5.1/include)
include_directories(/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/getfem/5.1/include)
include_directories(${mkSuperluInc})
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DGMM_USES_BLAS -std=c++14" )
#-I-DGMM_USES_BLAS -I-I
#Can manually add the sources using the set command as follows:
#set(SOURCES src/mainapp.cpp src/Student.cpp)
#However, the file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES ${CMAKE_SOURCE_DIR}/*.cpp)
MESSAGE( " STATUS CMAKE SOURCE DIR:         " ${CMAKE_SOURCE_DIR} )
MESSAGE( " STATUS SOURCES:         " ${SOURCES} )
link_directories(/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/getfem/5.1/lib)
LINK_DIRECTORIES(/u/sw/pkgs/toolchains/gcc-glibc/5/pkgs/openblas/0.2.17/lib )
SET( SUPERLU_LIB  "$ENV{mkSuperluLib}" )
LINK_DIRECTORIES(${SUPERLU_LIB})
SET( QHULL_LIB  "$ENV{mkQhullLib}" )
LINK_DIRECTORIES(${QHULL_LIB})
add_executable(biot ${SOURCES})
TARGET_LINK_LIBRARIES(biot getfem qhull)
TARGET_LINK_LIBRARIES(biot openblas)
TARGET_LINK_LIBRARIES(biot superlu)

