cmake_minimum_required(VERSION 3.7)
project(itebd)

set(CMAKE_CXX_STANDARD 11)
set(ITENSOR_BASE_DIR /Users/hyhui/ITensor)

include_directories(${ITENSOR_BASE_DIR})
add_library(itensor STATIC IMPORTED)
set_property(TARGET itensor PROPERTY IMPORTED_LOCATION ${ITENSOR_BASE_DIR}/lib/libitensor.a)
set(ITENSOR_FLAGS "-DPLATFORM_lapack -D__ASSERT_MACROS_DEFINE_VERSIONS_WITHOUT_UNDERSCORES=0")
set(ITENSOR_LINK_FLAGS "-DPLATFORM_lapack -L${ITENSOR_BASE_DIR}/lib -litensor -lpthread -L/usr/lib -lblas -llapack")

set(SOURCE_FILES cond-mat.0605597.cpp)
add_executable(itebd ${SOURCE_FILES})
target_link_libraries(itebd itensor "${ITENSOR_LINK_FLAGS}")