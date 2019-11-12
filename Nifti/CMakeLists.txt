project(Nifti)
cmake_minimum_required(VERSION 2.7)

set(OSX_SDK_ROOT /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs)

set(SDK_VERSIONS_TO_CHECK 10.10)
set(CMAKE_OSX_ARCHITECTURES "x86_64" CACHE STRING "Force build for 64-bit" FORCE)
set(CMAKE_OSX_DEPLOYMENT_TARGET "${SDK_VERSION}" CACHE STRING "Force build for 64-bit" FORCE)
set(CMAKE_OSX_SYSROOT "${TEST_OSX_SYSROOT}" CACHE PATH "Force build for 64-bit" FORCE)

set(NIFTI_SOURCES ${PROJECT_SOURCE_DIR}/include/nifti1_io.c
${PROJECT_SOURCE_DIR}/include/znzlib.c)

set( NIFTI_HEADERS ${PROJECT_SOURCE_DIR}/include/nifti/nifti1.h
${PROJECT_SOURCE_DIR}/include/nifti1_io.h
${PROJECT_SOURCE_DIR}/include/znzlib.h)

add_library(Nifti STATIC ${NIFTI_SOURCES})

install( FILES ${NIFTI_HEADERS} DESTINATION include/Nifti)
install( TARGETS Nifti DESTINATION ${CMAKE_BINARY_DIR}/lib)
