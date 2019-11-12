PROJECT(CGRAPH-LIBTIM)
CMAKE_MINIMUM_REQUIRED(VERSION 2.7)

set(OSX_SDK_ROOT /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs)

set(TEST_OSX_SYSROOT /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk)
set(CMAKE_OSX_ARCHITECTURES "x86_64" CACHE STRING "Force build for 64-bit" FORCE)
set(CMAKE_OSX_DEPLOYMENT_TARGET "${SDK_VERSION}" CACHE STRING "Force build for 64-bit" FORCE)

set(CMAKE_OSX_SYSROOT "${TEST_OSX_SYSROOT}" CACHE PATH "Force build for 64-bit" FORCE)

set(LD_LIBRARY_PATH ./)

### LIBRARIES ###
## LIBRARY STATIC NIFTI
include_directories(Nifti/include)
add_subdirectory(Nifti)

### FILES ###
file(GLOB SRC_FILES sources/src/*.cpp)
include_directories(sources/inc)

### EXECUTABLE ###
# configure an executable
ADD_EXECUTABLE(shaping_elongM ${SRC_FILES})
TARGET_LINK_LIBRARIES(shaping_elongM Nifti)
