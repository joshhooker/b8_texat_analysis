cmake_minimum_required(VERSION 3.9)
project(b8_texat_analysis)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

#Find ROOT
find_package(ROOT REQUIRED)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
include_directories(${ROOT_INCLUDE_DIR})
link_directories(${ROOT_LIBRARY_DIR})

add_executable(b8_texat_analysis
        calcExcitationEnergy/calcExcitationEnergy.C
        calcExcitationEnergy/testCS.out
        b8_methane.dat
        CubicSpline.h
        cuts.root
        EnergyLoss.dat
        EnergyLoss.h
        gainFile.dat
        MakeChain.C
        proton_methane.dat
        runSpectra.sh
        siCalibration.dat
        Spectra.C
        Spectra.h)

target_link_libraries(b8_texat_analysis ${ROOT_LIBRARIES})
