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

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")

set(SCRIPTS
    averageBeamEnergy.out
    b8_methane.dat
    cuts.root
    EnergyLoss.dat
    gainFile.dat
    proton_methane.dat
    siCalibration.dat
    )

foreach(_script ${SCRIPTS})
    configure_file(
        ${PROJECT_SOURCE_DIR}/${_script}
        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${_script}
        COPYONLY
    )
endforeach()

add_executable(b8_texat_analysis
        Spectra.C Hough2D.C FitTrack.C)

target_link_libraries(b8_texat_analysis ${ROOT_LIBRARIES})