cmake_minimum_required(VERSION 3.9)
project(b8_texat_analysis)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

#Find ROOT
find_package(ROOT REQUIRED)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
include_directories(${ROOT_INCLUDE_DIR})
link_directories(${ROOT_LIBRARY_DIR})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")

set(SCRIPTS
    averageBeamEnergy.dat
    b8_methane.dat
    centerGain.dat
    csiCalibration.dat
    csReg3.out
    cuts.root
    EnergyLoss.dat
    leftChainHits.out
    leftStripHits.out
    proton_methane.dat
    rightChainHits.out
    rightStripHits.out
    siCalibration.dat
    )

foreach(_script ${SCRIPTS})
    configure_file(
        ${PROJECT_SOURCE_DIR}/${_script}
        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${_script}
        COPYONLY
    )
endforeach()

IF(EXISTS "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/compile_commands.json")
    EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/compile_commands.json
        ${PROJECT_SOURCE_DIR}/compile_commands.json
    )
ENDIF()

add_executable(b8_texat_analysis
  Spectra.cpp Hough2D.cpp HoughTrack.cpp EnergyLoss.cpp)

target_link_libraries(b8_texat_analysis ${ROOT_LIBRARIES})
