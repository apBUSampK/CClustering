#  Copyright (c) 2026 Savva Savenkov, Artemii Novikov

#  Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
#  associated documentation files (the "Software"), to deal in the Software without restriction,
#  including without limitation the rights to use, copy, modify, merge, publish, distribute,
#  sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:

#  The above copyright notice and this permission notice shall be included in all copies or
#  substantial portions of the Software.

#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
#  NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
#  DAMAGES OR OTHEROUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#  SOFTWARE.

cmake_minimum_required(VERSION 3.22)
project(
    CPreEquilibrium
    VERSION 1.1
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)

find_package(COLA REQUIRED)
find_package(Geant4 REQUIRED)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(${Geant4_USE_FILE})

set(CMAKE_INSTALL_PREFIX ${COLA_DIR})
set(COLA_MODULE_NAME CPreEquilibrium)

# Add library target
file(
    GLOB
    srcs
    CONFIGURE_DEPENDS
    Wrapper/*.cc
    RepulsionStage/*.cc
    MSTClustering/*.cc
    ExcitationEnergy/*.cc
)

add_library(CPreEquilibrium SHARED ${srcs})

target_compile_features(CPreEquilibrium PRIVATE cxx_std_17)

target_include_directories(
    CPreEquilibrium
    PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/MSTClustering>"
           "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/RepulsionStage>"
           "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/ExcitationEnergy>"
           "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/Wrapper>"
           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${COLA_MODULE_NAME}>
)

target_link_libraries(CPreEquilibrium PUBLIC COLA ${Geant4_LIBRARIES})

# Specify public headers for installation
set(PROJECT_PUBLIC_HEADERS
    "Wrapper/CPreEquilibriumFactory.hh;Wrapper/CPreEquilibriumModule.hh;RepulsionStage/Repulsion.hh;MSTClustering/MSTClustering.hh;MSTClustering/CoordinateMSTClustering.hh;\
ExcitationEnergy/ExcitationEnergy.hh"
)

set_target_properties(
    CPreEquilibrium
    PROPERTIES PUBLIC_HEADER "${PROJECT_PUBLIC_HEADERS}"
               VERSION "${PROJECT_VERSION}"
               SOVERSION "${PROJECT_VERSION_MAJOR}"
               INSTALL_RPATH_USE_LINK_PATH TRUE
               BUILD_WITH_INSTALL_RPATH TRUE
)

# Check for COLA_FETCH_MODE to prevent install instructions from generating and change data path.
if (NOT DEFINED ENV{COLA_FETCH_MODE})
    target_compile_definitions(
        CPreEquilibrium PRIVATE DATA_INSTALL="${CMAKE_INSTALL_PREFIX}/data/${COLA_MODULE_NAME}/"
    )

    # Configure config files
    configure_package_config_file(
        "${PROJECT_SOURCE_DIR}/cmake/${COLA_MODULE_NAME}Config.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/${COLA_MODULE_NAME}Config.cmake"
        INSTALL_DESTINATION lib/cmake/${COLA_MODULE_NAME}
    )

    write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/${COLA_MODULE_NAME}ConfigVersion.cmake
        VERSION ${CPreEquilibrium_Version}
        COMPATIBILITY AnyNewerVersion
    )

    # Install library
    install(
        TARGETS CPreEquilibrium
        EXPORT CPreEquilibriumExport
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/${COLA_MODULE_NAME}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${COLA_MODULE_NAME}
        INCLUDES
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${COLA_MODULE_NAME}
    )

    # Install export file and config files
    install(EXPORT CPreEquilibriumExport
            DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${COLA_MODULE_NAME}
    )

    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${COLA_MODULE_NAME}Config.cmake"
                  "${CMAKE_CURRENT_BINARY_DIR}/${COLA_MODULE_NAME}ConfigVersion.cmake"
            DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${COLA_MODULE_NAME}
    )

    #install .dat files for excitation energy calculation
    install(FILES ${PROJECT_SOURCE_DIR}/ExcitationEnergy/CorrectedALADINParameters.dat
            DESTINATION data/${COLA_MODULE_NAME}
    )

    setup_cola_quality_tools()
else ()
    target_compile_definitions(
        CPreEquilibrium PRIVATE DATA_INSTALL="${PROJECT_SOURCE_DIR}/ExcitationEnergy/"
    )
endif ()
