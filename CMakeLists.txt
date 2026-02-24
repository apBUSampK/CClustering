cmake_minimum_required(VERSION 3.22)
project(CPreEquilibrium VERSION 1.1)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_STANDARD 17)

find_package(COLA REQUIRED)
find_package(Geant4 REQUIRED)

include(CMakePackageConfigHelpers)
include(${Geant4_USE_FILE})

set(CMAKE_INSTALL_PREFIX ${COLA_DIR})
set(COLA_MODULE_NAME CPreEquilibrium)

# Add library target
set(srcs Wrapper/CPreEquilibriumFactory.cc MSTClustering/MSTClustering.cc MSTClustering/CoordinateMSTClustering.cc
        ExcitationEnergy/Ericson.cc ExcitationEnergy/GaimardSchmidt.cc ExcitationEnergy/ExcitationEnergy.cc)

add_library(CPreEquilibrium SHARED ${srcs})

target_include_directories(CPreEquilibrium PUBLIC
        "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/MSTClustering>"
        "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/RepulsionStage>"
        "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/ExcitationEnergy>"
        "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/Wrapper>"
        $<INSTALL_INTERFACE:include/${COLA_MODULE_NAME}>)

target_link_libraries(CPreEquilibrium PUBLIC COLA ${Geant4_LIBRARIES})

# Specify public headers for installation
set(PROJECT_PUBLIC_HEADERS "Wrapper/CPreEquilibriumFactory.hh;Wrapper/CPreEquilibrium.hh;RepulsionStage/Repulsion.hh;MSTClustering/MSTClustering.hh;MSTClustering/CoordinateMSTClustering.hh;\
ExcitationEnergy/Ericson.hh;ExcitationEnergy/ExcitationEnergy.hh;ExcitationEnergy/GaimardSchmidt.hh")

set_target_properties(CPreEquilibrium PROPERTIES PUBLIC_HEADER "${PROJECT_PUBLIC_HEADERS}")

# Add installation path for library to find .dat files:
add_compile_definitions(DATA_INSTALL="${CMAKE_INSTALL_PREFIX}/data/${COLA_MODULE_NAME}/")

# Configure config files
configure_package_config_file(
        "${PROJECT_SOURCE_DIR}/cmake/${COLA_MODULE_NAME}Config.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/${COLA_MODULE_NAME}Config.cmake"
        INSTALL_DESTINATION lib/cmake/${COLA_MODULE_NAME}
)

write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/${COLA_MODULE_NAME}ConfigVersion.cmake
        COMPATIBILITY AnyNewerVersion)

# Install library
install(TARGETS CPreEquilibrium
        EXPORT CPreEquilibriumExport
        LIBRARY DESTINATION lib/${COLA_MODULE_NAME}
        PUBLIC_HEADER DESTINATION include/${COLA_MODULE_NAME}
        INCLUDES DESTINATION include/${COLA_MODULE_NAME})

# Install export file and config files
install(EXPORT CPreEquilibriumExport
        DESTINATION lib/cmake/${COLA_MODULE_NAME})

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${COLA_MODULE_NAME}Config.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/${COLA_MODULE_NAME}ConfigVersion.cmake"
        DESTINATION lib/cmake/${COLA_MODULE_NAME})

#install .dat files for excitation energy calculation
install(FILES ${PROJECT_SOURCE_DIR}/ExcitationEnergy/CorrectedALADINParameters.dat
        DESTINATION data/${COLA_MODULE_NAME})