cmake_minimum_required(VERSION 3.22)
project(CClustering VERSION 1.1)

set(CMAKE_CXX_STANDARD 17)

find_package(COLA REQUIRED)
find_package(Geant4 REQUIRED)

include(CMakePackageConfigHelpers)
include(${Geant4_USE_FILE})

set(CMAKE_INSTALL_PREFIX ${COLA_DIR})
set(COLA_MODULE_NAME CClustering)

# Add library target
set(srcs Wrapper/CClusteringFactory.cc MST/MSTClustering.cc MST/CoordinateMSTClustering.cc RepulsionStage/Repulsion.cc
        ExcitationEnergy/Ericson.cc ExcitationEnergy/GaimardSchmidt.cc ExcitationEnergy/ExcitationEnergy.cc )
add_library(CClustering SHARED ${srcs})

target_include_directories(CClustering PUBLIC
        "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/MST>"
        "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/RepulsionStage>"
        "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/ExcitationEnergy>"
        "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/Wrapper>"
        $<INSTALL_INTERFACE:include/${COLA_MODULE_NAME}>)

target_link_libraries(CClustering PUBLIC COLA ${Geant4_LIBRARIES})

# Specify public headers for installation
set(PROJECT_PUBLIC_HEADERS "Wrapper/CClusteringFactory.hh;Wrapper/CClustering.hh;RepulsionStage/Repulsion.hh;MST/MSTClustering.hh;MST/CoordinateMSTClustering.hh;\
ExcitationEnergy/Ericson.hh;ExcitationEnergy/ExcitationEnergy.hh;ExcitationEnergy/GaimardSchmidt.hh")

set_target_properties(CClustering PROPERTIES PUBLIC_HEADER "${PROJECT_PUBLIC_HEADERS}")

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
install(TARGETS CClustering
        EXPORT CClusteringExport
        LIBRARY DESTINATION lib/${COLA_MODULE_NAME}
        PUBLIC_HEADER DESTINATION include/${COLA_MODULE_NAME}
        INCLUDES DESTINATION include/${COLA_MODULE_NAME})

# Install export file and config files
install(EXPORT CClusteringExport
        DESTINATION lib/cmake/${COLA_MODULE_NAME})

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${COLA_MODULE_NAME}Config.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/${COLA_MODULE_NAME}ConfigVersion.cmake"
        DESTINATION lib/cmake/${COLA_MODULE_NAME})
