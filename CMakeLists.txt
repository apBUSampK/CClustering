cmake_minimum_required(VERSION 3.22)
project(CClustering VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)

find_package(COLA REQUIRED)

include(CMakePackageConfigHelpers)

set(CMAKE_INSTALL_PREFIX ${COLA_DIR})
set(COLA_MODULE_NAME CClustering)

# Add library target
set(srcs MST/MSTClustering.h MST/MSTClustering.cpp)
add_library(CClustering SHARED ${srcs})

target_include_directories(CClustering PUBLIC
        "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/MST>"
        $<INSTALL_INTERFACE:include/${COLA_MODULE_NAME}>)

target_link_libraries(CClustering PUBLIC COLA)

# Specify public headers for installation
set(PROJECT_PUBLIC_HEADERS "MST/MSTClustering.h")

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
