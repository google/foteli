# Copyright 2024 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

cmake_minimum_required(VERSION 3.10)

project(Foteli
  LANGUAGES CXX
  VERSION 1.0.0
  DESCRIPTION "Local tone-mapping for HDR images")

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

find_package(PkgConfig)
if (PKG_CONFIG_FOUND)
  pkg_check_modules(OIIO IMPORTED_TARGET OpenImageIO)
endif ()

add_library(foteli foteli.h foteli.cc)
set_target_properties(foteli
  PROPERTIES
    IMPORTED_LINK_INTERFACE_LANGUAGES CXX
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED YES)

include(GenerateExportHeader)
generate_export_header(foteli)
target_include_directories(foteli PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/foteli>")

if (OIIO_FOUND)
  add_executable(foteli_main main.cc)
  target_link_libraries(foteli_main foteli PkgConfig::OIIO)
else ()
  message(WARNING "OpenImageIO not found; not building the command-line executable")
endif ()

find_package(Doxygen)
if (DOXYGEN_FOUND)
  doxygen_add_docs(doc foteli.h)
endif ()

install(FILES foteli.h "${CMAKE_CURRENT_BINARY_DIR}/foteli_export.h"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/foteli")

install(TARGETS foteli EXPORT FoteliTargets)
install(EXPORT FoteliTargets DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Foteli")
configure_package_config_file(
  cmake/FoteliConfig.cmake.in
  FoteliConfig.cmake
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Foteli")
write_basic_package_version_file(
  FoteliConfigVersion.cmake
  COMPATIBILITY AnyNewerVersion)
install(FILES
          "${CMAKE_CURRENT_BINARY_DIR}/FoteliConfig.cmake"
          "${CMAKE_CURRENT_BINARY_DIR}/FoteliConfigVersion.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Foteli")

if (TARGET foteli_main)
  install(PROGRAMS "$<TARGET_FILE:foteli_main>" TYPE BIN RENAME foteli)
endif ()
