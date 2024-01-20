# ======================================================================== #
# Copyright 2019-2023 Ingo Wald                                            #
#                                                                          #
# Licensed under the Apache License, Version 2.0 (the "License");          #
# you may not use this file except in compliance with the License.         #
# You may obtain a copy of the License at                                  #
#                                                                          #
#     http://www.apache.org/licenses/LICENSE-2.0                           #
#                                                                          #
# Unless required by applicable law or agreed to in writing, software      #
# distributed under the License is distributed on an "AS IS" BASIS,        #
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. #
# See the License for the specific language governing permissions and      #
# limitations under the License.                                           #
# ======================================================================== #

cmake_minimum_required(VERSION 2.8.12)
cmake_policy(SET CMP0048 NEW)

set(CUDA_USE_STATIC_CUDA_RUNTIME ON)
set(CMAKE_CUDA_USE_STATIC_CUDA_RUNTIME ON)
set(CMAKE_USE_STATIC_CUDA_RUNTIME ON)
set(CMAKE_CUDA_USE_STATIC_RUNTIME ON)

project(Prime-OWL LANGUAGES C CXX CUDA VERSION 0.2.0)


set(CMAKE_CXX_STANDARD 14 CACHE STRING "")

if(NOT SET_UP_CONFIGURATIONS_DONE)
  set(SET_UP_CONFIGURATIONS_DONE 1)

  # No reason to set CMAKE_CONFIGURATION_TYPES if it's not a multiconfig generator
  # Also no reason mess with CMAKE_BUILD_TYPE if it's a multiconfig generator.
  if(CMAKE_CONFIGURATION_TYPES) # multiconfig generator?
    set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "" FORCE) 
  else()
    if(NOT CMAKE_BUILD_TYPE)
      #            message("Defaulting to release build.")
      set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
    endif()
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY HELPSTRING "Choose the type of build")
    # set the valid options for cmake-gui drop-down list
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release")
  endif()
  
  SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
  SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
  SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
endif()

# re policty cmp0104
if (NOT (DEFINED CMAKE_CUDA_ARCHITECTURES))
  set(CMAKE_CUDA_ARCHITECTURES OFF)
endif()


if (NOT (TARGET OWL))
  add_subdirectory(${PROJECT_SOURCE_DIR}/submodules/owl owl EXCLUDE_FROM_ALL)
endif()
 
enable_testing()

# ------------------------------------------------------------------
# the primer library itself, with both GPU and CPU backends
# ------------------------------------------------------------------
add_subdirectory(primer)

# ------------------------------------------------------------------
# python bindings for this library
# ------------------------------------------------------------------
option(PRIMER_ENABLE_PYTHON "Enable 'py-primer' python bindings?" ON)
if (PRIMER_ENABLE_PYTHON)
  add_subdirectory(pyPrimer)
endif()

# ------------------------------------------------------------------
# tutorial/samples
# ------------------------------------------------------------------
add_subdirectory(samples)
