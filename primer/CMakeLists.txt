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


# ------------------------------------------------------------------
# include owl-based implementation of primer API
# ------------------------------------------------------------------
add_library(primer-interface INTERFACE)
target_include_directories(primer-interface INTERFACE
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}/primer/include
  ${CMAKE_CURRENT_LIST_DIR}
  )
target_link_libraries(primer-interface INTERFACE owl::owl)

# ------------------------------------------------------------------
# build abstract primer abstraction layer that implements the API.
# ------------------------------------------------------------------
set(SOURCES
  primer-common/Geom.h
  primer-common/Group.h
  primer-common/Context.h
  primer-common/primer.cpp
  )

add_library(primer_static STATIC ${SOURCES})
target_link_libraries(primer_static PUBLIC primer-interface)
set_target_properties(primer_static
  PROPERTIES
  CUDA_USE_STATIC_RUNTIME ON
  CUDA_SEPARABLE_COMPILATION ON
  POSITION_INDEPENDENT_CODE ON
  CUDA_RESOLVE_DEVICE_SYMBOLS ON
  CUDA_USE_STATIC_CUDA_RUNTIME ON
#  CUDA_ARCHITECTURES OFF
  )
target_link_libraries(primer_static PUBLIC CUDA::cuda_driver CUDA::cudart_static)

add_library(primer_dynamic SHARED ${SOURCES})
target_link_libraries(primer_dynamic PUBLIC primer-interface)
set_target_properties(primer_dynamic
  PROPERTIES
  CUDA_USE_STATIC_RUNTIME ON
  CUDA_SEPARABLE_COMPILATION ON
  POSITION_INDEPENDENT_CODE ON
  CUDA_RESOLVE_DEVICE_SYMBOLS ON
  CUDA_USE_STATIC_CUDA_RUNTIME ON
#  CUDA_ARCHITECTURES OFF
  )
find_package(CUDAToolkit REQUIRED)
target_link_libraries(primer_dynamic PUBLIC CUDA::cuda_driver CUDA::cudart_static)

# allow actual backend(s) in subdirector(ies) to link to this library
# in current directory
cmake_policy(SET CMP0079 NEW)

# ------------------------------------------------------------------
# include owl-based implementation of primer API. this will build a
# owl-based context that links itself to primer and creates an actual,
# owl-based implementation of this API
# ------------------------------------------------------------------
add_subdirectory(owl-prime)

# ------------------------------------------------------------------
# include embree-based cpu fallback for primer API. this will build a
# embree-based host device; this will _not_ work with device-side
# arrays.
# ------------------------------------------------------------------
add_subdirectory(embree-prime)

