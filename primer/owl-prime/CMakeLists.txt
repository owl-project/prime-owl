# ======================================================================== #
# Copyright 2019-2021 Ingo Wald                                            #
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

# builds device programs that the owl-prime device will use on the device
embed_ptx(
  OUTPUT_TARGET
  prime-owl-ptx
  
  PTX_LINK_LIBRARIES
  primer-interface
  
  SOURCES
  deviceCode.cu
  )

add_library(owl-prime STATIC
  Triangles.h
  Triangles.cpp
  Spheres.h
  Spheres.cpp
  Model.h
  Model.cpp
  Group.h
  Group.cpp
  Context.h
  Context.cu
  )
target_link_libraries(owl-prime
  prime-owl-ptx
  primer-interface
  )
set_target_properties(owl-prime
  PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  CUDA_USE_STATIC_CUDA_RUNTIME ON
  CUDA_USE_STATIC_RUNTIME ON
  )

# aaand: link ourselves into primer
target_link_libraries(primer_dynamic PUBLIC owl-prime)
target_link_libraries(primer_static  PUBLIC owl-prime)


