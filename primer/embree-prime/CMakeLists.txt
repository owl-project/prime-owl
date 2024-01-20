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

option(PRIMER_BUILD_EMBREE "Build Embree-backend" OFF)
if (NOT PRIMER_BUILD_EMBREE)
  return()
endif()

find_package(Embree 3.13.5 REQUIRED)
if (NOT Embree_FOUND)
  return()
endif()

add_library(embree-prime STATIC
  Geom.h
  Triangles.h
  Triangles.cpp
  Group.h
  Group.cpp
  Model.h
  Model.cpp
  Context.h
  Context.cpp
  )

target_compile_definitions(embree-prime PUBLIC
  -DPRIMER_HAVE_EMBREE=1
  )
target_link_libraries(embree-prime
  PRIVATE
  primer-interface
  embree
  )
set_property(TARGET embree-prime
  PROPERTY POSITION_INDEPENDENT_CODE ON
  )

# aaand: link ourselves into primer
target_link_libraries(primer_dynamic PUBLIC embree-prime)
target_link_libraries(primer_static  PUBLIC embree-prime)
