# ======================================================================== #
# Copyright 2020-2023 Ingo Wald                                            #
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

find_package(PythonLibs 3)
find_package(Python 3 COMPONENTS Interpreter Development)

if (Python_FOUND)
  add_subdirectory(../submodules/pybind11 pybind11_build EXCLUDE_FROM_ALL)
  
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION pybind11::thin_lto)

  pybind11_add_module(pypr
    Context.cpp
    Mesh.cpp
    Model.cpp
    bindings.cpp
    )
  target_link_libraries(pypr
    PRIVATE
    primer_static
    )

  add_executable(fromArrays_generateReference
    fromArrays_generateReference.cpp
    )
  target_link_libraries(fromArrays_generateReference
    primer_static
    )
else()
  message("#primer: could not find python developer packages, disabling py-primer python bindings")
endif()
