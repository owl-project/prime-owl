// ======================================================================== //
// Copyright 2020-2021 Ingo Wald                                            //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#include "pyPrimer/Mesh.h"
#include "pyPrimer/Context.h"

namespace pypr {
  
  Mesh::Mesh(Context::SP context,
             int userID,
             const pybind11::buffer &vertex_buffer,
             const pybind11::buffer &index_buffer)
    : context(context)
  {
    assert(context);
    std::cout << "#pypr: creating mesh..." << std::endl;

    // 'parse' vertex buffer
    pybind11::buffer_info vertex_info = vertex_buffer.request();
    if (vertex_info.format != pybind11::format_descriptor<float>::format())
      throw std::runtime_error("vertex buffer is not in floats");

    pybind11::array_t<float> vertex_array
      = pybind11::cast<pybind11::array_t<float>>(vertex_buffer);
    int vertex_count = vertex_array.size() / 3;
    const vec3f *vertex = (const vec3f *)vertex_array.request().ptr;
    std::cout << "#pypr:  -> found " << vertex_count << " vertices" << std::endl;

    // 'parse' index buffer
    pybind11::buffer_info index_info = index_buffer.request();
    if (index_info.format != pybind11::format_descriptor<int>::format())
      throw std::runtime_error("index buffer is not in ints");

    pybind11::array_t<int> index_array
      = pybind11::cast<pybind11::array_t<int>>(index_buffer);
    int index_count = index_array.size() / 3;
    const vec3i *index = (vec3i *)index_array.request().ptr;
    std::cout << "#pypr:  -> found " << index_count << " indices" << std::endl;

    primer = opMeshCreate(context->primer,
                          userID,
                          (const float*)vertex,vertex_count,sizeof(vec3f),
                          (const int *)index,index_count,sizeof(vec3ui));
  }

  Mesh::~Mesh()
  {}
  
}
