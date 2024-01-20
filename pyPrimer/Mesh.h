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

#pragma once

#include "pyPrimer/common.h"

namespace pypr {
  
  struct Context;

  struct Mesh : public std::enable_shared_from_this<Mesh>
  {
    typedef std::shared_ptr<Mesh> SP;

    Mesh(std::shared_ptr<Context> context,
         int userID,
         const pybind11::buffer &vertex_buffer,
         const pybind11::buffer &index_buffer);
    
    virtual ~Mesh();

    std::shared_ptr<Context> context;
    OPGeom primer;
  };
  
}
