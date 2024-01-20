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
#include <owl/common/math/box.h>
#include "pyPrimer/Mesh.h"

// #include "pyOWL/Module.h"
// #include "pyOWL/GeomType.h"
// #include "pyOWL/Buffer.h"
// #include "pyOWL/Geom.h"
// #include "pyOWL/Group.h"
// #include "pyOWL/MissProg.h"
// #include "pyOWL/RayGen.h"

namespace pypr {

  struct Ray {
    float org_x;
    float org_y;
    float org_z;
    float t_min;
    float dir_x;
    float dir_y;
    float dir_z;
    float t_max;
  };

  struct Hit {
    int primID;
    int geomID;
    int instID;
    float t;
    float u;
    float v;
  };
  
  
  struct Context : public std::enable_shared_from_this<Context> {
    typedef std::shared_ptr<Context> SP;

    Context();

    virtual ~Context();

    OPContext primer = 0;
  };


  /* just a helper function for the samples */
  // void save_png_rgba8(Buffer::SP buffer,
  //                     const std::vector<int> &fbSize,
  //                     const std::string &fileName);
}
