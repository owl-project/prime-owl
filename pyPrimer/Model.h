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
#include "Mesh.h"

namespace pypr {

  struct Context;

  struct Model : public std::enable_shared_from_this<Model>
  {
    typedef std::shared_ptr<Model> SP;

    Model(std::shared_ptr<Context> context,
          const std::vector<Mesh::SP> &meshes);

    /*! this is the "SOA"-like version, where origin, direction, tmin, mtax are all separate arrays */
    void trace(const pybind11::buffer &orgs,
               const pybind11::buffer &dirs,
               const pybind11::buffer &tmins,
               const pybind11::buffer &tmaxs,
               const pybind11::buffer &hitIDs,
               const pybind11::buffer &hitCoords
               );
    /*! this is the more "native"-version AoS variant of calling
      trace, where all rays are specified in a single array with 8
      floats per ray, and hits are returned with two arrays (one with
      three ints per hit for the hit IDs, the other with three floats
      per hit for (t,u,v) coordinates. In the rays array, each ray is
      represented as 8 consecutive floats, in the format of
      [org_x,org_y,org_z,tMin,dir_x,dir_y,dir_z,tMax]. In the hitIDs
      array it's three ints meaning [primID,geomID,instID] */
    void traceRays(const pybind11::buffer &rays,
                   const pybind11::buffer &hits);
      
    std::shared_ptr<Context> context;
    OPModel primer = 0;
  };
  
}
