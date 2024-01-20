// ======================================================================== //
// Copyright 2019-2023 Ingo Wald                                            //
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

#include "embree-prime/Context.h"
#include "embree-prime/Group.h"
#include "embree-prime/Triangles.h"
#include "embree-prime/Spheres.h"

namespace ep {

  Group::Group(Context *context,
               std::vector<OPGeom> &geoms)
    : primer::Group(geoms)
  {
    PING;
    this->handle = rtcNewScene(context->device);
    PING;
    PRINT(geoms.size());
    for (auto _geom : geoms) {
      ep::Geom *geom = (ep::Geom *)_geom;
      PING;
      PRINT(geom->geom);
      PRINT(geom->userID);
      rtcAttachGeometryByID(this->handle,geom->geom,geom->userID);
      PING;
      rtcEnableGeometry(geom->geom);
      PING;
      rtcReleaseGeometry(geom->geom);
      PING;
    }
    PING;
    rtcCommitScene(this->handle);
    PING;
  }
  
  Group::~Group()
  {}
  
} // ::op
