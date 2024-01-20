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

#include "embree-prime/Triangles.h"
#include "embree-prime/Context.h"

namespace ep {

  Triangles::Triangles(Context *context,
             uint32_t userID,
             const std::vector<vec3f> &vertices,
             const std::vector<vec3i> &indices)
    : Geom(context,userID),
      vertices(vertices),
      indices(indices)
  {
    this->geom  = rtcNewGeometry(context->device,RTC_GEOMETRY_TYPE_TRIANGLE);
    this->indicesBuffer
      = rtcNewSharedBuffer(context->device,
                           this->indices.data(),
                           this->indices.size()*sizeof(vec3i));
    rtcSetGeometryBuffer(this->geom, RTC_BUFFER_TYPE_INDEX, 0,
                         RTC_FORMAT_UINT3, this->indicesBuffer, 0,
                         sizeof(vec3i), this->indices.size());
    this->verticesBuffer
      = rtcNewSharedBuffer(context->device,
                           this->vertices.data(),
                           this->vertices.size()*sizeof(vec3f));
    rtcSetGeometryBuffer(this->geom, RTC_BUFFER_TYPE_VERTEX, 0,
                         RTC_FORMAT_FLOAT3, this->verticesBuffer, 0,
                         sizeof(vec3f), this->vertices.size());
    rtcSetGeometryUserData(this->geom,(void*)this);
    rtcCommitGeometry(this->geom);
  }

}
