// ======================================================================== //
// Copyright 2019-2021 Ingo Wald                                            //
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

#include "primer-common/Context.h"
#include "embree-prime/Triangles.h"
#include "embree-prime/Spheres.h"
#include "embree-prime/Model.h"

namespace ep {
  using namespace owl::common;
  
  struct Context : public primer::Context {

    Context();
    
    /*! create a mesh from 'fat' triangles (where each triangles is
      made up of 9 floats for the three vertex coordinates; no
      indexing) */
    primer::Geom *createTriangles(uint32_t userID,
                             size_t numTriangles,
                             const float *v0x, 
                             const float *v0y, 
                             const float *v0z, 
                             const float *v1x, 
                             const float *v1y, 
                             const float *v1z, 
                             const float *v2x, 
                             const float *v2y, 
                             const float *v2z,
                             size_t strideInBytes) override;

    /*! create a mesh from vertex array and index array */
    primer::Geom *createTriangles(uint32_t userID,
                             /* vertex array */
                             const vec3f *vertices,
                             size_t numVertices,
                             size_t sizeOfVertexInBytes,
                             /* index array */
                             const vec3i *indices,
                             size_t numIndices,
                             size_t sizeOfIndesStructInBytes) override;

    /*! creates a new spheres geometry from the given centers and radii */
    primer::Geom *createSpheres(uint32_t userID,
                                   int numSpheres,
                                   const float *xArray, 
                                   const float *yArray, 
                                   const float *zArray,
                                   size_t strideCenters,
                                   const float *rArray,
                                   size_t strideRadii) override;
                             
    
    primer::Model *createModel(std::vector<OPInstance> &instances) override;

    primer::Group *createGroup(std::vector<OPGeom> &geoms) override;
    
    RTCDevice device = 0;
  };
  
}
