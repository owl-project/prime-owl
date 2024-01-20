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
#include "embree-prime/Triangles.h"
#include "embree-prime/Group.h"

primer::Context *primer::Context::createEmbreeContext()
{ return new ep::Context; }

namespace ep {
  
  Context::Context()
  {
    device = rtcNewDevice("");
    LOG("created new embree device " << (int*)device);
  }
  
  template<typename T>
  const T &getWithOffset(const T *base, int idx, size_t strideInBytes)
  {
    unsigned char *ptr = (unsigned char *)base;
    ptr += idx * strideInBytes;
    return *(T*)ptr;
  }
  
  /*! creates a new spheres geometry from the given centers and radii */
  primer::Geom *Context::createSpheres(uint32_t userID,
                                        int numSpheres,
                                        const float *xArray, 
                                        const float *yArray, 
                                        const float *zArray,
                                        size_t strideCenters,
                                        const float *rArray,
                                        size_t strideRadii)
  {
    throw std::runtime_error("spheres not implemented in embree backend");
  }
                             
  primer::Geom *Context::createTriangles(uint32_t userID,
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
                                    size_t strideInBytes)
  {
    std::vector<vec3f> vertices(3*numTriangles);
    for (int i=0;i<numTriangles;i++) {
      vertices[3*i+0] = {
                         getWithOffset(v0x,i,strideInBytes),
                         getWithOffset(v0y,i,strideInBytes),
                         getWithOffset(v0z,i,strideInBytes)
      };
      vertices[3*i+1] = {
                         getWithOffset(v1x,i,strideInBytes),
                         getWithOffset(v1y,i,strideInBytes),
                         getWithOffset(v1z,i,strideInBytes)
      };
      vertices[3*i+2] = {
                         getWithOffset(v2x,i,strideInBytes),
                         getWithOffset(v2y,i,strideInBytes),
                         getWithOffset(v2z,i,strideInBytes)
      };
    }
    
    std::vector<vec3i> indices(numTriangles);
    for (int i=0;i<numTriangles;i++) {
      indices[i] = vec3i(3*i)+vec3i{0,1,2};
    }

    return new Triangles(this,userID,vertices,indices);
  }

  /*! create a mesh from vertex array and index array */
  primer::Geom *Context::createTriangles(uint32_t userID,
                                    /* vertex array */
                                    const vec3f *vertices,
                                    size_t numVertices,
                                    size_t vertexStrideInBytes,
                                    /* index array */
                                    const vec3i *indices,
                                    size_t numIndices,
                                    size_t indexStrideInBytes)
  {
    std::vector<vec3f> d_vertices(numVertices);
    std::vector<vec3i> d_indices(numIndices);
    
    for (int i=0;i<numVertices;i++) 
      d_vertices[i] = getWithOffset(vertices,i,vertexStrideInBytes);
    for (int i=0;i<numIndices;i++) 
      d_indices[i] = getWithOffset(indices,i,indexStrideInBytes);
    
    return new Triangles(this,userID,d_vertices,d_indices);
  }

  primer::Group *Context::createGroup(std::vector<OPGeom> &geoms) 
  {
    return new ep::Group(this,geoms);
  }

  
  primer::Model *Context::createModel(std::vector<OPInstance> &instances)
  {
    return new ep::Model(this,instances);
  }
  
}
