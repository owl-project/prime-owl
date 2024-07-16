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

#include "owl-prime/Context.h"
#include "owl-prime/Triangles.h"
#include "owl-prime/Spheres.h"
#include "owl-prime/Group.h"

primer::Context *primer::Context::createOffloadContext(int gpuID)
{ return new op::Context(gpuID); }

extern "C" char deviceCode_ptx[];

namespace op {
  
  OWLVarDecl Context::lpVariables[]
  = {
     { "rays",      OWL_RAW_POINTER, OWL_OFFSETOF(Context::LPData,rays) },
     { "hits",      OWL_RAW_POINTER, OWL_OFFSETOF(Context::LPData,hits) },
     { "activeIDs", OWL_RAW_POINTER, OWL_OFFSETOF(Context::LPData,activeIDs) },
     { "model",     OWL_GROUP,       OWL_OFFSETOF(Context::LPData,model) },
     { "numRays",   OWL_INT,         OWL_OFFSETOF(Context::LPData,numRays) },
     { "isSOA",     OWL_INT,         OWL_OFFSETOF(Context::LPData,isSOA) },
     { "numHitsPerRay",   OWL_INT,   OWL_OFFSETOF(Context::LPData,numHitsPerRay) },
     { "flags",     OWL_ULONG,       OWL_OFFSETOF(Context::LPData,flags) },
     { nullptr /* end of list sentinel */ }
  };

  Context::Context(int gpuID)
  {
    if (gpuID < 0) gpuID = 0;
    
    owl = owlContextCreate(&gpuID,1);
    module = owlModuleCreate(owl,deviceCode_ptx);
    rayGen = owlRayGenCreate(owl,module,"traceRays",sizeof(int),nullptr,0);
    launchParams = owlParamsCreate(owl,sizeof(LPData),lpVariables,-1);
    
    meshGeomType = owlGeomTypeCreate(owl,OWL_TRIANGLES,
                                     sizeof(Triangles::SBTData),
                                     Triangles::variables,-1);
    owlGeomTypeSetClosestHit(meshGeomType,0,module,"Triangles");
    owlGeomTypeSetAnyHit(meshGeomType,0,module,"MultiHit");

    spheresGeomType = owlGeomTypeCreate(owl,OWL_GEOM_USER,
                                        sizeof(Spheres::SBTData),
                                        Spheres::variables,-1);
    owlGeomTypeSetBoundsProg(spheresGeomType,module,"Spheres");
    owlGeomTypeSetIntersectProg(spheresGeomType,0,module,"Spheres");
    owlGeomTypeSetClosestHit(spheresGeomType,0,module,"Spheres");
    owlGeomTypeSetAnyHit(spheresGeomType,0,module,"MultiHit");

    owlBuildPrograms(owl);
    owlBuildPipeline(owl);
  }
  
  void Context::checkSBT()
  {
    if (!sbtDirty) return;

    owlBuildSBT(owl);
    owlBuildPipeline(owl);
    sbtDirty = false;
  }

  template<typename T>
  inline __both__
  const T &getWithOffset(const T *base, int idx, size_t strideInBytes)
  {
    unsigned char *ptr = (unsigned char *)base;
    ptr += idx * strideInBytes;
    return *(T*)ptr;
  }

  __global__ void copySpheres(float4 *spheres,
                              int numSpheres,
                              const float *x,
                              const float *y,
                              const float *z,
                              int centerStride,
                              const float *r,
                              int rStride)
  {
    int tid = threadIdx.x+blockIdx.x*blockDim.x;
    if (tid >= numSpheres) return;
    spheres[tid]
      = {
         getWithOffset(x,tid,centerStride),
         getWithOffset(y,tid,centerStride),
         getWithOffset(z,tid,centerStride),
         getWithOffset(r,tid,rStride)
    };
  }

  /*! creates a new spheres geometry from the given centers and radii */
  primer::Geom *Context::createSpheres(uint32_t userID,
                                       int numSpheres,
                                       /* vertex array */
                                       const float *xArray, 
                                       const float *yArray, 
                                       const float *zArray,
                                       size_t strideCenters,
                                       const float *rArray,
                                       size_t strideRadii)
  {
    OWLBuffer spheresBuffer = 0;
    if (isDeviceAccessible(xArray)) {
      spheresBuffer = owlDeviceBufferCreate(owl,OWL_FLOAT4,numSpheres,0);
      copySpheres<<<divRoundUp(numSpheres,1024),1024>>>
        ((float4*)owlBufferGetPointer(spheresBuffer,0),
         numSpheres,
         xArray,
         yArray,
         zArray,strideCenters,
         rArray,strideRadii);
    } else {
      std::vector<float4> aos(numSpheres);
      for (int i=0;i<numSpheres;i++) {
        aos[i] = {
                  getWithOffset(xArray,i,strideCenters),
                  getWithOffset(yArray,i,strideCenters),
                  getWithOffset(zArray,i,strideCenters),
                  getWithOffset(rArray,i,strideRadii)
        };
      }
      spheresBuffer = owlDeviceBufferCreate(owl,OWL_FLOAT4,numSpheres,aos.data());
    }
    return new Spheres(this,userID,
                       spheresBuffer,numSpheres);
    
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
    OWLBuffer vertexBuffer
      = owlManagedMemoryBufferCreate(owl,OWL_FLOAT3,3*numTriangles,0);
    vec3f *vertices = (vec3f*)owlBufferGetPointer(vertexBuffer,0);
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
    
    OWLBuffer indexBuffer
      = owlManagedMemoryBufferCreate(owl,OWL_INT3,numTriangles,0);
    
    vec3i *indices = (vec3i*)owlBufferGetPointer(indexBuffer,0);
    // TODO: do that in CUDA (so would also work on device pointers)
    for (int i=0;i<numTriangles;i++) {
      indices[i] = vec3i(3*i)+vec3i{0,1,2};
    }

    return new Triangles(this,userID,
                         vertexBuffer,3*numTriangles,
                         indexBuffer,numTriangles);
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
    // TODO: do all this without copies if these are already device pointers
    OWLBuffer vertexBuffer
      = owlManagedMemoryBufferCreate(owl,OWL_FLOAT3,numVertices,0);
    vec3f *d_vertices = (vec3f*)owlBufferGetPointer(vertexBuffer,0);
    
    OWLBuffer indexBuffer
      = owlManagedMemoryBufferCreate(owl,OWL_INT3,numIndices,0);
    vec3i *d_indices = (vec3i*)owlBufferGetPointer(indexBuffer,0);

    for (int i=0;i<numVertices;i++) 
      d_vertices[i] = getWithOffset(vertices,i,vertexStrideInBytes);
    for (int i=0;i<numIndices;i++) 
      d_indices[i] = getWithOffset(indices,i,indexStrideInBytes);
    
    return new Triangles(this,userID,
                         vertexBuffer,numVertices,
                         indexBuffer,numIndices);
  }

  primer::Group *Context::createGroup(std::vector<OPGeom> &geoms) 
  {
    return new op::Group(this,geoms);
  }

  primer::Model *Context::createModel(const std::vector<OPGroup>  &groups,
                                      const std::vector<affine3f> &xfms,
                                      const std::vector<int>      &userIDs)
  {
    return new op::Model(this,groups,xfms,userIDs);
  }
  
} // ::op
