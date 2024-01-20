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

#pragma once

#include "primer-common/Context.h"
#include "owl-prime/Triangles.h"
#include "owl-prime/Spheres.h"
#include "owl-prime/Model.h"
#include "owl-prime/Group.h"

#define CUDA_CALL(a) OWL_CUDA_CALL(a)
#define CUDA_SYNC_CHECK OWL_CUDA_SYNC_CHECK

namespace op {
  using namespace owl::common;
  
  inline bool isDeviceAccessible(const void *ptr)
  {
    cudaPointerAttributes attributes = {};
    // do NOT check for error: in CUDA<10, passing a host pointer will
    // actually create a cuda error, so let's just call and then
    // ignore/clear the error
    cudaPointerGetAttributes(&attributes,ptr);
    cudaGetLastError();

    return attributes.devicePointer != 0;
  }
  
  struct Context : public primer::Context {
    /*! launch parameters */
    struct LPData
    {
      OPTraceFlags flags;
      OptixTraversableHandle model;
      int          numRays;

      // ------------------------------------------------------------------
      // array-of-structs variants; when tracing soa rays there may be
      // null or invalid (do NOT put into a union - unions aren't
      // supported by owl variabels)
      // ------------------------------------------------------------------
      primer::Ray *rays;
      primer::Hit *hits;
      int         *activeIDs;
      /*! only used for shadow rays */
      // int        *isOccluded;

      // ------------------------------------------------------------------
      // soa variants; when tracing aos these may be null or invalid
      // (do NOT put into a union - unions aren't supported by owl
      // variabels). Note some of these _may_ be null even in soa mode
      // (eg, shadow ray will not set u or v, anyway, so may be null,
      // and even for regular rays the user may not want to know u or
      // v, etc)
      // ------------------------------------------------------------------
      struct {
        struct {
          float *org_x;
          float *org_y;
          float *org_z;
          float *dir_x;
          float *dir_y;
          float *dir_z;
          float *t_min;
          float *t_max;
        } ray;
        struct {
          int *primID;
          int *geomID;
          int *instID;
          int *t;
          int *u;
          int *v;
        } hit;
      } soa;
      
      int         isShadow;
      int         isSOA;
      /*! if 0, this will be a regular 'find closest hit' call;
        otherwise, this specifies how many hits we have reserved for
        each ray, and will return the N closest hits */
      int         numHitsPerRay;
    };
    
    Context(int gpuID);

    static OWLVarDecl lpVariables[];

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
                                   /* centers */
                                   const float *xArray, 
                                   const float *yArray, 
                                   const float *zArray,
                                   size_t strideCenters,
                                   /* radii */
                                   const float *rArray,
                                   size_t strideRadii) override;
    
    primer::Model *createModel(std::vector<OPInstance> &instances) override;
    
    primer::Group *createGroup(std::vector<OPGeom> &geoms) override;

    void checkSBT();

    bool        sbtDirty        = true;
    OWLContext  owl             = 0;
    OWLGeomType meshGeomType    = 0;
    OWLGeomType spheresGeomType = 0;
    OWLParams   launchParams    = 0;
    OWLRayGen   rayGen          = 0;
    OWLModule   module          = 0;
  };
  
}
