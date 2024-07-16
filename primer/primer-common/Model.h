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

#include "primer-common/Ray.h"
#include "primer-common/Context.h"
#include "primer-common/StridedArray.h"

namespace primer {
  using namespace owl::common;
  
  struct Model {

    // Model(std::vector<OPInstance> &instances)
    //   : instances(instances)
    // {}
    
    virtual void trace(Ray *rays,
                       Hit *hits,
                       int  numRaysAndHits,
                       int *activeRayIndices,
                       int  numActice,
                       OPTraceFlags flags) = 0;

    /*! finds the (up to) N closest hits along a ray. For each ray the
        hits array stores N "slots" for hits; teh trace will fill in
        the N closest hits. If there are less than N hits along a ray,
        the next hit slot after the last found one will have its
        primID set to -1 */
    virtual void traceN(Ray *rays,
                        Hit *hits,
                        int  numRays,
                        int  numHitsPerRay,
                        OPTraceFlags flags) = 0;
    
    /*! the SoA variant that operates on arbitrary strided arrays -
        need to have this in the backend (not just rearrange on
        front-end) because the arrays might not even be accessible on
        the host, so re-arrangement has to be done whatever device the
        user created. Any of the arrays other than org and dir may be
        null (in which case their contents get set to defaluts (for
        rays) or simply not written to (for hit data)); but arrays
        that re not null must have the given number of elemnets in
        them. */
    virtual void trace(size_t numRays,
                       const float *ray_org_x, size_t stride_ray_org_x,
                       const float *ray_org_y, size_t stride_ray_org_y,
                       const float *ray_org_z, size_t stride_ray_org_z,
                       const float *ray_dir_x, size_t stride_ray_dir_x,
                       const float *ray_dir_y, size_t stride_ray_dir_y,
                       const float *ray_dir_z, size_t stride_ray_dir_z,
                       const float *ray_tmin, size_t stride_ray_tmin,
                       const float *ray_tmax, size_t stride_ray_tmax,
                       int *hit_instID, size_t stride_hit_instID,
                       int *hit_geomID, size_t stride_hit_geomID,
                       int *hit_primID, size_t stride_hit_primID,
                       float *hit_t, size_t stride_hit_t,
                       float *hit_u, size_t stride_hit_u,
                       float *hit_v, size_t stride_hit_v,
                       OPTraceFlags flags) = 0;
    
    virtual void build() {}

    // const std::vector<OPInstance> instances;    
    Context *context;
  };
  
}
