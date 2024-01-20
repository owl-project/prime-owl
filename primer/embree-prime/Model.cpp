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

#include "embree-prime/Model.h"
#include "embree-prime/Group.h"
#include "embree-prime/Context.h"
#include "primer-common/StridedArray.h"

#define EP_NOTIMPLEMENTED throw std::runtime_error("#primer: this function isn't implemented yet for the embree backend")

namespace ep {
  using primer::StridedArray;
  
  Model::Model(Context *context,
               std::vector<OPInstance> &instances)
    : primer::Model(instances),
      context(context)
  {
    this->handle = rtcNewScene(context->device);
    for (auto inst : instances) {
      RTCGeometry ig = rtcNewGeometry(context->device,RTC_GEOMETRY_TYPE_INSTANCE);
      Group *group = (Group *)inst.group;
      rtcSetGeometryInstancedScene(ig, group->handle);
      // rtcSetGeometryTransform(ig, 0, RTC_FORMAT_FLOAT4X3_ROW_MAJOR, &inst.xfm);
      rtcSetGeometryTransform(ig, 0, RTC_FORMAT_FLOAT4X3_COLUMN_MAJOR, &inst.xfm);
      rtcCommitGeometry(ig);

      rtcAttachGeometryByID(this->handle,ig,inst.userID);
      rtcEnableGeometry(ig);
      rtcReleaseGeometry(ig);
    }
    rtcCommitScene(this->handle);
  }
  
  void Model::trace(Ray *rays,
                    Hit *hits,
                    int  numRays,
                    int *activeIDs,
                    int  numActive,
                    OPTraceFlags flags)
  {
    RTCIntersectContext context;
    for (int it=0;it<(activeIDs?numActive:numRays);it++) {
      int rayID = activeIDs?activeIDs[it]:it;
      if (rayID < 0) continue;

      rtcInitIntersectContext(&context);
      Ray &ray = rays[rayID];
      Hit &hit = hits[rayID];
      
      RTCRayHit er;
      er.ray.org_x = ray.origin.x;    // x coordinate of ray origin
      er.ray.org_y = ray.origin.y;    // y coordinate of ray origin
      er.ray.org_z = ray.origin.z;    // z coordinate of ray origin
      er.ray.tnear = ray.tMin;        // start of ray segment
      
      er.ray.dir_x = ray.direction.x; // x coordinate of ray direction
      er.ray.dir_y = ray.direction.y; // y coordinate of ray direction
      er.ray.dir_z = ray.direction.z; // z coordinate of ray direction
      er.ray.time  = 0.f;             // time of this ray for motion blur

      er.ray.tfar  = ray.tMax;
      er.ray.mask  = -1;  // ray mask
      er.ray.id    =  0;  // ray ID
      er.ray.flags =  0;  // ray flags

      er.hit.primID    = -1;
      er.hit.instID[0] = -1;
      er.hit.geomID    = RTC_INVALID_GEOMETRY_ID;

      rtcIntersect1(handle,&context,&er);

      if (er.hit.geomID == RTC_INVALID_GEOMETRY_ID) {
        hit.primID = -1;
        hit.instID = -1;
        hit.geomID = -1;
      } else {
        hit.primID = er.hit.primID;
        hit.instID = er.hit.instID[0];
        hit.geomID = er.hit.geomID;
        hit.u = er.hit.u;
        hit.v = er.hit.v;
        hit.t = er.ray.tfar;
      }
      
    }
  }

  void Model::build() 
  {
    /* nothing to be done, we currently build in model creation :-/ */
  }
    
  /*! the SoA variant that operates on arbitrary strided arrays -
    need to have this in the backend (not just rearrange on
    front-end) because the arrays might not even be accessible on
    the host, so re-arrangement has to be done whatever device the
    user created. Any of the arrays other than org and dir may be
    null (in which case their contents get set to defaluts (for
    rays) or simply not written to (for hit data)); but arrays
    that re not null must have the given number of elemnets in
    them. */
  void Model::trace(size_t numRays,
                    const float *ray_org_x, size_t stride_ray_org_x,
                    const float *ray_org_y, size_t stride_ray_org_y,
                    const float *ray_org_z, size_t stride_ray_org_z,
                    const float *ray_dir_x, size_t stride_ray_dir_x,
                    const float *ray_dir_y, size_t stride_ray_dir_y,
                    const float *ray_dir_z, size_t stride_ray_dir_z,
                    const float *ray_tmin, size_t stride_ray_tmin,
                    const float *ray_tmax, size_t stride_ray_tmax,
                    int *_hit_instID, size_t stride_hit_instID,
                    int *_hit_geomID, size_t stride_hit_geomID,
                    int *_hit_primID, size_t stride_hit_primID,
                    float *_hit_t, size_t stride_hit_t,
                    float *_hit_u, size_t stride_hit_u,
                    float *_hit_v, size_t stride_hit_v,
                    OPTraceFlags flags) 
  {
    std::vector<Ray> rays(numRays);
    std::vector<Hit> hits(numRays);
    
    for (int i=0;i<numRays;i++) {
      rays[i].origin[0] = primer::strided(ray_org_x,stride_ray_org_x,i);
      rays[i].origin[1] = primer::strided(ray_org_y,stride_ray_org_y,i);
      rays[i].origin[2] = primer::strided(ray_org_z,stride_ray_org_z,i);
      rays[i].direction[0] = primer::strided(ray_dir_x,stride_ray_dir_x,i);
      rays[i].direction[1] = primer::strided(ray_dir_y,stride_ray_dir_y,i);
      rays[i].direction[2] = primer::strided(ray_dir_z,stride_ray_dir_z,i);
    }
    if (ray_tmin) 
      for (int i=0;i<numRays;i++) 
        rays[i].tMin = primer::strided(ray_tmin,stride_ray_tmin,i);
    else
      for (int i=0;i<numRays;i++)
        rays[i].tMin = 1e-6f;
    if (ray_tmax) 
      for (int i=0;i<numRays;i++) 
        rays[i].tMax = primer::strided(ray_tmax,stride_ray_tmax,i);
    else
      for (int i=0;i<numRays;i++)
        rays[i].tMax = std::numeric_limits<float>::infinity();
        
    this->trace(rays.data(),hits.data(),
                numRays,(int*)nullptr,
                0,flags);
      
    StridedArray<int> hit_primID(_hit_primID,stride_hit_primID);
    StridedArray<int> hit_geomID(_hit_geomID,stride_hit_geomID);
    StridedArray<int> hit_instID(_hit_instID,stride_hit_instID);
    StridedArray<float> hit_t(_hit_t,stride_hit_t);
    StridedArray<float> hit_u(_hit_u,stride_hit_u);
    StridedArray<float> hit_v(_hit_v,stride_hit_v);

    if (hit_primID) 
      for (int i=0;i<numRays;i++)
        hit_primID[i] = hits[i].primID;
    if (hit_geomID) 
      for (int i=0;i<numRays;i++)
        hit_geomID[i] = hits[i].geomID;
    if (hit_instID) 
      for (int i=0;i<numRays;i++)
        hit_instID[i] = hits[i].instID;
    if (hit_t) 
      for (int i=0;i<numRays;i++)
        hit_t[i] = hits[i].t;
    if (hit_u) 
      for (int i=0;i<numRays;i++)
        hit_u[i] = hits[i].u;
    if (hit_v) 
      for (int i=0;i<numRays;i++)
        hit_v[i] = hits[i].v;
  }
  
  /*! finds the (up to) N closest hits along a ray. For each ray the
    hits array stores N "slots" for hits; teh trace will fill in
    the N closest hits. If there are less than N hits along a ray,
    the next hit slot after the last found one will have its
    primID set to -1 */
  void Model::traceN(Ray *rays,
                     Hit *hits,
                     int  numRays,
                     int  numHitsPerRay,
                     OPTraceFlags flags) 
  {
    EP_NOTIMPLEMENTED;    
  }
    
} // ::ep
