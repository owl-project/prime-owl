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

#include "owl-prime/Model.h"
#include "owl-prime/Context.h"
#include "owl-prime/Group.h"
#include "primer-common/StridedArray.h"

namespace op {
  using primer::StridedArray;
  
  Model::Model(Context *context,
               const std::vector<OPGroup>  &groups,
               const std::vector<affine3f> &xfms,
               const std::vector<int>      &userIDs)
    : // primer::Model(context),
      context(context)
  {
    // std::vector<affine3f> owlXfms;
    std::vector<OWLGroup> owlGroups;
    for (auto &g : groups) {
      op::Group *group = (op::Group *)g;
      if (group->trianglesGroup) {
        owlGroups.push_back(group->trianglesGroup);
        // owlXfms.push_back((const affine3f&)inst.xfm);
        // owlInstIDs.push_back(inst.userID);
      }
      if (group->userGroup) {
        owlGroups.push_back(group->userGroup);
        // owlXfms.push_back((const affine3f&)inst.xfm);
        // owlInstIDs.push_back(inst.userID);
      }
    }
    handle
      = owlInstanceGroupCreate(context->owl,
                               owlGroups.size(),
                               owlGroups.data(),
                               (const uint32_t*)userIDs.data(),
                               (const float *)xfms.data());
    owlGroupBuildAccel(handle);
  }

  void Model::trace(Ray *rays,
                    Hit *hits,
                    int  numRays,
                    int *activeIDs,
                    int  numActive,
                    OPTraceFlags flags)
  {
    static std::mutex mutex;
    std::lock_guard<std::mutex> lock(mutex);

    context->checkSBT();
    
    auto &lp = context->launchParams;
    const int numHits = numRays;
    
    Ray *d_rays = 0;
    Hit *d_hits = 0;
    int *d_activeIDs = 0;
    
    if (isDeviceAccessible(rays)) {
      d_rays = rays;
    } else {
      CUDA_CALL(Malloc(&d_rays,numRays*sizeof(Ray)));
      CUDA_CALL(Memcpy(d_rays,rays,numRays*sizeof(Ray),cudaMemcpyDefault));
    }
    if (isDeviceAccessible(hits)) {
      d_hits = hits;
    } else {
      CUDA_CALL(Malloc(&d_hits,numRays*sizeof(Hit)));
      CUDA_CALL(Memcpy(d_hits,hits,numRays*sizeof(Hit),cudaMemcpyDefault));
    }
    if (isDeviceAccessible(activeIDs)) {
      d_activeIDs = activeIDs;
    } else {
      CUDA_CALL(Malloc(&d_activeIDs,numActive*sizeof(int)));
      CUDA_CALL(Memcpy(d_activeIDs,activeIDs,numActive*sizeof(int),cudaMemcpyDefault));
    }
    owlParamsSetPointer(lp,"rays",d_rays);
    owlParamsSetPointer(lp,"hits",d_hits);
    owlParamsSetPointer(lp,"activeIDs",d_activeIDs);
    owlParamsSet1i(lp,"numRays",numActive?numActive:numRays);
    owlParamsSet1i(lp,"numHitsPerRay",0);
    owlParamsSet1i(lp,"isSOA",0);
    owlParamsSet1ul(lp,"flags",flags);
    owlParamsSetGroup(lp,"model",handle);

    owlLaunch2D(context->rayGen,numRays,1,lp);
    
    if (d_hits != hits) {
      CUDA_CALL(Memcpy(hits,d_hits,numRays*sizeof(Hit),cudaMemcpyDefault));
    }
    
    if (d_rays != rays)
      CUDA_CALL(Free(d_rays));
    if (d_hits != hits)
      CUDA_CALL(Free(d_hits));
    if (d_activeIDs != activeIDs)
      CUDA_CALL(Free(d_activeIDs));
    
    CUDA_SYNC_CHECK();
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
    static std::mutex mutex;
    std::lock_guard<std::mutex> lock(mutex);

    context->checkSBT();
    
    auto &lp = context->launchParams;
    const int numHits = numRays;
    
    Ray *d_rays = 0;
    Hit *d_hits = 0;
    int *d_activeIDs = 0;
    
    if (isDeviceAccessible(rays)) {
      d_rays = rays;
    } else {
      CUDA_CALL(Malloc(&d_rays,numRays*sizeof(Ray)));
      CUDA_CALL(Memcpy(d_rays,rays,numRays*sizeof(Ray),cudaMemcpyDefault));
    }
    if (isDeviceAccessible(hits)) {
      d_hits = hits;
    } else {
      CUDA_CALL(Malloc(&d_hits,numRays*numHitsPerRay*sizeof(Hit)));
      CUDA_CALL(Memcpy(d_hits,hits,numRays*numHitsPerRay*sizeof(Hit),cudaMemcpyDefault));
    }
    owlParamsSetPointer(lp,"rays",d_rays);
    owlParamsSetPointer(lp,"hits",d_hits);
    owlParamsSetPointer(lp,"activeIDs",nullptr);
    owlParamsSet1i(lp,"numRays",numRays);
    owlParamsSet1i(lp,"numHitsPerRay",numHitsPerRay);
    owlParamsSet1i(lp,"isSOA",0);
    owlParamsSet1ul(lp,"flags",flags);
    owlParamsSetGroup(lp,"model",handle);

    owlLaunch2D(context->rayGen,numRays,1,lp);
    
    if (d_hits != hits) {
      CUDA_CALL(Memcpy(hits,d_hits,numRays*numHitsPerRay*sizeof(Hit),cudaMemcpyDefault));
    }
    
    if (d_rays != rays)
      CUDA_CALL(Free(d_rays));
    if (d_hits != hits)
      CUDA_CALL(Free(d_hits));
    
    CUDA_SYNC_CHECK();
  }
    

  void Model::build() 
  {
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
    if (!isDeviceAccessible(ray_org_x)) {
      /* arrays are not accessible on the device; need to first
         upload. to do that, first rearrange locally on the host ino
         AoS layout, then call aos trace to do the upload, trace, and
         download */
      Ray *rays = 0;
      Hit *hits = 0;
      cudaMallocManaged(&rays,numRays*sizeof(Ray));
      cudaMallocManaged(&hits,numRays*sizeof(Hit));

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
        
      this->trace(rays,hits,
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

      cudaFree(rays);
      cudaFree(hits);

      return;
    }

    throw std::runtime_error("trace() on SoA data is currently implemented only for host-side data...");
  }
  
} // ::op
