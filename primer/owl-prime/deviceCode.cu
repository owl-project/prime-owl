// ======================================================================== //
// Copyright 2019-2022 Ingo Wald                                            //
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

#include "owl-prime/deviceCode.h"
#include "owl-prime/Context.h"
#include "owl-prime/Triangles.h"

using namespace owl::common;
using op::Hit;

extern "C" __constant__ op::Context::LPData optixLaunchParams;

struct MultiHitPRD
{
  /*! the list of hits to be used for this trace */
  Hit  *hits;
  /*! number of hits already stored in the list */
  int   numHitsFound;
  /*! the (current) maximum t after which we *know* no more closer
      hits can be found. this is either the last entry in the list (if
      it is full), or ray.tmax (if not) */
  float cutOff;
};

  
OPTIX_BOUNDS_PROGRAM(Spheres)(const void  *geomData,
                              box3f       &primBounds,
                              const int    primID)
{
  const auto &self = *(const op::Spheres::SBTData*)geomData;
  float4 sphere = self.spheres[primID];
  primBounds.lower.x = sphere.x - sphere.w;
  primBounds.lower.y = sphere.y - sphere.w;
  primBounds.lower.z = sphere.z - sphere.w;
  primBounds.upper.x = sphere.x + sphere.w;
  primBounds.upper.y = sphere.y + sphere.w;
  primBounds.upper.z = sphere.z + sphere.w;
}

OPTIX_INTERSECT_PROGRAM(Spheres)()
{
  const auto &self
    = owl::getProgramData<op::Spheres::SBTData>();
  const int primID = optixGetPrimitiveIndex();
  
  const vec3f org  = optixGetWorldRayOrigin();
  const vec3f dir  = optixGetWorldRayDirection();
  float hit_t      = optixGetRayTmax();
  const float tmin = optixGetRayTmin();

  const float4 sphere = self.spheres[primID];
  const vec3f center = {sphere.x,sphere.y,sphere.z};
  const float radius = sphere.w;
  const vec3f oc = org - center;
         
  const float a = dot(dir,dir);
  const float b = dot(oc, dir);
  const float c = dot(oc, oc) - radius * radius;
  const float discriminant = b * b - a * c;
  
  if (discriminant < 0.f) return;

  {
    float temp = (-b - sqrtf(discriminant)) / a;
    if (temp < hit_t && temp > tmin) 
      hit_t = temp;
  }
      
  {
    float temp = (-b + sqrtf(discriminant)) / a;
    if (temp < hit_t && temp > tmin) 
      hit_t = temp;
  }
  if (hit_t < optixGetRayTmax()) {
    optixReportIntersection(hit_t, 0);
  }
}

/*! closest hit program: fills in a hit structure for the current ray */
OPTIX_CLOSEST_HIT_PROGRAM(Spheres)()
{
  const auto &self = owl::getProgramData<op::Spheres::SBTData>();
  Hit &hit = owl::getPRD<Hit>();
  hit.t = optixGetRayTmax();
  hit.primID = optixGetPrimitiveIndex();
  hit.geomID = self.userID;
  hit.instID = optixGetInstanceIndex();
  hit.u      = 0.f;
  hit.v      = 0.f;
}

/*! closest hit program: fills in a hit structure for the current ray */
OPTIX_CLOSEST_HIT_PROGRAM(Triangles)()
{
  const auto &self  = owl::getProgramData<op::Triangles::SBTData>();
  Hit &hit = owl::getPRD<Hit>();
  hit.t = optixGetRayTmax();
  hit.primID = optixGetPrimitiveIndex();
  hit.geomID = self.userID;
  hit.instID = optixGetInstanceIndex();
  hit.u      = optixGetTriangleBarycentrics().x;
  hit.v      = optixGetTriangleBarycentrics().y;
}

/*! closest hit program: used for the multi-hit kernel with up to N
    hits per ray. this should never get called if not in a multi-hit
    trace (we use the disable_ah flag to turn it off for other
    queries), so inside this code we can safely assume that there's a
    multi-hit PRD */
OPTIX_ANY_HIT_PROGRAM(MultiHit)()
{
  const auto &self = owl::getProgramData<op::Triangles::SBTData>();
  // if (lp.numHitsPerRay == 0) {
  //   // NOT multi-hit, but TRACE_CONTINUE
  //   Hit &hit = owl::getPRD<Hit>();
  // } else {
    MultiHitPRD &prd = owl::getPRD<MultiHitPRD>();
    Hit *hitList = prd.hits;

    Hit thisHit;
    thisHit.t      = optixGetRayTmax();
    thisHit.primID = optixGetPrimitiveIndex();
    thisHit.geomID = self.userID;
    thisHit.instID = optixGetInstanceIndex();
    thisHit.u      = optixGetTriangleBarycentrics().x;
    thisHit.v      = optixGetTriangleBarycentrics().y;
    if (thisHit.t >= prd.cutOff) {
      // ACCEPT hit to optix (this will move new ray.tmax to thisHit.t),
      // but don't save in list.
      return;
    } else {
      // we KNOW that the list is not yet full, OR that this hit is
      // closer than the furthest so-far found one.
      auto &lp = optixLaunchParams;
      int insertPos
        = (prd.numHitsFound < lp.numHitsPerRay)
        ? /* list not yet full: append*/prd.numHitsFound
        : /* list full, overwrite furthest*/(lp.numHitsPerRay-1);
      while (insertPos > 0 && hitList[insertPos-1].t > thisHit.t) {
        hitList[insertPos] = hitList[insertPos-1];
        --insertPos;
      }
      hitList[insertPos] = thisHit;
      prd.numHitsFound++;
      if (prd.numHitsFound >= lp.numHitsPerRay) {
        prd.cutOff = min(prd.cutOff,hitList[lp.numHitsPerRay-1].t);
        prd.numHitsFound = lp.numHitsPerRay;
      }
      if (thisHit.t < prd.cutOff)
        optixIgnoreIntersection();
    }
  // }
}

OPTIX_RAYGEN_PROGRAM(traceRays)()
{
  Hit hit;
  hit.clear();
  int rayID
    = owl::getLaunchIndex().x
    + owl::getLaunchDims().x
    * owl::getLaunchIndex().y;

  auto &lp = optixLaunchParams;
  if (rayID >= lp.numRays) return;

  if (lp.activeIDs) {
    rayID = lp.activeIDs[rayID];
    if (rayID < 0) return;
  }
  
  owl::Ray ray;
  if (lp.isSOA) {
    ray = owl::Ray(vec3f(lp.soa.ray.org_x[rayID],
                         lp.soa.ray.org_y[rayID],
                         lp.soa.ray.org_z[rayID]),
                   vec3f(lp.soa.ray.dir_x[rayID],
                         lp.soa.ray.dir_y[rayID],
                         lp.soa.ray.dir_z[rayID]),
                   lp.soa.ray.t_min[rayID],
                   lp.soa.ray.t_max[rayID]);
  } else {
    ray = owl::Ray((const vec3f&)lp.rays[rayID].origin,
                   (const vec3f&)lp.rays[rayID].direction,
                   lp.rays[rayID].tMin,
                   lp.rays[rayID].tMax);
  }
  if (ray.tmin < ray.tmax) {
    if (lp.numHitsPerRay == 0) {
      if (lp.flags & OP_TRACE_CONTINUE) {
        uint32_t rayFlags = 0;
        owl::traceRay(lp.model,ray,hit,rayFlags);
      } else {
        // this is the "defualt" path where each ray finds _the_ closest hit
        uint32_t rayFlags = OPTIX_RAY_FLAG_DISABLE_ANYHIT;
        if (lp.flags & OP_TRACE_FLAGS_TERMINATE_ON_FIRST_HIT)
          rayFlags |= OPTIX_RAY_FLAG_TERMINATE_ON_FIRST_HIT;
        owl::traceRay(lp.model,ray,hit,rayFlags);
      }
    } else {
      // in this path each ray can find and store up to N different
      // hits; we do this via a anyhit program (and disable closest hit)
      uint32_t rayFlags = OPTIX_RAY_FLAG_DISABLE_CLOSESTHIT;
      MultiHitPRD multiHitPrd;
      multiHitPrd.numHitsFound = 0;
      multiHitPrd.cutOff       = ray.tmax;
      multiHitPrd.hits         = lp.hits+lp.numHitsPerRay*rayID;
      owl::traceRay(lp.model,ray,multiHitPrd,rayFlags);
      for (int i=multiHitPrd.numHitsFound; i<lp.numHitsPerRay; i++)
        multiHitPrd.hits[i].primID = -1;
    }
  }
  if (lp.isSOA) {
    if (lp.soa.hit.primID)
      lp.soa.hit.primID[rayID] = hit.primID;
    if (lp.soa.hit.geomID)
      lp.soa.hit.geomID[rayID] = hit.geomID;
    if (lp.soa.hit.instID)
      lp.soa.hit.instID[rayID] = hit.instID;
  } else {
    lp.hits[rayID] = hit;
  }
}


