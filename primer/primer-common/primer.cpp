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

/*! \file primer/common/primer.cpp Implements the 'primer' API */

// the public API that we're implementing in this file ...
#include <primer/primer.h>
// for our implementation
#include "primer-common/Context.h"
#include "primer-common/Model.h"

#define WARN(a) std::cout << OWL_TERMINAL_RED << a << OWL_TERMINAL_DEFAULT << std::endl;

#undef OP_API
#define OP_API extern "C"

namespace primer {

  /*! enables logging/debugging messages if true, suppresses them if
    false */
#ifdef NDEBUG
  bool Context::logging = false;
#else
  bool Context::logging = true;
#endif
  
#if PRIMER_HAVE_EMBREE
  /* the "real" Context *Context::createEmbreeContext() will be
     defined in embree/Context.cpp */
#else
  Context *Context::createEmbreeContext()
  { 
    throw std::runtime_error("user explicitly requested CPU context, but this library was compiled without embree support");
  }
#endif

  Geom *castCheck(OPGeom _geom)
  {
    assert(_geom);
    return (Geom *)_geom;
  }

  Group *castCheck(OPGroup _group)
  {
    assert(_group);
    return (Group *)_group;
  }

  Model *castCheck(OPModel _model)
  {
    assert(_model);
    return (Model *)_model;
  }

  Context *castCheck(OPContext _context)
  {
    assert(_context);
    return (Context *)_context;
  }
  
}

#undef OP_API
#define OP_API extern "C"

using namespace primer;

/*! trace all the rays in the given array of input rays, and write
 *  results into the given array of output rays. Unlike
 *  opTraceAsync, this call is synchronous, so will block the
 *  calling thread until all rays are traced and all hits have been
 *  written. Input and output arrays have to remain valid and
 *  unmodified while teh call is active, but can be modified or
 *  released the moment this call returns */
OP_API void opTrace(/*! the model to trace into */
                    OPModel model,
                    /*! number of rays to trace - both arrays must have as
                     *  many entries */
                    size_t numRays,
                    /*! array of rays; must be (at least) as many as
                     *  'numRays' */
                    OPRay *arrayOfRays,
                    /*! array of where to write the results; as many as
                     *  'numRays' */
                    OPHit *arrayOfHits,
                    OPTraceFlags flags) 
{
  assert(model);
  castCheck(model)->trace((primer::Ray *)arrayOfRays,
                          (primer::Hit *)arrayOfHits,
                          (int)numRays,
                          0,0,
                          flags);
}


/*! trace all the rays in the given array of input rays, and write
 *  results into the given array of output rays. For each ray,
 *  return (up to) N hits, with the this for ray nuber 0 stored in
 *  hits[0..N-1], those for ray1 stored in hits[N..2N-1], etc. For
 *  each ray, the hits are stored in ascending order, and are the N
 *  respectively closest hits alon gthis ray. If a ray has less than
 *  N hits then the other hits in that ray's section of the hit
 *  array will be set to primID=-1 */
OP_API void opTraceNHits(/*! the model to trace into */
                         OPModel model,
                         /*! number of rays to trace - both arrays must have as
                          *  many entries */
                         size_t numRays,
                         /*! array of rays; must be (at least) as many as
                          *  'numRays' */
                         OPRay *arrayOfRays,
                         /*! the number of maximum "N" hits to be
                           stored per ray */
                         int maxHitsPerRay,
                         /*! array of where to write the results; as many as
                          *  'numRays' */
                         OPHit *arrayOfHits,
                         /*! trace flags that can fine-tune how the trace
                          *  executes. */
                         OPTraceFlags traceFlags)
{
  assert(model);
  castCheck(model)->traceN((primer::Ray *)arrayOfRays,
                           (primer::Hit *)arrayOfHits,
                           (int)numRays,
                           (int)maxHitsPerRay,
                           traceFlags);
}


OP_API void opTraceIndexed(/*! the model to trace into */
                           OPModel model,
                           /*! number of ray indices - this is the
                             number of rays to be traced. */
                           size_t   numActiveRayIndices,
                           /*! optionally, a array of ray IDs that specify which
                            *  of the rays in 'arrayOfRays' to trace. if null, we
                            *  impliictly assume the list of active rays/hits is
                            *  [0,1,2,...numRays) */
                           int32_t *arrayOfActiveRayIndices,

                           /*! size of ray and hit arrays, required
                             for being able to offload those
                             arrays if required */
                           size_t sizeOfRayAndHitArray,

                           /*! array of rays; must be (at least) as many as
                            *  'numRays' */
                           OPRay *arrayOfRays,
                           /*! array of where to write the results; as many as
                            *  'numRays' */
                           OPHit *arrayOfHits,
                           /*! trace flags that can fine-tune how the trace
                            *  executes. */
                           OPTraceFlags flags) 
{
  assert(model);
  castCheck(model)->trace((primer::Ray *)arrayOfRays,
                          (primer::Hit *)arrayOfHits,
                          (int)sizeOfRayAndHitArray,
                          arrayOfActiveRayIndices,
                          (int)numActiveRayIndices,
                          flags);
}

OP_API OPContext opContextCreate(OPContextType contextType,
                                 /*! which GPU to use, '0' being the first
                                  *  GPU, '1' the second, etc. '-1' means
                                  *  'use host CPU only' */
                                 int32_t gpuToUse)
{
  if (contextType == OP_CONTEXT_HOST_FORCE_CPU) {
    return (OPContext)primer::Context::createEmbreeContext();
  }

  try {
    return (OPContext)primer::Context::createOffloadContext(gpuToUse);
  } catch (...) {
    if (contextType == OP_CONTEXT_REQUIRE_GPU) {
      WARN("could not create required GPU context");
      //throw std::runtime_error("could not create required GPU context");
      return (OPContext)0;
    }
    WARN("could not create a GPU context; falling back to CPU context");
    return (OPContext)primer::Context::createEmbreeContext();
  }
}

OP_API
OPGeom opMeshCreate(OPContext _context,
                    uint32_t userGeomIDtoUseInHits,
                    /* vertex array */
                    const float *vertices,
                    size_t numVertices,
                    size_t sizeOfVertexInBytes,
                    /* index array */
                    const int   *indices,
                    size_t numIndices,
                    size_t sizeOfIndexStructInBytes)
{
  Context *context = castCheck(_context);
  Geom *mesh = context->createTriangles
    (/* ID */userGeomIDtoUseInHits,
     /* vertex array */
     (const vec3f*)vertices,numVertices,sizeOfVertexInBytes,
     /* index array */
     (const vec3i*)indices,numIndices,sizeOfIndexStructInBytes);
  assert(mesh);
  return (OPGeom)mesh;
}
  
/*! creates a geometry that contains one or more spheres; each
  sphere specified by origin and radius */
OP_API OPGeom opSpheresFromArrays(OPContext _context,
                                  uint32_t userGeomIDtoUseInHits,
                                  size_t   numSpheres,
                                  /*! sphere centers, each assuming
                                    x,y,z being stored
                                    consecutively */
                                  const float *centers_x,
                                  const float *centers_y,
                                  const float *centers_z,
                                  size_t centersStrideInBytes,
                                  /* array of sphere radii */
                                  const float *radii,
                                  size_t radiiStrideInBytes)
{
  Context *context = castCheck(_context);
  Geom *spheres = context->createSpheres
    (/* ID */userGeomIDtoUseInHits,
     numSpheres,
     /* centers */
     centers_x, 
     centers_y, 
     centers_z, 
     centersStrideInBytes,
     /* radii */
     radii,
     radiiStrideInBytes);
  assert(spheres);
  return (OPGeom)spheres;
}

/*! creates a geometry that contains one or more spheres; each
  sphere is specified by a float4 whose xyz coordinate specify the
  position, and w coordinate specifies radius */
OP_API OPGeom opSpheres4f(OPContext context,
                          uint32_t userGeomIDtoUseInHits,
#ifdef __CUDACC__
                          const float4 *spheres,
#else
                          const primer_float4 *spheres,
#endif
                          size_t   numSpheres)
{
  return opSpheresFromArrays(context,
                             userGeomIDtoUseInHits,
                             numSpheres,
                             &spheres->x,
                             &spheres->y,
                             &spheres->z,
                             sizeof(primer_float4),
                             &spheres->w,
                             sizeof(primer_float4));
}
  


/*! creates a mesh from an array of "fat" triangles (ie, one struct
  with three vertices each, w/o indices). Can be used for both SoA
  and AoS layouts */
OP_API OPGeom opMeshFromTriangles(OPContext _context,
                                  uint32_t userMeshIDtoUseInHits,
                                  size_t numTriangles,
                                  /* vertex 0 */
                                  const float *v0x,
                                  const float *v0y,
                                  const float *v0z,
                                  /* vertex 1 */
                                  const float *v1x,
                                  const float *v1y,
                                  const float *v1z,
                                  /* vertex 2 */
                                  const float *v2x,
                                  const float *v2y,
                                  const float *v2z,
                                  size_t stride)
{
  Context *context = castCheck(_context);
  Geom *mesh = context->createTriangles(userMeshIDtoUseInHits,
                                        numTriangles,
                                        v0x,v0y,v0z,
                                        v1x,v1y,v1z,
                                        v2x,v2y,v2z,
                                        stride);
  assert(mesh);
  return (OPGeom)mesh;
}

/*! creates a "trivial" model from a single triangle mesh that is
 *  specified through a single array of "fat" triangles; ie, where
 *  each tringle is its own struct that contains all three vertices
 *  and whatever else it might need (ie, NO separate vertex or index
 *  arrays). The model will automatically be built w/ opModelBuild() */
OP_API OPModel opModelFromTriangles(OPContext _context,
                                    size_t numTriangles,
                                    const float *triangle0vertex0,
                                    const float *triangle0vertex1,
                                    const float *triangle0vertex2,
                                    size_t strideInBytes)
{
  OPGeom mesh = opMeshFromTriangles(_context,
                                    /* ID */0,
                                    numTriangles,
                                    /* vertex 0 */
                                    triangle0vertex0+0,
                                    triangle0vertex0+1,
                                    triangle0vertex0+2,
                                    /* vertex 1 */
                                    triangle0vertex1+0,
                                    triangle0vertex1+1,
                                    triangle0vertex1+2,
                                    /* vertex 2 */
                                    triangle0vertex2+0,
                                    triangle0vertex2+1,
                                    triangle0vertex2+2,
                                    /* stride */
                                    strideInBytes);//sizeof(vec3f));
  return opModelFromGeoms(_context,&mesh,1);
}
  
OP_API OPModel opModelFromGeoms(OPContext _context,
                                OPGeom *geoms,
                                size_t numGeoms)
{
  OPInstance inst;
  inst.userID = 0;
  (affine3f&)inst.xfm = affine3f();
  inst.group = opGroupCreate(_context,geoms,1);
  return opModelCreate(_context,&inst,1);
}

OP_API OPModel opModelFromGeom(OPContext context,
                               OPGeom mesh)
{
  return opModelFromGeoms(context,&mesh,1);
}


/*! SOA/Fortran verion of opTrace, but allowing for arbitrary
 *  strides: traces all the rays in the given array of input rays,
 *  and write results into the given array of output rays. Unlike
 *  opTraceAsync, this call is synchronous, so will block the
 *  calling thread until all rays are traced and all hits have been
 *  written.
 *
 *  Input and output arrays have to remain valid and unmodified
 *  while teh call is active, but can be modified or released the
 *  moment this call returns. 
 *
 *  Note that all IDs returned in the hit point will use _C_
 *  indexing convention (starting with 0), not fortran indexing */
OP_API
void opTraceStridedSOA(/*! the model to trace into */
                       OPModel model,
                       /*! number of rays to trace - both arrays must have as
                        *  many entries */
                       size_t numRays,
                       /*! @{ arrays that describe the rays to trace: */
                       const float *ray_org_xs,
                       size_t stride_ray_org_xs,
                       const float *ray_org_ys, 
                       size_t stride_ray_org_ys,
                       const float *ray_org_zs, 
                       size_t stride_ray_org_zs,
                       const float *ray_dir_xs,
                       size_t stride_ray_dir_xs,
                       const float *ray_dir_ys,
                       size_t stride_ray_dir_ys,
                       const float *ray_dir_zs,
                       size_t stride_ray_dir_zs,
                       const float *ray_tmins,
                       size_t stride_ray_tmins,
                       const float *ray_tmaxs,
                       size_t stride_ray_tmaxs,
                       /*! @} */
                       /*! @{ arrays that will hold all the elemnets of the
                        *  computed hits; those that are null will get
                        *  ignored */
                       int32_t *hit_instIDs,
                       size_t stride_hit_instIDs,
                       int32_t *hit_geomIDs,
                       size_t stride_hit_geomIDs,
                       int32_t *hit_primIDs,
                       size_t stride_hit_primIDs,
                       float *hit_ts,
                       size_t stride_hit_ts,
                       float *hit_us,
                       size_t stride_hit_us,
                       float *hit_vs,
                       size_t stride_hit_vs,
                       /*! trace flags that can fine-tune how the trace
                        *  executes. */
                       OPTraceFlags flags)
{
  castCheck(model)->trace(numRays,
                          ray_org_xs, stride_ray_org_xs,
                          ray_org_ys, stride_ray_org_ys,
                          ray_org_zs, stride_ray_org_zs,
                          ray_dir_xs, stride_ray_dir_xs,
                          ray_dir_ys, stride_ray_dir_ys,
                          ray_dir_zs, stride_ray_dir_zs,
                          ray_tmins, stride_ray_tmins,
                          ray_tmaxs, stride_ray_tmaxs,
                          hit_instIDs, stride_hit_instIDs,
                          hit_geomIDs, stride_hit_geomIDs,
                          hit_primIDs, stride_hit_primIDs,
                          hit_ts, stride_hit_ts,
                          hit_us, stride_hit_us,
                          hit_vs, stride_hit_vs,
                          flags);
}


OP_API
OPGroup opGroupCreate(OPContext _context,
                      OPGeom *_geoms,
                      int numGeoms)
{
  Context *context = castCheck(_context);
  std::vector<OPGeom> geoms;
  for (int i=0;i<numGeoms;i++) {
    assert(_geoms[i]);
    geoms.push_back(_geoms[i]);
  }
  Group *group = context->createGroup(geoms);
  assert(group);
  return (OPGroup)group;
}

OP_API
OPModel opModelCreate(OPContext    _context,
                      OPInstance  *instances,
                      int          numInstances)
{
  Context *context = castCheck(_context);
  std::vector<OPInstance> _instances(numInstances);
  std::copy(instances,instances+numInstances,_instances.data());
  Model *model = context->createModel(_instances);
  assert(model);
  return (OPModel)model;
}
