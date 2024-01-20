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

#ifndef PRIMER_PRIMER_H
#define PRIMER_PRIMER_H

#include <sys/types.h>
#include <stdint.h>

# define OP_API extern "C"

#ifdef __cplusplus
# define OP_IF_CPP(a) a
#else
# define OP_IF_CPP(a) /* nothing */
#endif

struct primer_float4 { float x,y,z,w; };

/*! defines data layout for a input ray for those functions that
  operate on arrays of rays in AoS layout */
struct OPRay {
  float origin[3];
  float tMin;
  float direction[3];
  float tMax;
};

/*! sruct describing a ray's hit; when user supplied device-side
    arrays this array has to be 16-byte aligned; for host-side data
    arrays that get automatically uploaded this is not required. */
struct OPHit {
  /*! the 'primitive' within a given mesh. a value of -1 here will
    indicate that nothing was hit */
  int primID;
  
  /*! the (user-provided) geometry ID that is associated with the
    given mesh that the intersected primitmive was in */
  int geomID;
  
  /*! the (user provided) instance ID of the instance that the
    intersected primitive was in */
  int instID;
  
  /*! distance to, and barycentric surface coordinates of the hit
     point */
  float t, u, v;
  
  /*! internal variables to be used by back-end */
  int internal[2];
};

/*! the OPContext is the root context object for all of primer; it
  manages creation of new groups, models, meshes, etc; baesd on the
  *type* of context it may do this in various different ways */
typedef struct _OPContext  *OPContext;

/*! a OPModel is a object that contains one or more geometries
  (optionally using groups and/or instances), and that rays can be
  traced against */
typedef struct _OPModel    *OPModel;

/*! groups can be used to group one or more geometries together for
  later instantiation. All geometries in a group have to be of the
  same type (ie, you can have a group with multiple meshes, and
  another for multiple spheres geometrys; but you can not have one
  group that contains both meshes and sphere geometries. */
typedef struct _OPGroup    *OPGroup;

/*! a "geometry" that contains one of more geometric primitives like
    spheres, triangles, etc. Note geomtries are "virtual" in the sense
    that this geom could be either a "spheres" geometry, a "boxes"
    geometry, depending on how it was created */
typedef struct _OPGeom     *OPGeom;

#define OP_TRACE_FLAGS_DEFAULT 0ull
/*! if specified, will terminate traversal on the first intersection
 *  encountered */
#define OP_TRACE_FLAGS_TERMINATE_ON_FIRST_HIT  (1ull<<0)
#if 0
// not yet supported
#define OP_TRACE_FLAGS_CULL_BACK_FACING_TRIANGLES  (1ull<<1)
#define OP_TRACE_FLAGS_CULL_FRONT_FACING_TRIANGLES (1ull<<2)
#endif

/*! if specified, traversal will start at the latest hit-point stored
    in the hits array (i.e., the resut from the last trace call), and
    will find the respectively next hit along that ray. */
#define OP_TRACE_CONTINUE  (1ull<<3)

/*! specifies which kind of context shuld be created (ie, whether the
 *  context is supposed to be created; in paritulcar, where the model
 *  data and ray/hit arrays will _natively_ be sitting */
typedef enum {
              /*! tells primer that all data arrays (meshes, rays, and
                hits) will all be accessible on the host, and that
                it is free to choose whatever backend (CPU or GPU)
                it wants to create for that. Typically this will
                mean that it will try to create a OWL/OptiX
                "offload" backend, but that if for some reaon it
                cannot do this it can still fall back to a host-side
                trace. This is the expected default behavior. (Note
                that this would still be a valid target even for a
                CUDA application that uses CUDA managed memory -
                though if the app can do CUDA, anyway, then it might
                just as well requirest a GPU backend). */
              OP_CONTEXT_HOST_DEFAULT,
              OP_CONTEXT_DEFAULT=OP_CONTEXT_HOST_DEFAULT,
              
              /*! tells primer that all data arrays (meshes, rays, and
                hits reside in GPU memory, and thus *require* a GPU
                backend. If none such can be created for some reason
                (typically because no suitable GPU is present),
                context creation will fail and return an invalid
                context. */
              OP_CONTEXT_REQUIRE_GPU,
              OP_CONTEXT_GPU=OP_CONTEXT_REQUIRE_GPU,
              
              /*! explicitly requests a CPU back-end where all tracing
                happens on the host CPU. If the input data (meshes,
                rays, hits) are only accessible on the GPU, then
                this will lead to runtime errors (crashes) */
              OP_CONTEXT_HOST_FORCE_CPU,
              
              /*! explicitly requests a 'offload' device where all
                data resides on the host, but all trace will happen
                on the GPU. If no such backend can be created
                (typically because the system does not have a
                suitable GPU) then context creation will fail and
                return a null context. */
              OP_CONTEXT_HOST_FORCE_OFFLOAD
} OPContextType;

typedef uint64_t OPTraceFlags;


#ifndef __cplusplus
extern "C" {
#endif

  OP_API OPContext opContextCreate(OPContextType contextType,
                                   /*! which GPU to use, '0' being the first
                                    *  GPU, '1' the second, etc. '-1' means
                                    *  'use host CPU only' */
                                   int32_t gpuToUse);
  
  OP_API void  opContextDestroy(OPContext context);


  // ==================================================================
  // (Trangle-)Meshes: triangle meshes are made up of one or more
  // triangles, specified through a set of vertices and a set of
  // vertex indices. For convenience we also allow for creating meshes
  // from individual triangles, and/or various different memory
  // layouts for specifying these vertices and/or indices.
  // ==================================================================
  OP_API OPGeom opMeshCreate(OPContext context,
                             uint32_t userMeshIDtoUseInHits,
                             /* vertex array */
                             const float *vertices,
                             size_t numVertices,
                             size_t sizeOfVertexInBytes,
                             /* index array */
                             const int   *indices,
                             size_t numTriangles,
                             size_t sizeOfIndexStructInBytes);
  
  OP_API OPGeom opMeshCreateSOA(OPContext context,
                                uint32_t userMeshIDtoUseInHits,
                                /* index array */
                                /*! to help bindings for fortran: if set to
                                 *  1, we assume that a vertex of index of
                                 *  'i' would actually refer to (C-style)
                                 *  position 'i-1' in the array; should be 0
                                 *  for regular C-style indexing where the
                                 *  first array element is A[0]. */
                                int32_t indicesAssumeFortranIndexing,
                                const int   *indices_x,
                                const int   *indices_y,
                                const int   *indices_z,
                                /* vertex array */
                                size_t numVertices,
                                const float *vertices_x,
                                const float *vertices_y,
                                const float *vertices_z);

  /*! creates a mesh from an array of "fat" triangles (ie, one struct
      with three vertices each, w/o indices). Can be used for both SoA
      and AoS layouts */
  OP_API OPGeom opMeshFromTriangles(OPContext context,
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
                                    size_t stride);

  // ==================================================================
  // Spheres: a spheres geometry is made up of one or more spheres,
  // each specified by a center and a radius.
  // ==================================================================

  /*! creates a geometry that contains one or more spheres; each
      sphere specified by origin and radius */
  OP_API OPGeom opSpheresFromArrays(OPContext context,
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
                                    size_t radiiStrideInBytes);

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
                            size_t   numSpheres);
  

  // ==================================================================
  // model creation. Models can contain anything from a single
  // geometry to multiple geometries in multiple groups, and possibly
  // intances thereof. Internally all models will be organized with
  // instances and groups; but to allow users that do not use
  // instances or groups (or even want to know what those are!) there
  // are also convenience functions to create a model from a single
  // geometry, or from a "plain" list of geometries.
  // ==================================================================
  
  /*! a affine transform where the first three vectors describe the
      linear transform, and 'p' is the translation. This is
      binary-compatible to both embree's AffineSpace<3> as well as
      owl's owl::common::affine3f. */
  struct OPInstance {
    OPGroup     group OP_IF_CPP(=0);
    /*! the user-provided instance ID that will be returned in
        Hit::instID */
    uint32_t    userID OP_IF_CPP(=0);
    float       xfm[3][4];
  };
  
  
  OP_API
  OPGroup opGroupCreate(OPContext _context,
                        OPGeom *geoms,
                        int numGeoms);

  OP_API
  OPModel opModelCreate(OPContext   context,
                        OPInstance *instances,
                        int         numInstances);






  
  
  // ==================================================================
  // some helper functions for creating complete scenes from a single
  // mesh, or a single set of "fat" triangles, etc.
  // ==================================================================
  
  OP_API OPModel opModelFromGeoms(OPContext context,
                                  OPGeom *geoms,
                                  size_t  numGeoms);
  
  /*! creates a "trivial" model from a single triangle mesh (ie,
   *  creates a model that contains a single object with a single
   *  mesh, and a single instnace of thta object, with a 'dummy'
   *  identify transform) The model will automatically be built w/
   *  opModelBuild() */
  OP_API OPModel opModelFromMeshSOA(OPContext context,
                                    size_t numVertices,
                                    float *vertices_x,
                                    float *vertices_y,
                                    float *vertices_z,
                                    size_t numIndices,
                                    int32_t *indices_x,
                                    int32_t *indices_y,
                                    int32_t *indices_z,
                                    int32_t theseAreFortranIndices);
  
  
  /*! creates a "trivial" model from a single triangle mesh that is
   *  specified through a single array of "fat" triangles; ie, where
   *  each tringle is its own struct that contains all three vertices
   *  and whatever else it might need (ie, NO separate vertex or index
   *  arrays). The model will automatically be built w/ opModelBuild() */
  OP_API OPModel opModelFromTriangles(OPContext context,
                                      size_t numTriangles,
                                      const float *triangle0vertex0,
                                      const float *triangle0vertex1,
                                      const float *triangle0vertex2,
                                      size_t strideInBytes);
  
  // ==================================================================
  // how to actually trace (arrays) of rays and get (arrays) of hits
  // ==================================================================
  
  
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
                      /*! trace flags that can fine-tune how the trace
                       *  executes. */
                      OPTraceFlags traceFlags OP_IF_CPP(= OP_TRACE_FLAGS_DEFAULT));

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
                           OPHit *arrayOfHitsWithNHitsPerRay,
                           /*! trace flags that can fine-tune how the trace
                            *  executes. */
                           OPTraceFlags traceFlags OP_IF_CPP(= OP_TRACE_FLAGS_DEFAULT));

  /*! trace all the rays in the given array of input rays, and write
   *  results into the given array of output rays. Unlike
   *  opTraceAsync, this call is synchronous, so will block the
   *  calling thread until all rays are traced and all hits have been
   *  written. Input and output arrays have to remain active and
   *  unmodified while teh call is active, but can be modified or
   *  released the moment this call returns */
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
                             size_t sizeOfRayAndHitArray,
                             /*! trace flags that can fine-tune how the trace
                              *  executes. */
                             /*! array of rays; must be (at least) as many as
                              *  'numRays' */
                             OPRay *arrayOfRays,
                             /*! array of where to write the results; as many as
                              *  'numRays' */
                             OPHit *arrayOfHits,
                             /*! size of ray and hit arrays, required
                               for being able to offload those
                               arrays if required */
                             OPTraceFlags traceFlags OP_IF_CPP(= OP_TRACE_FLAGS_DEFAULT));
  
  /*! SOA/Fortran verion of opTrace: traces all the rays in the given
   *  array of input rays, and write results into the given array of
   *  output rays. Unlike opTraceAsync, this call is synchronous, so
   *  will block the calling thread until all rays are traced and all
   *  hits have been written. 
   *
   *  Input and output arrays have to remain valid and unmodified
   *  while teh call is active, but can be modified or released the
   *  moment this call returns. 
   *
   *  Note that all IDs returned in the hit point will use _C_
   *  indexing convention (starting with 0), not fortran indexing */
  OP_API void opTraceSOA(/*! the model to trace into */
                         OPModel model,
                         /*! number of rays to trace - both arrays must have as
                          *  many entries */
                         size_t numRays,
                         /*! @{ arrays that describe the rays to trace: */
                         float *ray_origin_xs, 
                         float *ray_origin_ys, 
                         float *ray_origin_zs, 
                         float *ray_direction_xs,
                         float *ray_direction_ys,
                         float *ray_direction_zs,
                         float *ray_t_mins,
                         float *ray_t_maxs,
                         /*! @} */
                         /*! @{ arrays that will hold all the elemnets of the
                          *  computed hits; those that are null will get
                          *  ignored */
                         int32_t *hit_instIDs,
                         int32_t *hit_geomIDs,
                         int32_t *hit_primIDs,
                         float *hit_ts,
                         float *hit_us,
                         float *hit_vs,
                         /*! trace flags that can fine-tune how the trace
                          *  executes. */
                         OPTraceFlags traceFlags OP_IF_CPP(= OP_TRACE_FLAGS_DEFAULT));

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
                         const float *ray_origin_xs,
                         size_t stride_ray_origin_xs,
                         const float *ray_origin_ys, 
                         size_t stride_ray_origin_ys,
                         const float *ray_origin_zs, 
                         size_t stride_ray_origin_zs,
                         const float *ray_direction_xs,
                         size_t stride_ray_direction_xs,
                         const float *ray_direction_ys,
                         size_t stride_ray_direction_ys,
                         const float *ray_direction_zs,
                         size_t stride_ray_direction_zs,
                         const float *ray_t_mins,
                         size_t stride_ray_tmins,
                         const float *ray_t_maxs,
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
                         size_t stride_hit_uIDs,
                         float *hit_vs,
                         size_t stride_hit_vs,
                         /*! trace flags that can fine-tune how the trace
                          *  executes. */
                         OPTraceFlags traceFlags OP_IF_CPP(= OP_TRACE_FLAGS_DEFAULT));
  
#ifndef __cplusplus
}
#endif

#endif
