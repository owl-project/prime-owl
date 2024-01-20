# (py)Primer: A "Optix Prime" like Abstraction for "Ray Tracing beyond Rendering"

<!-- ![](png/samplePicture.png) -->

This project aims at providing an as-simple-as-possible abstraction to
the problem of "tracing and intersecting rays with geometric objects",
with the explicit goal of making enabling this technology also for
applications *beyond* rendering.

To do this, this project explicitly aims at the following key design criteria:

- ability to use this library from non-C++ sources: in particular,
  support for interacting through C, Fortran, and Python

- portability across both GPU (via OptiX) and CPU (via Embree
  fallback)

- fast BVH construction and ray traversal through use of
  (hardware-)accelerated ray tracing libraries in the backend

- a "optix prime" like "ray traversal and intersection" API to tracing
  rays and returning intersections: tracing takes an array of rays,
  and return a buffer of hits, nothing else. In paricular, *no*
  support for any sort of "rendering" (i.e., no "materials", "paths",
  "colors", or any of that)

- an intentionally simple API that is explicitly aimed at those users
  that are *not* domain experts in graphics.
  
# API Overview

## C/C++/CUDA Interface

Generally speaking the primer API covers two kinds of operations:
Building a primer-"model" of the geometry that one wants to trace rays
against; and shooting (wavefronts of) rays to get the resulting
ray/model intersections (or "hits").

### Tracing (Wavefronts of) Rays (in C/C++/CUDA)

Tracing rays in C/C++/CUDA is very simple: 

1) include `"primer/primer.h"`

2) use primer API functions to build a `OPModel` representation of the
geometry to be traced again (see below)

3) create one array of rays (of type `OPRay`) to be traced, and a
second array of hits (`OPHit`) to store the results in.

4) call `opTrace(model,numRays,rays,hits)` to trace those rays.

Primer supports both host- and device-side arrays for rays and hits;
so you can _either_ use a CUDA-renderer that operates on device-side
arrays, or a host-side renderer that operates on host-side arrays.

The primer API is strictly C99, so can be used from C, C++, and CUDA
without changes. Rays and hits are simple structs, so can also be read
and written from either of these languages, on both host and device.

### Building an `OPModel` of a given scene/model

Before rays can be traced, one first needs to tell primer what
geometry to trace against. In primer, this is called a "model" (of
type `OPModel`). Primer supports two-level geometries with triangle
meshes, 'groups' of multiple meshes that logically go together, and
instances of such groups (future versions may also support spheres, or
other geometry types than just triangles).

Before any model can be created, one first needs to create a primer
"context" (`OPContext`), which basically managed how exactly the
primer operations are to be executed. In particular, primer supports
the following types of context:

- `OP_CONTEXT_GPU`: tracing happens on a GPU, ray and hit arrays must
  be device-accessible arrays (ie, either CUDA managed or CUDA device
  memory). This requires a suitable GPU or will cause an error.
  
- `OP_CONTEXT_HOST_FORCE_CPU`: tracing is executed on the CPU even if
  a suitable GPU is available; ray and hit arrays must be host
  readable. This requires an installed version of embree for host-side tracing.
  
- `OP_CONTEXT_HOST_FORCE_OFFLOAD`: this will trace on the GPU (where
  available), but will automatically "offload" any potential host-side
  ray and hit arrays to the GPU as required. 
  
- `OP_CONTEXT_DEFAULT`: this will try to create an offload device
  (where trace happens on the GPU, and which can consume both host-
  and device-side ray and hit arrays), but will fall back to a host
  backend if no suitable GPU can be found.
  
One a context has been created (via `opContextCreate(...)`), there are
different ways for creating a model. Under the hood primer uses a
two-level scene representation with instances, and there are a variety
of functions to create triangle meshes (`opMeshCreate(...)`,
`opMeshCreateSOA(...)`, opMeshFromTriangles(...)`), Groups
(`opGroupCreate(...)`), or a model over instances of such groups
(`opModelFromInstances(...)`).

In addition, there are also several helper functions for simpler cases
where theinput model is a single triangle mesh
(`opModelFromGeom(...)`), etc.

## Fortran Interface

Using primer from Fortran works through the same `extern C` linkage
shared library/DLL mechanism as for C/C++. However, to make the library more
friendly to native fortran applications the API explicitly also supports:

- support for describing geometry with various fortran-style
  structure-of-array data layouts instead of array-of-structure
  layouts

- support for describing geometry indexing "fortran style" (i.e.,
  starting with index 1 instead of the C/C++ way of starting with 0)
  
- (NOT YET IMPLEMENTED:) support for using models, rays, and hits with
  fortran-style int64 and double-precision coordinates (internally all
  data is single-precision, but will be converted on the fly as
  required).

## Python Interface

To enable the use of this library through/from python primer also
comes with a `pypr` python module that offers `pybind11` based python
bindings to the C/C++ backend. Some high-level features of this:

- All "bulk" data (arrays of vertices and indices for scene creation,
  and arrays of rays/hits during trace) are passed and returned
  through `numpy` arrays.
  
- A pythonic way of accessing the library through a loadable module, such as

``` python
#!/usr/bin/python3
import pypr
import numpy as np
...
##
## create model from mesh(es) of triangles
##
mesh_vertices= np. ...
mesh_indices = np. ...
mesh         = pypr.create_mesh(0,mesh_vertices,mesh_indices)
model        = pypr.create_model([mesh])
##
## allocate array of rays and generate them using python/numpy
##
ray_orgs = np.ndarray((num_rays,3),dtype=np.float32)
ray_dirs = np.ndarray((num_rays,3),dtype=np.float32)
...
##
## allocate array(s) of hits:
##
hit_ints   = np.ndarray((num_rays,3),dtype=np.int32)
hit_floats = np.ndarray((num_rays,3),dtype=np.float32)
##
## and finally, trace the rays:
##
model.trace(ray_orgs,ray_dirs,ray_tmins,ray_tmaxs,hit_ints,hit_floats)
```

There are a few simply python examples in the `pyPrimer/` directory

To use this python extension module: make sure to do a

``` bash
export PYTHONPATH=.
```

before trying to run your code using this module. E.g.:

``` bash
export PYTHONPATH=.
python3 pyPrimer/simplestPossibleExample.py 
```

# Building this library

This projects builds on top of OWL;
http://owl-project.github.io. Prime has pretty much the same
dependencies and build process as owl, so maybe a good way to start is
to get owl to build first, using the instructions in that project -
then prime should build easily, and in exactly the same way.

To build, you need both projects checked out as sisters in the
directory tree (we'll switch to git submodule later on), then prime
should pick up owl automatically, and will build it as part of its own
build process (ie, you do _not_ have to have a version of owl
"installed" on your system).

## Building the CPU-fallback


### Install Dependencies

``` bash
sudo apt install cmake g++ build-essential libtbb-dev
```

For the pyton part, you probably also want to do:

``` bash
sudo apt install sudo apt install python3-pip
pip3 install numpy
```



### Build Embree

(alternatively, install from binary distribution, if you know how to)

``` bash
git clone --recursive https://github.com/embree/embree -b v3.9.0
mkdir embree/bin
cd embree/bin
cmake .. -DEMBREE_ISPC_SUPPORT=OFF -DEMBREE_TUTORIALS=OFF
make 
sudo make install
```

### Building Primer

``` bash
git clone --recursive https://gitlab.com/ingowald/prime-owl
mkdir prime-owl/bin
cd prime-owl/bin
cmake .. 
make 
```

Install CUDA, e.g.:
``` bash
wget https://developer.download.nvidia.com/compute/cuda/11.4.2/local_installers/cuda_11.4.2_470.57.02_linux.runsudo sh cuda_11.4.2_470.57.02_linux.run
```


# Known Bugs/Limitations

- currently supporting only triangles and triangle meshes

- double and int64 format conversions not yet implemented



