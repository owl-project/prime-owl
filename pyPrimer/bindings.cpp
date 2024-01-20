// ======================================================================== //
// Copyright 2020-2021 Ingo Wald                                            //
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

#include "pyPrimer/Context.h"
#include "pyPrimer/Model.h"

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

using namespace pypr;

namespace pypr {

  struct CurrentContext {
    static Context::SP create()
    {
      handle = std::make_shared<Context>();
      return handle;
    }

    static Context::SP get()
    {
      if (!handle)
        throw std::runtime_error("#primer: no context active"
                                 " - did you forget pypr.context_create()?");
      return handle;
    }
    
    static Context::SP handle;
  };
  Context::SP CurrentContext::handle = {};

  /*! creates a new context, and makes it current */
  std::shared_ptr<Context> createContext()
  {
    return CurrentContext::create();
  }

  /*! create a mesh in the current default context */
  Mesh::SP createMesh(int userID,
                      const pybind11::buffer &vertices,
                      const pybind11::buffer &indices)
  { return std::make_shared<Mesh>(CurrentContext::get(),userID,vertices,indices); }
                      
  /*! create a model in the current default context */
  Model::SP createModel(const std::vector<Mesh::SP> &meshes)
  { return std::make_shared<Model>(CurrentContext::get(),meshes); }
                      
}

PYBIND11_MODULE(pypr, m) {

  PYBIND11_NUMPY_DTYPE(pypr::Ray,org_x,org_y,org_z,t_min,dir_x,dir_y,dir_z,t_max );
  PYBIND11_NUMPY_DTYPE(pypr::Hit,primID,geomID,instID,t,u,v);
    
  // optional module docstring
  m.doc() = "pyPr: prime-OWL python wrappers";

  // create one context; almost all functions tare then per context
  m.def("context_create", &pypr::createContext,
        "Creates a 'primer' Context object");
  m.def("create_mesh",
        &pypr::createMesh);
  m.def("create_model",
        &pypr::createModel);

  // create one context; almost all functions tare then per context
  // m.def("mesh_create", &pypr::createMesh,
  //       "Creates a 'primer' Mesh object");
  
  // m.def("mesh_create", &pypr::defaultContext_createMesh,
  //       "Creates a 'primer' Context object");

  // define add function
  // m.def("save_png_rgba8", &pypr::save_png_rgba8,
  //       "Saves a OWLBuffer of rgba8 values to a given file, in PNG format");
  
  // -------------------------------------------------------
  auto model
    = pybind11::class_<pypr::Model,
                 std::shared_ptr<Model>>(m, "model");
  model.def("trace",
            &Model::trace, "traces array of rays and produces array of hits");
  model.def("trace",
            &Model::traceRays, "traces an array of rays, and fills in an array of hits. Each ray in the rays array is of the format [org_x,org_y,org_z,tMin,dir_x,dir_y,dir_z,tMax]; each hit in the hit arry is of the format [primID,geomID,instID,t,u,v]");
  
  // -------------------------------------------------------
  auto mesh
    = pybind11::class_<pypr::Mesh,
                 std::shared_ptr<Mesh>>(m, "mesh");

#if 0  
  // -------------------------------------------------------
  auto geomType
    = pybind11::class_<pypr::GeomType,
                 std::shared_ptr<GeomType>>(m, "GeomType");
  geomType.def("set_closest_hit", &pypr::GeomType::setClosestHit);

  // -------------------------------------------------------
  auto geom
    = pybind11::class_<pypr::Geom,
                 std::shared_ptr<Geom>>(m, "Geom");

  geom.def("set_vertices",&pypr::Geom::setVertices);
  geom.def("set_indices", &pypr::Geom::setIndices);
  geom.def("set_buffer", &pypr::Geom::setBuffer);
  geom.def("set_3f", &pypr::Geom::set3f);
  
  // -------------------------------------------------------
  auto missProg
    = pybind11::class_<pypr::MissProg,
                 std::shared_ptr<MissProg>>(m, "MissProg");

  missProg.def("set_3f", &pypr::MissProg::set3f);
  
  // -------------------------------------------------------
  auto rayGen
    = pybind11::class_<pypr::RayGen,
                 std::shared_ptr<RayGen>>(m, "RayGen");

  rayGen.def("set_2i", &pypr::RayGen::set2i);
  rayGen.def("set_3f", &pypr::RayGen::set3f);
  rayGen.def("set_buffer", &pypr::RayGen::setBuffer);
  rayGen.def("set_group", &pypr::RayGen::setGroup);
  rayGen.def("launch",
              &pypr::RayGen::launch2D);
  
  // -------------------------------------------------------
  auto group
    = pybind11::class_<pypr::Group,
                 std::shared_ptr<Group>>(m, "Group");
  group.def("build_accel", &pypr::Group::buildAccel);

  // -------------------------------------------------------
  auto buffer
    = pybind11::class_<pypr::Buffer,
                 std::shared_ptr<Buffer>>(m, "Buffer");
#endif
  
  // -------------------------------------------------------
  auto context
    = pybind11::class_<pypr::Context,
                 std::shared_ptr<Context>>(m, "Context");

  context.def(pybind11::init<>());
#if 0
  context.def("create_mesh",
               &pypr::Context::createMesh);
  context.def("module_from_file",
              &pypr::Context::createModuleFromFile);
  context.def("module_create",
              &pypr::Context::createModuleFromString);
  context.def("geom_type_create", &pypr::Context::createGeomType);
  context.def("geom_create", &pypr::Context::createGeom);
  context.def("miss_prog_create", &pypr::Context::createMissProg);
  context.def("ray_gen_create", &pypr::Context::createRayGen);
  context.def("triangles_geom_group_create", &pypr::Context::createTrianglesGeomGroup);
  context.def("instance_group_create", &pypr::Context::createInstanceGroup);
  context.def("device_buffer_create", &pypr::Context::createDeviceBuffer);
  context.def("host_pinned_buffer_create", &pypr::Context::createHostPinnedBuffer);

  context.def("build_programs",&pypr::Context::buildPrograms);
  context.def("build_pipeline",&pypr::Context::buildPipeline);
  context.def("build_SBT",&pypr::Context::buildSBT);
  context.def("context_destroy",&pypr::Context::destroy);
  
  // geom kinds 
  context.attr("GEOM_TRIANGLES") = pybind11::int_((int)OWL_GEOM_TRIANGLES);
  context.attr("GEOM_USER")      = pybind11::int_((int)OWL_GEOM_USER);
  // variable/data types
  context.attr("BUFPTR") = pybind11::int_((int)OWL_BUFPTR);
  context.attr("INT")    = pybind11::int_((int)OWL_INT);
  context.attr("INT2")   = pybind11::int_((int)OWL_INT2);
  context.attr("INT3")   = pybind11::int_((int)OWL_INT3);
  context.attr("INT4")   = pybind11::int_((int)OWL_INT4);
  context.attr("FLOAT")  = pybind11::int_((int)OWL_FLOAT);
  context.attr("FLOAT2") = pybind11::int_((int)OWL_FLOAT2);
  context.attr("FLOAT3") = pybind11::int_((int)OWL_FLOAT3);
  context.attr("FLOAT4") = pybind11::int_((int)OWL_FLOAT4);
#endif
}
