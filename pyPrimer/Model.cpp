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

#include "pyPrimer/Model.h"
#include "pyPrimer/Context.h"

namespace pypr {

  Model::Model(std::shared_ptr<Context> context,
               const std::vector<Mesh::SP> &_meshes)
    : context(context)
  {
    std::vector<OPGeom> meshes;
    for (auto m : _meshes) meshes.push_back(m->primer);
    this->primer = opModelFromGeoms(context->primer,meshes.data(),meshes.size());
  }

  void Model::trace(const pybind11::buffer &orgs_buffer,
                    const pybind11::buffer &dirs_buffer,
                    const pybind11::buffer &tmins_buffer,
                    const pybind11::buffer &tmaxs_buffer,
                    const pybind11::buffer &hitIDs_buffer,
                    const pybind11::buffer &hitCoords_buffer
                    )
  {
    // ------------------------------------------------------------------
    pybind11::buffer_info orgs_info = orgs_buffer.request();
    if (orgs_info.format != pybind11::format_descriptor<float>::format())
      throw std::runtime_error("orgs buffer is not in floats");
    pybind11::array_t<float> orgs_array
      = pybind11::cast<pybind11::array_t<float, pybind11::array::c_style |
	  pybind11::array::forcecast>>(orgs_buffer);
    int orgs_count = orgs_array.size() / 3;
    const vec3f *orgs = (const vec3f *)orgs_array.request().ptr;

    // ------------------------------------------------------------------
    pybind11::buffer_info dirs_info = dirs_buffer.request();
    if (dirs_info.format != pybind11::format_descriptor<float>::format())
      throw std::runtime_error("dirs buffer is not in floats");
    pybind11::array_t<float> dirs_array
      = pybind11::cast<pybind11::array_t<float, pybind11::array::c_style |
	  pybind11::array::forcecast>>(dirs_buffer);
    int dirs_count = dirs_array.size() / 3;
    const vec3f *dirs = (const vec3f *)dirs_array.request().ptr;
    if (dirs_count != orgs_count)
      throw std::runtime_error("number of ray directions does not match number of ray origins");

    // ------------------------------------------------------------------
    pybind11::buffer_info tmins_info = tmins_buffer.request();
    if (tmins_info.format != pybind11::format_descriptor<float>::format())
      throw std::runtime_error("tmins buffer is not in floats");
    pybind11::array_t<float> tmins_array
      = pybind11::cast<pybind11::array_t<float, pybind11::array::c_style |
	  pybind11::array::forcecast>>(tmins_buffer);
    int tmins_count = tmins_array.size();
    const float *tmins = (const float *)tmins_array.request().ptr;
    if (tmins_count != orgs_count)
      throw std::runtime_error("number of ray tmins does not match number of ray origins");

    // ------------------------------------------------------------------
    pybind11::buffer_info tmaxs_info = tmaxs_buffer.request();
    if (tmaxs_info.format != pybind11::format_descriptor<float>::format())
      throw std::runtime_error("tmaxs buffer is not in floats");
    pybind11::array_t<float> tmaxs_array
      = pybind11::cast<pybind11::array_t<float, pybind11::array::c_style |
	  pybind11::array::forcecast>>(tmaxs_buffer);
    int tmaxs_count = tmaxs_array.size();
    const float *tmaxs = (const float *)tmaxs_array.request().ptr;
    if (tmaxs_count != orgs_count)
      throw std::runtime_error("number of ray tmaxes does not match number of ray origins");

    // ------------------------------------------------------------------
    pybind11::buffer_info hitCoords_info = hitCoords_buffer.request();
    if (hitCoords_info.format != pybind11::format_descriptor<float>::format())
      throw std::runtime_error("hitCoords buffer is not in floats");
    pybind11::array_t<float> hitCoords_array
      = pybind11::cast<pybind11::array_t<float>>(hitCoords_buffer);
    int hitCoords_count = hitCoords_array.size() / 3;
    vec3f *hitCoords = (vec3f *)hitCoords_array.request().ptr;
    if (hitCoords_count != orgs_count)
      throw std::runtime_error("number of hit coords does not match number of ray origins");

    // ------------------------------------------------------------------
    pybind11::buffer_info hitIDs_info = hitIDs_buffer.request();
    if (hitIDs_info.format != pybind11::format_descriptor<int>::format())
      throw std::runtime_error("hitIDs buffer is not in ints");
    pybind11::array_t<int> hitIDs_array
      = pybind11::cast<pybind11::array_t<int>>(hitIDs_buffer);
    int hitIDs_count = hitIDs_array.size() / 3;
    vec3i *hitIDs = (vec3i *)hitIDs_array.request().ptr;
    if (hitIDs_count != orgs_count)
      throw std::runtime_error("number of hit coords does not match number of ray origins");

    opTraceStridedSOA(this->primer,
                      orgs_count,
                      // origins
                      &orgs->x,
                      sizeof(vec3f),
                      &orgs->y,
                      sizeof(vec3f),
                      &orgs->z,
                      sizeof(vec3f),
                      // directions
                      &dirs->x,
                      sizeof(vec3f),
                      &dirs->y,
                      sizeof(vec3f),
                      &dirs->z,
                      sizeof(vec3f),
                      // tmin/max
                      tmins,
                      sizeof(float),
                      tmaxs,
                      sizeof(float),
                      // hit IDs
                      &hitIDs->x,
                      sizeof(vec3i),
                      &hitIDs->y,
                      sizeof(vec3i),
                      &hitIDs->z,
                      sizeof(vec3i),
                      // hit coords
                      &hitCoords->x,
                      sizeof(vec3f),
                      &hitCoords->y,
                      sizeof(vec3f),
                      &hitCoords->z,
                      sizeof(vec3f));
  }

  void Model::traceRays(const pybind11::buffer &rays_buffer,
                        const pybind11::buffer &hits_buffer)
  {
    // ------------------------------------------------------------------
    pybind11::buffer_info rays_info = rays_buffer.request();
    // if (rays_info.format != pybind11::format_descriptor<pypr::Ray>::format())
    //   throw std::runtime_error("rays buffer is not in <OPRay> format");
    pybind11::array_t<pypr::Ray> rays_array
      = pybind11::cast<pybind11::array_t<pypr::Ray, pybind11::array::c_style |
	  pybind11::array::forcecast>>(rays_buffer);
    int rays_count = rays_array.size();
    OPRay *rays = (OPRay *)rays_array.request().ptr;

    // ------------------------------------------------------------------
    pybind11::buffer_info hits_info = hits_buffer.request();
    // if (hits_info.format != pybind11::format_descriptor<PyOPHit>::format())
    //   throw std::runtime_error("hits buffer is not in <OPHit> format");
    pybind11::array_t<pypr::Hit> hits_array
      = pybind11::cast<pybind11::array_t<pypr::Hit, pybind11::array::c_style |
	  pybind11::array::forcecast>>(hits_buffer);
    int hits_count = hits_array.size();
    OPHit *hits = (OPHit *)hits_array.request().ptr;
    if (hits_count != rays_count)
      throw std::runtime_error("dimensions of rays and hits array must be the same");

    opTrace(this->primer,
            rays_count,
            // origins
            rays,
            hits);
  }
}
