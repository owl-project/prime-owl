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

// public API
#include <primer/primer.h>
// internal stuff
#include <owl/common/math/vec.h>
#include <primer-common/Ray.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <random>
#include <cmath>
#include <fstream>

using namespace owl::common;

typedef std::vector<box3f> Boxes;
/*! array of patches; each patch is (vec3f anchorPoint, vec3f du, vec3f dv) */
struct Patch { vec3f anchor, du, dv; float value=0.f; };
typedef std::vector<Patch> Patches;
struct WorldConfig {
#if 1
  int   numBoxes = 100;
  vec2f boxSizeRange = { .2f,10.f };
  float worldSize = 50.f;
#else
  int   numBoxes = 1000;
  vec2f boxSizeRange = { .2,10 };
  float worldSize = 100.f;
#endif
};

struct Model {
  std::vector<vec3f> vertices;
  std::vector<vec3i> indices;
};


Boxes   createBoxes(const WorldConfig &config);
Patches createPatches(const Boxes &boxes, int faceRes = 16);
Model   createModel(const Boxes &boxes);
std::vector<OPRay> createRays(const Patches &patches, int numSamplesPerPatch=64);
std::vector<OPHit> computeHits(const Model &model,
                               const std::vector<OPRay> &rays);
void integratePatches(Patches &patches, const std::vector<OPHit> &hits);
void savePatchesToOBJ(const std::string &fileName,
                      const Patches &patches);


Boxes   createBoxes(const WorldConfig &config)
{
  std::random_device r;
  std::mt19937 rng(r());
 
  std::uniform_real_distribution<float>
    size_sample(config.boxSizeRange.x, config.boxSizeRange.y);
  std::uniform_real_distribution<float>
    pos_sample(0.f, config.worldSize);

  
  Boxes boxes(config.numBoxes);
  for (auto &box : boxes) {
    vec3f pos = vec3f(pos_sample(rng),
                      pos_sample(rng),
                      pos_sample(rng));
    vec3f size = vec3f(size_sample(rng),
                       size_sample(rng),
                       size_sample(rng));
    box = box3f{pos-0.5f*size,pos+0.5f*size};
  }
  return boxes;
}

Model   createModel(const Boxes &boxes)
{
  Model result;
  const int NUM_VERTICES = 8;
  static const vec3f unitBoxVertices[NUM_VERTICES] =
    {
      { 0.f,  0.f,  0.f},
      {+1.f,  0.f,  0.f},
      {+1.f, +1.f,  0.f},
      { 0.f, +1.f,  0.f},
      { 0.f, +1.f, +1.f},
      {+1.f, +1.f, +1.f},
      {+1.f,  0.f, +1.f},
      { 0.f,  0.f, +1.f},
    };

  const int NUM_INDICES = 12;
  static const vec3i unitBoxIndices[NUM_INDICES] =
    {
      {0, 2, 1}, //face front
      {0, 3, 2},
      {2, 3, 4}, //face top
      {2, 4, 5},
      {1, 2, 5}, //face right
      {1, 5, 6},
      {0, 7, 4}, //face left
      {0, 4, 3},
      {5, 4, 7}, //face back
      {5, 7, 6},
      {0, 6, 7}, //face bottom
      {0, 1, 6}
    };

  for (auto box : boxes) {
    int beginIndex = (int)result.vertices.size();
    for (auto vtx : unitBoxVertices) 
      result.vertices.push_back(box.lower + vtx*box.size());
    for (auto idx : unitBoxIndices)
      result.indices.push_back(beginIndex+idx);
  }
  return result;
}

void addFace(Patches &patches,
             const vec3f anchor,
             const vec3f du,
             const vec3f dv,
             int faceRes)
{
  for (int i=0;i<faceRes;i++)
    for (int j=0;j<faceRes;j++) {
      Patch p;
      p.anchor = anchor+(i+0.f)*du/faceRes+(j+0.f)*dv/faceRes;
      p.du = du * (1.f/faceRes);
      p.dv = dv * (1.f/faceRes);
      patches.push_back(p);
    }
}

Patches createPatches(const Boxes &boxes, int faceRes)
{
  Patches result;
  for (auto box : boxes) {
    addFace(result,box.lower,box.size()*vec3f(0,0,1),box.size()*vec3f(0,1,0),faceRes);
    addFace(result,box.lower,box.size()*vec3f(1,0,0),box.size()*vec3f(0,0,1),faceRes);
    addFace(result,box.lower,box.size()*vec3f(0,1,0),box.size()*vec3f(1,0,0),faceRes);

    addFace(result,box.upper,box.size()*vec3f(0,-1,0),box.size()*vec3f(0,0,-1),faceRes);
    addFace(result,box.upper,box.size()*vec3f(0,0,-1),box.size()*vec3f(-1,0,0),faceRes);
    addFace(result,box.upper,box.size()*vec3f(-1,0,0),box.size()*vec3f(0,-1,0),faceRes);
  }
  return result;
}

std::vector<OPRay> createRays(const Patches &patches, int numSamplesPerPatch)
{
  std::random_device r;
  std::mt19937 rng(r());
  
  std::uniform_real_distribution<float> uniform(0.f,1.f);
  
  std::vector<OPRay> results;
  for (auto patch : patches) {
    for (int i=0;i<numSamplesPerPatch;i++) {
      vec3f N = normalize(cross(patch.du,patch.dv));
      vec3f org
        = patch.anchor
        + uniform(rng)*patch.du
        + uniform(rng)*patch.dv;
      // org = org+1e-3f*N;
    
      float r1 = uniform(rng);
      float r2 = uniform(rng);
      vec3f dirSample;
      dirSample.z = r1;
      dirSample.x = cosf(2.f*(float)M_PI*r2)*sqrtf(1.f-r1*r1);
      dirSample.y = sinf(2.f*(float)M_PI*r2)*sqrtf(1.f-r1*r1);
      vec3f dir = xfmVector(owl::common::frame(N),dirSample);
      // dir = N;
      
      float tmin = 1e-3f;
      float tmax = 1e20f;

      OPRay ray;
      ray.origin[0] = org.x;
      ray.origin[1] = org.y;
      ray.origin[2] = org.z;
      ray.direction[0] = dir.x;
      ray.direction[1] = dir.y;
      ray.direction[2] = dir.z;
      ray.tMin = tmin;
      ray.tMax = tmax;
      results.push_back(ray);
    }
  }
  return results;
}

  inline vec3f hue_to_rgb(float hue)
  {
    float s = saturate( hue ) * 6.0f;
    float r = saturate( fabsf(s - 3.f) - 1.0f );
    float g = saturate( 2.0f - fabsf(s - 2.0f) );
    float b = saturate( 2.0f - fabsf(s - 4.0f) );
    return vec3f(r, g, b); 
  }
    
  inline vec3f temperature_to_rgb(float t)
  {
    float K = 4.0f / 6.0f;
    float h = K - K * t;
    float v = .5f + 0.5f * t;    return v * hue_to_rgb(h);
  }
    
                                    
  inline
  vec3f heatMap(float t)
  {
#if 1
    return temperature_to_rgb(t);
#else
    if (t < .25f) return lerp(vec3f(0.f,1.f,0.f),vec3f(0.f,1.f,1.f),(t-0.f)/.25f);
    if (t < .5f)  return lerp(vec3f(0.f,1.f,1.f),vec3f(0.f,0.f,1.f),(t-.25f)/.25f);
    if (t < .75f) return lerp(vec3f(0.f,0.f,1.f),vec3f(1.f,1.f,1.f),(t-.5f)/.25f);
    if (t < 1.f)  return lerp(vec3f(1.f,1.f,1.f),vec3f(1.f,0.f,0.f),(t-.75f)/.25f);
    return vec3f(1.f,0.f,0.f);
#endif
  }

void savePatchesToOBJ(const std::string &fileName,
                      const Patches &patches)
{
  std::ofstream obj(fileName+".obj");
  std::ofstream mtl(fileName+".mtl");
  obj << "mtllib "<<fileName<<".mtl" << std::endl;
  for (auto patch : patches) {
    vec3f v0 = patch.anchor;
    vec3f v1 = patch.anchor+patch.du;
    vec3f v2 = patch.anchor+patch.du+patch.dv;
    vec3f v3 = patch.anchor+patch.dv;
    static int g_primID = 0;
    int primID = g_primID++;
    vec3f c = heatMap(patch.value);//owl::randomColor(primID);
    mtl << "newmtl m"  << primID << std::endl;
    mtl << "Kd " << c.x << " " << c.y << " " << c.z << std::endl;
    mtl << "Ka " << c.x << " " << c.y << " " << c.z << std::endl;
    mtl << std::endl;
    obj << "usemtl m" << primID << std::endl;
    obj << "v " << v0.x << " " << v0.y << " " << v0.z << std::endl;
    obj << "v " << v1.x << " " << v1.y << " " << v1.z << std::endl;
    obj << "v " << v2.x << " " << v2.y << " " << v2.z << std::endl;
    obj << "v " << v3.x << " " << v3.y << " " << v3.z << std::endl;
    obj << "f -1 -2 -3 -4" << std::endl;
  }
}


std::vector<OPHit> computeHits(const Model &model,
                               const std::vector<OPRay> &rays)
{
  std::vector<OPHit> hits(rays.size());
  OPContext context = opContextCreate(OP_CONTEXT_HOST_DEFAULT,0);
  OPGeom mesh = opMeshCreate(context,0,
                             (const float *)model.vertices.data(),
                             model.vertices.size(),
                             sizeof(vec3f),
                             (const int *)model.indices.data(),
                             model.indices.size(),
                             sizeof(vec3i));
  OPModel _model = opModelFromGeoms(context,&mesh,1);

  opTrace(_model,rays.size(),
          (OPRay*)rays.data(),
          (OPHit*)hits.data());
  return hits;
}

void integratePatches(Patches &patches, const std::vector<OPHit> &hits)
{
  int numRaysPerPatch = int(hits.size()/patches.size());
  for (int patchID=0;patchID<patches.size();patchID++) {
    int unoccluded = 0;
    for (int i=0;i<numRaysPerPatch;i++)
      if (hits[patchID*numRaysPerPatch+i].primID < 0)
        unoccluded++;
    patches[patchID].value = unoccluded / float(numRaysPerPatch);
  }
}

template<typename T>
void write(std::ofstream &out, const T &v)
{ out.write((char*)&v,sizeof(v)); }

template<typename T>
void writeVec(std::ofstream &out, const std::vector<T> &vec)
{ for (auto &v : vec) write(out,v); }

int main(int ac, char **av)
{
  std::string resultFileBase = "fromArrays-testData";
  WorldConfig config;
  Boxes boxes = createBoxes(config);
  Model model = createModel(boxes);
  
  Patches patches = createPatches(boxes); 
  std::vector<OPRay> rays = createRays(patches);

  std::vector<OPHit> hits = computeHits(model,rays);
  integratePatches(patches, hits);

  std::ofstream vertices(resultFileBase+".vertices",std::ios::binary);
  writeVec(vertices,model.vertices);
  
  std::ofstream indices(resultFileBase+".indices",std::ios::binary);
  writeVec(indices,model.indices);

  std::ofstream ray_orgs(resultFileBase+".ray_orgs",std::ios::binary);
  std::ofstream ray_dirs(resultFileBase+".ray_dirs",std::ios::binary);
  std::ofstream ray_tmins(resultFileBase+".ray_tmins",std::ios::binary);
  std::ofstream ray_tmaxs(resultFileBase+".ray_tmaxs",std::ios::binary);
  for (auto ray : rays) {
    write(ray_orgs,ray.origin);
    write(ray_dirs,ray.direction);
    write(ray_tmins,ray.tMin);
    write(ray_tmaxs,ray.tMax);
  }

  std::ofstream hit_ids(resultFileBase+".hit_ids",std::ios::binary);
  std::ofstream hit_coords(resultFileBase+".hit_coords",std::ios::binary);
  for (auto hit : hits) {
    write(hit_ids,hit.instID);
    write(hit_ids,hit.geomID);
    write(hit_ids,hit.primID);
    write(hit_coords,hit.t);
    write(hit_coords,hit.u);
    write(hit_coords,hit.v);
  }

  savePatchesToOBJ(resultFileBase+"-result",patches);
}
