// ======================================================================== //
// Copyright 2018-2023 Ingo Wald                                            //
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

// helper stuff to more easily create test geometry
#include "owl/common/math/AffineSpace.h"
// primer itself
#include "primer/primer.h"
// std stuff
#include <vector>
#include <iostream>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"
#include <random>

using namespace owl::common;

unsigned char to255(float f)
{
  if (f <= 0.f) return 0;
  if (f >= 1.f) return 255;
  return (unsigned char)(f*255.9f);
}

void savePNG(const std::string &fileName,
             vec2i fbSize,
             const vec3f *pixels)
{
  std::vector<unsigned char> rgba;
  for (int iy=fbSize.y-1;iy>=0;--iy) 
    for (int ix=0;ix<fbSize.x;ix++) {
      vec3f pixel = pixels[ix+fbSize.x*iy];
      rgba.push_back(to255(pixel.x));
      rgba.push_back(to255(pixel.y));
      rgba.push_back(to255(pixel.z));
      rgba.push_back(255);
    }
  std::cout << "#owl-prime.rtWeekend: writing image " << fileName << std::endl;
  stbi_write_png(fileName.c_str(),fbSize.x,fbSize.y,4,
                 rgba.data(),fbSize.x*sizeof(uint32_t));
}

inline float randomFloat()
{
  static std::random_device rd;  // Will be used to obtain a seed for the random number engine
  static std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  static std::uniform_real_distribution<> dis(0.f, 1.f);
  return dis(gen);
}

inline vec3f randomPoint()
{
  return vec3f(randomFloat(),randomFloat(),randomFloat());
}

inline vec3f randomDirection()
{
  while (true) {
    vec3f v = 2.f*randomPoint()-vec3f(1.f);
    if (dot(v,v) <= 1.f)
      return normalize(v);
  }
}

int main(int, char **)
{
  std::cout << "creating primer context" << std::endl;
  OPContext context
    = opContextCreate(OP_CONTEXT_DEFAULT,0);

  std::vector<float4> spheres;
  for (int i=0;i<100;i++) {
    vec3f P = 8.f*(randomPoint()-vec3f(.2f));
    float r = .2f + sqrtf(.3f*randomFloat());
    spheres.push_back(make_float4(P.x,P.y,P.z,r));
  }

  
  OPGeom geom
    = opSpheres4f(context,
                  /* user-supplied geometry ID */ 0,
                  spheres.data(),spheres.size());

  OPModel model
    = opModelFromGeoms(context,
                       &geom,1);
  
  std::cout << "generating rays..." << std::endl;
  vec3f up(0,1,0);
  vec3f at(0,0,0);
  vec3f from(-3,-2,-1);

  vec2i fbSize(800,600);
  vec3f dir = normalize(at-from);
  float imagePlaneHeight = 10.f;
  vec3f horiz = (imagePlaneHeight * fbSize.x/float(fbSize.y)) * normalize(cross(dir,up));
  vec3f vert  = imagePlaneHeight * normalize(cross(horiz,dir));

  std::vector<OPRay> rays;
  for (int iy=0;iy<fbSize.y;iy++)
    for (int ix=0;ix<fbSize.x;ix++) {
      float du = (ix+.5f)/fbSize.x;
      float dv = (iy+.5f)/fbSize.y;
      OPRay ray;
      (vec3f&)ray.origin = from + (du-.5f) * horiz + (dv-.5f) * vert;
      (vec3f&)ray.direction = dir;
      ray.tMin = 0.f;
      ray.tMax = INFINITY;
      rays.push_back(ray);
    }
  std::cout << "tracing rays..." << std::endl;
  std::vector<OPHit> hits(rays.size());
  opTrace(model,rays.size(),rays.data(),hits.data(),OP_TRACE_FLAGS_DEFAULT);

  std::vector<vec3f> pixels(rays.size());
  // ------------------------------------------------------------------
  std::cout << "creating 'primID' image" << std::endl;
  for (int i=0;i<rays.size();i++)
    pixels[i] = randomColor(hits[i].primID);
  savePNG("opSampleSpheres_primIDs.png",fbSize,pixels.data());
}
