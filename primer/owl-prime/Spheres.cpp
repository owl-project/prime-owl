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

#include "Spheres.h"
#include "Context.h"

namespace op {

  OWLVarDecl Spheres::variables[]
  = {
     { "userID", OWL_INT, OWL_OFFSETOF(Spheres::SBTData,userID) },
     { "spheres", OWL_BUFPTR, OWL_OFFSETOF(Spheres::SBTData,spheres) },
     { nullptr }
  };
  
  Spheres::Spheres(Context *context,
                   uint32_t userID,
                   OWLBuffer spheresBuffer,
                   int numSpheres)
    : UserGeom(context,userID),
      spheresBuffer(spheresBuffer),
      numSpheres(numSpheres)
  {
    geom = owlGeomCreate(context->owl,
                         context->spheresGeomType);
    owlGeomSetPrimCount(geom,numSpheres);
    owlGeomSet1i(geom,"userID",userID);
    owlGeomSetBuffer(geom,"spheres",spheresBuffer);
  }

} // ::op
