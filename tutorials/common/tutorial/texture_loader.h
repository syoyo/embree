// ======================================================================== //
// Copyright 2009-2015 Intel Corporation                                    //
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

#pragma once

#include "sys/platform.h"
#include "sys/filename.h"
#include "sys/vector.h"
#include "math/vec2.h"
#include "math/vec3.h"
#include "math/affinespace.h"
#include "scene.h"

#include <vector>
#include <memory>

namespace embree
{
  /*! read texture from disk */
  OBJScene::Texture *loadTexture(const FileName& fileName);
  OBJScene::Texture *loadPtexTexture(const FileName& fileName);   
}
