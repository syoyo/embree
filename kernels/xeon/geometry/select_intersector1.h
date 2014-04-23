// ======================================================================== //
// Copyright 2009-2013 Intel Corporation                                    //
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

#include "common/ray.h"

namespace embree
{
  template<typename I0, typename I1>
  struct Select2Intersector1
  {
    typedef void Primitive;

    struct Precalculations 
    {
      __forceinline Precalculations (const Ray& ray) 
        : pre0(ray), pre1(ray) {}

    public:
      typename I0::Precalculations pre0;
      typename I1::Precalculations pre1;
    };

    static __forceinline void intersect(const Precalculations& pre, Ray& ray, size_t ty, void* prims, size_t num, void* geom)
    {
      if (likely(ty == 0)) 
        I0::intersect(pre.pre0,ray,ty,prims,num,geom);
      else if (unlikely(ty == 1))
        I1::intersect(pre.pre1,ray,ty,prims,num,geom);
    }

    static __forceinline bool occluded(const Precalculations& pre, Ray& ray, size_t ty, void* prims, size_t num, void* geom) 
    {
      if (likely(ty == 0)) 
        return I0::occluded(pre.pre0,ray,ty,prims,num,geom);
      else if (unlikely(ty == 1))
        return I1::occluded(pre.pre1,ray,ty,prims,num,geom);
      else 
        return false;
    }
  };

  template<typename I0, typename I1, typename I2>
  struct Select3Intersector1
  {
    typedef void Primitive;

    struct Precalculations 
    {
      __forceinline Precalculations (const Ray& ray) 
        : pre0(ray), pre1(ray), pre2(ray) {}

    public:
      typename I0::Precalculations pre0;
      typename I1::Precalculations pre1;
      typename I2::Precalculations pre2;
    };

    static __forceinline void intersect(const Precalculations& pre, Ray& ray, size_t ty, void* prims, size_t num, void* geom)
    {
      if (likely(ty == 0)) 
        I0::intersect(pre.pre0,ray,ty,prims,num,geom);

      else switch (ty) {
        case 1: I1::intersect(pre.pre1,ray,ty,prims,num,geom); break;
        case 2: I2::intersect(pre.pre2,ray,ty,prims,num,geom); break;
        }
    }

    static __forceinline bool occluded(const Precalculations& pre, Ray& ray, size_t ty, void* prims, size_t num, void* geom) 
    {
      if (likely(ty == 0)) 
        return I0::occluded(pre.pre0,ray,ty,prims,num,geom);

      else switch (ty) {
        case 1: return I1::occluded(pre.pre1,ray,ty,prims,num,geom);
        case 2: return I2::occluded(pre.pre2,ray,ty,prims,num,geom);
        }
      return false;
    }
  };
}
