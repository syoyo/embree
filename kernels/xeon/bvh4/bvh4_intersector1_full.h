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

#include "bvh4.h"
#include "common/ray.h"
#include "common/stack_item.h"

namespace embree
{
  namespace isa
  {
    /*! BVH4 single ray traversal implementation. */
    template<typename PrimitiveIntersector>
      class BVH4Intersector1Full 
    {
      static const size_t stackSize = 1+3*BVH4::maxDepth;

      struct __aligned(16) StackItem 
      {
      public:
        __forceinline static void swap2(StackItem& a, StackItem& b) 
        { 
#if defined(__AVX__)
          ssef sse_a = load4f(&a);
          ssef sse_b = load4f(&b);
          store4f(&a,sse_b);
          store4f(&b,sse_a);
#else
          StackItem t = b; b = a; a = t;
#endif
        }

        __forceinline friend bool operator<(const StackItem& s1, const StackItem& s2) {
          return s1.tNear > s2.tNear;
        }

        /*! Sort 2 stack items. */
        __forceinline friend void sort(StackItem& s1, StackItem& s2) {
          if (s2.tNear < s1.tNear) swap2(s2,s1);
        }
        
        /*! Sort 3 stack items. */
        __forceinline friend void sort(StackItem& s1, StackItem& s2, StackItem& s3)
        {
          if (s2.tNear < s1.tNear) swap2(s2,s1);
          if (s3.tNear < s2.tNear) swap2(s3,s2);
          if (s2.tNear < s1.tNear) swap2(s2,s1);
        }
    
        /*! Sort 4 stack items. */
        __forceinline friend void sort(StackItem& s1, StackItem& s2, StackItem& s3, StackItem& s4)
        {
          if (s2.tNear < s1.tNear) swap2(s2,s1);
          if (s4.tNear < s3.tNear) swap2(s4,s3);
          if (s3.tNear < s1.tNear) swap2(s3,s1);
          if (s4.tNear < s2.tNear) swap2(s4,s2);
          if (s3.tNear < s2.tNear) swap2(s3,s2);
        }

      public:
        size_t ref;
        float tNear,tFar;
      };

    private:
      static __forceinline size_t intersectBox(const BVH4::UANode* node, const sse3f& org, const sse3f& rdir, const sse3f& org_rdir, 
                                               const size_t nearX, const size_t nearY, const size_t nearZ, ssef& tNear, ssef& tFar);
      static __forceinline size_t intersectBox(const BVH4::CANode* node, const sse3f& org, const sse3f& rdir, const sse3f& org_rdir, 
                                               const size_t nearX, const size_t nearY, const size_t nearZ, ssef& tNear, ssef& tFar);
      static size_t intersectBox(const BVH4::UUNode* node, Ray& ray, const sse3f& org, const sse3f& dir, ssef& tNear, ssef& tFar);
      static size_t intersectBox(const BVH4::CUNode* node, Ray& ray, const sse3f& org, const sse3f& dir, ssef& tNear, ssef& tFar);

#if defined(__AVX__)
    static size_t intersectBox(const BVH4::UUNode* node, Ray& ray, 
                                      const avx3f& ray_org_dir, const sse3f& ray_org, const sse3f& ray_dir, 
                                      ssef& tNear, ssef& tFar);
#endif
      
    public:
      static void intersect(const BVH4* This, Ray& ray);
      static void occluded (const BVH4* This, Ray& ray);
    };
  }
}
