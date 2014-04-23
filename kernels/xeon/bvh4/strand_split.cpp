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
// distributed under the License is dist ributed on an "AS IS" BASIS,       //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#include "strand_split.h"

namespace embree
{
  void StrandSplit::find (BezierRefList& beziers, float bezierCost)
  {
    /* first try to split two hair strands */
    BezierRefList::block_iterator_unsafe i = beziers;
    Vec3fa axis0 = normalize(i->p3 - i->p0);
    float bestCos = 1.0f;
    Bezier1 bestI = *i;

    for (i; i; i++) {
      Vec3fa axisi = i->p3 - i->p0;
      float leni = length(axisi);
      if (leni == 0.0f) continue;
      axisi /= leni;
      float cos = abs(dot(axisi,axis0));
      if (cos < bestCos) { bestCos = cos; bestI = *i; }
    }
    Vec3fa axis1 = normalize(bestI.p3-bestI.p0);

    split.axis0 = axis0;
    split.axis1 = axis1;
    
    size_t size0 = 0;
    size_t size1 = 0;
    BBox3fa bounds0 = empty;
    BBox3fa bounds1 = empty;
    const LinearSpace3fa space0 = frame(axis0).transposed();
    const LinearSpace3fa space1 = frame(axis1).transposed();
    
    for (BezierRefList::block_iterator_unsafe i=beziers; i; i++)
    {
      const Bezier1& prim = *i; 
      const Vec3fa axisi = normalize(prim.p3-prim.p0);
      const float cos0 = abs(dot(axisi,axis0));
      const float cos1 = abs(dot(axisi,axis1));
      if (cos0 > cos1) {
        bounds0.extend(prim.bounds(space0));
        size0++;
      }
      else {
        bounds1.extend(prim.bounds(space1));
        size1++;
      }
    }

    split.cost = bezierCost*(size0*safeHalfArea(bounds0) + size1*safeHalfArea(bounds1));
  }
    
  void StrandSplit::Split::split(size_t threadIndex, PrimRefBlockAlloc<Bezier1>* alloc, BezierRefList& prims, BezierRefList& lprims, BezierRefList& rprims)
  {
    BezierRefList::item* lblock = lprims.insert(alloc->malloc(threadIndex));
    BezierRefList::item* rblock = rprims.insert(alloc->malloc(threadIndex));
    
    while (BezierRefList::item* block = prims.take()) 
    {
      for (size_t i=0; i<block->size(); i++) 
      {
        const Bezier1& prim = block->at(i); 
        const Vec3fa axisi = normalize(prim.p3-prim.p0);
        const float cos0 = abs(dot(axisi,axis0));
        const float cos1 = abs(dot(axisi,axis1));
        
        if (cos0 > cos1) 
        {
          if (likely(lblock->insert(prim))) continue; 
          lblock = lprims.insert(alloc->malloc(threadIndex));
          lblock->insert(prim);
        } 
        else 
        {
          if (likely(rblock->insert(prim))) continue;
          rblock = rprims.insert(alloc->malloc(threadIndex));
          rblock->insert(prim);
        }
      }
      alloc->free(threadIndex,block);
    }
  }
}
