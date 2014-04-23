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

#include "object_type_partition.h"

namespace embree
{
  ObjectTypePartitioning::ObjectTypePartitioning (TriRefList& tris, float triCost, BezierRefList& beziers, float bezierCost)
  {
    PrimInfo pinfo;
    pinfo.num = 0;
    pinfo.numTriangles = 0;
    pinfo.numBeziers = 0;
    pinfo.geomBounds = empty;
    pinfo.centBounds = empty;
    
    BBox3fa triGeomBounds = empty;
    BBox3fa triCentBounds = empty;
    TriRefList::iterator t=tris;
    while (TriRefBlock* block = t.next()) 
    {
      pinfo.num += block->size();
      pinfo.numTriangles += block->size();
      for (size_t i=0; i<block->size(); i++)
      {
        const BBox3fa bounds = block->at(i).bounds(); 
        triGeomBounds.extend(bounds);
        triCentBounds.extend(center2(bounds));
      }
    }
    pinfo.geomBounds.extend(triGeomBounds);
    pinfo.centBounds.extend(triCentBounds);

    BBox3fa bezierGeomBounds = empty;
    BBox3fa bezierCentBounds = empty;
    BezierRefList::iterator b=beziers;
    while (BezierRefBlock* block = b.next()) 
    {
      pinfo.num += block->size();
      pinfo.numBeziers += block->size();
      for (size_t i=0; i<block->size(); i++)
      {
        const BBox3fa bounds = block->at(i).bounds(); 
        bezierGeomBounds.extend(bounds);
        bezierCentBounds.extend(center2(bounds));
      }
    }
    pinfo.geomBounds.extend(bezierGeomBounds);
    pinfo.centBounds.extend(bezierCentBounds);

    if (pinfo.numTriangles != 0 && pinfo.numBeziers != 0)
      split.cost = blocks(pinfo.numTriangles)*triCost*safeHalfArea(triGeomBounds) + pinfo.numBeziers*bezierCost*safeHalfArea(bezierGeomBounds);
    else 
      split.cost = inf;
  }
    
  void ObjectTypePartitioning::Split::split(size_t threadIndex, PrimRefBlockAlloc<PrimRef>* alloc, TriRefList& prims, TriRefList& lprims, TriRefList& rprims)
  {
    lprims = prims;
  }

  void ObjectTypePartitioning::Split::split(size_t threadIndex, PrimRefBlockAlloc<Bezier1>* alloc, BezierRefList& prims, BezierRefList& lprims, BezierRefList& rprims)
  {
    rprims = prims;
  }
}
