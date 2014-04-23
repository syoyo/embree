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

#include "spatial_center_split.h"
#include "common/scene.h"

namespace embree
{
  extern Scene* g_scene; // FIXME: remove me

  const float leftSplitPos = 0.51f;
  const float rightSplitPos = 0.49f;

  SpatialCenterSplit::SpatialCenterSplit(TriRefList& tris, float triCost)
    : triCost(triCost), mapping(tris)
  {
    lcounts = rcounts = 0;
    lgeomBounds[0] = lgeomBounds[1] = lgeomBounds[2] = empty;
    rgeomBounds[0] = rgeomBounds[1] = rgeomBounds[2] = empty;
    lcentBounds[0] = lcentBounds[1] = lcentBounds[2] = empty;
    rcentBounds[0] = rcentBounds[1] = rcentBounds[2] = empty;

    bin(tris);
    best();
  }

  __forceinline SpatialCenterSplit::Mapping::Mapping (TriRefList& tris)
  {
    num = 0;
    BBox3fa geomBounds = empty;
    for (TriRefList::block_iterator_unsafe i=tris; i; i++) {
      geomBounds.extend(i->bounds()); num++;
    }

    center = 0.5f*center2(geomBounds);
    bleft  = geomBounds.lower + leftSplitPos *geomBounds.size();
    bright = geomBounds.lower + rightSplitPos*geomBounds.size();
  }

  void SpatialCenterSplit::bin(TriRefList& tris)
  {
    for (TriRefList::block_iterator_unsafe i=tris; i; i++)
    {
      PrimRef primi = *i;
      TriangleMesh* mesh = (TriangleMesh*) g_scene->get(i->geomID());
      TriangleMesh::Triangle tri = mesh->triangle(i->primID());
      const Vec3fa v0 = mesh->vertex(tri.v[0]);
      const Vec3fa v1 = mesh->vertex(tri.v[1]);
      const Vec3fa v2 = mesh->vertex(tri.v[2]);
      
      /*! map even and odd primitive to bin */
      const BBox3fa prim = primi.bounds(); const Vec3fa center = Vec3fa(center2(prim));
      
      /*! calculate for each dimension if primitive goes to the left and right for spatial splits */
      const Vec3ba left  = mapping.left (prim);
      const Vec3ba right = mapping.right(prim);

      /*! Test spatial split in center of each dimension. */
      for (int maxDim=0; maxDim<3; maxDim++)
      {
        /*! Choose better side for primitives that can be put left or right. */
        if (left[maxDim] && right[maxDim]) {
          if (prim.upper[maxDim]-mapping.bright[maxDim] > mapping.bleft[maxDim]-prim.lower[maxDim]) {
            rgeomBounds[maxDim].extend(prim); rcentBounds[maxDim].extend(center); rcounts[maxDim]++;
          } else {
            lgeomBounds[maxDim].extend(prim); lcentBounds[maxDim].extend(center); lcounts[maxDim]++;
          }
        }
        /*! These definitely go to the left. */
        else if (left[maxDim]) {
          lgeomBounds[maxDim].extend(prim); lcentBounds[maxDim].extend(center); lcounts[maxDim]++;
        }
        /*! These definitely go to the right. */
        else if (right[maxDim]) {
          rgeomBounds[maxDim].extend(prim); rcentBounds[maxDim].extend(center); rcounts[maxDim]++;
        }
        /*! These have to get split and put to left and right. */
        else {
          //PrimRef left,right; geom->split(primi,maxDim,mapping.center[maxDim],left,right);
          PrimRef left,right; splitTriangle(primi,maxDim,mapping.center[maxDim],v0,v1,v2,left,right);
          lgeomBounds[maxDim].extend(left. bounds()); lcentBounds[maxDim].extend(center2(left. bounds())); lcounts[maxDim]++;
          rgeomBounds[maxDim].extend(right.bounds()); rcentBounds[maxDim].extend(center2(right.bounds())); rcounts[maxDim]++;
        }
      }
    }
  }

  void SpatialCenterSplit::best()
  {
    split.mapping = mapping;

    /* find best spatial split */
    for (int i=0; i<3; i++) 
    {
      //size_t numFailed = pinfo.numFailed + (size_t(lcounts[i]) == pinfo.num && size_t(rcounts[i]) == pinfo.num);
      //if (unlikely(numFailed > maxFailed)) continue;
      //float u = float(ulp)*max(abs(pinfo.geomBounds.lower[i]),abs(pinfo.geomBounds.upper[i]));
      //bool flat = pinfo.geomBounds.upper[i] - pinfo.geomBounds.lower[i] <= 128.0f*u;
      //if (unlikely(flat)) continue;
      float sah = halfArea(lgeomBounds[i])*blocks(lcounts[i]) + halfArea(rgeomBounds[i])*blocks(rcounts[i]);
      sah *= triCost;
      if (unlikely(sah >= split.spatialSAH)) continue;
      split.spatialSAH = sah;
      split.sdim = i;
      //split.numFailed = (int)numFailed;
      split.numSpatialSplits += lcounts[i] + rcounts[i] - mapping.num;
    }
  }
      
  void SpatialCenterSplit::Split::split(size_t threadIndex, PrimRefBlockAlloc<PrimRef>* alloc, TriRefList& prims, TriRefList& lprims_o, TriRefList& rprims_o) const
  {
    TriRefList::item* lblock = lprims_o.insert(alloc->malloc(threadIndex));
    TriRefList::item* rblock = rprims_o.insert(alloc->malloc(threadIndex));
    
    /* sort each primitive to left, right, or left and right */
    while (atomic_set<TriRefBlock>::item* block = prims.take()) 
    {
      for (size_t i=0; i<block->size(); i++) 
      {
        const PrimRef& prim = block->at(i); 
        const BBox3fa bounds = prim.bounds();
        const bool left  = mapping.left(bounds )[sdim];
        const bool right = mapping.right(bounds)[sdim];

        /*! Choose better side for primitives that can be put left or right. */
        if (left && right) {
          if (prim.upper[sdim]-mapping.bright[sdim] > mapping.bleft[sdim]-prim.lower[sdim]) {
            if (likely(rblock->insert(prim))) continue;
            rblock = rprims_o.insert(alloc->malloc(threadIndex));
            rblock->insert(prim);          
          } else {
            if (likely(lblock->insert(prim))) continue; 
            lblock = lprims_o.insert(alloc->malloc(threadIndex));
            lblock->insert(prim);
          }
        }
        /*! These definitely go to the left. */
        else if (left) {
            if (likely(lblock->insert(prim))) continue; 
            lblock = lprims_o.insert(alloc->malloc(threadIndex));
            lblock->insert(prim);
        }
        /*! These definitely go to the right. */
        else if (right) {
          if (likely(rblock->insert(prim))) continue;
          rblock = rprims_o.insert(alloc->malloc(threadIndex));
          rblock->insert(prim);          
        }
        /*! These have to get split and put to left and right. */
        else 
        {
          /* split and sort to left and right */
          TriangleMesh* mesh = (TriangleMesh*) g_scene->get(prim.geomID());
          TriangleMesh::Triangle tri = mesh->triangle(prim.primID());
          const Vec3fa v0 = mesh->vertex(tri.v[0]);
          const Vec3fa v1 = mesh->vertex(tri.v[1]);
          const Vec3fa v2 = mesh->vertex(tri.v[2]);

          PrimRef left,right;
          splitTriangle(prim,sdim,mapping.center[sdim],v0,v1,v2,left,right);

          if (!lblock->insert(left)) {
            lblock = lprims_o.insert(alloc->malloc(threadIndex));
            lblock->insert(left);
          }
        
          if (!rblock->insert(right)) {
            rblock = rprims_o.insert(alloc->malloc(threadIndex));
            rblock->insert(right);
          }
        }
      }
      alloc->free(threadIndex,block);
    }
  }
}
