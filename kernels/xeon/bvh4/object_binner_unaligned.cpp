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

#include "object_binner_unaligned.h"
#include "common/scene.h"

namespace embree
{
  const size_t ObjectSplitBinnerUnaligned::maxBins;

  extern Scene* g_scene; // FIXME: remove me

  static BBox3fa primSpaceBounds(const LinearSpace3fa& space, const PrimRef& prim)
  {
    TriangleMesh* mesh = (TriangleMesh*) g_scene->get(prim.geomID());
    TriangleMesh::Triangle tri = mesh->triangle(prim.primID());
    BBox3fa bounds = empty;
    bounds.extend(xfmPoint(space,mesh->vertex(tri.v[0])));
    bounds.extend(xfmPoint(space,mesh->vertex(tri.v[1])));
    bounds.extend(xfmPoint(space,mesh->vertex(tri.v[2])));
    return bounds;
  }

  ObjectSplitBinnerUnaligned::ObjectSplitBinnerUnaligned(const LinearSpace3fa& space, TriRefList& triangles, float triCost) 
    : space(space), triCost(triCost), bezierCost(0.0f)
  {
    add(triangles);
    setup_binning();
    bin(triangles);
    best();
  }

  ObjectSplitBinnerUnaligned::ObjectSplitBinnerUnaligned(const LinearSpace3fa& space, 
                                                                       TriRefList& triangles, float triCost, BezierRefList& beziers, float bezierCost) 
    : space(space), triCost(triCost), bezierCost(bezierCost)
  {
    add(triangles);
    add(beziers);
    setup_binning();
    bin(triangles);
    bin(beziers);
    best();
  }

  void ObjectSplitBinnerUnaligned::add(TriRefList& prims)
  {
    TriRefList::iterator i=prims;
    while (TriRefBlock* block = i.next()) {
      info(block->base(),block->size());
    }
  }

  void ObjectSplitBinnerUnaligned::add(BezierRefList& prims)
  {
    BezierRefList::iterator i=prims;
    while (BezierRefBlock* block = i.next()) {
      info(block->base(),block->size());
    }
  }

  void ObjectSplitBinnerUnaligned::setup_binning()
  {
    new (&mapping) Mapping(pinfo);
    for (size_t i=0; i<mapping.size(); i++) {
      triCounts[i] = bezierCounts[i] = 0;
      geomBounds[i][0] = geomBounds[i][1] = geomBounds[i][2] = empty;
    }
  }

  void ObjectSplitBinnerUnaligned::bin(TriRefList& prims)
  {
    TriRefList::iterator j=prims;
    while (TriRefBlock* block = j.next())
      bin(block->base(),block->size());
  }

  void ObjectSplitBinnerUnaligned::bin(BezierRefList& prims)
  {
    BezierRefList::iterator j=prims;
    while (BezierRefBlock* block = j.next())
      bin(block->base(),block->size());
  }

  void ObjectSplitBinnerUnaligned::info(const PrimRef* prims, size_t num)
  {
    BBox3fa geomBounds = empty;
    BBox3fa centBounds = empty;
    for (size_t i=0; i<num; i++)
    {
      const BBox3fa bounds = primSpaceBounds(space,prims[i]); 
      geomBounds.extend(bounds);
      centBounds.extend(center2(bounds));
    }
    pinfo.num += num;
    pinfo.numTriangles += num;
    pinfo.geomBounds.extend(geomBounds);
    pinfo.centBounds.extend(centBounds);
  }

  void ObjectSplitBinnerUnaligned::info(const Bezier1* prims, size_t num)
  {
    BBox3fa geomBounds = empty;
    BBox3fa centBounds = empty;
    for (size_t i=0; i<num; i++)
    {
      const BBox3fa bounds = prims[i].bounds(space); 
      geomBounds.extend(bounds);
      centBounds.extend(center2(bounds));
    }
    pinfo.num += num;
    pinfo.numBeziers += num;
    pinfo.geomBounds.extend(geomBounds);
    pinfo.centBounds.extend(centBounds);
  }

  void ObjectSplitBinnerUnaligned::bin(const PrimRef* prims, size_t num)
  {
    if (num == 0) return;
    
    size_t i; for (i=0; i<num-1; i+=2)
    {
      /*! map even and odd primitive to bin */
      const BBox3fa prim0 = primSpaceBounds(space,prims[i+0]); const Vec3ia bin0 = mapping.bin(prim0); const Vec3fa center0 = Vec3fa(center2(prim0));
      const BBox3fa prim1 = primSpaceBounds(space,prims[i+1]); const Vec3ia bin1 = mapping.bin(prim1); const Vec3fa center1 = Vec3fa(center2(prim1));
      
      /*! increase bounds for bins for even primitive */
      const int b00 = bin0.x; triCounts[b00][0]++; geomBounds[b00][0].extend(prim0); 
      const int b01 = bin0.y; triCounts[b01][1]++; geomBounds[b01][1].extend(prim0); 
      const int b02 = bin0.z; triCounts[b02][2]++; geomBounds[b02][2].extend(prim0); 
      
      /*! increase bounds of bins for odd primitive */
      const int b10 = bin1.x; triCounts[b10][0]++; geomBounds[b10][0].extend(prim1); 
      const int b11 = bin1.y; triCounts[b11][1]++; geomBounds[b11][1].extend(prim1); 
      const int b12 = bin1.z; triCounts[b12][2]++; geomBounds[b12][2].extend(prim1); 
    }
    
    /*! for uneven number of primitives */
    if (i < num)
    {
      /*! map primitive to bin */
      const BBox3fa prim0 = primSpaceBounds(space,prims[i]); const Vec3ia bin0 = mapping.bin(prim0); const Vec3fa center0 = Vec3fa(center2(prim0));
      
      /*! increase bounds of bins */
      const int b00 = bin0.x; triCounts[b00][0]++; geomBounds[b00][0].extend(prim0); 
      const int b01 = bin0.y; triCounts[b01][1]++; geomBounds[b01][1].extend(prim0); 
      const int b02 = bin0.z; triCounts[b02][2]++; geomBounds[b02][2].extend(prim0); 
    }
  }

  void ObjectSplitBinnerUnaligned::bin(const Bezier1* prims, size_t num)
  {
    for (size_t i=0; i<num; i++)
    {
      /*! map even and odd primitive to bin */
      const BBox3fa prim0 = prims[i+0].bounds(space); const Vec3ia bin0 = mapping.bin(prim0); const Vec3fa center0 = Vec3fa(center2(prim0));
      
      /*! increase bounds for bins for even primitive */
      const int b00 = bin0.x; bezierCounts[b00][0]++; geomBounds[b00][0].extend(prim0); 
      const int b01 = bin0.y; bezierCounts[b01][1]++; geomBounds[b01][1].extend(prim0); 
      const int b02 = bin0.z; bezierCounts[b02][2]++; geomBounds[b02][2].extend(prim0); 
    }
  }
  
  void ObjectSplitBinnerUnaligned::best()
  {
    Vec3fa rAreas [maxBins];      //!< area of bounds of primitives on the right
    Vec3ia rTriCounts[maxBins];      //!< blocks of primitives on the right
    Vec3ia rBezierCounts[maxBins];      //!< blocks of primitives on the right
    
    /* sweep from right to left and compute parallel prefix of merged bounds */
    assert(mapping.size() > 0);
    Vec3ia triCount = 0, bezierCount = 0; 
    BBox3fa bx = empty; BBox3fa by = empty; BBox3fa bz = empty;
    for (size_t i=mapping.size()-1; i>0; i--)
    {
      triCount += triCounts[i];
      bezierCount += bezierCounts[i];
      rTriCounts[i] = triCount;
      rBezierCounts[i] = bezierCount;
      bx = merge(bx,geomBounds[i][0]); rAreas[i][0] = halfArea(bx);
      by = merge(by,geomBounds[i][1]); rAreas[i][1] = halfArea(by);
      bz = merge(bz,geomBounds[i][2]); rAreas[i][2] = halfArea(bz);
    }
    
    /* sweep from left to right and compute SAH */
    Vec3ia ii = 1; Vec3fa bestSAH = pos_inf; Vec3ia bestPos = 0;
    triCount = 0; bezierCount = 0; bx = empty; by = empty; bz = empty;
    for (size_t i=1; i<mapping.size(); i++, ii+=1)
    {
      triCount += triCounts[i-1];
      bezierCount += bezierCounts[i-1];
      bx = merge(bx,geomBounds[i-1][0]); float Ax = halfArea(bx);
      by = merge(by,geomBounds[i-1][1]); float Ay = halfArea(by);
      bz = merge(bz,geomBounds[i-1][2]); float Az = halfArea(bz);
      const Vec3fa lArea = Vec3fa(Ax,Ay,Az);
      const Vec3fa rArea = rAreas[i];
      const Vec3fa triSAH    = lArea*Vec3fa(blocks(triCount)) + rArea*Vec3fa(blocks(rTriCounts[i]));
      const Vec3fa bezierSAH = lArea*Vec3fa(bezierCount     ) + rArea*Vec3fa(rBezierCounts[i]     );
      const Vec3fa sah = triCost*triSAH + bezierCost*bezierSAH;
      bestPos = select(lt_mask(sah,bestSAH),ii ,bestPos);
      bestSAH = select(lt_mask(sah,bestSAH),sah,bestSAH);
    }

    /* find best dimension */
    for (int i=0; i<3; i++) 
    {
      if (unlikely(pinfo.centBounds.lower[i] >= pinfo.centBounds.upper[i])) 
        continue;
      
      if (bestSAH[i] < split.cost && bestPos[i] != 0) {
        split.dim = i;
        split.pos = bestPos[i];
        split.cost = bestSAH[i];
      }
    }

    split.mapping = mapping;
    split.space = space;
  }

  void ObjectSplitBinnerUnaligned::Split::split(size_t thread, PrimRefBlockAlloc<PrimRef>* alloc, TriRefList& prims, TriRefList& lprims, TriRefList& rprims)
  {
    TriRefList::item* lblock = lprims.insert(alloc->malloc(thread));
    TriRefList::item* rblock = rprims.insert(alloc->malloc(thread));
    
    while (TriRefList::item* block = prims.take()) 
    {
      for (size_t i=0; i<block->size(); i++) 
      {
        const PrimRef& prim = block->at(i); 
                
        if (mapping.bin_unsafe(primSpaceBounds(space,prim))[dim] < pos) 
        {
          if (likely(lblock->insert(prim))) continue; 
          lblock = lprims.insert(alloc->malloc(thread));
          lblock->insert(prim);
        } 
        else 
        {
          if (likely(rblock->insert(prim))) continue;
          rblock = rprims.insert(alloc->malloc(thread));
          rblock->insert(prim);
        }
      }
      alloc->free(thread,block);
    }
  }

  void ObjectSplitBinnerUnaligned::Split::split(size_t thread, PrimRefBlockAlloc<Bezier1>* alloc, BezierRefList& prims, BezierRefList& lprims, BezierRefList& rprims)
  {
    BezierRefList::item* lblock = lprims.insert(alloc->malloc(thread));
    BezierRefList::item* rblock = rprims.insert(alloc->malloc(thread));
    
    while (BezierRefList::item* block = prims.take()) 
    {
      for (size_t i=0; i<block->size(); i++) 
      {
        const Bezier1& prim = block->at(i); 
                
        if (mapping.bin_unsafe(prim.bounds(space))[dim] < pos) 
        {
          if (likely(lblock->insert(prim))) continue; 
          lblock = lprims.insert(alloc->malloc(thread));
          lblock->insert(prim);
        } 
        else 
        {
          if (likely(rblock->insert(prim))) continue;
          rblock = rprims.insert(alloc->malloc(thread));
          rblock->insert(prim);
        }
      }
      alloc->free(thread,block);
    }
  }
}
