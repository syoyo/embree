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

#include "spatial_split_binner.h"
#include "common/scene.h"

namespace embree
{
  extern Scene* g_scene; // FIXME: remove me

  void SpatialSplit::compute(TriRefList& tris, float triCost)
  {
    this->triCost = triCost;
    this->bezierCost = 0.0f;

    /* calculate geometry bounds */
    geomBounds = empty;
    for (TriRefList::block_iterator_unsafe i=tris; i; i++)
      geomBounds.extend(i->bounds());

    /* calculate binning function */
    ofs  = (ssef) geomBounds.lower;
    diag = (ssef) geomBounds.size();
    scale = select(diag != 0.0f,rcp(diag) * ssef(BINS * 0.99f),ssef(0.0f));

    /* initialize bins */
    for (size_t i=0; i<BINS; i++) {
      bounds[i][0] = bounds[i][1] = bounds[i][2] = bounds[i][3] = empty;
      numTriBegin[i] = numTriEnd[i] = 0;
      numBezierBegin[i] = numBezierEnd[i] = 0;
    }

    bin(tris);
    best();
  }

  void SpatialSplit::compute(BezierRefList& beziers, float bezierCost)
  {
    this->triCost = 0.0f;
    this->bezierCost = bezierCost;

    /* calculate geometry bounds */
    geomBounds = empty;
    for (BezierRefList::block_iterator_unsafe i=beziers; i; i++)
      geomBounds.extend(i->bounds());

    /* calculate binning function */
    ofs  = (ssef) geomBounds.lower;
    diag = (ssef) geomBounds.size();
    scale = select(diag != 0.0f,rcp(diag) * ssef(BINS * 0.99f),ssef(0.0f));

    /* initialize bins */
    for (size_t i=0; i<BINS; i++) {
      bounds[i][0] = bounds[i][1] = bounds[i][2] = bounds[i][3] = empty;
      numTriBegin[i] = numTriEnd[i] = 0;
      numBezierBegin[i] = numBezierEnd[i] = 0;
    }

    bin(beziers);
    best();
  }

  void SpatialSplit::compute(TriRefList& tris, float triCost, BezierRefList& beziers, float bezierCost)
  {
    this->triCost = triCost;
    this->bezierCost = bezierCost;

    /* calculate geometry bounds */
    geomBounds = empty;
    for (TriRefList::block_iterator_unsafe i=tris; i; i++)
      geomBounds.extend(i->bounds());

    for (BezierRefList::block_iterator_unsafe i=beziers; i; i++)
      geomBounds.extend(i->bounds());

    /* calculate binning function */
    ofs  = (ssef) geomBounds.lower;
    diag = (ssef) geomBounds.size();
    scale = select(diag != 0.0f,rcp(diag) * ssef(BINS * 0.99f),ssef(0.0f));

    /* initialize bins */
    for (size_t i=0; i<BINS; i++) {
      bounds[i][0] = bounds[i][1] = bounds[i][2] = bounds[i][3] = empty;
      numTriBegin[i] = numTriEnd[i] = 0;
      numBezierBegin[i] = numBezierEnd[i] = 0;
    }

    bin(tris);
    bin(beziers);
    best();
  }

  void SpatialSplit::bin(TriRefList& tris)
  {
    for (TriRefList::block_iterator_unsafe i=tris; i; i++)
    {
      PrimRef prim = *i;
      TriangleMesh* mesh = (TriangleMesh*) g_scene->get(i->geomID());
      TriangleMesh::Triangle tri = mesh->triangle(i->primID());

      const BBox3fa primBounds = prim.bounds();
      const ssei startbin = clamp(floori((ssef(primBounds.lower)-ofs)*scale),ssei(0),ssei(BINS-1));
      const ssei endbin   = clamp(floori((ssef(primBounds.upper)-ofs)*scale),ssei(0),ssei(BINS-1));

      const Vec3fa v0 = mesh->vertex(tri.v[0]);
      const Vec3fa v1 = mesh->vertex(tri.v[1]);
      const Vec3fa v2 = mesh->vertex(tri.v[2]);

      for (size_t dim=0; dim<3; dim++) 
      {
        size_t bin;
        PrimRef rest = prim;
        for (bin=startbin[dim]; bin<endbin[dim]; bin++) 
        {
          const float pos = float(bin+1)/scale[dim]+ofs[dim];
          
          PrimRef left,right;
          splitTriangle(prim,dim,pos,v0,v1,v2,left,right);

          bounds[bin][dim].extend(left.bounds());
          rest = right;
        }
        numTriBegin[startbin[dim]][dim]++;
        numTriEnd  [endbin  [dim]][dim]++;
        bounds[bin][dim].extend(rest.bounds());
      }
    }
  }

  void SpatialSplit::bin(BezierRefList& beziers)
  {
    for (BezierRefList::block_iterator_unsafe i=beziers; i; i++)
    {
      const Vec3fa v0 = i->p0;
      const Vec3fa v1 = i->p3;
      const ssei bin0 = clamp(floori((ssef(v0)-ofs)*scale),ssei(0),ssei(BINS-1));
      const ssei bin1 = clamp(floori((ssef(v1)-ofs)*scale),ssei(0),ssei(BINS-1));
      const ssei startbin = min(bin0,bin1);
      const ssei endbin   = max(bin0,bin1);

      for (size_t dim=0; dim<3; dim++) 
      {
        size_t bin;
        Bezier1 curve = *i;
        for (bin=startbin[dim]; bin<endbin[dim]; bin++) // FIXME: one can prevent many transformations in this loop here !!!
        {
          const float pos = float(bin+1)/scale[dim]+ofs[dim];
          //const Vec3fa plane(space.vx[dim],space.vy[dim],space.vz[dim],-pos);
          const Vec3fa plane(dim == 0,dim == 1,dim == 2,-pos); // FIXME: optimize
          Bezier1 bincurve,restcurve; 
          if (curve.split(plane,bincurve,restcurve)) {
            const BBox3fa cbounds = bincurve.bounds();
            bounds[bin][dim].extend(cbounds);
            curve = restcurve;
          }
        }
        numBezierBegin[startbin[dim]][dim]++;
        numBezierEnd  [endbin  [dim]][dim]++;
        const BBox3fa cbounds = curve.bounds();
        bounds[bin][dim].extend(cbounds);
      }
    }
  }

    void SpatialSplit::best()
    {
      /* sweep from right to left and compute parallel prefix of merged bounds */
      Vec3fa rAreas[BINS];
      Vec3ia rTriCounts[BINS];
      Vec3ia rBezierCounts[BINS];
      Vec3ia triCount = 0, bezierCount = 0; BBox3fa bx = empty; BBox3fa by = empty; BBox3fa bz = empty;
      for (size_t i=BINS-1; i>0; i--)
    {
      triCount += numTriEnd[i];
      bezierCount += numBezierEnd[i];
      rTriCounts[i] = triCount;
      rBezierCounts[i] = bezierCount;
      bx.extend(bounds[i][0]); rAreas[i][0] = halfArea(bx);
      by.extend(bounds[i][1]); rAreas[i][1] = halfArea(by);
      bz.extend(bounds[i][2]); rAreas[i][2] = halfArea(bz);
    }
    
    /* sweep from left to right and compute SAH */
    Vec3ia ii = 1; Vec3fa bestSAH = pos_inf; Vec3ia bestPos = 0;
    triCount = 0; bezierCount = 0; bx = empty; by = empty; bz = empty;
    for (size_t i=1; i<BINS; i++, ii+=1)
    {
      triCount += numTriBegin[i-1];
      bezierCount += numBezierBegin[i-1];
      bx.extend(bounds[i-1][0]); float Ax = halfArea(bx);
      by.extend(bounds[i-1][1]); float Ay = halfArea(by);
      bz.extend(bounds[i-1][2]); float Az = halfArea(bz);
      const Vec3fa lArea = Vec3fa(Ax,Ay,Az);
      const Vec3fa rArea = rAreas[i];
      const Vec3fa triSAH    = lArea*Vec3fa(blocks(triCount)) + rArea*Vec3fa(blocks(rTriCounts[i]));
      const Vec3fa bezierSAH = lArea*Vec3fa(bezierCount     ) + rArea*Vec3fa(rBezierCounts[i]     );
      const Vec3fa sah = triCost*triSAH + bezierCost*bezierSAH;
      bestPos  = select(lt_mask(sah,bestSAH),ii ,bestPos);
      bestSAH  = select(lt_mask(sah,bestSAH),sah,bestSAH);
    }
    
    /* find best dimension */
    split.ofs = ofs;
    split.scale = scale;
    split.cost = inf;
    split.dim = -1;
    split.pos = 0.0f;
    split.ipos = 0;
    split.numSpatialSplits = 0;
    
    float bestCost = inf;
    for (size_t dim=0; dim<3; dim++) 
    {
      /* ignore zero sized dimensions */
      if (unlikely(scale[dim] == 0.0f)) 
        continue;
      
      /* test if this is a better dimension */
      if (bestSAH[dim] < bestCost && bestPos[dim] != 0) {
        split.dim = dim;
        split.pos = bestPos[dim]/scale[dim]+ofs[dim];
        split.ipos = bestPos[dim];
        split.cost = bestSAH[dim];
        bestCost = bestSAH[dim];
      }
    }

    if (split.dim == -1)
      return;

    for (size_t i=0; i<split.ipos; i++) {
      split.numSpatialSplits += numTriBegin[i][split.dim]-numTriEnd[i][split.dim];
      split.numSpatialSplits += numBezierBegin[i][split.dim]-numBezierEnd[i][split.dim];
    }
    assert(((ssize_t)split.numSpatialSplits) >= 0);
    }
      
  void SpatialSplit::Split::split(size_t threadIndex, PrimRefBlockAlloc<PrimRef>* alloc, TriRefList& prims, 
				  TriRefList& lprims_o, PrimInfo& linfo, TriRefList& rprims_o, PrimInfo& rinfo) const
  {
    size_t lnum = 0, rnum = 0;
    BBox3fa lgeomBounds = empty, rgeomBounds = empty;
    BBox3fa lcentBounds = empty, rcentBounds = empty;

    TriRefList::item* lblock = lprims_o.insert(alloc->malloc(threadIndex));
    TriRefList::item* rblock = rprims_o.insert(alloc->malloc(threadIndex));
    
    /* sort each primitive to left, right, or left and right */
    while (atomic_set<TriRefBlock>::item* block = prims.take()) 
    {
      for (size_t i=0; i<block->size(); i++) 
      {
        const PrimRef& prim = block->at(i); 
        const BBox3fa bounds = prim.bounds();

        /* sort to the left side */
        if (bounds.lower[dim] <= pos && bounds.upper[dim] <= pos)
        {
	  lgeomBounds.extend(bounds);
	  lcentBounds.extend(center2(bounds));
	  lnum++;

          if (likely(lblock->insert(prim))) continue; 
          lblock = lprims_o.insert(alloc->malloc(threadIndex));
          lblock->insert(prim);
          continue;
        }

        /* sort to the right side */
        if (bounds.lower[dim] >= pos && bounds.upper[dim] >= pos)
        {
	  rgeomBounds.extend(bounds);
	  rcentBounds.extend(center2(bounds));
	  rnum++;

          if (likely(rblock->insert(prim))) continue;
          rblock = rprims_o.insert(alloc->malloc(threadIndex));
          rblock->insert(prim);
          continue;
        }

        /* split and sort to left and right */
        TriangleMesh* mesh = (TriangleMesh*) g_scene->get(prim.geomID());
        TriangleMesh::Triangle tri = mesh->triangle(prim.primID());
        const Vec3fa v0 = mesh->vertex(tri.v[0]);
        const Vec3fa v1 = mesh->vertex(tri.v[1]);
        const Vec3fa v2 = mesh->vertex(tri.v[2]);

        PrimRef left,right;
        splitTriangle(prim,dim,pos,v0,v1,v2,left,right);

	lgeomBounds.extend(bounds);
	lcentBounds.extend(center2(bounds));
	lnum++;

        if (!lblock->insert(left)) {
          lblock = lprims_o.insert(alloc->malloc(threadIndex));
          lblock->insert(left);
        }

	rgeomBounds.extend(bounds);
	rcentBounds.extend(center2(bounds));
	rnum++;
        
        if (!rblock->insert(right)) {
          rblock = rprims_o.insert(alloc->malloc(threadIndex));
          rblock->insert(right);
        }
      }
      alloc->free(threadIndex,block);
    }

    linfo.geomBounds.extend(lgeomBounds);
    linfo.centBounds.extend(lcentBounds);
    linfo.numTriangles += lnum;
    linfo.num += lnum;

    rinfo.geomBounds.extend(rgeomBounds);
    rinfo.centBounds.extend(rcentBounds);
    rinfo.numTriangles += rnum;
    rinfo.num += rnum;
  }

  void SpatialSplit::Split::split(size_t threadIndex, PrimRefBlockAlloc<Bezier1>* alloc, BezierRefList& prims, 
				  BezierRefList& lprims_o, PrimInfo& linfo, BezierRefList& rprims_o, PrimInfo& rinfo) const
  {
    size_t lnum = 0, rnum = 0;
    BBox3fa lgeomBounds = empty, rgeomBounds = empty;
    BBox3fa lcentBounds = empty, rcentBounds = empty;

    /* calculate splitting plane */
    //const Vec3fa plane(space.vx[dim],space.vy[dim],space.vz[dim],-pos);
    const Vec3fa plane(dim == 0,dim == 1,dim == 2,-pos); // FIXME: optimize
    
    BezierRefList::item* lblock = lprims_o.insert(alloc->malloc(threadIndex));
    BezierRefList::item* rblock = rprims_o.insert(alloc->malloc(threadIndex));
    
    /* sort each primitive to left, right, or left and right */
    while (BezierRefList::item* block = prims.take()) 
    {
      for (size_t i=0; i<block->size(); i++) 
      {
        const Bezier1& prim = block->at(i); 
	const BBox3fa bounds = prim.bounds();
        const float p0p = prim.p0[dim]-pos; //dot(prim.p0,plane)+plane.w;
        const float p3p = prim.p3[dim]-pos; //dot(prim.p3,plane)+plane.w;

        /* sort to the left side */
        if (p0p <= 0.0f && p3p <= 0.0f)
        {
	  lgeomBounds.extend(bounds);
	  lcentBounds.extend(center2(bounds));
	  lnum++;

          if (likely(lblock->insert(prim))) continue; 
          lblock = lprims_o.insert(alloc->malloc(threadIndex));
          lblock->insert(prim);
          continue;
        }

        /* sort to the right side */
        if (p0p >= 0.0f && p3p >= 0.0f)
        {
	  rgeomBounds.extend(bounds);
	  rcentBounds.extend(center2(bounds));
	  rnum++;

          if (likely(rblock->insert(prim))) continue;
          rblock = rprims_o.insert(alloc->malloc(threadIndex));
          rblock->insert(prim);
          continue;
        }

        /* split and sort to left and right */
        Bezier1 left,right;
        if (prim.split(plane,left,right)) 
        {
	  lgeomBounds.extend(bounds);
	  lcentBounds.extend(center2(bounds));
	  lnum++;

          if (!lblock->insert(left)) {
            lblock = lprims_o.insert(alloc->malloc(threadIndex));
            lblock->insert(left);
          }

	  rgeomBounds.extend(bounds);
	  rcentBounds.extend(center2(bounds));
	  rnum++;

          if (!rblock->insert(right)) {
            rblock = rprims_o.insert(alloc->malloc(threadIndex));
            rblock->insert(right);
          }
          continue;
        }

        /* insert to left side as fallback */
	lgeomBounds.extend(bounds);
	lcentBounds.extend(center2(bounds));
	lnum++;

        if (!lblock->insert(prim)) {
          lblock = lprims_o.insert(alloc->malloc(threadIndex));
          lblock->insert(prim);
        }
      }
      alloc->free(threadIndex,block);
    }

    linfo.geomBounds.extend(lgeomBounds);
    linfo.centBounds.extend(lcentBounds);
    linfo.numBeziers += lnum;
    linfo.num += lnum;

    rinfo.geomBounds.extend(rgeomBounds);
    rinfo.centBounds.extend(rcentBounds);
    rinfo.numBeziers += rnum;
    rinfo.num += rnum;
  }
}
