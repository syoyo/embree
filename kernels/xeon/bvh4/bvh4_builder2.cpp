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
// distributed under the License is dist ributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#include "bvh4.h"
#include "bvh4_builder2.h"
#include "bvh4_statistics.h"
#include "common/scene_bezier_curves.h"
#include "common/scene_triangle_mesh.h"
#include "geometry/bezier1.h"
#include "geometry/triangle4.h"
#include "bvh4_rotate.h"

#include "builders/splitter_fallback.h"

#include "common/scene_triangle_mesh.h"

namespace embree
{
  Scene* g_scene = NULL; // FIXME: remove me
  BVH4Builder2* g_builder2 = NULL;

  BBox3fa primSpaceBounds(const LinearSpace3fa& space, const PrimRef& prim)
  {
    TriangleMesh* mesh = (TriangleMesh*) g_scene->get(prim.geomID());
    TriangleMesh::Triangle tri = mesh->triangle(prim.primID());
    BBox3fa bounds = empty;
    bounds.extend(xfmPoint(space,mesh->vertex(tri.v[0])));
    bounds.extend(xfmPoint(space,mesh->vertex(tri.v[1])));
    bounds.extend(xfmPoint(space,mesh->vertex(tri.v[2])));
    return bounds;
  }

  BVH4Builder2::BVH4Builder2 (BVH4* bvh, BuildSource* source, void* geometry, const size_t minLeafSize, const size_t maxLeafSize)
    : source(source), geometry(geometry), minLeafSize(minLeafSize), maxLeafSize(maxLeafSize), bvh(bvh)
  {
    size_t maxLeafTris    = BVH4::maxLeafBlocks*bvh->primTys[0]->blockSize;
    size_t maxLeafBeziers = BVH4::maxLeafBlocks*bvh->primTys[1]->blockSize;
    if (maxLeafTris < this->maxLeafSize) this->maxLeafSize = maxLeafTris;
    //if (maxLeafBeziers < this->maxLeafSize) this->maxLeafSize = maxLeafBeziers; // FIXME: keep separate for tris and beziers
  }

  void BVH4Builder2::split_fallback(size_t threadIndex, PrimRefBlockAlloc<PrimRef>* alloc, TriRefList& prims, TriRefList& lprims, TriRefList& rprims)
  {  
    size_t num = 0;
    atomic_set<TriRefBlock>::item* lblock = lprims.insert(alloc->malloc(threadIndex));
    atomic_set<TriRefBlock>::item* rblock = rprims.insert(alloc->malloc(threadIndex));
    
    while (atomic_set<TriRefBlock>::item* block = prims.take()) 
    {
      for (size_t i=0; i<block->size(); i++) 
      {
        const PrimRef& prim = block->at(i); 
        
        if (num&1) 
        {
          num++;
          if (likely(lblock->insert(prim))) continue; 
          lblock = lprims.insert(alloc->malloc(threadIndex));
          lblock->insert(prim);
        } else {
          num++;
          if (likely(rblock->insert(prim))) continue;
          rblock = rprims.insert(alloc->malloc(threadIndex));
          rblock->insert(prim);
        }
      }
      alloc->free(threadIndex,block);
    }
  }

  void BVH4Builder2::split_fallback(size_t threadIndex, PrimRefBlockAlloc<Bezier1>* alloc, BezierRefList& prims, BezierRefList& lprims, BezierRefList& rprims)
  {  
    size_t num = 0;
    atomic_set<BezierRefBlock>::item* lblock = lprims.insert(alloc->malloc(threadIndex));
    atomic_set<BezierRefBlock>::item* rblock = rprims.insert(alloc->malloc(threadIndex));
    
    while (atomic_set<BezierRefBlock>::item* block = prims.take()) 
    {
      for (size_t i=0; i<block->size(); i++) 
      {
        const Bezier1& prim = block->at(i); 
        
        if (num&1) 
        {
          num++;
          if (likely(lblock->insert(prim))) continue; 
          lblock = lprims.insert(alloc->malloc(threadIndex));
          lblock->insert(prim);
        } else {
          num++;
          if (likely(rblock->insert(prim))) continue;
          rblock = rprims.insert(alloc->malloc(threadIndex));
          rblock->insert(prim);
        }
      }
      alloc->free(threadIndex,block);
    }
  }

  const PrimInfo BVH4Builder2::computePrimInfo(TriRefList& tris, BezierRefList& beziers)
  {
    PrimInfo pinfo;
    BBox3fa geomBounds = empty;
    BBox3fa centBounds = empty;
    TriRefList::iterator t=tris;
    while (TriRefBlock* block = t.next()) 
    {
      pinfo.num += block->size();
      pinfo.numTriangles += block->size();
      for (size_t i=0; i<block->size(); i++)
      {
        const BBox3fa bounds = block->at(i).bounds(); 
        geomBounds.extend(bounds);
        centBounds.extend(center2(bounds));
      }
    }

    BezierRefList::iterator b=beziers;
    while (BezierRefBlock* block = b.next()) 
    {
      pinfo.num += block->size();
      pinfo.numBeziers += block->size();
      for (size_t i=0; i<block->size(); i++)
      {
        const BBox3fa bounds = block->at(i).bounds(); 
        geomBounds.extend(bounds);
        centBounds.extend(center2(bounds));
      }
    }
    pinfo.geomBounds = geomBounds;
    pinfo.centBounds = centBounds;
    return pinfo;
  }

  const PrimInfo BVH4Builder2::computePrimInfo(BezierRefList& beziers)
  {
    PrimInfo pinfo;
    BBox3fa geomBounds = empty;
    BBox3fa centBounds = empty;
    BezierRefList::iterator b=beziers;
    while (BezierRefBlock* block = b.next()) 
    {
      pinfo.num += block->size();
      pinfo.numBeziers += block->size();
      for (size_t i=0; i<block->size(); i++)
      {
        const BBox3fa bounds = block->at(i).bounds(); 
        geomBounds.extend(bounds);
        centBounds.extend(center2(bounds));
      }
    }
    pinfo.geomBounds = geomBounds;
    pinfo.centBounds = centBounds;
    return pinfo;
  }

  const PrimInfo BVH4Builder2::computePrimInfo(TriRefList& tris)
  {
    PrimInfo pinfo;
    BBox3fa geomBounds = empty;
    BBox3fa centBounds = empty;
    TriRefList::iterator t=tris;
    while (TriRefBlock* block = t.next()) 
    {
      pinfo.num += block->size();
      pinfo.numTriangles += block->size();
      for (size_t i=0; i<block->size(); i++)
      {
        const BBox3fa bounds = block->at(i).bounds(); 
        geomBounds.extend(bounds);
        centBounds.extend(center2(bounds));
      }
    }
    pinfo.geomBounds = geomBounds;
    pinfo.centBounds = centBounds;
    return pinfo;
  }

  const BBox3fa BVH4Builder2::computeAlignedBounds(TriRefList& tris)
  {
    BBox3fa bounds = empty;
    for (TriRefList::block_iterator_unsafe i=tris; i; i++)
      bounds.extend(i->bounds());
    return bounds;
  }

  const BBox3fa BVH4Builder2::computeAlignedBounds(BezierRefList& beziers)
  {
    BBox3fa bounds = empty;
    for (BezierRefList::block_iterator_unsafe i=beziers; i; i++)
      bounds.extend(i->bounds());
    return bounds;
  }

  const BBox3fa BVH4Builder2::computeAlignedBounds(TriRefList& tris, BezierRefList& beziers) {
    return merge(computeAlignedBounds(tris),computeAlignedBounds(beziers));
  }

  const NAABBox3fa BVH4Builder2::computeAlignedBounds(TriRefList& tris, const LinearSpace3fa& space)
  {
    BBox3fa bounds = empty;
    for (TriRefList::block_iterator_unsafe i=tris; i; i++)
      bounds.extend(primSpaceBounds(space,*i));
    return NAABBox3fa(space,bounds);
  }

  const NAABBox3fa BVH4Builder2::computeAlignedBounds(BezierRefList& beziers, const LinearSpace3fa& space)
  {
    BBox3fa bounds = empty;
    for (BezierRefList::block_iterator_unsafe i=beziers; i; i++)
      bounds.extend(i->bounds(space));
    return NAABBox3fa(space,bounds);
  }

  const NAABBox3fa BVH4Builder2::computeAlignedBounds(TriRefList& tris, BezierRefList& beziers, const LinearSpace3fa& space) {
    return NAABBox3fa(space,merge(computeAlignedBounds(tris,space).bounds,computeAlignedBounds(beziers,space).bounds));
  }

  const LinearSpace3fa BVH4Builder2::computeHairSpace(BezierRefList& prims)
  {
    size_t N = BezierRefList::block_iterator_unsafe(prims).size();
    if (N == 0)
      return one; // FIXME: can cause problems with compression

    float bestArea = inf;
    LinearSpace3fa bestSpace = one;
    BBox3fa bestBounds = empty;

    size_t k=0;
    for (BezierRefList::block_iterator_unsafe i = prims; i; i++)
    {
      if ((k++) % ((N+3)/4)) continue;
      //size_t k = begin + rand() % (end-begin);
      const Vec3fa axis = normalize(i->p3 - i->p0);
      if (length(i->p3 - i->p0) < 1E-9) continue;
      const LinearSpace3fa space0 = frame(axis).transposed();
      //const LinearSpace3fa space0 = LinearSpace3fa::rotate(Vec3fa(0,0,1),2.0f*float(pi)*drand48())*frame(axis).transposed();
      const LinearSpace3fa space = clamp(space0);
      BBox3fa bounds = empty;
      float area = 0.0f;
      for (BezierRefList::block_iterator_unsafe j = prims; j; j++) {
        const BBox3fa cbounds = j->bounds(space);
        area += halfArea(cbounds);
        bounds.extend(cbounds);
      }

      if (area <= bestArea) {
        bestBounds = bounds;
        bestSpace = space;
        bestArea = area;
      }
    }
    //assert(bestArea != (float)inf); // FIXME: can get raised if all selected curves are points
/*#ifdef DEBUG
    if (bestArea == (float)inf)
      {
        std::cout << "WARNING: bestArea == (float)inf" << std::endl; 
      }
      #endif*/

    return bestSpace;
  }

  void BVH4Builder2::build(size_t threadIndex, size_t threadCount) 
  {
    size_t numBeziers = 0;
    size_t numTriangles = 0;
    Scene* scene = (Scene*) geometry;
    g_scene = scene; // FIXME: hack
    g_builder2 = this; // FIXME: hack
    for (size_t i=0; i<scene->size(); i++) 
    {
      Geometry* geom = scene->get(i);
      if (!geom->isEnabled()) continue;

      if (geom->type == BEZIER_CURVES) {
        BezierCurves* set = (BezierCurves*) geom;
        numBeziers  += set->numCurves;
      }

      if (geom->type == TRIANGLE_MESH) {
        TriangleMesh* set = (TriangleMesh*) geom;
        if (set->numTimeSteps != 1) continue;
        numTriangles += set->numTriangles;
      }
    }
    size_t numPrimitives = numBeziers + numTriangles;

    remainingSpatialSplits = 4.0f*numPrimitives; // FIXME: hardcoded constant
    bvh->init(numPrimitives+remainingSpatialSplits);
    if (numPrimitives == 0) return;

    if (g_verbose >= 2) 
      std::cout << "building " + bvh->name() + " with SAH builder ... " << std::flush;

    double t0 = 0.0, t1 = 0.0f;
    if (g_verbose >= 2 || g_benchmark)
      t0 = getSeconds();
    
    /* first generate primrefs */
    //size_t numTriangles = 0;
    size_t numVertices = 0;
    TriRefList tris;
    BezierRefList beziers;

    //size_t numTris = 0;
    //size_t numBeziers = 0;
    BBox3fa geomBounds = empty;
    BBox3fa centBounds = empty;

    for (size_t i=0; i<scene->size(); i++) 
    {
      Geometry* geom = scene->get(i);
      if (!geom->isEnabled()) continue;

      if (geom->type == TRIANGLE_MESH) 
      {
        TriangleMesh* set = (TriangleMesh*) geom;
        if (set->numTimeSteps != 1) continue;
        //numTriangles += set->numTriangles;
        //numVertices  += set->numVertices;
        for (size_t j=0; j<set->numTriangles; j++) {
          const BBox3fa bounds = set->bounds(j);
          geomBounds.extend(bounds);
          centBounds.extend(center2(bounds));
          const PrimRef prim = PrimRef(bounds,i,j);
           TriRefList::item* block = tris.head();
          if (block == NULL || !block->insert(prim)) {
            block = tris.insert(allocTriRefs.malloc(threadIndex));
            block->insert(prim);
          }
        }
      }

      if (geom->type == BEZIER_CURVES) 
      {
        BezierCurves* set = (BezierCurves*) geom;
        //numBeziers  += set->numCurves;
        //numVertices += set->numVertices;
        for (size_t j=0; j<set->numCurves; j++) {
          const int ofs = set->curve(j);
          const Vec3fa& p0 = set->vertex(ofs+0);
          const Vec3fa& p1 = set->vertex(ofs+1);
          const Vec3fa& p2 = set->vertex(ofs+2);
          const Vec3fa& p3 = set->vertex(ofs+3);
          const Bezier1 bezier(p0,p1,p2,p3,0,1,i,j);
          //bounds.extend(subdivideAndAdd(threadIndex,prims,bezier,enablePreSubdivision));
          const BBox3fa bounds = bezier.bounds();
          geomBounds.extend(bounds);
          centBounds.extend(center2(bounds));
          BezierRefList::item* block = beziers.head();
          if (block == NULL || !block->insert(bezier)) {
            block = beziers.insert(allocBezierRefs.malloc(threadIndex));
            block->insert(bezier);
          }
        }
      }
    }
    PrimInfo pinfo = computePrimInfo(tris,beziers);
    GeneralSplit split; computeSplit(pinfo,tris,beziers,split,pinfo.geomBounds);

    /* perform binning */
    bvh->numPrimitives = numPrimitives;
    bvh->bounds = geomBounds;
    //if (primTy.needVertices) bvh->numVertices = numVertices; // FIXME
    //else                     bvh->numVertices = 0;


#if 1
    BuildTask task(&bvh->root,1,tris,beziers,pinfo,split,geomBounds);
    recurseTask(threadIndex,task);
    //for (int i=0; i<5; i++) BVH4Rotate::rotate(bvh,*task.dst); 
#else
    BuildTask task(&bvh->root,1,tris,beziers,pinfo,split,geomBounds);
    numActiveTasks = 1;
    tasks.push_back(task);
    push_heap(tasks.begin(),tasks.end());
    TaskScheduler::executeTask(threadIndex,threadCount,_task_build_parallel,this,threadCount,"BVH4Builder2::build_parallel");
    //TaskScheduler::executeTask(threadIndex,threadCount,_task_build_parallel,this,1,"BVH4Builder2::build_parallel");
#endif
    //bvh->root = recurse(threadIndex,1,tris,beziers,pinfo,split);

    PRINT(remainingSpatialSplits);

    /* free all temporary blocks */
    Alloc::global.clear();

    if (g_verbose >= 2 || g_benchmark) 
      t1 = getSeconds();
    
    if (g_verbose >= 2) {
      std::cout << "[DONE]" << std::endl;
      std::cout << "  dt = " << 1000.0f*(t1-t0) << "ms, perf = " << 1E-6*double(source->size())/(t1-t0) << " Mprim/s" << std::endl;
      std::cout << BVH4Statistics(bvh).str();
    }
  }

  void BVH4Builder2::computeSplit(PrimInfo& pinfo, TriRefList& tris, BezierRefList& beziers, GeneralSplit& split, const NAABBox3fa& nodeBounds)
  {
    float bestSAH = inf;
    const float triCost = 1.0f; // FIXME:
    const float bezierCost = BVH4::intCost; // FIXME:

    ObjectTypePartitioning object_type(tris,triCost,beziers,bezierCost);
    float object_type_split_sah = object_type.split.splitSAH() + BVH4::travCostAligned*halfArea(nodeBounds.bounds);
    bestSAH = min(bestSAH,object_type_split_sah);

    StrandSplit strand_split;
    float strand_split_sah = inf;
    /*if (pinfo.numTriangles == 0 && pinfo.numBeziers > 0) {
      strand_split.find(beziers,bezierCost);
      strand_split_sah = strand_split.split.splitSAH() + BVH4::travCostUnaligned*halfArea(nodeBounds.bounds);
      bestSAH = min(bestSAH,strand_split_sah);
      }*/
    
    ObjectSplitBinner object_binning_aligned(tris,triCost,beziers,bezierCost);
    float object_binning_aligned_sah = object_binning_aligned.split.splitSAH() + BVH4::travCostAligned*halfArea(nodeBounds.bounds);;
    bestSAH = min(bestSAH,object_binning_aligned_sah);

    bool enableSpatialSplits = false;
    //bool enableSpatialSplits = remainingSpatialSplits > 0;
    SpatialSplit spatial_binning_aligned(tris,triCost,beziers,bezierCost);
    float spatial_binning_aligned_sah = spatial_binning_aligned.split.splitSAH() + BVH4::travCostAligned*halfArea(nodeBounds.bounds);;
    if (enableSpatialSplits) 
      bestSAH = min(bestSAH,spatial_binning_aligned_sah );
    
    const LinearSpace3fa hairspace = computeHairSpace(beziers);
    
    ObjectSplitBinnerUnaligned object_binning_unaligned;
    object_binning_unaligned.compute(hairspace,tris,triCost,beziers,bezierCost);
    float object_binning_unaligned_sah = object_binning_unaligned.split.splitSAH() + BVH4::travCostUnaligned*halfArea(nodeBounds.bounds);;
    bestSAH = min(bestSAH,object_binning_unaligned_sah);
    
    if (bestSAH == float(inf))
      new (&split) GeneralSplit(object_binning_aligned.pinfo.size());
    else if (bestSAH == strand_split_sah)
      new (&split) GeneralSplit(strand_split.split);
    else if (bestSAH == object_binning_aligned_sah)
      new (&split) GeneralSplit(object_binning_aligned.split,true);
    else if (enableSpatialSplits && bestSAH == spatial_binning_aligned_sah) {
      new (&split) GeneralSplit(spatial_binning_aligned.split,true);
      atomic_add(&remainingSpatialSplits,-spatial_binning_aligned.split.numSpatialSplits);
    }
    else if (bestSAH == object_binning_unaligned_sah)
      new (&split) GeneralSplit(object_binning_unaligned.split);
    else if (bestSAH == object_type_split_sah) {
      new (&split) GeneralSplit(object_type.split);
    }
    else
      throw std::runtime_error("internal error");
  }

  void BVH4Builder2::computeSplit(PrimInfo& pinfo, TriRefList& tris, GeneralSplit& split)
  {
    float bestSAH = inf;
    const float triCost = 1.0f; // FIXME:

    ObjectSplitBinner object_binning_aligned(tris,triCost);
    float object_binning_aligned_sah = object_binning_aligned.split.splitSAH() + BVH4::travCostAligned*halfArea(pinfo.geomBounds);;
    bestSAH = min(bestSAH,object_binning_aligned_sah);

    bool enableSpatialSplits = false;
    //bool enableSpatialSplits = remainingSpatialSplits > 0;
    //SpatialSplit spatial_binning_aligned(tris,triCost);
    SpatialCenterSplit spatial_binning_aligned(tris,triCost);
    float spatial_binning_aligned_sah = spatial_binning_aligned.split.splitSAH() + BVH4::travCostAligned*halfArea(pinfo.geomBounds);
    if (enableSpatialSplits) bestSAH = min(bestSAH,spatial_binning_aligned_sah);
    
    if (bestSAH == float(inf))
      new (&split) GeneralSplit(object_binning_aligned.pinfo.size());
    else if (bestSAH == object_binning_aligned_sah)
      new (&split) GeneralSplit(object_binning_aligned.split,true);
    else if (enableSpatialSplits && bestSAH == spatial_binning_aligned_sah) {
      new (&split) GeneralSplit(spatial_binning_aligned.split,true);
      atomic_add(&remainingSpatialSplits,-spatial_binning_aligned.split.numSpatialSplits);
    }
    else {
      throw std::runtime_error("internal error");
    }
  }

  BVH4::NodeRef BVH4Builder2::leaf(size_t threadIndex, size_t depth, TriRefList& prims, const PrimInfo& pinfo)
  {
    size_t N = bvh->primTys[0]->blocks(TriRefList::block_iterator_unsafe(prims).size());

    if (N > (size_t)BVH4::maxLeafBlocks) {
      //std::cout << "WARNING: Loosing " << N-BVH4::maxLeafBlocks << " primitives during build!" << std::endl;
      std::cout << "!" << std::flush;
      N = (size_t)BVH4::maxLeafBlocks;
    }
    Scene* scene = (Scene*) geometry;

    if (bvh->primTys[0] == &SceneTriangle4::type) 
    {
      char* leaf = bvh->allocPrimitiveBlocks(threadIndex,0,N);

      TriRefList::block_iterator_unsafe iter(prims);
      for (size_t j=0; j<N; j++) 
      {
        void* This = leaf+j*bvh->primTys[0]->bytes;
        ssei geomID = -1, primID = -1, mask = -1;
        sse3f v0 = zero, v1 = zero, v2 = zero;
    
        for (size_t i=0; i<4 && iter; i++, iter++)
        {
          const PrimRef& prim = *iter;
          const TriangleMesh* mesh = scene->getTriangleMesh(prim.geomID());
          const TriangleMesh::Triangle& tri = mesh->triangle(prim.primID());
          const Vec3fa& p0 = mesh->vertex(tri.v[0]);
          const Vec3fa& p1 = mesh->vertex(tri.v[1]);
          const Vec3fa& p2 = mesh->vertex(tri.v[2]);
          geomID [i] = prim.geomID();
          primID [i] = prim.primID();
          mask   [i] = mesh->mask;
          v0.x[i] = p0.x; v0.y[i] = p0.y; v0.z[i] = p0.z;
          v1.x[i] = p1.x; v1.y[i] = p1.y; v1.z[i] = p1.z;
          v2.x[i] = p2.x; v2.y[i] = p2.y; v2.z[i] = p2.z;
        }
        new (This) Triangle4(v0,v1,v2,geomID,primID,mask);
      }
      //assert(!iter);
    
      /* free all primitive blocks */
      while (TriRefList::item* block = prims.take())
        allocTriRefs.free(threadIndex,block);
        
      return bvh->encodeLeaf(leaf,N,0);
    } 
    else 
      throw std::runtime_error("unknown primitive type");
  }

  BVH4::NodeRef BVH4Builder2::leaf(size_t threadIndex, size_t depth, BezierRefList& prims, const PrimInfo& pinfo)
  {
    size_t N = BezierRefList::block_iterator_unsafe(prims).size();

    if (N > (size_t)BVH4::maxLeafBlocks) {
      //std::cout << "WARNING: Loosing " << N-BVH4::maxLeafBlocks << " primitives during build!" << std::endl;
      std::cout << "!" << std::flush;
      N = (size_t)BVH4::maxLeafBlocks;
    }
    Scene* scene = (Scene*) geometry;

    if (bvh->primTys[1] == &Bezier1Type::type) 
    { 
      Bezier1* leaf = (Bezier1*) bvh->allocPrimitiveBlocks(threadIndex,1,N);
      BezierRefList::block_iterator_unsafe iter(prims);
      for (size_t i=0; i<N; i++) { leaf[i] = *iter; iter++; }
      //assert(!iter);

      /* free all primitive blocks */
      while (BezierRefList::item* block = prims.take())
        allocBezierRefs.free(threadIndex,block);

      return bvh->encodeLeaf((char*)leaf,N,1);
    } 
    else if (bvh->primTys[1] == &SceneBezier1i::type) 
    {
      Bezier1i* leaf = (Bezier1i*) bvh->allocPrimitiveBlocks(threadIndex,1,N);
      BezierRefList::block_iterator_unsafe iter(prims);
      for (size_t i=0; i<N; i++) {
        const Bezier1& curve = *iter; iter++;
        const BezierCurves* in = (BezierCurves*) scene->get(curve.geomID);
        const Vec3fa& p0 = in->vertex(in->curve(curve.primID));
        leaf[i] = Bezier1i(&p0,curve.geomID,curve.primID,-1); // FIXME: support mask
      }

      /* free all primitive blocks */
      while (BezierRefList::item* block = prims.take())
        allocBezierRefs.free(threadIndex,block);

      return bvh->encodeLeaf((char*)leaf,N,1);
    }
    else 
      throw std::runtime_error("unknown primitive type");
  }

  BVH4::NodeRef BVH4Builder2::leaf(size_t threadIndex, size_t depth, TriRefList& tris, BezierRefList& beziers, const PrimInfo& pinfo)
  {
    size_t Ntris    = TriRefList   ::block_iterator_unsafe(tris   ).size(); // FIXME:
    size_t Nbeziers = BezierRefList::block_iterator_unsafe(beziers).size();
    if (Ntris == 0) 
      return leaf(threadIndex,depth,beziers,pinfo);
    else if (Nbeziers == 0)
      return leaf(threadIndex,depth,tris,pinfo);
    
    BVH4::UANode* node = bvh->allocUANode(threadIndex);
    node->set(0,computeAlignedBounds(tris   ),leaf(threadIndex,depth+1,tris   ,pinfo));
    node->set(1,computeAlignedBounds(beziers),leaf(threadIndex,depth+1,beziers,pinfo));
    return bvh->encodeNode(node);
  }

  void BVH4Builder2::processTrianglesAndBeziers(size_t threadIndex, BuildTask& task, BuildTask task_o[BVH4::N], size_t& N)
  {
    /*! compute leaf and split cost */
    const float leafSAH  = task.pinfo.leafSAH (bvh->primTys[0]->intCost,BVH4::intCost);
    const float splitSAH = task.split.splitSAH() + BVH4::travCostAligned*halfArea(task.nodeBounds.bounds);

    /*! create a leaf node when threshold reached or SAH tells us to stop */
    if (task.pinfo.size() <= minLeafSize || task.depth > BVH4::maxBuildDepth || (task.pinfo.size() <= maxLeafSize && leafSAH <= splitSAH)) {
      *task.dst = leaf(threadIndex,task.depth,task.tris,task.beziers,task.pinfo); N = 0; return;
    }
    
    /*! initialize child list */
    TriRefList    ctris   [BVH4::N]; ctris   [0] = task.tris;
    BezierRefList cbeziers[BVH4::N]; cbeziers[0] = task.beziers;
    GeneralSplit  csplit[BVH4::N];   csplit  [0] = task.split;
    PrimInfo      cpinfo[BVH4::N]; cpinfo[0] = task.pinfo;
    size_t numChildren = 1;
    bool aligned = true;

    /*! split until node is full or SAH tells us to stop */
    do {
      
      /*! find best child to split */
      float bestSAH = 0; 
      ssize_t bestChild = -1;
      for (size_t i=0; i<numChildren; i++) 
      {
        float dSAH = csplit[i].splitSAH()-cpinfo[i].leafSAH(1.0f,BVH4::intCost);
        if (cpinfo[i].size() <= minLeafSize) continue; 
        if (cpinfo[i].size() > maxLeafSize) dSAH = min(0.0f,dSAH); //< force split for large jobs
        if (dSAH <= bestSAH) { bestChild = i; bestSAH = dSAH; }
        //if (area(cpinfo[i].geomBounds) > bestSAH) { bestChild = i; bestSAH = area(cpinfo[i].geomBounds); }
      }
      if (bestChild == -1) break;
      
      /*! perform best found split and find new splits */
      aligned &= csplit[bestChild].aligned;
      TriRefList    ltris,   rtris;    csplit[bestChild].split(threadIndex,&allocTriRefs,   ctris[bestChild],   ltris,rtris);
      BezierRefList lbeziers,rbeziers; csplit[bestChild].split(threadIndex,&allocBezierRefs,cbeziers[bestChild],lbeziers,rbeziers);
      PrimInfo linfo = computePrimInfo(ltris,lbeziers);
      PrimInfo rinfo = computePrimInfo(rtris,rbeziers);
      GeneralSplit lsplit; computeSplit(linfo,ltris,lbeziers,lsplit,linfo.geomBounds); // FIXME: ,linfo.geomBounds not correct
      GeneralSplit rsplit; computeSplit(rinfo,rtris,rbeziers,rsplit,rinfo.geomBounds);
      ctris[bestChild  ] = ltris; cbeziers[bestChild  ] = lbeziers; cpinfo[bestChild] = linfo; csplit[bestChild  ] = lsplit;
      ctris[numChildren] = rtris; cbeziers[numChildren] = rbeziers; cpinfo[numChildren] = rinfo; csplit[numChildren] = rsplit;
      numChildren++;
      
    } while (numChildren < BVH4::N);
    
    /*! create an aligned node */
    if (aligned)
    {
      BVH4::UANode* node = bvh->allocUANode(threadIndex);
      for (size_t i=0; i<numChildren; i++) {
        const BBox3fa bounds = computeAlignedBounds(ctris[i],cbeziers[i]);
        node->set(i,bounds);
        new (&task_o[i]) BuildTask(&node->child(i),task.depth+1,ctris[i],cbeziers[i],cpinfo[i],csplit[i],bounds);
      }
      *task.dst = bvh->encodeNode(node);
      N = numChildren;
    } 

    /*! create an unaligned node */
    else
    {
      BVH4::UUNode* node = bvh->allocUUNode(threadIndex);
      for (size_t i=0; i<numChildren; i++) {
        const LinearSpace3fa hairspace = computeHairSpace(cbeziers[i]);
        const NAABBox3fa bounds = computeAlignedBounds(ctris[i],cbeziers[i],hairspace);
        node->set(i,bounds);
        new (&task_o[i]) BuildTask(&node->child(i),task.depth+1,ctris[i],cbeziers[i],cpinfo[i],csplit[i],bounds);
      }
      *task.dst = bvh->encodeNode(node);
      N = numChildren;
    }
  }

  void BVH4Builder2::processTriangles(size_t threadIndex, BuildTask& task, BuildTask task_o[BVH4::N], size_t& N)
  {
    /*! compute leaf and split cost */
    const float leafSAH  = task.pinfo.triSAH (bvh->primTys[0]->intCost);
    const float splitSAH = task.split.splitSAH() + BVH4::travCost*halfArea(task.nodeBounds.bounds);

    /*! create a leaf node when threshold reached or SAH tells us to stop */
    if (task.pinfo.size() <= minLeafSize || task.depth > BVH4::maxBuildDepth || (task.pinfo.size() <= maxLeafSize && leafSAH <= splitSAH)) {
      *task.dst = leaf(threadIndex,task.depth,task.tris,task.pinfo); N = 0; return;
    }
    
    /*! initialize child list */
    TriRefList    ctris   [BVH4::N]; ctris   [0] = task.tris;
    GeneralSplit  csplit[BVH4::N];   csplit  [0] = task.split;
    PrimInfo      cpinfo[BVH4::N]; cpinfo[0] = task.pinfo;
    size_t numChildren = 1;

    /*! split until node is full or SAH tells us to stop */
    do {
      
      /*! find best child to split */
      float bestSAH = 0; 
      ssize_t bestChild = -1;
      for (size_t i=0; i<numChildren; i++) 
      {
        float dSAH = csplit[i].splitSAH()-cpinfo[i].triSAH(bvh->primTys[0]->intCost);
        if (cpinfo[i].size() <= minLeafSize) continue; 
        if (cpinfo[i].size() > maxLeafSize) dSAH = min(0.0f,dSAH); //< force split for large jobs
        if (dSAH <= bestSAH) { bestChild = i; bestSAH = dSAH; }
      }
      if (bestChild == -1) break;
      
      /*! perform best found split and find new splits */
      TriRefList    ltris,   rtris;    csplit[bestChild].split(threadIndex,&allocTriRefs,   ctris[bestChild],   ltris,rtris);
      PrimInfo linfo = computePrimInfo(ltris);
      PrimInfo rinfo = computePrimInfo(rtris);
      GeneralSplit lsplit; computeSplit(linfo,ltris,lsplit);
      GeneralSplit rsplit; computeSplit(rinfo,rtris,rsplit);
      ctris[bestChild  ] = ltris; cpinfo[bestChild] = linfo; csplit[bestChild  ] = lsplit;
      ctris[numChildren] = rtris; cpinfo[numChildren] = rinfo; csplit[numChildren] = rsplit;
      numChildren++;
      
    } while (numChildren < BVH4::N);
    
    /*! create an aligned node */
    BVH4::UANode* node = bvh->allocUANode(threadIndex);
    for (size_t i=0; i<numChildren; i++) {
      node->set(i,cpinfo[i].geomBounds);
      new (&task_o[i]) BuildTask(&node->child(i),task.depth+1,ctris[i],cpinfo[i],csplit[i],cpinfo[i].geomBounds);
    }
    *task.dst = bvh->encodeNode(node);
    N = numChildren;
  }

  void BVH4Builder2::recurseTask(size_t threadIndex, BuildTask& task)
  {
    size_t numChildren;
    BuildTask tasks[BVH4::N];
    //processTriangles(threadIndex,task,tasks,numChildren);
    processTrianglesAndBeziers(threadIndex,task,tasks,numChildren);
    for (size_t i=0; i<numChildren; i++) 
      recurseTask(threadIndex,tasks[i]);
  }

  void BVH4Builder2::task_build_parallel(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, TaskScheduler::Event* event) 
  {
    while (numActiveTasks) 
    {
      taskMutex.lock();
      if (tasks.size() == 0) {
        taskMutex.unlock();
        continue;
      }

      /* take next task from heap */
      BuildTask task = tasks.front();
      pop_heap(tasks.begin(),tasks.end());
      tasks.pop_back();
      taskMutex.unlock();

      /* recursively finish task */
      if (task.pinfo.size() < 512) {
        atomic_add(&numActiveTasks,-1);
        recurseTask(threadIndex,task);
      }
      
      /* execute task and add child tasks */
      else 
      {
        size_t numChildren;
        BuildTask ctasks[BVH4::N];
        processTrianglesAndBeziers(threadIndex,task,ctasks,numChildren);
        taskMutex.lock();
        for (size_t i=0; i<numChildren; i++) {
          atomic_add(&numActiveTasks,+1);
          tasks.push_back(ctasks[i]);
          push_heap(tasks.begin(),tasks.end());
        }
        atomic_add(&numActiveTasks,-1);
        taskMutex.unlock();
      }
    }
  }
  
  Builder* BVH4Builder2ObjectSplit4 (void* accel, BuildSource* source, void* geometry, const size_t minLeafSize, const size_t maxLeafSize) {
    return new BVH4Builder2((BVH4*)accel,source,geometry,minLeafSize,maxLeafSize);
  }
}
