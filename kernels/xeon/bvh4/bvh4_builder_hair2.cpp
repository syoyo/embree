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
#include "bvh4_builder_hair2.h"
#include "bvh4_statistics.h"
#include "common/scene_bezier_curves.h"
#include "geometry/bezier1.h"

namespace embree
{
  BVH4BuilderHair2::BVH4BuilderHair2 (BVH4* bvh, Scene* scene)
    : scene(scene), minLeafSize(1), maxLeafSize(inf), bvh(bvh)
  {
    size_t maxLeafTris    = BVH4::maxLeafBlocks*bvh->primTys[0]->blockSize;
    size_t maxLeafBeziers = BVH4::maxLeafBlocks*bvh->primTys[1]->blockSize;
    if (maxLeafBeziers < this->maxLeafSize) this->maxLeafSize = maxLeafBeziers; // FIXME: keep separate for tris and beziers
  }

  const PrimInfo BVH4BuilderHair2::computePrimInfo(BezierRefList& beziers)
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

  const BBox3fa BVH4BuilderHair2::computeAlignedBounds(BezierRefList& beziers)
  {
    BBox3fa bounds = empty;
    for (BezierRefList::block_iterator_unsafe i=beziers; i; i++)
      bounds.extend(i->bounds());
    return bounds;
  }

  const NAABBox3fa BVH4BuilderHair2::computeAlignedBounds(BezierRefList& beziers, const LinearSpace3fa& space)
  {
    BBox3fa bounds = empty;
    for (BezierRefList::block_iterator_unsafe i=beziers; i; i++)
      bounds.extend(i->bounds(space));
    return NAABBox3fa(space,bounds);
  }

  const NAABBox3fa BVH4BuilderHair2::computeHairSpace(BezierRefList& prims)
  {
    size_t N = BezierRefList::block_iterator_unsafe(prims).size();
    if (N == 0)
      return empty; // FIXME: can cause problems with compression

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
    return NAABBox3fa(bestSpace,bestBounds);
  }

  void BVH4BuilderHair2::build(size_t threadIndex, size_t threadCount) 
  {
    size_t numBeziers = 0;
    size_t numTriangles = 0;
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
    PrimInfo pinfo = computePrimInfo(beziers);
    GeneralSplit split; computeSplit(pinfo,beziers,split,pinfo.geomBounds);

    /* perform binning */
    bvh->numPrimitives = numPrimitives;
    bvh->bounds = geomBounds;
    //if (primTy.needVertices) bvh->numVertices = numVertices; // FIXME
    //else                     bvh->numVertices = 0;


#if 0
    BuildTask task(&bvh->root,1,beziers,pinfo,split,geomBounds);
    recurseTask(threadIndex,task);
    //for (int i=0; i<5; i++) BVH4Rotate::rotate(bvh,*task.dst); 
#else
    BuildTask task(&bvh->root,1,beziers,pinfo,split,geomBounds);
    numActiveTasks = 1;
    tasks.push_back(task);
    push_heap(tasks.begin(),tasks.end());
    TaskScheduler::executeTask(threadIndex,threadCount,_task_build_parallel,this,threadCount,"BVH4BuilderHair2::build_parallel");
    //TaskScheduler::executeTask(threadIndex,threadCount,_task_build_parallel,this,1,"BVH4BuilderHair2::build_parallel");
#endif
    //bvh->root = recurse(threadIndex,1,tris,beziers,pinfo,split);

    PRINT(remainingSpatialSplits);
    
    /* free all temporary blocks */
    Alloc::global.clear();

    if (g_verbose >= 2 || g_benchmark) 
      t1 = getSeconds();
    
    if (g_verbose >= 2) {
      std::cout << "[DONE]" << std::endl;
      std::cout << "  dt = " << 1000.0f*(t1-t0) << "ms, perf = " << 1E-6*double(numPrimitives)/(t1-t0) << " Mprim/s" << std::endl;
      std::cout << BVH4Statistics(bvh).str();
    }
  }

  BVH4::NodeRef BVH4BuilderHair2::leaf(size_t threadIndex, size_t depth, BezierRefList& prims, const PrimInfo& pinfo)
  {
    size_t N = BezierRefList::block_iterator_unsafe(prims).size();

    if (N > (size_t)BVH4::maxLeafBlocks) {
      //std::cout << "WARNING: Loosing " << N-BVH4::maxLeafBlocks << " primitives during build!" << std::endl;
      std::cout << "!" << std::flush;
      N = (size_t)BVH4::maxLeafBlocks;
    }

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

  void BVH4BuilderHair2::computeSplit(PrimInfo& pinfo, BezierRefList& beziers, GeneralSplit& split, const NAABBox3fa& nodeBounds) //,  bool isAligned)
  {
    float bestSAH = inf;
    const float bezierCost = BVH4::intCost;
    //const int travCostAligned = isAligned ? BVH4::travCostAligned : BVH4::travCostUnaligned;
    float travCostAligned = BVH4::travCostAligned;

    ObjectSplitBinner object_binning_aligned(beziers,bezierCost);
    float object_binning_aligned_sah = object_binning_aligned.split.splitSAH() + travCostAligned*halfArea(nodeBounds.bounds);
    bestSAH = min(bestSAH,object_binning_aligned_sah);

    /*bool enableSpatialSplits = false;
    //bool enableSpatialSplits = remainingSpatialSplits > 0;
    SpatialSplit spatial_binning_aligned(beziers,bezierCost);
    float spatial_binning_aligned_sah = spatial_binning_aligned.split.splitSAH() + travCostAligned*halfArea(nodeBounds.bounds);
    if (enableSpatialSplits) 
    bestSAH = min(bestSAH,spatial_binning_aligned_sah );*/
    
    const NAABBox3fa hairspace = computeHairSpace(beziers);
    
    ObjectSplitBinnerUnaligned object_binning_unaligned(hairspace.space,beziers,bezierCost);
    float object_binning_unaligned_sah = object_binning_unaligned.split.splitSAH() + BVH4::travCostUnaligned*halfArea(nodeBounds.bounds);
    bestSAH = min(bestSAH,object_binning_unaligned_sah);
    
    if (bestSAH == float(inf))
      new (&split) GeneralSplit(object_binning_aligned.pinfo.size());

    else if (bestSAH == object_binning_aligned_sah)
      new (&split) GeneralSplit(object_binning_aligned.split,true);

    /*else if (enableSpatialSplits && bestSAH == spatial_binning_aligned_sah) {
      new (&split) GeneralSplit(spatial_binning_aligned.split,true);
      atomic_add(&remainingSpatialSplits,-spatial_binning_aligned.split.numSpatialSplits);
      }*/
    else if (bestSAH == object_binning_unaligned_sah)
      new (&split) GeneralSplit(object_binning_unaligned.split);

    else
      throw std::runtime_error("internal error");
  }

  void BVH4BuilderHair2::processBeziers(size_t threadIndex, BuildTask& task, BuildTask task_o[BVH4::N], size_t& N)
  {
    /*! compute leaf and split cost */
    const float leafSAH  = task.pinfo.bezierSAH (BVH4::intCost);
    const float splitSAH = task.split.splitSAH() + BVH4::travCostAligned*halfArea(task.nodeBounds.bounds);

    /*! create a leaf node when threshold reached or SAH tells us to stop */
    if (task.pinfo.size() <= minLeafSize || task.depth > BVH4::maxBuildDepth || (task.pinfo.size() <= maxLeafSize && leafSAH <= splitSAH)) {
      *task.dst = leaf(threadIndex,task.depth,task.beziers,task.pinfo); N = 0; return;
    }
    
    /*! initialize child list */
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
        float dSAH = csplit[i].splitSAH()-cpinfo[i].bezierSAH(BVH4::intCost);
        if (cpinfo[i].size() <= minLeafSize) continue; 
        if (cpinfo[i].size() > maxLeafSize) dSAH = min(0.0f,dSAH); //< force split for large jobs
        if (dSAH <= bestSAH) { bestChild = i; bestSAH = dSAH; }
        //if (area(cpinfo[i].geomBounds) > bestSAH) { bestChild = i; bestSAH = area(cpinfo[i].geomBounds); }
      }
      if (bestChild == -1) break;
      
      /*! perform best found split and find new splits */
      aligned &= csplit[bestChild].aligned;
      PrimInfo linfo,rinfo;
      BezierRefList lbeziers,rbeziers; csplit[bestChild].split(threadIndex,&allocBezierRefs,cbeziers[bestChild],lbeziers,linfo,rbeziers,rinfo);
      GeneralSplit lsplit; computeSplit(linfo,lbeziers,lsplit,linfo.geomBounds); // FIXME: ,linfo.geomBounds not correct
      GeneralSplit rsplit; computeSplit(rinfo,rbeziers,rsplit,rinfo.geomBounds);
      cbeziers[bestChild  ] = lbeziers; cpinfo[bestChild] = linfo; csplit[bestChild  ] = lsplit;
      cbeziers[numChildren] = rbeziers; cpinfo[numChildren] = rinfo; csplit[numChildren] = rsplit;
      numChildren++;
      
    } while (numChildren < BVH4::N);
    
    /*! create an aligned node */
    if (aligned)
    {
      BVH4::UANode* node = bvh->allocUANode(threadIndex);
      for (size_t i=0; i<numChildren; i++) {
        node->set(i,cpinfo[i].geomBounds);
        new (&task_o[i]) BuildTask(&node->child(i),task.depth+1,cbeziers[i],cpinfo[i],csplit[i],cpinfo[i].geomBounds);
      }
      *task.dst = bvh->encodeNode(node);
      N = numChildren;
    } 

    /*! create an unaligned node */
    else
    {
      BVH4::UUNode* node = bvh->allocUUNode(threadIndex);
      for (size_t i=0; i<numChildren; i++) {
        const NAABBox3fa bounds = computeHairSpace(cbeziers[i]);
        node->set(i,bounds);
        new (&task_o[i]) BuildTask(&node->child(i),task.depth+1,cbeziers[i],cpinfo[i],csplit[i],bounds);
      }
      *task.dst = bvh->encodeNode(node);
      N = numChildren;
    }
  }

  void BVH4BuilderHair2::recurseTask(size_t threadIndex, BuildTask& task)
  {
    size_t numChildren;
    BuildTask tasks[BVH4::N];
    processBeziers(threadIndex,task,tasks,numChildren);
    for (size_t i=0; i<numChildren; i++) 
      recurseTask(threadIndex,tasks[i]);
  }

  void BVH4BuilderHair2::task_build_parallel(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, TaskScheduler::Event* event) 
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
        processBeziers(threadIndex,task,ctasks,numChildren);
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
  
  Builder* BVH4BuilderHair2_ (BVH4* accel, Scene* scene) {
    return new BVH4BuilderHair2(accel,scene);
  }
}
