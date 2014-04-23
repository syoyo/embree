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

#include "general_split.h"

namespace embree
{
  /* BVH builder. The builder is multi-threaded and implements 3
   * different build strategies: 1) Small tasks are finished in a
   * single thread (BuildTask) 2) Medium sized tasks are split into
   * two tasks using a single thread (SplitTask) and 3) Large tasks are
   * split using multiple threads on one processor. */

  class BVH4BuilderHair2 : public Builder
  {
    ALIGNED_CLASS;
  public:

    /*! Type shortcuts */
    typedef typename BVH4::BaseNode    Node;
    typedef typename BVH4::NodeRef NodeRef;

    typedef PrimRefBlockT<PrimRef> TriRefBlock;
    typedef atomic_set<TriRefBlock> TriRefList;

    typedef PrimRefBlockT<Bezier1> BezierRefBlock;
    typedef atomic_set<BezierRefBlock> BezierRefList;
      
    /*! stores all info to build a subtree */
    struct BuildTask
    {
      __forceinline BuildTask () {}

      __forceinline BuildTask (BVH4::NodeRef* dst, size_t depth, BezierRefList& beziers, const PrimInfo& pinfo, const GeneralSplit& split, const NAABBox3fa& bounds)
        : dst(dst), depth(depth), beziers(beziers), pinfo(pinfo), split(split), nodeBounds(bounds) {}

    public:
      __forceinline friend bool operator< (const BuildTask& a, const BuildTask& b) {
        return area(a.pinfo.geomBounds) < area(b.pinfo.geomBounds);
      //  return a.pinfo.size() < b.pinfo.size();
        }
      
    public:
      BVH4::NodeRef* dst;
      size_t depth;
      BezierRefList beziers;
      PrimInfo pinfo;
      GeneralSplit split;
      NAABBox3fa nodeBounds;
    };

  public:

    static const PrimInfo computePrimInfo(BezierRefList& beziers);

    /*! calculate bounds for range of primitives */
    static const BBox3fa computeAlignedBounds(BezierRefList& beziers);
    
    /*! calculate bounds for range of primitives */
    static const NAABBox3fa computeAlignedBounds(BezierRefList& beziers, const LinearSpace3fa& space);

    static const LinearSpace3fa computeHairSpace(BezierRefList& prims);

    /*! builder entry point */
    void build(size_t threadIndex, size_t threadCount);

    /*! Constructor. */
    BVH4BuilderHair2 (BVH4* bvh, Scene* scene);

    /*! creates a leaf node */
    NodeRef leaf   (size_t threadIndex, size_t depth, BezierRefList& prims, const PrimInfo& pinfo);

    void computeSplit(PrimInfo& pinfo, BezierRefList& beziers, GeneralSplit& split, const NAABBox3fa& nodeBounds);
    
    TASK_RUN_FUNCTION(BVH4BuilderHair2,task_build_parallel);
    
    /*! execute single task and create subtasks */
    void processBeziers(size_t threadIndex, BuildTask& task, BuildTask task_o[BVH4::N], size_t& N);

    /*! recursive build function for aligned and non-aligned bounds */
    void recurseTask(size_t threadIndex, BuildTask& task);

  private:
    //BuildSource* source;      //!< build source interface
    //void* geometry;           //!< input geometry
    Scene* scene;
    

  public:
    size_t minLeafSize;                 //!< minimal size of a leaf
    size_t maxLeafSize;                 //!< maximal size of a leaf
    PrimRefBlockAlloc<Bezier1> allocBezierRefs;                 //!< Allocator for primitive blocks

  public:
    BVH4* bvh;                      //!< Output BVH4

    MutexSys taskMutex;
    volatile atomic_t numActiveTasks;
    volatile atomic_t remainingSpatialSplits;
    std::vector<BuildTask> tasks;

  };
}
