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

#include "object_binner.h"
#include "object_binner_unaligned.h"
#include "spatial_split_binner.h"
#include "spatial_center_split.h"
#include "object_type_partition.h"
#include "strand_split.h"

namespace embree
{
  /* BVH builder. The builder is multi-threaded and implements 3
   * different build strategies: 1) Small tasks are finished in a
   * single thread (BuildTask) 2) Medium sized tasks are split into
   * two tasks using a single thread (SplitTask) and 3) Large tasks are
   * split using multiple threads on one processor. */

  class BVH4Builder2 : public Builder
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
    
    class __aligned(16) GeneralSplit
    {
      enum Type { OBJECT_SPLIT, OBJECT_SPLIT_UNALIGNED, SPATIAL_SPLIT, SPATIAL_CENTER_SPLIT, TYPE_SPLIT, FALLBACK_SPLIT, STRAND_SPLIT };

    public:

      __forceinline GeneralSplit () {}

      __forceinline GeneralSplit (size_t N) 
        : type(FALLBACK_SPLIT), aligned(true), split_sah(inf), num(N) {}

      __forceinline GeneralSplit(const StrandSplit::Split& split) {
        type = STRAND_SPLIT; aligned = false;
        split_sah = split.splitSAH();
        new (data) StrandSplit::Split(split);
      }

      __forceinline GeneralSplit(const ObjectSplitBinner::Split& split, bool aligned_in) {
        type = OBJECT_SPLIT; aligned = aligned_in;
        split_sah = split.splitSAH();
        new (data) ObjectSplitBinner::Split(split);
      }

      __forceinline GeneralSplit(const ObjectSplitBinnerUnaligned::Split& split) {
        type = OBJECT_SPLIT_UNALIGNED; aligned = false;
        split_sah = split.splitSAH();
        new (data) ObjectSplitBinnerUnaligned::Split(split);
      }

      __forceinline GeneralSplit(const SpatialSplit::Split& split, bool aligned_in) {
        type = SPATIAL_SPLIT; aligned = aligned_in;
        split_sah = split.splitSAH();
        new (data) SpatialSplit::Split(split);
      }

      __forceinline GeneralSplit(const SpatialCenterSplit::Split& split, bool aligned_in) {
        type = SPATIAL_CENTER_SPLIT; aligned = aligned_in;
        split_sah = split.splitSAH();
        new (data) SpatialCenterSplit::Split(split);
      }

      __forceinline GeneralSplit(const ObjectTypePartitioning::Split& split) {
        type = TYPE_SPLIT; aligned = true;
        split_sah = split.splitSAH();
        new (data) ObjectTypePartitioning::Split(split);
      }

      __forceinline void split(size_t threadIndex, PrimRefBlockAlloc<PrimRef>* alloc, TriRefList& prims, TriRefList& lprims, TriRefList& rprims) 
      {
        switch (type) {
        case STRAND_SPLIT   :  break;  // FIXME: should clear prims array?       
        case OBJECT_SPLIT : ((ObjectSplitBinner::Split* )data)->split(threadIndex,alloc,prims,lprims,rprims); break;          
        case OBJECT_SPLIT_UNALIGNED : ((ObjectSplitBinnerUnaligned::Split* )data)->split(threadIndex,alloc,prims,lprims,rprims); break;          
        case SPATIAL_SPLIT: ((SpatialSplit::Split*)data)->split(threadIndex,alloc,prims,lprims,rprims); break;          
        case SPATIAL_CENTER_SPLIT: ((SpatialCenterSplit::Split*)data)->split(threadIndex,alloc,prims,lprims,rprims); break;          
        case TYPE_SPLIT   : ((ObjectTypePartitioning::Split*   )data)->split(threadIndex,alloc,prims,lprims,rprims); break;         
        default           : split_fallback(threadIndex,alloc,prims,lprims,rprims); break;
        }
      }

      void split(size_t threadIndex, PrimRefBlockAlloc<Bezier1>* alloc, BezierRefList& prims, BezierRefList& lprims, BezierRefList& rprims)
      {
        switch (type) {
        case STRAND_SPLIT   : ((StrandSplit::Split*   )data)->split(threadIndex,alloc,prims,lprims,rprims); break;         
        case OBJECT_SPLIT : ((ObjectSplitBinner::Split* )data)->split(threadIndex,alloc,prims,lprims,rprims); break;          
        case OBJECT_SPLIT_UNALIGNED : ((ObjectSplitBinnerUnaligned::Split* )data)->split(threadIndex,alloc,prims,lprims,rprims); break;          
        case SPATIAL_SPLIT: ((SpatialSplit::Split*)data)->split(threadIndex,alloc,prims,lprims,rprims); break;          
        case SPATIAL_CENTER_SPLIT: PING; throw std::runtime_error("implement me"); break; // FIXME: //((SpatialCenterSplit::Split*)data)->split(threadIndex,alloc,prims,lprims,rprims); break;          
        case TYPE_SPLIT   : ((ObjectTypePartitioning::Split*   )data)->split(threadIndex,alloc,prims,lprims,rprims); break;         
        default           : split_fallback(threadIndex,alloc,prims,lprims,rprims); break;
        }
      }

      __forceinline float splitSAH() const { return split_sah; }
      
      bool aligned;
    private:
      Type type;
      float split_sah;
      size_t num;
      __aligned(16) char data[256];
    };
  
    /*! stores all info to build a subtree */
    struct BuildTask
    {
      __forceinline BuildTask () {}

      __forceinline BuildTask (BVH4::NodeRef* dst, size_t depth, TriRefList& tris, BezierRefList& beziers, const PrimInfo& pinfo, const GeneralSplit& split, const NAABBox3fa& bounds)
        : dst(dst), depth(depth), tris(tris), beziers(beziers), pinfo(pinfo), split(split), nodeBounds(bounds) {}

      __forceinline BuildTask (BVH4::NodeRef* dst, size_t depth, TriRefList& tris, const PrimInfo& pinfo, const GeneralSplit& split, const NAABBox3fa& bounds)
        : dst(dst), depth(depth), tris(tris), pinfo(pinfo), split(split), nodeBounds(bounds) {}

    public:
      __forceinline friend bool operator< (const BuildTask& a, const BuildTask& b) {
        return area(a.pinfo.geomBounds) < area(b.pinfo.geomBounds);
      //  return a.pinfo.size() < b.pinfo.size();
        }
      
    public:
      BVH4::NodeRef* dst;
      size_t depth;
      TriRefList tris;
      BezierRefList beziers;
      PrimInfo pinfo;
      GeneralSplit split;
      NAABBox3fa nodeBounds;
    };

  public:

    static const PrimInfo computePrimInfo(TriRefList& triangles);
    static const PrimInfo computePrimInfo(BezierRefList& beziers);
    static const PrimInfo computePrimInfo(TriRefList& tris, BezierRefList& beziers);

    /*! calculate bounds for range of primitives */
    static const BBox3fa computeAlignedBounds(TriRefList& tris);
    static const BBox3fa computeAlignedBounds(BezierRefList& beziers);
    static const BBox3fa computeAlignedBounds(TriRefList& tris, BezierRefList& beziers);
    
    /*! calculate bounds for range of primitives */
    static const NAABBox3fa computeAlignedBounds(TriRefList& tris, const LinearSpace3fa& space);
    static const NAABBox3fa computeAlignedBounds(BezierRefList& beziers, const LinearSpace3fa& space);
    static const NAABBox3fa computeAlignedBounds(TriRefList& tris, BezierRefList& beziers, const LinearSpace3fa& space);

    static const LinearSpace3fa computeHairSpace(BezierRefList& prims);

    static void split_fallback(size_t threadIndex, PrimRefBlockAlloc<PrimRef>* alloc, TriRefList& prims, TriRefList& lprims, TriRefList& rprims);
    static void split_fallback(size_t threadIndex, PrimRefBlockAlloc<Bezier1>* alloc, BezierRefList& prims, BezierRefList& lprims, BezierRefList& rprims);

    /*! builder entry point */
    void build(size_t threadIndex, size_t threadCount);

    /*! Constructor. */
    BVH4Builder2 (BVH4* bvh, BuildSource* source, void* geometry, const size_t minLeafSize = 1, const size_t maxLeafSize = inf);

    /*! creates a leaf node */
    NodeRef leaf   (size_t threadIndex, size_t depth, TriRefList& prims   , const PrimInfo& pinfo);
    NodeRef leaf   (size_t threadIndex, size_t depth, BezierRefList& prims, const PrimInfo& pinfo);
    NodeRef leaf   (size_t threadIndex, size_t depth, TriRefList& tris, BezierRefList& beziers, const PrimInfo& pinfo);

    void computeSplit(PrimInfo& pinfo, TriRefList& tris, BezierRefList& beziers, GeneralSplit& split, const NAABBox3fa& nodeBounds);
    void computeSplit(PrimInfo& pinfo, TriRefList& tris, GeneralSplit& split);
    
    TASK_RUN_FUNCTION(BVH4Builder2,task_build_parallel);
    
    /*! execute single task and create subtasks */
    void processTrianglesAndBeziers(size_t threadIndex, BuildTask& task, BuildTask task_o[BVH4::N], size_t& N);
    void processTriangles(size_t threadIndex, BuildTask& task, BuildTask task_o[BVH4::N], size_t& N);

    /*! recursive build function for aligned and non-aligned bounds */
    void recurseTask(size_t threadIndex, BuildTask& task);

  private:
    BuildSource* source;      //!< build source interface
    void* geometry;           //!< input geometry
    

  public:
    size_t minLeafSize;                 //!< minimal size of a leaf
    size_t maxLeafSize;                 //!< maximal size of a leaf
    PrimRefBlockAlloc<Bezier1> allocBezierRefs;                 //!< Allocator for primitive blocks
    PrimRefBlockAlloc<PrimRef> allocTriRefs;                 //!< Allocator for primitive blocks

  public:
    BVH4* bvh;                      //!< Output BVH4

    MutexSys taskMutex;
    volatile atomic_t numActiveTasks;
    volatile atomic_t remainingSpatialSplits;
    std::vector<BuildTask> tasks;

  };
}
