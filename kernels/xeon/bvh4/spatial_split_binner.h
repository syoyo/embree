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

namespace embree
{
  struct SpatialSplit
  {
    /*! Maximal number of bins. */
    static const size_t BINS = 16;
    static const size_t logBlockSize = 2;

    public:

    typedef PrimRefBlockT<PrimRef> TriRefBlock;
    typedef atomic_set<TriRefBlock> TriRefList;

    typedef PrimRefBlockT<Bezier1> BezierRefBlock;
    typedef atomic_set<BezierRefBlock> BezierRefList;

    /*! Compute the number of blocks occupied for each dimension. */
    __forceinline static Vec3ia blocks(const Vec3ia& a) { return (a+Vec3ia((1 << logBlockSize)-1)) >> logBlockSize; }
    
    /*! Compute the number of blocks occupied in one dimension. */
    __forceinline static int  blocks(size_t a) { return (int)((a+((1LL << logBlockSize)-1)) >> logBlockSize); }

    SpatialSplit (TriRefList& tris, float triCost);
    SpatialSplit (BezierRefList& beziers, float bezierCost);
    SpatialSplit (TriRefList& tris, float triCost, BezierRefList& beziers, float bezierCost);
    
    class Split
    {
    public:
      
      /*! create an invalid split by default */
      __forceinline Split () : dim(0), pos(0), ipos(0), cost(inf), numSpatialSplits(0) {}
      
      /*! return SAH cost of performing the split */
      __forceinline float splitSAH() const { return cost; } 
      
      void split(size_t threadIndex, PrimRefBlockAlloc<PrimRef>* alloc, TriRefList& prims, TriRefList& lprims, TriRefList& rprims) const;
      
      void split(size_t threadIndex, PrimRefBlockAlloc<Bezier1>* alloc, BezierRefList& prims, BezierRefList& lprims, BezierRefList& rprims) const;
      
    public:
      //LinearSpace3fa space;
      float pos;
      int ipos;
      int dim;
      float cost;
      ssef ofs,scale;
      size_t numSpatialSplits;
    };
    
    void bin(TriRefList& tris);
    void bin(BezierRefList& beziers);
    void best();
    
    Split split;
    
    BBox3fa geomBounds;
    ssef ofs;
    ssef diag;
    ssef scale;
    
    BBox3fa bounds[BINS][4];
    Vec3ia    numTriBegin[BINS];
    Vec3ia    numTriEnd[BINS];
    Vec3ia    numBezierBegin[BINS];
    Vec3ia    numBezierEnd[BINS];
    
    float triCost;
    float bezierCost;
  };
}
