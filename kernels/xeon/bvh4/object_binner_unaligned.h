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
  class ObjectSplitBinnerUnaligned
  {
    /*! Maximal number of bins. */
    static const size_t maxBins = 32;
    static const size_t logBlockSize = 2;

  public:
    
    typedef ObjectSplitBinner::Mapping Mapping;

    typedef PrimRefBlockT<PrimRef> TriRefBlock;
    typedef atomic_set<TriRefBlock> TriRefList;

    typedef PrimRefBlockT<Bezier1> BezierRefBlock;
    typedef atomic_set<BezierRefBlock> BezierRefList;

    /*! Compute the number of blocks occupied for each dimension. */
    __forceinline static Vec3ia blocks(const Vec3ia& a) { return (a+Vec3ia((1 << logBlockSize)-1)) >> logBlockSize; }
    
    /*! Compute the number of blocks occupied in one dimension. */
    __forceinline static int  blocks(size_t a) { return (int)((a+((1LL << logBlockSize)-1)) >> logBlockSize); }

    /*! Stores information about an object split. */
    class Split
    {
    public:
      
      /*! create an invalid split by default */
      __forceinline Split () : dim(0), pos(0), cost(inf) {}
      
      /*! return SAH cost of performing the split */
      __forceinline float splitSAH() const { return cost; } 

      void split(size_t threadIndex, PrimRefBlockAlloc<PrimRef>* alloc, TriRefList& prims, TriRefList& lprims, PrimInfo& linfo, TriRefList& rprims, PrimInfo& rinfo);

      void split(size_t threadIndex, PrimRefBlockAlloc<Bezier1>* alloc, BezierRefList& prims, BezierRefList& lprims, PrimInfo& linfo, BezierRefList& rprims, PrimInfo& rinfo);

    public:
      Mapping mapping;    //!< Mapping to bins
      int dim;            //!< Best object split dimension
      int pos;            //!< Best object split position
      float cost;         //!< SAH cost of performing best object split
      LinearSpace3fa space;
    };
    
    /*! default constructor */
    __forceinline ObjectSplitBinnerUnaligned () {};
    void compute(const LinearSpace3fa& space, TriRefList& prims, float triCost);
    void compute(const LinearSpace3fa& space, BezierRefList& beziers, float bezierCost);
    void compute(const LinearSpace3fa& space, TriRefList& prims, float triCost, BezierRefList& beziers, float bezierCost);
  
  private:

    void add(TriRefList& prims);
    void add(BezierRefList& prims);

    void setup_binning();
    
    void bin(TriRefList& prims);

    void bin(BezierRefList& prims);

    /*! calculate the best possible split */
    void best();

    void info(const PrimRef* prims, size_t num);

    void info(const Bezier1* prims, size_t num);

    /*! bin an array of primitives */
    void bin(const PrimRef* prim, size_t num);

    /*! bin an array of beziers */
    void bin(const Bezier1* prim, size_t num);
      
  public:
    PrimInfo pinfo;                //!< bounding information of geometry
    Mapping mapping;               //!< mapping from geometry to the bins
    LinearSpace3fa space;
    
    /* initialize binning counter and bounds */
    Vec3ia  triCounts[maxBins];    
    Vec3ia  bezierCounts[maxBins];    
    BBox3fa geomBounds[maxBins][4]; 

    Split split;
    float triCost;
    float bezierCost;
  };
}
