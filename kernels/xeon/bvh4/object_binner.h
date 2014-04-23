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

#include "builders/primrefalloc.h"
#include "builders/primrefblock.h"
#include "geometry/primitive.h"
#include "geometry/bezier1.h"

namespace embree
{
  /*! stores bounding information for a set of primitives */
  class PrimInfo
  {
    static const size_t logBlockSize = 2;

    /*! Compute the number of blocks occupied in one dimension. */
    __forceinline static int  blocks(size_t a) { return (int)((a+((1LL << logBlockSize)-1)) >> logBlockSize); }

  public:
    __forceinline PrimInfo () 
      : num(0), numTriangles(0), numBeziers(0), geomBounds(empty), centBounds(empty) {}
    
    __forceinline PrimInfo (size_t numTriangles, size_t numBeziers, const BBox3fa& geomBounds, const BBox3fa& centBounds) 
      : num(numTriangles+numBeziers), numTriangles(numTriangles), numBeziers(numBeziers), geomBounds(geomBounds), centBounds(centBounds) {}

    /*! returns the number of primitives */
    __forceinline size_t size() const { 
      return num; 
    }

    __forceinline float triSAH(float triCost) const { 
      return halfArea(geomBounds)*triCost*blocks(numTriangles);
    }

    __forceinline float bezierSAH(float bezierCost) const { 
      return halfArea(geomBounds)*bezierCost*numBeziers; 
    }
    
    __forceinline float leafSAH(float triCost, float bezierCost) const { 
      return halfArea(geomBounds)*(triCost*blocks(numTriangles) + bezierCost*numBeziers); 
    }

  public:
    size_t num;          //!< number of primitives
    size_t numTriangles, numBeziers;
    BBox3fa geomBounds;   //!< geometry bounds of primitives
    BBox3fa centBounds;   //!< centroid bounds of primitives
  };
  
  class ObjectSplitBinner
  {
    /*! Maximal number of bins. */
    static const size_t maxBins = 32;
    static const size_t logBlockSize = 2;

    typedef PrimRefBlockT<PrimRef> TriRefBlock;
    typedef atomic_set<TriRefBlock> TriRefList;

    typedef PrimRefBlockT<Bezier1> BezierRefBlock;
    typedef atomic_set<BezierRefBlock> BezierRefList;

  public:
    
    /*! Compute the number of blocks occupied for each dimension. */
    __forceinline static Vec3ia blocks(const Vec3ia& a) { return (a+Vec3ia((1 << logBlockSize)-1)) >> logBlockSize; }
    
    /*! Compute the number of blocks occupied in one dimension. */
    __forceinline static int  blocks(size_t a) { return (int)((a+((1LL << logBlockSize)-1)) >> logBlockSize); }

    /*! mapping from bounding boxes to bins */
    class Mapping
    {
    public:
      
      /*! default constructor */
      __forceinline Mapping () {}
      
      /*! construct from primitive info */
      __forceinline Mapping (const PrimInfo& pinfo)
      {
        num   = min(maxBins,size_t(4.0f + 0.05f*pinfo.num));
        ofs   = pinfo.centBounds.lower;
        scale = rcp(pinfo.centBounds.upper - pinfo.centBounds.lower) * Vec3fa(float(num));
      }
      
      /*! returns number of bins */
      __forceinline size_t size() const { return num; }
      
      /*! Computes the bin numbers for each dimension for a box. */
      __forceinline Vec3ia bin_unsafe(const BBox3fa& box) const {
        return Vec3ia((center2(box) - ofs)*scale-Vec3fa(0.5f));
      }
      
      /*! Computes the bin numbers for each dimension for a box. */
      __forceinline Vec3ia bin(const BBox3fa& box) const {
        return clamp(bin_unsafe(box),Vec3ia(0),Vec3ia(int(num-1)));
      }
      
    private:
      size_t num;       //!< number of bins to use
      Vec3fa ofs;        //!< offset to compute bin
      Vec3fa scale;      //!< scaling factor to compute bin
    };
    
    /*! Stores information about an object split. */
    class Split
    {
    public:
      
      /*! create an invalid split by default */
      __forceinline Split () : dim(0), pos(0), cost(inf) {}
      
      /*! return SAH cost of performing the split */
      __forceinline float splitSAH() const { return cost; } 

      void split(size_t threadIndex, PrimRefBlockAlloc<PrimRef>* alloc, TriRefList& prims, TriRefList& lprims, TriRefList& rprims);

      void split(size_t threadIndex, PrimRefBlockAlloc<Bezier1>* alloc, BezierRefList& prims, BezierRefList& lprims, BezierRefList& rprims);

    public:
      Mapping mapping;    //!< Mapping to bins
      int dim;            //!< Best object split dimension
      int pos;            //!< Best object split position
      float cost;         //!< SAH cost of performing best object split
      LinearSpace3fa space;
    };
    
    /*! default constructor */
    ObjectSplitBinner (TriRefList& prims, float triCost);
    ObjectSplitBinner (BezierRefList& beziers, float bezierCost);
    ObjectSplitBinner (TriRefList& prims, float triCost, BezierRefList& beziers, float bezierCost);
  
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
    
    /* initialize binning counter and bounds */
    Vec3ia  triCounts[maxBins];    
    Vec3ia  bezierCounts[maxBins];    
    BBox3fa geomBounds[maxBins][4]; 

    Split split;
    float triCost;
    float bezierCost;
  };
}
