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
  struct SpatialCenterSplit
  {
  public:

    static const size_t logBlockSize = 2;

    typedef PrimRefBlockT<PrimRef> TriRefBlock;
    typedef atomic_set<TriRefBlock> TriRefList;

    /*! mapping from bounding boxes to bins and spatial mapping */
    class Mapping
    {
    public:
      
      /*! default constructor */
      __forceinline Mapping () {}
      
      /*! construct from centroid bounds and number of primitives */
      Mapping (TriRefList& tris);
/*      {
        center = 0.5f*center2(pinfo.geomBounds);
        bleft  = pinfo.geomBounds.lower + leftSplitPos *pinfo.geomBounds.size();
        bright = pinfo.geomBounds.lower + rightSplitPos*pinfo.geomBounds.size();
        }*/
      
      /*! returns number of bins */
      //__forceinline size_t size() const { return num; }

      /*! Computes if object is on the left of a spatial split. */
      __forceinline Vec3ba left(const BBox3fa& box) const {
        return le_mask(box.upper,bleft);
      }

      /*! Computes if object is on the right of a spatial split. */
      __forceinline Vec3ba right(const BBox3fa& box) const {
        return ge_mask(box.lower,bright);
      }

      /*! for spatial split */
    public:
      size_t num;
      Vec3fa center;     //!< center for spatial splits
      Vec3fa bleft;      //!< if upper bound of a primitive is smaller it goes to the left
      Vec3fa bright;     //!< if lower bound of a primitive is larger it goes to the left
    };
    
    /*! Stores information about an object split. */
    class Split
    {
    public:
      
      /*! create an invalid split by default */
      __forceinline Split () 
        : spatialSAH(inf), sdim(0), numFailed(0), numSpatialSplits(0) {}
      
      /*! return SAH cost of performing the split */
      __forceinline float splitSAH() const { return spatialSAH; } 
      
      void split(size_t threadIndex, PrimRefBlockAlloc<PrimRef>* alloc, TriRefList& prims, TriRefList& lprims, TriRefList& rprims) const;

      /*! best spatial split */
    public:
      float spatialSAH;    //!< SAH cost of performing best spatial split 
      int sdim;            //!< best spatial split dimension
      int numFailed;       //!< number of times a spatial split failed
      Mapping mapping;    //!< Mapping to bins
      int numSpatialSplits;
    };
    
    /*! default constructor */
    SpatialCenterSplit(TriRefList& tris, float triCost);

  private:
    
    /*! Compute the number of blocks occupied for each dimension. */
    __forceinline static Vec3ia blocks(const Vec3ia& a) { return (a+Vec3ia((1 << logBlockSize)-1)) >> logBlockSize; }
    
    /*! Compute the number of blocks occupied in one dimension. */
    __forceinline static int  blocks(size_t a) { return (int)((a+((1LL << logBlockSize)-1)) >> logBlockSize); }
    
    void bin(TriRefList& tris);
    void best();
    
  public:
    //PrimInfo pinfo;                 //!< bounding information of geometry
    Mapping mapping;                //!< mapping from geometry to the bins
    Split split;

    /* counter and bounds for spatial binning */
  public:
    Vec3ia lcounts, rcounts;        //< left and right count for spatial split
    BBox3fa lgeomBounds[4];        //< left geometry bounds for spatial split
    BBox3fa rgeomBounds[4];        //< right geometry bounds for spatial split
    BBox3fa lcentBounds[4];        //< left centroid bounds for spatial split
    BBox3fa rcentBounds[4];        //< right centroid bounds for spatial split

    float triCost;
  };
}
