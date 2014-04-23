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
  class __aligned(16) GeneralSplit
    {
      typedef PrimRefBlockT<PrimRef> TriRefBlock;
    typedef atomic_set<TriRefBlock> TriRefList;

    typedef PrimRefBlockT<Bezier1> BezierRefBlock;
    typedef atomic_set<BezierRefBlock> BezierRefList;

      enum Type { OBJECT_SPLIT, OBJECT_SPLIT_UNALIGNED, SPATIAL_SPLIT, SPATIAL_CENTER_SPLIT, TYPE_SPLIT, FALLBACK_SPLIT, STRAND_SPLIT };

    public:

      __forceinline GeneralSplit () {}

      __forceinline GeneralSplit (size_t N) 
        : type(FALLBACK_SPLIT), aligned(true), split_sah(inf), num(N) {}

      /*__forceinline GeneralSplit(const StrandSplit::Split& split) {
        type = STRAND_SPLIT; aligned = false;
        split_sah = split.splitSAH();
        new (data) StrandSplit::Split(split);
	}*/

      __forceinline GeneralSplit(const ObjectSplitBinner::Split& split) {
        type = OBJECT_SPLIT; aligned = true;
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

      /*__forceinline GeneralSplit(const SpatialCenterSplit::Split& split, bool aligned_in) {
        type = SPATIAL_CENTER_SPLIT; aligned = aligned_in;
        split_sah = split.splitSAH();
        new (data) SpatialCenterSplit::Split(split);
      }

      __forceinline GeneralSplit(const ObjectTypePartitioning::Split& split) {
        type = TYPE_SPLIT; aligned = true;
        split_sah = split.splitSAH();
        new (data) ObjectTypePartitioning::Split(split);
	}*/

      __forceinline void split(size_t threadIndex, PrimRefBlockAlloc<PrimRef>* alloc, TriRefList& prims, TriRefList& lprims, PrimInfo& linfo, TriRefList& rprims, PrimInfo& rinfo) 
      {
        switch (type) {
	  //case STRAND_SPLIT   :  break;  // FIXME: should clear prims array?       
        case OBJECT_SPLIT : ((ObjectSplitBinner::Split* )data)->split(threadIndex,alloc,prims,lprims,linfo,rprims,rinfo); break;          
        case OBJECT_SPLIT_UNALIGNED : ((ObjectSplitBinnerUnaligned::Split* )data)->split(threadIndex,alloc,prims,lprims,linfo,rprims,rinfo); break;          
        case SPATIAL_SPLIT: ((SpatialSplit::Split*)data)->split(threadIndex,alloc,prims,lprims,linfo,rprims,rinfo); break;          
	  //case SPATIAL_CENTER_SPLIT: ((SpatialCenterSplit::Split*)data)->split(threadIndex,alloc,prims,lprims,rprims); break;          
	  //case TYPE_SPLIT   : ((ObjectTypePartitioning::Split*   )data)->split(threadIndex,alloc,prims,lprims,rprims); break;         
        default           : ObjectSplitBinner::split_fallback(threadIndex,alloc,prims,lprims,linfo,rprims,rinfo); break;
        }
      }

      void split(size_t threadIndex, PrimRefBlockAlloc<Bezier1>* alloc, BezierRefList& prims, BezierRefList& lprims, PrimInfo& linfo, BezierRefList& rprims, PrimInfo& rinfo)
      {
        switch (type) {
	  //case STRAND_SPLIT   : ((StrandSplit::Split*   )data)->split(threadIndex,alloc,prims,lprims,rprims); break;         
        case OBJECT_SPLIT : ((ObjectSplitBinner::Split* )data)->split(threadIndex,alloc,prims,lprims,linfo,rprims,rinfo); break;          
        case OBJECT_SPLIT_UNALIGNED : ((ObjectSplitBinnerUnaligned::Split* )data)->split(threadIndex,alloc,prims,lprims,linfo,rprims,rinfo); break;          
        case SPATIAL_SPLIT: ((SpatialSplit::Split*)data)->split(threadIndex,alloc,prims,lprims,linfo,rprims,rinfo); break;          
	  //case SPATIAL_CENTER_SPLIT: PING; throw std::runtime_error("implement me"); break; // FIXME: //((SpatialCenterSplit::Split*)data)->split(threadIndex,alloc,prims,lprims,rprims); break;          
	  //case TYPE_SPLIT   : ((ObjectTypePartitioning::Split*   )data)->split(threadIndex,alloc,prims,lprims,rprims); break;         
        default           : ObjectSplitBinner::split_fallback(threadIndex,alloc,prims,lprims,linfo,rprims,rinfo); break;
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
}
