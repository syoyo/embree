// ======================================================================== //
// Copyright 2009-2015 Intel Corporation                                    //
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

#include "bvh4.h"

namespace embree
{
  namespace isa
  {
    class BVH4Refit : public Builder
    {
      ALIGNED_CLASS;
    public:
      
      /*! Type shortcuts */
      typedef BVH4::Node    Node;
      typedef BVH4::NodeRef NodeRef;
      
    public:
      
      void build(size_t threadIndex, size_t threadCount);

      void clear();

      /*! Constructor. */
      BVH4Refit (BVH4* bvh, Builder* builder, TriangleMesh* mesh, size_t mode);

      ~BVH4Refit();

      void refit_sequential(size_t threadIndex, size_t threadCount);

      virtual BBox3fa update(void* prim, size_t N, TriangleMesh* mesh) const = 0;
      
    private:
      size_t annotate_tree_sizes(NodeRef& ref);
      void calculate_refit_roots ();
      
      BBox3fa leaf_bounds(NodeRef& ref);
      BBox3fa node_bounds(NodeRef& ref);
      BBox3fa recurse_bottom(NodeRef& ref);
      BBox3fa recurse_top(NodeRef& ref);
      
    private:
      TriangleMesh* mesh;
      bool listMode;
      
    public:
      Builder* builder;
      BVH4* bvh;                      //!< BVH to refit
      std::vector<NodeRef*> roots;    //!< List of equal sized subtrees for bvh refit
    };

    template<typename Primitive>
      class BVH4RefitT : public BVH4Refit
    {
    public:
      BVH4RefitT (BVH4* bvh, Builder* builder, TriangleMesh* mesh, size_t mode)
        : BVH4Refit(bvh,builder,mesh,mode) {}

    private:      
      BBox3fa update(void* prim, size_t N, TriangleMesh* mesh) const 
      {
        BBox3fa bounds = empty;
        for (size_t i=0; i<N; i++)
          bounds.extend(((Primitive*)prim)[i].update(mesh));
        return bounds;
      }
    };
  }
}
