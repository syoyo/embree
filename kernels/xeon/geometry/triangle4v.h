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

#include "primitive.h"

namespace embree
{
  /*! Stores the vertices of 4 triangles in struct of array layout. */
  struct Triangle4v
  { 
    typedef sseb simdb;
    typedef ssef simdf;
    typedef ssei simdi;

  public:
    struct Type : public PrimitiveType 
    {
      Type ();
      size_t size(const char* This) const;
    };
    static Type type;

  public:

    /*! returns maximal number of stored triangles */
    static __forceinline size_t max_size() { return 4; }
    
     /*! returns required number of primitive blocks for N primitives */
    static __forceinline size_t blocks(size_t N) { return (N+max_size()-1)/max_size(); }
   
  public:

    /*! Default constructor. */
    __forceinline Triangle4v () {}

    /*! Construction from vertices and IDs. */
    __forceinline Triangle4v (const sse3f& v0, const sse3f& v1, const sse3f& v2, const ssei& geomIDs, const ssei& primIDs)
      : v0(v0), v1(v1), v2(v2), geomIDs(geomIDs), primIDs(primIDs) {}
    
    /*! Returns a mask that tells which triangles are valid. */
    __forceinline sseb valid() const { return geomIDs != ssei(-1); }

    /*! Returns true if the specified triangle is valid. */
    __forceinline bool valid(const size_t i) const { assert(i<4); return geomIDs[i] != -1; }

    /*! Returns the number of stored triangles. */
    __forceinline size_t size() const { return __bsf(~movemask(valid())); }

    /*! returns the geometry IDs */
    __forceinline ssei geomID() const { return geomIDs; }
    __forceinline int geomID(const size_t i) const { assert(i<4); return geomIDs[i]; }

    /*! returns the primitive IDs */
    __forceinline ssei primID() const { return primIDs; }
    __forceinline int  primID(const size_t i) const { assert(i<4); return primIDs[i]; }

    /*! calculate the bounds of the triangles */
    __forceinline BBox3fa bounds() const 
    {
      sse3f lower = min(v0,v1,v2);
      sse3f upper = max(v0,v1,v2);
      sseb mask = valid();
      lower.x = select(mask,lower.x,ssef(pos_inf));
      lower.y = select(mask,lower.y,ssef(pos_inf));
      lower.z = select(mask,lower.z,ssef(pos_inf));
      upper.x = select(mask,upper.x,ssef(neg_inf));
      upper.y = select(mask,upper.y,ssef(neg_inf));
      upper.z = select(mask,upper.z,ssef(neg_inf));
      return BBox3fa(Vec3fa(reduce_min(lower.x),reduce_min(lower.y),reduce_min(lower.z)),
                     Vec3fa(reduce_max(upper.x),reduce_max(upper.y),reduce_max(upper.z)));
    }
    
    /*! non temporal store */
    __forceinline static void store_nt(Triangle4v* dst, const Triangle4v& src)
    {
      store4f_nt(&dst->v0.x,src.v0.x);
      store4f_nt(&dst->v0.y,src.v0.y);
      store4f_nt(&dst->v0.z,src.v0.z);
      store4f_nt(&dst->v1.x,src.v1.x);
      store4f_nt(&dst->v1.y,src.v1.y);
      store4f_nt(&dst->v1.z,src.v1.z);
      store4f_nt(&dst->v2.x,src.v2.x);
      store4f_nt(&dst->v2.y,src.v2.y);
      store4f_nt(&dst->v2.z,src.v2.z);
      store4i_nt(&dst->geomIDs,src.geomIDs);
      store4i_nt(&dst->primIDs,src.primIDs);
    }

    /*! fill triangle from triangle list */
    __forceinline void fill(atomic_set<PrimRefBlock>::block_iterator_unsafe& prims, Scene* scene, const bool list)
    {
      ssei vgeomID = -1, vprimID = -1;
      sse3f v0 = zero, v1 = zero, v2 = zero;
      
      for (size_t i=0; i<4 && prims; i++, prims++)
      {
	const PrimRef& prim = *prims;
	const size_t geomID = prim.geomID();
        const size_t primID = prim.primID();
        const TriangleMesh* __restrict__ const mesh = scene->getTriangleMesh(geomID);
        const TriangleMesh::Triangle& tri = mesh->triangle(primID);
        const Vec3fa& p0 = mesh->vertex(tri.v[0]);
        const Vec3fa& p1 = mesh->vertex(tri.v[1]);
        const Vec3fa& p2 = mesh->vertex(tri.v[2]);
        vgeomID [i] = geomID;
        vprimID [i] = primID;
        v0.x[i] = p0.x; v0.y[i] = p0.y; v0.z[i] = p0.z;
        v1.x[i] = p1.x; v1.y[i] = p1.y; v1.z[i] = p1.z;
        v2.x[i] = p2.x; v2.y[i] = p2.y; v2.z[i] = p2.z;
      }
      Triangle4v::store_nt(this,Triangle4v(v0,v1,v2,vgeomID,vprimID));
    }

    /*! fill triangle from triangle list */
    __forceinline void fill(const PrimRef* prims, size_t& begin, size_t end, Scene* scene, const bool list)
    {
      ssei vgeomID = -1, vprimID = -1;
      sse3f v0 = zero, v1 = zero, v2 = zero;
      
      for (size_t i=0; i<4 && begin<end; i++, begin++)
      {
	const PrimRef& prim = prims[begin];
        const size_t geomID = prim.geomID();
        const size_t primID = prim.primID();
        const TriangleMesh* __restrict__ const mesh = scene->getTriangleMesh(geomID);
        const TriangleMesh::Triangle& tri = mesh->triangle(primID);
        const Vec3fa& p0 = mesh->vertex(tri.v[0]);
        const Vec3fa& p1 = mesh->vertex(tri.v[1]);
        const Vec3fa& p2 = mesh->vertex(tri.v[2]);
        vgeomID [i] = geomID;
        vprimID [i] = primID;
        v0.x[i] = p0.x; v0.y[i] = p0.y; v0.z[i] = p0.z;
        v1.x[i] = p1.x; v1.y[i] = p1.y; v1.z[i] = p1.z;
        v2.x[i] = p2.x; v2.y[i] = p2.y; v2.z[i] = p2.z;
      }
      Triangle4v::store_nt(this,Triangle4v(v0,v1,v2,vgeomID,vprimID));
    }

    /*! updates the primitive */
    __forceinline BBox3fa update(TriangleMesh* mesh)
    {
      BBox3fa bounds = empty;
      ssei vgeomID = -1, vprimID = -1;
      sse3f v0 = zero, v1 = zero, v2 = zero;
	
      for (size_t i=0; i<4; i++)
      {
        if (primID(i) == -1) break;
        const unsigned geomId = geomID(i);
        const unsigned primId = primID(i);
        const TriangleMesh::Triangle& tri = mesh->triangle(primId);
        const Vec3fa p0 = mesh->vertex(tri.v[0]);
        const Vec3fa p1 = mesh->vertex(tri.v[1]);
        const Vec3fa p2 = mesh->vertex(tri.v[2]);
        bounds.extend(merge(BBox3fa(p0),BBox3fa(p1),BBox3fa(p2)));
        vgeomID [i] = geomId;
        vprimID [i] = primId;
        v0.x[i] = p0.x; v0.y[i] = p0.y; v0.z[i] = p0.z;
        v1.x[i] = p1.x; v1.y[i] = p1.y; v1.z[i] = p1.z;
        v2.x[i] = p2.x; v2.y[i] = p2.y; v2.z[i] = p2.z;
      }
      new (this) Triangle4v(v0,v1,v2,vgeomID,vprimID);
      return bounds;
    }
   
  public:
    sse3f v0;       //!< 1st vertex of the triangles
    sse3f v1;       //!< 2nd vertex of the triangles
    sse3f v2;       //!< 3rd vertex of the triangle.
    ssei geomIDs;   //!< geometry ID
    ssei primIDs;   //!< primitive ID
  };
}
