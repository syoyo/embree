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

#include "embree2/rtcore.h"
#include "common/alloc.h"
#include "common/accel.h"
#include "common/scene.h"
#include "geometry/primitive.h"

//#define HIGH_BIT_ENCODING 1

namespace embree
{
  /*! Multi BVH with 4 children. Each node stores the bounding box of
   * it's 4 children as well as 4 child pointers. */
  class BVH4 : public Bounded
  {
    ALIGNED_CLASS;
  public:
    
    /*! forward declaration of node types */
    struct BaseNode;
    struct UANode;
    struct CANode;
    struct UUNode;
    struct CUNode;

    /*! branching width of the tree */
    static const size_t N = 4;

    /*! some pointer decoration masks */
    static const size_t deco_mask    = 0xF00000000000000FLL;  //!< these bits are used to decorate the pointer
    static const size_t barrier_mask = 0x0000000000000008LL;  //!< barrier bit to mark nodes during build
    static const size_t type_mask    = 0x0000000000000007LL;  //!< node and leaf type bits
    static const size_t alignment    = 0x0000000000000010LL;  //!< required pointer alignment

    /*! node types */
    static const size_t tyUA = 0;  //!< uncompressed axis aligned node
    static const size_t tyUU = 1;  //!< uncompressed unaligned node
    static const size_t tyCA = 2;  //!< compressed axis aligned node
    static const size_t tyCU = 3;  //!< compressed unaligned node
    static const size_t tyL0 = 4;  //!< first leaf node type
    static const size_t maxLeafTypes = 4;             //!< maximal number of leaf types

    /*! Empty node */
    static const size_t emptyNode = 0x0000000000000004LL;

    /*! Invalid node, used as marker in traversal */
    static const size_t invalidNode = 0x0FFFFFFFFFFFFFF4LL;

    /*! Maximal depth of the BVH. */
    static const size_t maxBuildDepth = 32;
    static const size_t maxBuildDepthLeaf = maxBuildDepth+16;
    static const size_t maxDepth = maxBuildDepthLeaf+maxBuildDepthLeaf+maxBuildDepth;
    
    /*! Maximal number of primitive blocks in a leaf. */
    //static const size_t maxLeafBlocks = 15;
    static const size_t maxLeafBlocks = 6;

    /*! Cost of one traversal step. */
    static const int travCost = 1; // FIXME: remove
    static const int travCostAligned = 1;
    static const int travCostUnaligned = 3;
    static const int intCost = 6;

    /*! Pointer that points to a node or a list of primitives */
    struct NodeRef
    {
      /*! Default constructor */
      __forceinline NodeRef () {}

      /*! Construction from integer */
      __forceinline NodeRef (size_t ptr) : ptr(ptr) { }

      /*! Cast to size_t */
      __forceinline operator size_t() const { return ptr; }

       /*! Prefetches the node this reference points to */
      __forceinline void prefetch() const {
	prefetchL1(((char*)(ptr & ~deco_mask))+0*64);
	prefetchL1(((char*)(ptr & ~deco_mask))+1*64);
	//prefetchL1(((char*)(ptr & ~deco_mask))+2*64);
	//prefetchL1(((char*)(ptr & ~deco_mask))+3*64);
      }

      /*! Sets the barrier bit. */
      __forceinline void setBarrier() { ptr |= barrier_mask; }
      
      /*! Clears the barrier bit. */
      __forceinline void clearBarrier() { ptr &= ~barrier_mask; }

      /*! Checks if this is an barrier. A barrier tells the top level tree rotations how deep to enter the tree. */
      __forceinline bool isBarrier() const { return ptr & barrier_mask; }

      /*! checks if this is a leaf */
      __forceinline bool isLeaf() const { return (ptr & type_mask) >= tyL0; }

      /*! checks if this is a node */
      __forceinline bool isNode() const { return (ptr & type_mask) < tyL0; }
      
      /*! checks for different node types */
      __forceinline bool isUANode() const { return (ptr & type_mask) == tyUA; }
      __forceinline bool isCANode() const { return (ptr & type_mask) == tyCA; }
      __forceinline bool isUUNode() const { return (ptr & type_mask) == tyUU; }
      __forceinline bool isCUNode() const { return (ptr & type_mask) == tyCU; }

      /*! returns node pointer */
      __forceinline       BaseNode* node()       { assert(isNode()); return (      BaseNode*)(ptr & ~deco_mask); }
      __forceinline const BaseNode* node() const { assert(isNode()); return (const BaseNode*)(ptr & ~deco_mask); }

      /*! returns node pointer */
      __forceinline       BaseNode* getNode()       { assert(isNode()); return (      BaseNode*)(ptr & ~deco_mask); }
      __forceinline const BaseNode* getNode() const { assert(isNode()); return (const BaseNode*)(ptr & ~deco_mask); }

      /*! returns uncompressed axis aligned node pointer */
      __forceinline       UANode* getUANode()       { assert(isUANode()); return (      UANode*)ptr; }
      __forceinline const UANode* getUANode() const { assert(isUANode()); return (const UANode*)ptr; }

      /*! returns compressed axis aligned node pointer */
      __forceinline       CANode* getCANode()       { assert(isCANode()); return (      CANode*)(ptr & ~deco_mask); }
      __forceinline const CANode* getCANode() const { assert(isCANode()); return (const CANode*)(ptr & ~deco_mask); }

      /*! returns uncompressed unaligned node pointer */
      __forceinline       UUNode* getUUNode()       { assert(isUUNode()); return (      UUNode*)(ptr & ~deco_mask); }
      __forceinline const UUNode* getUUNode() const { assert(isUUNode()); return (const UUNode*)(ptr & ~deco_mask); }

      /*! returns compressed unaligned node pointer */
      __forceinline       CUNode* getCUNode()       { assert(isCUNode()); return (      CUNode*)(ptr & ~deco_mask); }
      __forceinline const CUNode* getCUNode() const { assert(isCUNode()); return (const CUNode*)(ptr & ~deco_mask); }
      
      /*! returns leaf pointer */
      __forceinline char* getLeaf(size_t& num, size_t& type) const {
        assert(isLeaf());
        num  = ptr >> 60;
        type = (ptr & type_mask) - tyL0;
        return (char*)(ptr & ~deco_mask);
        //return (char*)(((ptr >> 4) << 8) >> 4); //~deco_mask);
      }

      /*! returns leaf pointer */
      /*__forceinline char* leaf(size_t& num) const {
        assert(isLeaf());
        num  = ptr >> 60;
        return (char*)(ptr & ~deco_mask);
        //return (char*)(((ptr >> 4) << 8) >> 4); //~deco_mask);
	}*/

    private:
      size_t ptr;
    };

    /*! Base node structure */
    struct BaseNode
    {
      /*! Clears the node. */
      __forceinline void clear() {
        for (size_t i=0; i<N; i++) children[i] = emptyNode;
      }

      /*! Sets bounding box and ID of child. */
      __forceinline void set(size_t i, const NodeRef& childID) {
        assert(i < N);
        children[i] = childID;
      }

      /*! Returns reference to specified child */
      __forceinline       NodeRef& child(size_t i)       { assert(i<N); return children[i]; }
      __forceinline const NodeRef& child(size_t i) const { assert(i<N); return children[i]; }

    public:
      NodeRef children[N];   //!< Pointer to the children (can be a node or leaf)
    };

    /*! Uncompressed node with axis aligned bounds */
    struct UANode : public BaseNode
    {
      enum { stride = sizeof(ssef) };

      /*! Clears the node. */
      __forceinline void clear() {
        lower_x = lower_y = lower_z = pos_inf; 
        upper_x = upper_y = upper_z = neg_inf;
        BaseNode::clear();
      }

      /*! Sets bounding box. */
      __forceinline void set(const size_t i, const BBox3fa& bounds) 
      {
        assert(i < N);
        lower_x[i] = bounds.lower.x; lower_y[i] = bounds.lower.y; lower_z[i] = bounds.lower.z;
        upper_x[i] = bounds.upper.x; upper_y[i] = bounds.upper.y; upper_z[i] = bounds.upper.z;
      }

      /*! Sets ID of child. */
      __forceinline void set(size_t i, const NodeRef& childID) {
        BaseNode::set(i,childID);
      }

      /*! Sets bounding box and ID of child. */
      __forceinline void set(const size_t i, const BBox3fa& bounds, const NodeRef& childID) {
        assert(i < N);
        lower_x[i] = bounds.lower.x; lower_y[i] = bounds.lower.y; lower_z[i] = bounds.lower.z;
        upper_x[i] = bounds.upper.x; upper_y[i] = bounds.upper.y; upper_z[i] = bounds.upper.z;
        BaseNode::set(i,childID);
      }

      /*! Returns bounds of node. */
      __forceinline BBox3fa bounds() const {
        const Vec3fa lower(reduce_min(lower_x),reduce_min(lower_y),reduce_min(lower_z));
        const Vec3fa upper(reduce_max(upper_x),reduce_max(upper_y),reduce_max(upper_z));
        return BBox3fa(lower,upper);
      }

      /*! Returns bounds of specified child. */
      __forceinline const BBox3fa bounds(const size_t i) const { 
        assert(i < N);
        const Vec3fa lower(lower_x[i],lower_y[i],lower_z[i]);
        const Vec3fa upper(upper_x[i],upper_y[i],upper_z[i]);
        return BBox3fa(lower,upper);
      }

      /*! Returns bounds of all children */
      __forceinline void bounds(BBox<ssef>& bounds0, BBox<ssef>& bounds1, BBox<ssef>& bounds2, BBox<ssef>& bounds3) const {
        transpose(lower_x,lower_y,lower_z,ssef(zero),bounds0.lower,bounds1.lower,bounds2.lower,bounds3.lower);
        transpose(upper_x,upper_y,upper_z,ssef(zero),bounds0.upper,bounds1.upper,bounds2.upper,bounds3.upper);
      }

      /*! returns 4 bounding boxes */
      __forceinline const BBoxSSE3f getBounds(const size_t nearX, const size_t nearY, const size_t nearZ) const 
      {
        const size_t farX  = nearX ^ sizeof(ssef), farY  = nearY ^ sizeof(ssef), farZ  = nearZ ^ sizeof(ssef);
        const ssef nearx = load4f((const char*)&lower_x+nearX);
        const ssef neary = load4f((const char*)&lower_y+nearY);
        const ssef nearz = load4f((const char*)&lower_z+nearZ);
        const ssef farx  = load4f((const char*)&lower_x+farX );
        const ssef fary  = load4f((const char*)&lower_y+farY );
        const ssef farz  = load4f((const char*)&lower_z+farZ );
        return BBoxSSE3f(Vec3<ssef>(nearx,neary,nearz),Vec3<ssef>(farx,fary,farz));
      }

      /*! Returns the extend of the bounds of the ith child */
      __forceinline Vec3fa extend(const size_t i) const {
        assert(i < N);
        return bounds(i).size();
      }

    public:
      ssef lower_x;           //!< X dimension of lower bounds of all 4 children.
      ssef upper_x;           //!< X dimension of upper bounds of all 4 children.
      ssef lower_y;           //!< Y dimension of lower bounds of all 4 children.
      ssef upper_y;           //!< Y dimension of upper bounds of all 4 children.
      ssef lower_z;           //!< Z dimension of lower bounds of all 4 children.
      ssef upper_z;           //!< Z dimension of upper bounds of all 4 children.
    };

    /*! Uncompressed node with unaligned bounds */
    struct UUNode : public BaseNode
    {
      /*! Clears the node. */
      __forceinline void clear() 
      {
	AffineSpace3fa empty = AffineSpace3fa::scale(Vec3fa(1E+19));
	naabb.l.vx = empty.l.vx;
	naabb.l.vy = empty.l.vy;
	naabb.l.vz = empty.l.vz;
	naabb.p    = empty.p;
        BaseNode::clear();
      }

      /*! Sets bounding box. */
      __forceinline void set(size_t i, const NAABBox3fa& b) 
      {
        assert(i < N);

        AffineSpace3fa space = b.space;
        space.p -= b.bounds.lower;
        space = AffineSpace3fa::scale(1.0f/max(Vec3fa(1E-19),b.bounds.upper-b.bounds.lower))*space;
        
        naabb.l.vx.x[i] = space.l.vx.x;
        naabb.l.vx.y[i] = space.l.vx.y;
        naabb.l.vx.z[i] = space.l.vx.z;

        naabb.l.vy.x[i] = space.l.vy.x;
        naabb.l.vy.y[i] = space.l.vy.y;
        naabb.l.vy.z[i] = space.l.vy.z;

        naabb.l.vz.x[i] = space.l.vz.x;
        naabb.l.vz.y[i] = space.l.vz.y;
        naabb.l.vz.z[i] = space.l.vz.z;

        naabb.p.x[i] = space.p.x;
        naabb.p.y[i] = space.p.y;
        naabb.p.z[i] = space.p.z;
      }

      /*! Sets ID of child. */
      __forceinline void set(size_t i, const NodeRef& childID) {
        BaseNode::set(i,childID);
      }

      /*! Sets bounding box and child. */
      __forceinline void set(size_t i, const NAABBox3fa& b, const NodeRef& childID) {
        set(i,b);
        set(i,childID);
      }

      /*! Returns the extend of the bounds of the ith child */
      __forceinline Vec3fa extend(size_t i) const {
        assert(i<N);
        const Vec3fa vx(naabb.l.vx.x[i],naabb.l.vx.y[i],naabb.l.vx.z[i]);
        const Vec3fa vy(naabb.l.vy.x[i],naabb.l.vy.y[i],naabb.l.vy.z[i]);
        const Vec3fa vz(naabb.l.vz.x[i],naabb.l.vz.y[i],naabb.l.vz.z[i]);
        const Vec3fa p (naabb.p   .x[i],naabb.p   .y[i],naabb.p   .z[i]);
        return rsqrt(vx*vx + vy*vy + vz*vz);
      }

    public:
      AffineSpaceSSE3f naabb;   //!< non-axis aligned bounding boxes (bounds are [0,1] in specified space)
    };

    /*! Compressed node with axis aligned bounds */
    struct CANode : public BaseNode
    {
      enum { stride = 4 };

      /*! Clears the node. */
      __forceinline void clear() 
      {
        offset = 0.0f; scale = 0.0f;
        lower_x[0] = lower_x[1] = lower_x[2] = lower_x[3] = 0;
        lower_y[0] = lower_y[1] = lower_y[2] = lower_y[3] = 0;
        lower_z[0] = lower_z[1] = lower_z[2] = lower_z[3] = 0;
        upper_x[0] = upper_x[1] = upper_x[2] = upper_x[3] = 0;
        upper_y[0] = upper_y[1] = upper_y[2] = upper_y[3] = 0;
        upper_z[0] = upper_z[1] = upper_z[2] = upper_z[3] = 0;
        BaseNode::clear();
      }

      /*! Sets non-axis aligned space of node and parent bounding box. */
      __forceinline void set(const NAABBox3fa& naabb) {
        offset = naabb.bounds.lower;
        scale  = naabb.bounds.size()/255.0f;
      }

      /*! Sets bounding box. */
      __forceinline void set(size_t i, const BBox3fa& bounds) 
      {
        assert(i < N);
        const Vec3fa lower = select(eq_mask(scale,Vec3fa(0.0f)),0.0f,(bounds.lower-Vec3fa(offset))/Vec3fa(scale));
        assert(lower.x >= 0.0f && lower.x <= 255.01f);
        assert(lower.y >= 0.0f && lower.y <= 255.01f);
        assert(lower.z >= 0.0f && lower.z <= 255.01f);
        lower_x[i] = (unsigned char) clamp(floorf(lower.x),0.0f,255.0f);
        lower_y[i] = (unsigned char) clamp(floorf(lower.y),0.0f,255.0f);
        lower_z[i] = (unsigned char) clamp(floorf(lower.z),0.0f,255.0f);
        const Vec3fa upper = select(eq_mask(scale,Vec3fa(0.0f)),0.0f,(bounds.upper-Vec3fa(offset))/Vec3fa(scale));
        assert(upper.x >= 0.0f && upper.x <= 255.01f);
        assert(upper.y >= 0.0f && upper.y <= 255.01f);
        assert(upper.z >= 0.0f && upper.z <= 255.01f);
        upper_x[i] = (unsigned char) clamp(ceilf(upper.x),0.0f,255.0f);
        upper_y[i] = (unsigned char) clamp(ceilf(upper.y),0.0f,255.0f);
        upper_z[i] = (unsigned char) clamp(ceilf(upper.z),0.0f,255.0f);
      }

      /*! Sets ID of child. */
      __forceinline void set(size_t i, const NodeRef& childID) {
        BaseNode::set(i,childID);
      }

      /*! returns ith bounding boxes */
      __forceinline BBox3fa getBounds(size_t i) const 
      {
        assert(i < N);
        const Vec3f lower((float)lower_x[i],lower_y[i],lower_z[i]);
        const Vec3f upper((float)upper_x[i],upper_y[i],upper_z[i]);
        return BBox3fa(Vec3fa(offset+scale*lower),Vec3fa(offset+scale*upper));
      }

      /*! returns 4 bounding boxes */
      __forceinline const BBoxSSE3f getBounds(const size_t nearX, const size_t nearY, const size_t nearZ) const 
      {
        const size_t farX  = nearX ^ 4, farY  = nearY ^ 4, farZ  = nearZ ^ 4;
        const ssef near_x = ssef(_mm_cvtepu8_epi32(*(ssei*)((char*)&this->lower_x+nearX)));
        const ssef near_y = ssef(_mm_cvtepu8_epi32(*(ssei*)((char*)&this->lower_y+nearY)));
        const ssef near_z = ssef(_mm_cvtepu8_epi32(*(ssei*)((char*)&this->lower_z+nearZ)));
        const ssef far_x  = ssef(_mm_cvtepu8_epi32(*(ssei*)((char*)&this->lower_x+farX )));
        const ssef far_y  = ssef(_mm_cvtepu8_epi32(*(ssei*)((char*)&this->lower_y+farY )));
        const ssef far_z  = ssef(_mm_cvtepu8_epi32(*(ssei*)((char*)&this->lower_z+farZ )));
        const Vec3<ssef> offset = *(Vec3fa*)&this->offset;
        const Vec3<ssef> scale  = *(Vec3fa*)&this->scale;
        return BBoxSSE3f(scale*Vec3<ssef>(near_x,near_y,near_z)+offset,
                          scale*Vec3<ssef>(far_x, far_y, far_z)+offset);
      }

      /*! Returns the extend of the bounds of the ith child */
      __forceinline Vec3fa extend(size_t i) const {
        assert(i < N);
        return getBounds(i).size();
      }

    public:
      Vec3f offset;               //!< offset to decompress bounds
      Vec3f scale;                //!< scale  to decompress bounds
      unsigned char lower_x[N];   //!< X dimension of lower bounds of all 4 children.
      unsigned char upper_x[N];   //!< X dimension of upper bounds of all 4 children.
      unsigned char lower_y[N];   //!< Y dimension of lower bounds of all 4 children.
      unsigned char upper_y[N];   //!< Y dimension of upper bounds of all 4 children.
      unsigned char lower_z[N];   //!< Z dimension of lower bounds of all 4 children.
      unsigned char upper_z[N];   //!< Z dimension of upper bounds of all 4 children.
    };

    /*! Compressed node with unaligned bounds */
    struct CUNode : public BaseNode
    {
      /*! Clears the node. */
      __forceinline void clear() 
      {
        xfm_vx[0] = 1; xfm_vx[1] = 0; xfm_vx[2] = 0;
        xfm_vy[0] = 0; xfm_vy[1] = 1; xfm_vy[2] = 0;
        xfm_vz[0] = 0; xfm_vz[1] = 0; xfm_vz[2] = 1;
        offset = 0.0f; scale = 0.0f;
        lower_x[0] = lower_x[1] = lower_x[2] = lower_x[3] = 1;
        lower_y[0] = lower_y[1] = lower_y[2] = lower_y[3] = 1;
        lower_z[0] = lower_z[1] = lower_z[2] = lower_z[3] = 1;
        upper_x[0] = upper_x[1] = upper_x[2] = upper_x[3] = 1;
        upper_y[0] = upper_y[1] = upper_y[2] = upper_y[3] = 1;
        upper_z[0] = upper_z[1] = upper_z[2] = upper_z[3] = 1;
        align[0] = align[1] = align[2] = align[3] = 0;
        BaseNode::clear();
      }

      /*! Sets non-axis aligned space of node and parent bounding box. */
      __forceinline void set(const NAABBox3fa& naabb) 
      {
        const LinearSpace3fa& space = naabb.space;
        const BBox3fa& bounds = naabb.bounds;
        xfm_vx[0] = (char) (127.0f*space.vx.x); assert(127.0f*space.vx.x >= -127.0f && 127.0f*space.vx.x <= 127.0f && truncf(127.0f*space.vx.x) == 127.0f*space.vx.x);
        xfm_vx[1] = (char) (127.0f*space.vx.y); assert(127.0f*space.vx.y >= -127.0f && 127.0f*space.vx.y <= 127.0f && truncf(127.0f*space.vx.y) == 127.0f*space.vx.y);
        xfm_vx[2] = (char) (127.0f*space.vx.z); assert(127.0f*space.vx.z >= -127.0f && 127.0f*space.vx.z <= 127.0f && truncf(127.0f*space.vx.z) == 127.0f*space.vx.z);
        xfm_vx[3] = 0;
        xfm_vy[0] = (char) (127.0f*space.vy.x); assert(127.0f*space.vy.x >= -127.0f && 127.0f*space.vy.x <= 127.0f && truncf(127.0f*space.vy.x) == 127.0f*space.vy.x);
        xfm_vy[1] = (char) (127.0f*space.vy.y); assert(127.0f*space.vy.y >= -127.0f && 127.0f*space.vy.y <= 127.0f && truncf(127.0f*space.vy.y) == 127.0f*space.vy.y);
        xfm_vy[2] = (char) (127.0f*space.vy.z); assert(127.0f*space.vy.z >= -127.0f && 127.0f*space.vy.z <= 127.0f && truncf(127.0f*space.vy.z) == 127.0f*space.vy.z);
        xfm_vy[3] = 0;
        xfm_vz[0] = (char) (127.0f*space.vz.x); assert(127.0f*space.vz.x >= -127.0f && 127.0f*space.vz.x <= 127.0f && truncf(127.0f*space.vz.x) == 127.0f*space.vz.x);
        xfm_vz[1] = (char) (127.0f*space.vz.y); assert(127.0f*space.vz.y >= -127.0f && 127.0f*space.vz.y <= 127.0f && truncf(127.0f*space.vz.y) == 127.0f*space.vz.y);
        xfm_vz[2] = (char) (127.0f*space.vz.z); assert(127.0f*space.vz.z >= -127.0f && 127.0f*space.vz.z <= 127.0f && truncf(127.0f*space.vz.z) == 127.0f*space.vz.z);
        xfm_vz[3] = 0;
        offset = 127.0f*bounds.lower;
        scale  = (127.0f*bounds.upper-127.0f*bounds.lower)/255.0f;
      }

      /*! Sets bounding box. */
      __forceinline void set(size_t i, const BBox3fa& bounds) 
      {
        assert(i < N);
        const Vec3fa lower = select(eq_mask(scale,Vec3fa(0.0f)),0.0f,(127.0f*bounds.lower-Vec3fa(offset))/Vec3fa(scale));
        assert(lower.x >= 0.0f && lower.x <= 255.01f); // FIXME: should be smaller than 255.0f
        assert(lower.y >= 0.0f && lower.y <= 255.01f);
        assert(lower.z >= 0.0f && lower.z <= 255.01f);
        lower_x[i] = (unsigned char) clamp(floorf(lower.x),0.0f,255.0f);
        lower_y[i] = (unsigned char) clamp(floorf(lower.y),0.0f,255.0f);
        lower_z[i] = (unsigned char) clamp(floorf(lower.z),0.0f,255.0f);
        const Vec3fa upper = select(eq_mask(scale,Vec3fa(0.0f)),0.0f,(127.0f*bounds.upper-Vec3fa(offset))/Vec3fa(scale));
        assert(upper.x >= 0.0f && upper.x <= 255.01f); // FIXME: should be smaller than 255.0f
        assert(upper.y >= 0.0f && upper.y <= 255.01f);
        assert(upper.z >= 0.0f && upper.z <= 255.01f);
        upper_x[i] = (unsigned char) clamp(ceilf(upper.x),0.0f,255.0f);
        upper_y[i] = (unsigned char) clamp(ceilf(upper.y),0.0f,255.0f);
        upper_z[i] = (unsigned char) clamp(ceilf(upper.z),0.0f,255.0f);
      }

      /*! Sets ID of child. */
      __forceinline void set(size_t i, const NodeRef& childID) {
        BaseNode::set(i,childID);
      }

      /*! returns transformation */
      __forceinline const LinearSpace3fa getXfm() const 
      {
        const ssef vx = ssef(_mm_cvtepi8_epi32(*(ssei*)&xfm_vx));
        const ssef vy = ssef(_mm_cvtepi8_epi32(*(ssei*)&xfm_vy));
        const ssef vz = ssef(_mm_cvtepi8_epi32(*(ssei*)&xfm_vz));
        return LinearSpace3fa((Vec3fa)vx,(Vec3fa)vy,(Vec3fa)vz);
      }

      /*! returns ith bounding boxes */
      __forceinline NAABBox3fa getBounds(size_t i) const 
      {
        assert(i < N);
        const Vec3f lower((float)lower_x[i],lower_y[i],lower_z[i]);
        const Vec3f upper((float)upper_x[i],upper_y[i],upper_z[i]);
        return NAABBox3fa(getXfm(),BBox3fa(Vec3fa(offset+scale*lower),Vec3fa(offset+scale*upper)));
      }

      /*! returns 4 bounding boxes */
      __forceinline BBoxSSE3f getBounds() const 
      {
        const ssef lower_x = ssef(_mm_cvtepu8_epi32(*(ssei*)&this->lower_x));
        const ssef lower_y = ssef(_mm_cvtepu8_epi32(*(ssei*)&this->lower_y));
        const ssef lower_z = ssef(_mm_cvtepu8_epi32(*(ssei*)&this->lower_z));
        const ssef upper_x = ssef(_mm_cvtepu8_epi32(*(ssei*)&this->upper_x));
        const ssef upper_y = ssef(_mm_cvtepu8_epi32(*(ssei*)&this->upper_y));
        const ssef upper_z = ssef(_mm_cvtepu8_epi32(*(ssei*)&this->upper_z));
        const Vec3<ssef> offset = *(Vec3fa*)&this->offset;
        const Vec3<ssef> scale  = *(Vec3fa*)&this->scale;
        return BBoxSSE3f(scale*Vec3<ssef>(lower_x,lower_y,lower_z)+offset,
                          scale*Vec3<ssef>(upper_x,upper_y,upper_z)+offset);
      }

      /*! Returns the extend of the bounds of the ith child */
      __forceinline Vec3fa extend(size_t i) const {
        assert(i < N);
        return getBounds(i).bounds.size()/127.0f;
      }

    public:
      char xfm_vx[4];             //!< 1st column of transformation
      char xfm_vy[4];             //!< 2nd column of transformation
      char xfm_vz[4];             //!< 3rd column of transformation
      char align[4];              
      Vec3f offset;               //!< offset to decompress bounds
      Vec3f scale;                //!< scale  to decompress bounds
      unsigned char lower_x[N];   //!< X dimension of lower bounds of all 4 children.
      unsigned char lower_y[N];   //!< Y dimension of lower bounds of all 4 children.
      unsigned char lower_z[N];   //!< Z dimension of lower bounds of all 4 children.
      unsigned char upper_x[N];   //!< X dimension of upper bounds of all 4 children.
      unsigned char upper_y[N];   //!< Y dimension of upper bounds of all 4 children.
      unsigned char upper_z[N];   //!< Z dimension of upper bounds of all 4 children.
    };

    /*! swap the children of two nodes */
    __forceinline static void swap(UANode* a, size_t i, UANode* b, size_t j)
    {
      assert(i<4 && j<4);
      std::swap(a->children[i],b->children[j]);
      std::swap(a->lower_x[i],b->lower_x[j]);
      std::swap(a->lower_y[i],b->lower_y[j]);
      std::swap(a->lower_z[i],b->lower_z[j]);
      std::swap(a->upper_x[i],b->upper_x[j]);
      std::swap(a->upper_y[i],b->upper_y[j]);
      std::swap(a->upper_z[i],b->upper_z[j]);
    }

    /*! compacts a node (moves empty children to the end) */
    __forceinline static void compact(UANode* a)
    {
      /* find right most filled node */
      ssize_t j=4;
      for (j=j-1; j>=0; j--)
        if (a->child(j) != emptyNode)
          break;

      /* replace empty nodes with filled nodes */
      for (ssize_t i=0; i<j; i++) {
        if (a->child(i) == emptyNode) {
          swap(a,i,a,j);
          for (j=j-1; j>i; j--)
            if (a->child(j) != emptyNode)
              break;
        }
      }
    }
    
  public:

    /*! BVH4 default constructor. */
    BVH4 (const PrimitiveType& primTy0, void* geometry = NULL);
    BVH4 (const PrimitiveType& primTy0, const PrimitiveType& primTy1, void* geometry = NULL);
    BVH4 (const PrimitiveType& primTy0, const PrimitiveType& primTy1, const PrimitiveType& primTy2, void* geometry = NULL);

    /*! BVH4 destruction */
    ~BVH4 ();

    /*! returns name of BVH */
    std::string name() const {
      std::string str = "BVH4<";
      if (primTys[0]) str += primTys[0]->name;
      for (size_t i=1; i<4; i++) {
        str += ",";
        if (primTys[i]) str += primTys[i]->name;
      }
      str+= ">";
      return str;
    }

    /*! BVH4 instantiations */
    static Accel* BVH4Bezier1i(Scene* scene);
    static Accel* BVH4Triangle1(Scene* scene);
    static Accel* BVH4Triangle4(Scene* scene);
    static Accel* BVH4Triangle8(Scene* scene);
    static Accel* BVH4Triangle1v(Scene* scene);
    static Accel* BVH4Triangle4v(Scene* scene);
    static Accel* BVH4Triangle4i(Scene* scene);
    
    static Accel* BVH4BVH4Triangle1Morton(Scene* scene);
    static Accel* BVH4BVH4Triangle1ObjectSplit(Scene* scene);
    static Accel* BVH4BVH4Triangle4ObjectSplit(Scene* scene);
    static Accel* BVH4BVH4Triangle1vObjectSplit(Scene* scene);
    static Accel* BVH4BVH4Triangle4vObjectSplit(Scene* scene);

    static Accel* BVH4Triangle1SpatialSplit(Scene* scene);
    static Accel* BVH4Triangle4SpatialSplit(Scene* scene);
    static Accel* BVH4Triangle8SpatialSplit(Scene* scene);
    static Accel* BVH4Triangle1ObjectSplit(Scene* scene);
    static Accel* BVH4Triangle4ObjectSplit(Scene* scene);
    static Accel* BVH4Triangle8ObjectSplit(Scene* scene);
    static Accel* BVH4Triangle1vObjectSplit(Scene* scene);
    static Accel* BVH4Triangle4vObjectSplit(Scene* scene);
    static Accel* BVH4Triangle4iObjectSplit(Scene* scene);

    static Accel* BVH4Triangle1ObjectSplit(TriangleMesh* mesh);
    static Accel* BVH4Triangle4ObjectSplit(TriangleMesh* mesh);
    static Accel* BVH4Triangle1vObjectSplit(TriangleMesh* mesh);
    static Accel* BVH4Triangle4vObjectSplit(TriangleMesh* mesh);
    static Accel* BVH4Triangle4Refit(TriangleMesh* mesh);

    /*! initializes the acceleration structure */
    void init (size_t numPrimitives = 0);

    /*! Clears the barrier bits of a subtree. */
    void clearBarrier(NodeRef& node);

    Ref<LinearAllocatorPerThread> alloc; // FIXME: why using reference?


    __forceinline UANode* allocUANode(size_t thread) {
      UANode* node = (UANode*) alloc->malloc(thread,sizeof(UANode),1 << 7); node->clear(); return node;
    }

    __forceinline CANode* allocCANode(size_t thread) {
      CANode* node = (CANode*) alloc->malloc(thread,sizeof(CANode),1 << 7); node->clear(); return node;
    }

    __forceinline UUNode* allocUUNode(size_t thread) {
      UUNode* node = (UUNode*) alloc->malloc(thread,sizeof(UUNode),1 << 7); node->clear(); return node;
    }

    __forceinline CUNode* allocCUNode(size_t thread) {
      CUNode* node = (CUNode*) alloc->malloc(thread,sizeof(CUNode),1 << 7); node->clear(); return node;
    }

    /*__forceinline char* allocPrimitiveBlocks(size_t thread, size_t ty, size_t num) {
      return (char*) alloc->malloc(thread,num*primTys[ty]->bytes,alignment);
      }*/

    __forceinline char* allocPrimitiveBlocks(size_t thread, size_t num) {
      return (char*) alloc->malloc(thread,num*primTy.bytes,1 << 6);
    }

    /*! Encodes a node */
    __forceinline NodeRef encodeNode(UANode* node) { assert(((size_t)node & deco_mask) == 0); return NodeRef((size_t)node | tyUA); }
    __forceinline NodeRef encodeNode(CANode* node) { assert(((size_t)node & deco_mask) == 0); return NodeRef((size_t)node | tyCA); }
    __forceinline NodeRef encodeNode(UUNode* node) { assert(((size_t)node & deco_mask) == 0); return NodeRef((size_t)node | tyUU); }
    __forceinline NodeRef encodeNode(CUNode* node) { assert(((size_t)node & deco_mask) == 0); return NodeRef((size_t)node | tyCU); }
    
    /*! Encodes a leaf */
    __forceinline NodeRef encodeLeaf(char* tri, size_t num, size_t type = 0) {
      assert(((size_t)tri & deco_mask) == 0); 
      assert(type < maxLeafTypes);
      return NodeRef((size_t)tri | (tyL0 + type) | (min(num,maxLeafBlocks) << 60));
    }

  public:
    
    /*! calculates the amount of bytes allocated */
    size_t bytesAllocated() 
    {
      if (nodes || primitives)
        return bytesNodes+bytesPrimitives+numVertices*sizeof(Vec3fa);
      else
        return alloc->bytes()+numVertices*sizeof(Vec3fa);
    }

  public:
    const PrimitiveType& primTy;       //!< primitive type stored in the BVH
    const PrimitiveType* primTys[4];    //!< primitive types stored in the BVH
    
    void* geometry;                    //!< pointer to additional data for primitive intersector
    NodeRef root;                      //!< Root node
    size_t numPrimitives;
    size_t numVertices;

    /*! data arrays for fast builders */
  public:
    void* nodes;
    size_t bytesNodes;
    void* primitives;
    size_t bytesPrimitives;
    std::vector<BVH4*> objects;
  };

  // FIXME: move the below code to somewhere else
  typedef void (*createTriangleMeshAccelTy)(TriangleMesh* mesh, BVH4*& accel, Builder*& builder); 
  typedef Builder* (*BVH4BuilderTopLevelFunc)(BVH4* accel, Scene* scene, const createTriangleMeshAccelTy createTriangleMeshAccel);

#define DECLARE_TOPLEVEL_BUILDER(symbol)                                         \
  namespace isa   { extern Builder* symbol(BVH4* accel, Scene* scene, const createTriangleMeshAccelTy createTriangleMeshAccel); } \
  namespace sse41 { extern Builder* symbol(BVH4* accel, Scene* scene, const createTriangleMeshAccelTy createTriangleMeshAccel); } \
  namespace avx   { extern Builder* symbol(BVH4* accel, Scene* scene, const createTriangleMeshAccelTy createTriangleMeshAccel); } \
  namespace avx2  { extern Builder* symbol(BVH4* accel, Scene* scene, const createTriangleMeshAccelTy createTriangleMeshAccel); } \
  BVH4BuilderTopLevelFunc symbol;
}
