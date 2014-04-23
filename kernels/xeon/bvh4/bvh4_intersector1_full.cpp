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

#include "bvh4_intersector1_full.h"

#include "geometry/select_intersector1.h"
#include "geometry/triangle4_intersector1_moeller.h"

#if defined(__AVX__)
#include "geometry/bezier1_intersector1.h"
#include "geometry/bezier1i_intersector1.h"
#endif

#define TFAR(x) x
#define NTFAR(x)   

namespace embree
{ 
  namespace isa
  {
#if 0

#if defined(__AVX__)

    template<typename PrimitiveIntersector>
    __forceinline size_t BVH4Intersector1Full<PrimitiveIntersector>::intersectBox(const BVH4::UANode* node, 
                                                                                  const sse3f& org, const sse3f& rdir, const sse3f& org_rdir, 
                                                                                  const size_t nearX, const size_t nearY, const size_t nearZ,
                                                                                  ssef& tNear, ssef& tFar)
    {
      const BBoxSSE3f bounds = node->getBounds(nearX,nearY,nearZ);

#if defined (__AVX2__)
      const ssef tNearX = msub(bounds.lower.x, rdir.x, org_rdir.x);
      const ssef tNearY = msub(bounds.lower.y, rdir.y, org_rdir.y);
      const ssef tNearZ = msub(bounds.lower.z, rdir.z, org_rdir.z);
      const ssef tFarX  = msub(bounds.upper.x, rdir.x, org_rdir.x);
      const ssef tFarY  = msub(bounds.upper.y, rdir.y, org_rdir.y);
      const ssef tFarZ  = msub(bounds.upper.z, rdir.z, org_rdir.z);
#else
      const ssef tNearX = (bounds.lower.x - org.x) * rdir.x;
      const ssef tNearY = (bounds.lower.y - org.y) * rdir.y;
      const ssef tNearZ = (bounds.lower.z - org.z) * rdir.z;
      const ssef tFarX  = (bounds.upper.x - org.x) * rdir.x;
      const ssef tFarY  = (bounds.upper.y - org.y) * rdir.y;
      const ssef tFarZ  = (bounds.upper.z - org.z) * rdir.z;
#endif
      
#if defined(__SSE4_1__)
      tNear = maxi(maxi(tNearX,tNearY),maxi(tNearZ,tNear));
      tFar  = mini(mini(tFarX ,tFarY ),mini(tFarZ ,tFar));
      const sseb vmask = cast(tNear) > cast(tFar);
      return movemask(vmask)^((1<<BVH4::N)-1);
#else
      tNear = max(tNearX,tNearY,tNearZ,tNear);
      tFar  = min(tFarX ,tFarY ,tFarZ ,tFar);
      const sseb vmask = tNear <= tFar;
      return movemask(vmask);
#endif
    }

    template<typename PrimitiveIntersector>
    __forceinline size_t BVH4Intersector1Full<PrimitiveIntersector>::intersectBox(const BVH4::UUNode* node, Ray& ray, 
                                                                                  const avx3f& ray_org_dir, const sse3f& ray_org, const sse3f& ray_dir, 
                                                                                  ssef& tNear, ssef& tFar)
    {
#if BVH4HAIR_COMPRESS_UNALIGNED_NODES
      const LinearSpace3fa xfm = node->getXfm();
      //const Vec3fa dir = xfmVector(xfm,ray.dir);
      const Vec3fa dir = madd(xfm.vx,(Vec3fa)ray_dir.x,madd(xfm.vy,(Vec3fa)ray_dir.y,xfm.vz*(Vec3fa)ray_dir.z));
      //const sse3f rdir = Vec3fa(one)/dir; // FIXME: not 100% safe
      const sse3f rdir = rcp_safe(dir); 
      //const Vec3fa org = xfmPoint(xfm,ray.org);
      const Vec3fa org = madd(xfm.vx,(Vec3fa)ray_org.x,madd(xfm.vy,(Vec3fa)ray_org.y,xfm.vz*(Vec3fa)ray_org.z));
      const sse3f vorg  = sse3f(org);
      const sse3f vrdir = sse3f(rdir);
      const BBoxSSE3f bounds = node->getBounds();
      const sse3f tLowerXYZ = (bounds.lower - vorg) * vrdir;
      const sse3f tUpperXYZ = (bounds.upper - vorg) * vrdir;
#else
#if 1
      const AffineSpaceSSE3f xfm = node->naabb;
      const sse3f dir = xfmVector(xfm,ray_dir);
      //const sse3f nrdir = sse3f(ssef(-1.0f))/dir; // FIXME: not 100% safe
      const sse3f nrdir = sse3f(ssef(-1.0f))*rcp_safe(dir);
      const sse3f org = xfmPoint(xfm,ray_org);
      const sse3f tLowerXYZ = org * nrdir;     // (Vec3fa(zero) - org) * rdir;
      const sse3f tUpperXYZ = tLowerXYZ - nrdir; // (Vec3fa(one ) - org) * rdir;
#else
      const AffineSpaceSSE3f xfm = node->naabb;
      const avxf orgdirx = madd(avxf(xfm.l.vx.x),ray_org_dir.x,
                                madd(avxf(xfm.l.vy.x),ray_org_dir.y,
                                     madd(avxf(xfm.l.vz.x),ray_org_dir.z,
                                          avxf(xfm.p.x,zero))));
      const avxf orgdiry = madd(avxf(xfm.l.vx.y),ray_org_dir.x,
                                madd(avxf(xfm.l.vy.y),ray_org_dir.y,
                                     madd(avxf(xfm.l.vz.y),ray_org_dir.z,
                                          avxf(xfm.p.y,zero))));
      const avxf orgdirz = madd(avxf(xfm.l.vx.z),ray_org_dir.x,
                                madd(avxf(xfm.l.vy.z),ray_org_dir.y,
                                     madd(avxf(xfm.l.vz.z),ray_org_dir.z,
                                          avxf(xfm.p.z,zero))));
      const sse3f org(extract<0>(orgdirx),extract<0>(orgdiry),extract<0>(orgdirz));
      const sse3f dir(extract<1>(orgdirx),extract<1>(orgdiry),extract<1>(orgdirz));
      const sse3f nrdir = sse3f(ssef(-1.0f))/dir; //rcp_safe(dir); // FIXME: not 100% safe
      const sse3f tLowerXYZ = org * nrdir;     // (Vec3fa(zero) - org) * rdir;
      const sse3f tUpperXYZ = tLowerXYZ - nrdir; // (Vec3fa(one ) - org) * rdir;
#endif
#endif
      
#if defined(__SSE4_1__)
      const ssef tNearX = mini(tLowerXYZ.x,tUpperXYZ.x);
      const ssef tNearY = mini(tLowerXYZ.y,tUpperXYZ.y);
      const ssef tNearZ = mini(tLowerXYZ.z,tUpperXYZ.z);
      const ssef tFarX  = maxi(tLowerXYZ.x,tUpperXYZ.x);
      const ssef tFarY  = maxi(tLowerXYZ.y,tUpperXYZ.y);
      const ssef tFarZ  = maxi(tLowerXYZ.z,tUpperXYZ.z);
      tNear = maxi(maxi(tNearX,tNearY),maxi(tNearZ,tNear));
      tFar  = mini(mini(tFarX ,tFarY ),mini(tFarZ ,tFar));
      const sseb vmask = tNear <= tFar;
      return movemask(vmask);
#else
      const ssef tNearX = min(tLowerXYZ.x,tUpperXYZ.x);
      const ssef tNearY = min(tLowerXYZ.y,tUpperXYZ.y);
      const ssef tNearZ = min(tLowerXYZ.z,tUpperXYZ.z);
      const ssef tFarX  = max(tLowerXYZ.x,tUpperXYZ.x);
      const ssef tFarY  = max(tLowerXYZ.y,tUpperXYZ.y);
      const ssef tFarZ  = max(tLowerXYZ.z,tUpperXYZ.z);
      tNear = max(tNearX,tNearY,tNearZ,tNear);
      tFar  = min(tFarX ,tFarY ,tFarZ ,tFar);
      const sseb vmask = tNear <= tFar;
      return movemask(vmask);
#endif
    }

    template<typename PrimitiveIntersector>
    void BVH4Intersector1Full<PrimitiveIntersector>::intersect(const BVH4* bvh, Ray& ray)
    {
      /*! perform per ray precalculations required by the primitive intersector */
      const typename PrimitiveIntersector::Precalculations pre(ray);

      /*! stack state */
      StackItem stack[stackSize];  //!< stack of nodes 
      StackItem* stackPtr = stack+1;        //!< current stack pointer
      StackItem* stackEnd = stack+stackSize;
      stack[0].ref = bvh->root;
      stack[0].tNear = ray.tnear;
      stack[0].tFar = ray.tfar;
      
      /*! offsets to select the side that becomes the lower or upper bound */
      const size_t nearX = ray.dir.x >= 0.0f ? 0*BVH4::UANode::stride : 1*BVH4::UANode::stride;
      const size_t nearY = ray.dir.y >= 0.0f ? 0*BVH4::UANode::stride : 1*BVH4::UANode::stride;
      const size_t nearZ = ray.dir.z >= 0.0f ? 0*BVH4::UANode::stride : 1*BVH4::UANode::stride;
      
      /*! load the ray into SSE registers */
      const sse3f org(ray.org.x,ray.org.y,ray.org.z);
      const sse3f dir(ray.dir.x,ray.dir.y,ray.dir.z);
      const Vec3fa ray_rdir = rcp_safe(ray.dir);
      const sse3f rdir(ray_rdir.x,ray_rdir.y,ray_rdir.z);
      const Vec3fa ray_org_rdir = ray.org*ray_rdir;
      const sse3f org_rdir(ray_org_rdir.x,ray_org_rdir.y,ray_org_rdir.z);
      const avx3f org_dir(avxf(org.x,dir.x),avxf(org.y,dir.y),avxf(org.z,dir.z));

      /* pop loop */
      while (true) pop:
      {
        /*! pop next node */
        if (unlikely(stackPtr == stack)) break;
        /*for (size_t i=1; i<stackPtr-&stack[0]; i++)
          if (stack[i-1].tNear < stack[i+0].tNear)
            StackItem::swap2(stack[i-1],stack[i+0]);*/
        stackPtr--;
        BVH4::NodeRef cur = BVH4::NodeRef(stackPtr->ref);
        ssef tNear = stackPtr->tNear;
        ssef tFar = min(stackPtr->tFar,ray.tfar);
        
        /*! if popped node is too far, pop next one */
        if (unlikely(_mm_cvtss_f32(tNear) > _mm_cvtss_f32(tFar)))
          continue;

        /* downtraversal loop */
        while (true)
        {
          /*! process nodes with aligned bounds */
          size_t mask;
          if (likely(cur.isUANode()))
            mask = intersectBox(cur.getUANode(),org,rdir,org_rdir,nearX,nearY,nearZ,tNear,tFar);

          /*! process nodes with unaligned bounds */
          else if (unlikely(cur.isUUNode()))
            mask = intersectBox(cur.getUUNode(),ray,org_dir,org,dir,tNear,tFar);

          /*! otherwise this is a leaf */
          else break;

          /*! if no child is hit, pop next node */
          STAT3(normal.trav_nodes,1,1,1);
          const BVH4::Node* node = cur.getNode();
          if (unlikely(mask == 0))
            goto pop;
          
          /*! one child is hit, continue with that child */
          size_t r = __bscf(mask);
          BVH4::NodeRef c0 = node->child(r);
          c0.prefetch();

          if (likely(mask == 0)) {
            cur = c0;  tNear = tNear[r]; tFar = tFar[r];
            assert(cur != BVH4::emptyNode);
            continue;
          }

          /*! two children are hit, push far child, and continue with closer child */
          const float n0 = tNear[r]; const float f0 = tFar[r]; 
          r = __bscf(mask);
          BVH4::NodeRef c1 = node->child(r); c1.prefetch(); const float n1 = tNear[r]; const float f1 = tFar[r];
          //assert(c0 != BVH4::emptyNode); // FIXME: enable
          //assert(c1 != BVH4::emptyNode); // FIXME: enable
          if (likely(mask == 0)) {
            assert(stackPtr < stackEnd); 
            if (n0 < n1) { 
              stackPtr->ref = c1; stackPtr->tNear = n1; stackPtr->tFar = f1; stackPtr++; 
              cur = c0; tNear = n0; tFar = f0;
              continue; 
            }
            else { 
              stackPtr->ref = c0; stackPtr->tNear = n0; stackPtr->tFar = f0; stackPtr++; 
              cur = c1; tNear = n1; tFar = f1; 
              continue; 
            }
          }
          
          /*! Here starts the slow path for 3 or 4 hit children. We push
           *  all nodes onto the stack to sort them there. */
          assert(stackPtr < stackEnd); 
          stackPtr->ref = c0; stackPtr->tNear = n0; stackPtr->tFar = f0; stackPtr++;
          assert(stackPtr < stackEnd); 
          stackPtr->ref = c1; stackPtr->tNear = n1; stackPtr->tFar = f1; stackPtr++;
          
          /*! three children are hit, push all onto stack and sort 3 stack items, continue with closest child */
          assert(stackPtr < stackEnd); 
          r = __bscf(mask);
          BVH4::NodeRef c = node->child(r); c.prefetch(); float n2 = tNear[r]; float f2 = tFar[r]; stackPtr->ref = c; stackPtr->tNear = n2; stackPtr->tFar = f2; stackPtr++;
          //assert(c != BVH4::emptyNode); // FIXME: enable
          if (likely(mask == 0)) {
            sort(stackPtr[-1],stackPtr[-2],stackPtr[-3]);
            cur = (BVH4::NodeRef) stackPtr[-1].ref; tNear = stackPtr[-1].tNear; tFar = stackPtr[-1].tFar; stackPtr--;
            continue;
          }

          /*! four children are hit, push all onto stack and sort 4 stack items, continue with closest child */
          assert(stackPtr < stackEnd); 
          r = __bscf(mask);
          c = node->child(r); c.prefetch(); float n3 = tNear[r]; float f3 = tFar[r]; stackPtr->ref = c; stackPtr->tNear = n3; stackPtr->tFar = f3; stackPtr++;
          //assert(c != BVH4::emptyNode); // FIXME: enable
          sort(stackPtr[-1],stackPtr[-2],stackPtr[-3],stackPtr[-4]);
          cur = (BVH4::NodeRef) stackPtr[-1].ref; tNear = stackPtr[-1].tNear; tFar = stackPtr[-1].tFar; stackPtr--;
        }
        
        /*! this is a leaf node */
        STAT3(normal.trav_leaves,1,1,1);
        size_t num, ty; void* prim = cur.getLeaf(num,ty);
        PrimitiveIntersector::intersect(pre,ray,ty,prim,num,bvh->geometry);
      }
      AVX_ZERO_UPPER();
    }

    template<typename PrimitiveIntersector>
    void BVH4Intersector1Full<PrimitiveIntersector>::occluded(const BVH4* bvh, Ray& ray) 
    {
      /*! perform per ray precalculations required by the primitive intersector */
      const typename PrimitiveIntersector::Precalculations pre(ray);

      /*! stack state */
      StackItem stack[stackSize];  //!< stack of nodes 
      StackItem* stackPtr = stack+1;        //!< current stack pointer
      StackItem* stackEnd = stack+stackSize;
      stack[0].ref = bvh->root;
      stack[0].tNear = ray.tnear;
      stack[0].tFar = ray.tfar;
      
      /*! offsets to select the side that becomes the lower or upper bound */
      const size_t nearX = ray.dir.x >= 0.0f ? 0*BVH4::UANode::stride : 1*BVH4::UANode::stride;
      const size_t nearY = ray.dir.y >= 0.0f ? 0*BVH4::UANode::stride : 1*BVH4::UANode::stride;
      const size_t nearZ = ray.dir.z >= 0.0f ? 0*BVH4::UANode::stride : 1*BVH4::UANode::stride;
      
      /*! load the ray into SSE registers */
      const sse3f org(ray.org.x,ray.org.y,ray.org.z);
      const sse3f dir(ray.dir.x,ray.dir.y,ray.dir.z);
      const Vec3fa ray_rdir = rcp_safe(ray.dir);
      const sse3f rdir(ray_rdir.x,ray_rdir.y,ray_rdir.z);
      const Vec3fa ray_org_rdir = ray.org*ray_rdir;
      const sse3f org_rdir(ray_org_rdir.x,ray_org_rdir.y,ray_org_rdir.z);
      const avx3f org_dir(avxf(org.x,dir.x),avxf(org.y,dir.y),avxf(org.z,dir.z));

      /* pop loop */
      while (true) pop:
      {
        /*! pop next node */
        if (unlikely(stackPtr == stack)) break;
        stackPtr--;
        BVH4::NodeRef cur = BVH4::NodeRef(stackPtr->ref);
        ssef tNear = stackPtr->tNear;
        ssef tFar = min(stackPtr->tFar,ray.tfar);
        
        /*! if popped node is too far, pop next one */
        if (unlikely(_mm_cvtss_f32(tNear) > _mm_cvtss_f32(tFar)))
          continue;

        /* downtraversal loop */
        while (true)
        {
          /*! process nodes with aligned bounds */
          size_t mask;
          if (likely(cur.isUANode()))
            mask = intersectBox(cur.getUANode(),org,rdir,org_rdir,nearX,nearY,nearZ,tNear,tFar);

          /*! process nodes with unaligned bounds */
          else if (unlikely(cur.isUUNode()))
            mask = intersectBox(cur.getUUNode(),ray,org_dir,org,dir,tNear,tFar);

          /*! otherwise this is a leaf */
          else break;

          /*! if no child is hit, pop next node */
          STAT3(shadow.trav_nodes,1,1,1);
          const BVH4::Node* node = cur.getNode();
          if (unlikely(mask == 0))
            goto pop;
          
          /*! one child is hit, continue with that child */
          size_t r = __bscf(mask);
          BVH4::NodeRef c0 = node->child(r); c0.prefetch();

          if (likely(mask == 0)) {
            cur = c0; tNear = tNear[r]; tFar = tFar[r];
            assert(cur != BVH4::emptyNode);
            continue;
          }
     
          /*! two children are hit, push far child, and continue with closer child */
           const float n0 = tNear[r]; const float f0 = tFar[r]; 
          r = __bscf(mask);
          BVH4::NodeRef c1 = node->child(r); c1.prefetch(); const float n1 = tNear[r]; const float f1 = tFar[r];
          //assert(c0 != BVH4::emptyNode); // FIXME: enable
          //assert(c1 != BVH4::emptyNode); // FIXME: enable
          if (likely(mask == 0)) {
            assert(stackPtr < stackEnd); 
            if (n0 < n1) { stackPtr->ref = c1; stackPtr->tNear = n1; stackPtr->tFar = f1; stackPtr++; cur = c0; tNear = n0; tFar = f0; continue; }
            else         { stackPtr->ref = c0; stackPtr->tNear = n0; stackPtr->tFar = f0; stackPtr++; cur = c1; tNear = n1; tFar = f1; continue; }
          }
          
          /*! Here starts the slow path for 3 or 4 hit children. We push
           *  all nodes onto the stack to sort them there. */
          assert(stackPtr < stackEnd); 
          stackPtr->ref = c0; stackPtr->tNear = n0; stackPtr->tFar = f0; stackPtr++;
          assert(stackPtr < stackEnd); 
          stackPtr->ref = c1; stackPtr->tNear = n1; stackPtr->tFar = f1; stackPtr++;
          
          /*! three children are hit, push all onto stack and sort 3 stack items, continue with closest child */
          assert(stackPtr < stackEnd); 
          r = __bscf(mask);
          BVH4::NodeRef c = node->child(r); c.prefetch(); float n2 = tNear[r]; float f2 = tFar[r]; stackPtr->ref = c; stackPtr->tNear = n2; stackPtr->tFar = f2; stackPtr++;
          //assert(c != BVH4::emptyNode); // FIXME: enable
          if (likely(mask == 0)) {
            sort(stackPtr[-1],stackPtr[-2],stackPtr[-3]);
            cur = (BVH4::NodeRef) stackPtr[-1].ref; tNear = stackPtr[-1].tNear; tFar = stackPtr[-1].tFar; stackPtr--;
            continue;
          }

          /*! four children are hit, push all onto stack and sort 4 stack items, continue with closest child */
          assert(stackPtr < stackEnd); 
          r = __bscf(mask);
          c = node->child(r); c.prefetch(); float n3 = tNear[r]; float f3 = tFar[r]; stackPtr->ref = c; stackPtr->tNear = n3; stackPtr->tFar = f3; stackPtr++;
          //assert(c != BVH4::emptyNode); // FIXME: enable
          sort(stackPtr[-1],stackPtr[-2],stackPtr[-3],stackPtr[-4]);
          cur = (BVH4::NodeRef) stackPtr[-1].ref; tNear = stackPtr[-1].tNear; tFar = stackPtr[-1].tFar; stackPtr--;
        }
        
        /*! this is a leaf node */
        STAT3(shadow.trav_leaves,1,1,1);
        size_t num,ty; void* prim = cur.getLeaf(num,ty);
        if (PrimitiveIntersector::occluded(pre,ray,ty,prim,num,bvh->geometry)) {
          ray.geomID = 0;
          break;
        }
      }
      AVX_ZERO_UPPER();
    }

#endif

#else

    template<typename PrimitiveIntersector>
    __forceinline size_t BVH4Intersector1Full<PrimitiveIntersector>::intersectBox(const BVH4::UANode* node, 
                                                                         const sse3f& org, const sse3f& rdir, const sse3f& org_rdir, 
                                                                         const size_t nearX, const size_t nearY, const size_t nearZ,
                                                                         ssef& tNear, ssef& tFar)
    {
      const BBoxSSE3f bounds = node->getBounds(nearX,nearY,nearZ);

#if defined (__AVX2__)
      const ssef tNearX = msub(bounds.lower.x, rdir.x, org_rdir.x);
      const ssef tNearY = msub(bounds.lower.y, rdir.y, org_rdir.y);
      const ssef tNearZ = msub(bounds.lower.z, rdir.z, org_rdir.z);
      const ssef tFarX  = msub(bounds.upper.x, rdir.x, org_rdir.x);
      const ssef tFarY  = msub(bounds.upper.y, rdir.y, org_rdir.y);
      const ssef tFarZ  = msub(bounds.upper.z, rdir.z, org_rdir.z);
#else
      const ssef tNearX = (bounds.lower.x - org.x) * rdir.x;
      const ssef tNearY = (bounds.lower.y - org.y) * rdir.y;
      const ssef tNearZ = (bounds.lower.z - org.z) * rdir.z;
      const ssef tFarX  = (bounds.upper.x - org.x) * rdir.x;
      const ssef tFarY  = (bounds.upper.y - org.y) * rdir.y;
      const ssef tFarZ  = (bounds.upper.z - org.z) * rdir.z;
#endif
      
#if defined(__SSE4_1__)
      tNear = maxi(maxi(tNearX,tNearY),maxi(tNearZ,tNear));
      tFar  = mini(mini(tFarX ,tFarY ),mini(tFarZ ,tFar));
      const sseb vmask = cast(tNear) > cast(tFar);
      return movemask(vmask)^((1<<BVH4::N)-1);
#else
      tNear = max(tNearX,tNearY,tNearZ,tNear);
      tFar  = min(tFarX ,tFarY ,tFarZ ,tFar);
      const sseb vmask = tNear <= tFar;
      return movemask(vmask);
#endif
    }

    template<typename PrimitiveIntersector>
    __forceinline size_t BVH4Intersector1Full<PrimitiveIntersector>::intersectBox(const BVH4::CANode* node, 
                                                                         const sse3f& org, const sse3f& rdir, const sse3f& org_rdir, 
                                                                         const size_t nearX, const size_t nearY, const size_t nearZ,
                                                                         ssef& tNear, ssef& tFar)
    {
      const BBoxSSE3f bounds = node->getBounds(nearX,nearY,nearZ);

#if defined (__AVX2__)
      const ssef tNearX = msub(bounds.lower.x, rdir.x, org_rdir.x);
      const ssef tNearY = msub(bounds.lower.y, rdir.y, org_rdir.y);
      const ssef tNearZ = msub(bounds.lower.z, rdir.z, org_rdir.z);
      const ssef tFarX  = msub(bounds.upper.x, rdir.x, org_rdir.x);
      const ssef tFarY  = msub(bounds.upper.y, rdir.y, org_rdir.y);
      const ssef tFarZ  = msub(bounds.upper.z, rdir.z, org_rdir.z);
#else
      const ssef tNearX = (bounds.lower.x - org.x) * rdir.x;
      const ssef tNearY = (bounds.lower.y - org.y) * rdir.y;
      const ssef tNearZ = (bounds.lower.z - org.z) * rdir.z;
      const ssef tFarX  = (bounds.upper.x - org.x) * rdir.x;
      const ssef tFarY  = (bounds.upper.y - org.y) * rdir.y;
      const ssef tFarZ  = (bounds.upper.z - org.z) * rdir.z;
#endif
      
#if defined(__SSE4_1__)
      tNear = maxi(maxi(tNearX,tNearY),maxi(tNearZ,tNear));
      tFar  = mini(mini(tFarX ,tFarY ),mini(tFarZ ,tFar));
      const sseb vmask = cast(tNear) > cast(tFar);
      return movemask(vmask)^((1<<BVH4::N)-1);
#else
      tNear = max(tNearX,tNearY,tNearZ,tNear);
      tFar  = min(tFarX ,tFarY ,tFarZ ,tFar);
      const sseb vmask = tNear <= tFar;
      return movemask(vmask);
#endif
    }

    template<typename PrimitiveIntersector>
    __forceinline size_t BVH4Intersector1Full<PrimitiveIntersector>::intersectBox(const BVH4::UUNode* node, Ray& ray, 
                                                                         const sse3f& ray_org, const sse3f& ray_dir, 
                                                                         ssef& tNear, ssef& tFar)
    {
      const AffineSpaceSSE3f xfm = node->naabb;
      const sse3f dir = xfmVector(xfm,ray_dir);
      //const sse3f nrdir = sse3f(ssef(-1.0f))/dir; // FIXME: not 100% safe
      const sse3f nrdir = sse3f(ssef(-1.0f))*rcp_safe(dir);
      const sse3f org = xfmPoint(xfm,ray_org);
      const sse3f tLowerXYZ = org * nrdir;     // (Vec3fa(zero) - org) * rdir;
      const sse3f tUpperXYZ = tLowerXYZ - nrdir; // (Vec3fa(one ) - org) * rdir;
      
#if defined(__SSE4_1__)
      const ssef tNearX = mini(tLowerXYZ.x,tUpperXYZ.x);
      const ssef tNearY = mini(tLowerXYZ.y,tUpperXYZ.y);
      const ssef tNearZ = mini(tLowerXYZ.z,tUpperXYZ.z);
      const ssef tFarX  = maxi(tLowerXYZ.x,tUpperXYZ.x);
      const ssef tFarY  = maxi(tLowerXYZ.y,tUpperXYZ.y);
      const ssef tFarZ  = maxi(tLowerXYZ.z,tUpperXYZ.z);
      tNear = maxi(maxi(tNearX,tNearY),maxi(tNearZ,tNear));
      tFar  = mini(mini(tFarX ,tFarY ),mini(tFarZ ,tFar));
      const sseb vmask = tNear <= tFar;
      return movemask(vmask);
#else
      const ssef tNearX = min(tLowerXYZ.x,tUpperXYZ.x);
      const ssef tNearY = min(tLowerXYZ.y,tUpperXYZ.y);
      const ssef tNearZ = min(tLowerXYZ.z,tUpperXYZ.z);
      const ssef tFarX  = max(tLowerXYZ.x,tUpperXYZ.x);
      const ssef tFarY  = max(tLowerXYZ.y,tUpperXYZ.y);
      const ssef tFarZ  = max(tLowerXYZ.z,tUpperXYZ.z);
      tNear = max(tNearX,tNearY,tNearZ,tNear);
      tFar  = min(tFarX ,tFarY ,tFarZ ,tFar);
      const sseb vmask = tNear <= tFar;
      return movemask(vmask);
#endif
    }

    template<typename PrimitiveIntersector>
    __forceinline size_t BVH4Intersector1Full<PrimitiveIntersector>::intersectBox(const BVH4::CUNode* node, Ray& ray, 
                                                                         const sse3f& ray_org, const sse3f& ray_dir, 
                                                                         ssef& tNear, ssef& tFar)
    {
      const LinearSpace3fa xfm = node->getXfm();
      //const Vec3fa dir = xfmVector(xfm,ray.dir);
      const Vec3fa dir = madd(xfm.vx,(Vec3fa)ray_dir.x,madd(xfm.vy,(Vec3fa)ray_dir.y,xfm.vz*(Vec3fa)ray_dir.z));
      //const sse3f rdir = Vec3fa(one)/dir; // FIXME: not 100% safe
      const sse3f rdir = rcp_safe(dir); 
      //const Vec3fa org = xfmPoint(xfm,ray.org);
      const Vec3fa org = madd(xfm.vx,(Vec3fa)ray_org.x,madd(xfm.vy,(Vec3fa)ray_org.y,xfm.vz*(Vec3fa)ray_org.z));
      const sse3f vorg  = sse3f(org);
      const sse3f vrdir = sse3f(rdir);
      const BBoxSSE3f bounds = node->getBounds();
      const sse3f tLowerXYZ = (bounds.lower - vorg) * vrdir;
      const sse3f tUpperXYZ = (bounds.upper - vorg) * vrdir;
      
#if defined(__SSE4_1__)
      const ssef tNearX = mini(tLowerXYZ.x,tUpperXYZ.x);
      const ssef tNearY = mini(tLowerXYZ.y,tUpperXYZ.y);
      const ssef tNearZ = mini(tLowerXYZ.z,tUpperXYZ.z);
      const ssef tFarX  = maxi(tLowerXYZ.x,tUpperXYZ.x);
      const ssef tFarY  = maxi(tLowerXYZ.y,tUpperXYZ.y);
      const ssef tFarZ  = maxi(tLowerXYZ.z,tUpperXYZ.z);
      tNear = maxi(maxi(tNearX,tNearY),maxi(tNearZ,tNear));
      tFar  = mini(mini(tFarX ,tFarY ),mini(tFarZ ,tFar));
      const sseb vmask = tNear <= tFar;
      return movemask(vmask);
#else
      const ssef tNearX = min(tLowerXYZ.x,tUpperXYZ.x);
      const ssef tNearY = min(tLowerXYZ.y,tUpperXYZ.y);
      const ssef tNearZ = min(tLowerXYZ.z,tUpperXYZ.z);
      const ssef tFarX  = max(tLowerXYZ.x,tUpperXYZ.x);
      const ssef tFarY  = max(tLowerXYZ.y,tUpperXYZ.y);
      const ssef tFarZ  = max(tLowerXYZ.z,tUpperXYZ.z);
      tNear = max(tNearX,tNearY,tNearZ,tNear);
      tFar  = min(tFarX ,tFarY ,tFarZ ,tFar);
      const sseb vmask = tNear <= tFar;
      return movemask(vmask);
#endif
    }

    template<typename PrimitiveIntersector>
    void BVH4Intersector1Full<PrimitiveIntersector>::intersect(const BVH4* bvh, Ray& ray)
    {
      /*! perform per ray precalculations required by the primitive intersector */
      const typename PrimitiveIntersector::Precalculations pre(ray);

      /*! stack state */
      StackItem stack[stackSize];  //!< stack of nodes 
      StackItem* stackPtr = stack+1;        //!< current stack pointer
      StackItem* stackEnd = stack+stackSize;
      stack[0].ref = bvh->root;
      stack[0].tNear = ray.tnear;
      TFAR(stack[0].tFar = ray.tfar);
      
      /*! offsets to select the side that becomes the lower or upper bound */
      const size_t nearX = ray.dir.x >= 0.0f ? 0*BVH4::UANode::stride : 1*BVH4::UANode::stride;
      const size_t nearY = ray.dir.y >= 0.0f ? 0*BVH4::UANode::stride : 1*BVH4::UANode::stride;
      const size_t nearZ = ray.dir.z >= 0.0f ? 0*BVH4::UANode::stride : 1*BVH4::UANode::stride;
      
      /*! load the ray into SSE registers */
      const sse3f org(ray.org.x,ray.org.y,ray.org.z);
      const sse3f dir(ray.dir.x,ray.dir.y,ray.dir.z);
      const Vec3fa ray_rdir = rcp_safe(ray.dir);
      const sse3f rdir(ray_rdir.x,ray_rdir.y,ray_rdir.z);
      const Vec3fa ray_org_rdir = ray.org*ray_rdir;
      const sse3f org_rdir(ray_org_rdir.x,ray_org_rdir.y,ray_org_rdir.z);

      /* pop loop */
      while (true) pop:
      {
        /*! pop next node */
        if (unlikely(stackPtr == stack)) break;
        /*for (size_t i=1; i<stackPtr-&stack[0]; i++)
          if (stack[i-1].tNear < stack[i+0].tNear)
            StackItem::swap2(stack[i-1],stack[i+0]);*/
        stackPtr--;
        BVH4::NodeRef cur = BVH4::NodeRef(stackPtr->ref);
        ssef tNear = stackPtr->tNear;
        TFAR(ssef tFar = min(stackPtr->tFar,ray.tfar));
        NTFAR(ssef tFar = ray.tfar);
        
        /*! if popped node is too far, pop next one */
        if (unlikely(_mm_cvtss_f32(tNear) > _mm_cvtss_f32(tFar)))
          continue;

        /* downtraversal loop */
        while (true)
        {
          NTFAR(tNear = ray.tnear);
          NTFAR(tFar  = ray.tfar);

          /*! process nodes with aligned bounds */
          size_t mask;
          if (likely(cur.isUANode()))
            mask = intersectBox(cur.getUANode(),org,rdir,org_rdir,nearX,nearY,nearZ,tNear,tFar);

          /*! process nodes with unaligned bounds */
          else if (unlikely(cur.isUUNode()))
            mask = intersectBox(cur.getUUNode(),ray,org,dir,tNear,tFar);

          /*! otherwise this is a leaf */
          else break;

          /*! if no child is hit, pop next node */
          STAT3(normal.trav_nodes,1,1,1);
          const BVH4::BaseNode* node = cur.getNode();
          if (unlikely(mask == 0))
            goto pop;
          
          /*! one child is hit, continue with that child */
          size_t r = __bscf(mask);
          BVH4::NodeRef c0 = node->child(r);
          c0.prefetch();

          if (likely(mask == 0)) {
            cur = c0;  TFAR(tNear = tNear[r]; tFar = tFar[r]); 
            assert(cur != BVH4::emptyNode);
            continue;
          }

          /*! two children are hit, push far child, and continue with closer child */
          const float n0 = tNear[r]; TFAR(const float f0 = tFar[r]); 
          r = __bscf(mask);
          BVH4::NodeRef c1 = node->child(r); c1.prefetch(); const float n1 = tNear[r]; TFAR(const float f1 = tFar[r]);
          //assert(c0 != BVH4::emptyNode); // FIXME: enable
          //assert(c1 != BVH4::emptyNode); // FIXME: enable
          if (likely(mask == 0)) {
            assert(stackPtr < stackEnd); 
            if (n0 < n1) { 
              stackPtr->ref = c1; stackPtr->tNear = n1; TFAR(stackPtr->tFar = f1); stackPtr++; 
              cur = c0; TFAR(tNear = n0; tFar = f0); 
              continue; 
            }
            else { 
              stackPtr->ref = c0; stackPtr->tNear = n0; TFAR(stackPtr->tFar = f0); stackPtr++; 
              cur = c1; TFAR(tNear = n1; tFar = f1); 
              continue; 
            }
          }
          
          /*! Here starts the slow path for 3 or 4 hit children. We push
           *  all nodes onto the stack to sort them there. */
          assert(stackPtr < stackEnd); 
          stackPtr->ref = c0; stackPtr->tNear = n0; TFAR(stackPtr->tFar = f0); stackPtr++;
          assert(stackPtr < stackEnd); 
          stackPtr->ref = c1; stackPtr->tNear = n1; TFAR(stackPtr->tFar = f1); stackPtr++;
          
          /*! three children are hit, push all onto stack and sort 3 stack items, continue with closest child */
          assert(stackPtr < stackEnd); 
          r = __bscf(mask);
          BVH4::NodeRef c = node->child(r); c.prefetch(); float n2 = tNear[r]; TFAR(float f2 = tFar[r]); stackPtr->ref = c; stackPtr->tNear = n2; TFAR(stackPtr->tFar = f2); stackPtr++;
          //assert(c != BVH4::emptyNode); // FIXME: enable
          if (likely(mask == 0)) {
            sort(stackPtr[-1],stackPtr[-2],stackPtr[-3]);
            cur = (BVH4::NodeRef) stackPtr[-1].ref; TFAR(tNear = stackPtr[-1].tNear; tFar = stackPtr[-1].tFar); stackPtr--;
            continue;
          }

          /*! four children are hit, push all onto stack and sort 4 stack items, continue with closest child */
          assert(stackPtr < stackEnd); 
          r = __bscf(mask);
          c = node->child(r); c.prefetch(); float n3 = tNear[r]; TFAR(float f3 = tFar[r]); stackPtr->ref = c; stackPtr->tNear = n3; TFAR(stackPtr->tFar = f3); stackPtr++;
          //assert(c != BVH4::emptyNode); // FIXME: enable
          sort(stackPtr[-1],stackPtr[-2],stackPtr[-3],stackPtr[-4]);
          cur = (BVH4::NodeRef) stackPtr[-1].ref; TFAR(tNear = stackPtr[-1].tNear; tFar = stackPtr[-1].tFar); stackPtr--;
        }
        
        /*! this is a leaf node */
        STAT3(normal.trav_leaves,1,1,1);
        size_t num, ty; void* prim = cur.getLeaf(num,ty);
        PrimitiveIntersector::intersect(pre,ray,ty,prim,num,bvh->geometry);
      }
      AVX_ZERO_UPPER();
    }

    template<typename PrimitiveIntersector>
    void BVH4Intersector1Full<PrimitiveIntersector>::occluded(const BVH4* bvh, Ray& ray) 
    {
      /*! perform per ray precalculations required by the primitive intersector */
      const typename PrimitiveIntersector::Precalculations pre(ray);

      /*! stack state */
      StackItem stack[stackSize];  //!< stack of nodes 
      StackItem* stackPtr = stack+1;        //!< current stack pointer
      StackItem* stackEnd = stack+stackSize;
      stack[0].ref = bvh->root;
      stack[0].tNear = ray.tnear;
      TFAR(stack[0].tFar = ray.tfar);
      
      /*! offsets to select the side that becomes the lower or upper bound */
      const size_t nearX = ray.dir.x >= 0.0f ? 0*BVH4::UANode::stride : 1*BVH4::UANode::stride;
      const size_t nearY = ray.dir.y >= 0.0f ? 0*BVH4::UANode::stride : 1*BVH4::UANode::stride;
      const size_t nearZ = ray.dir.z >= 0.0f ? 0*BVH4::UANode::stride : 1*BVH4::UANode::stride;
      
      /*! load the ray into SSE registers */
      const sse3f org(ray.org.x,ray.org.y,ray.org.z);
      const sse3f dir(ray.dir.x,ray.dir.y,ray.dir.z);
      const Vec3fa ray_rdir = rcp_safe(ray.dir);
      const sse3f rdir(ray_rdir.x,ray_rdir.y,ray_rdir.z);
      const Vec3fa ray_org_rdir = ray.org*ray_rdir;
      const sse3f org_rdir(ray_org_rdir.x,ray_org_rdir.y,ray_org_rdir.z);

      /* pop loop */
      while (true) pop:
      {
        /*! pop next node */
        if (unlikely(stackPtr == stack)) break;
        stackPtr--;
        BVH4::NodeRef cur = BVH4::NodeRef(stackPtr->ref);
        ssef tNear = stackPtr->tNear;
        TFAR(ssef tFar = min(stackPtr->tFar,ray.tfar));
        NTFAR(ssef tFar = ray.tfar);

        /*! if popped node is too far, pop next one */
        if (unlikely(_mm_cvtss_f32(tNear) > _mm_cvtss_f32(tFar)))
          continue;

        /* downtraversal loop */
        while (true)
        {
          NTFAR(tNear = ray.tnear);
          NTFAR(tFar  = ray.tfar);

          /*! process nodes with aligned bounds */
          size_t mask;
          if (likely(cur.isUANode()))
            mask = intersectBox(cur.getUANode(),org,rdir,org_rdir,nearX,nearY,nearZ,tNear,tFar);

          /*! process nodes with unaligned bounds */
          else if (unlikely(cur.isUUNode()))
            mask = intersectBox(cur.getUUNode(),ray,org,dir,tNear,tFar);

          /*! otherwise this is a leaf */
          else break;

          /*! if no child is hit, pop next node */
          STAT3(shadow.trav_nodes,1,1,1);
          const BVH4::BaseNode* node = cur.getNode();
          if (unlikely(mask == 0))
            goto pop;
          
          /*! one child is hit, continue with that child */
          size_t r = __bscf(mask);
          BVH4::NodeRef c0 = node->child(r); c0.prefetch();

          if (likely(mask == 0)) {
            cur = c0; TFAR(tNear = tNear[r]; tFar = tFar[r]);
            assert(cur != BVH4::emptyNode);
            continue;
          }
     
          /*! two children are hit, push far child, and continue with closer child */
          const float n0 = tNear[r]; TFAR(const float f0 = tFar[r]); 
          r = __bscf(mask);
          BVH4::NodeRef c1 = node->child(r); c1.prefetch(); const float n1 = tNear[r]; TFAR(const float f1 = tFar[r]);
          //assert(c0 != BVH4::emptyNode); // FIXME: enable
          //assert(c1 != BVH4::emptyNode); // FIXME: enable
          if (likely(mask == 0)) {
            assert(stackPtr < stackEnd); 
            if (n0 < n1) { stackPtr->ref = c1; stackPtr->tNear = n1; TFAR(stackPtr->tFar = f1); stackPtr++; cur = c0; TFAR(tNear = n0; tFar = f0); continue; }
            else         { stackPtr->ref = c0; stackPtr->tNear = n0; TFAR(stackPtr->tFar = f0); stackPtr++; cur = c1; TFAR(tNear = n1; tFar = f1); continue; }
          }
          
          /*! Here starts the slow path for 3 or 4 hit children. We push
           *  all nodes onto the stack to sort them there. */
          assert(stackPtr < stackEnd); 
          stackPtr->ref = c0; stackPtr->tNear = n0; TFAR(stackPtr->tFar = f0); stackPtr++;
          assert(stackPtr < stackEnd); 
          stackPtr->ref = c1; stackPtr->tNear = n1; TFAR(stackPtr->tFar = f1); stackPtr++;
          
          /*! three children are hit, push all onto stack and sort 3 stack items, continue with closest child */
          assert(stackPtr < stackEnd); 
          r = __bscf(mask);
          BVH4::NodeRef c = node->child(r); c.prefetch(); float n2 = tNear[r]; TFAR(float f2 = tFar[r]); stackPtr->ref = c; stackPtr->tNear = n2; TFAR(stackPtr->tFar = f2); stackPtr++;
          //assert(c != BVH4::emptyNode); // FIXME: enable
          if (likely(mask == 0)) {
            sort(stackPtr[-1],stackPtr[-2],stackPtr[-3]);
            cur = (BVH4::NodeRef) stackPtr[-1].ref; TFAR(tNear = stackPtr[-1].tNear; tFar = stackPtr[-1].tFar); stackPtr--;
            continue;
          }

          /*! four children are hit, push all onto stack and sort 4 stack items, continue with closest child */
          assert(stackPtr < stackEnd); 
          r = __bscf(mask);
          c = node->child(r); c.prefetch(); float n3 = tNear[r]; TFAR(float f3 = tFar[r]); stackPtr->ref = c; stackPtr->tNear = n3; TFAR(stackPtr->tFar = f3); stackPtr++;
          //assert(c != BVH4::emptyNode); // FIXME: enable
          sort(stackPtr[-1],stackPtr[-2],stackPtr[-3],stackPtr[-4]);
          cur = (BVH4::NodeRef) stackPtr[-1].ref; TFAR(tNear = stackPtr[-1].tNear; tFar = stackPtr[-1].tFar); stackPtr--;
        }
        
        /*! this is a leaf node */
        STAT3(shadow.trav_leaves,1,1,1);
        size_t num,ty; void* prim = cur.getLeaf(num,ty);
        if (PrimitiveIntersector::occluded(pre,ray,ty,prim,num,bvh->geometry)) {
          ray.geomID = 0;
          break;
        }
      }
      AVX_ZERO_UPPER();
    }
#endif
    
    //DEFINE_INTERSECTOR1(BVH4Triangle4Intersector1Moeller,BVH4Intersector1Full<Triangle4Intersector1MoellerTrumbore>);

#if defined(__AVX__)
    typedef Select2Intersector1<Triangle4Intersector1MoellerTrumbore, Bezier1Intersector1> Triangle4Bezier1Intersector1;
    DEFINE_INTERSECTOR1(BVH4Triangle4Bezier1Intersector1,BVH4Intersector1Full<Triangle4Bezier1Intersector1>);
    //DEFINE_INTERSECTOR1(BVH4Triangle4Bezier1Intersector1,BVH4Intersector1Full<Bezier1Intersector1>);

    DEFINE_INTERSECTOR1(BVH4Bezier1Intersector1,BVH4Intersector1Full<Bezier1Intersector1>);
    DEFINE_INTERSECTOR1(BVH4Bezier1iIntersector1,BVH4Intersector1Full<Bezier1iIntersector1>);
#endif
  }
}
