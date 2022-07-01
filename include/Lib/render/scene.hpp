#pragma once

#include <immintrin.h>
#include <bitset>

#include "Lib/core/.hpp"

#include "core.hpp"

namespace Render
{

struct AABB
{
    f32x3 min, max;

    bool hit(Ray const & ray) const
    {
        // TODO(bekorn): benchmark each method
        // https://godbolt.org/z/1xTx6oa4f

        // Book 2 Chapter 3.4
        // for (int a = 0; a < 3; a++)
        // {
        //     auto t0 = (min[a] - ray.pos[a]) / ray.dir[a];
        //     auto t1 = (max[a] - ray.pos[a]) / ray.dir[a];
        //     auto t_min = fmin(t0, t1);
        //     auto t_max = fmax(t0, t1);
        //     range.min = fmax(t_min, range.min);
        //     range.max = fmin(t_max, range.max);
        // }
        // return (range.max > range.min);


        // Book 2 Chapter 3.5
        // for (int a = 0; a < 3; a++) {
        //     auto invD = 1.0f / ray.dir[a];
        //     auto t0 = (min[a] - ray.pos[a]) * invD;
        //     auto t1 = (max[a] - ray.pos[a]) * invD;
        //     if (invD < 0.0f)
        //         std::swap(t0, t1);
        //     range.min = fmax(t0, range.min);
        //     range.max = fmin(t1, range.max);
        // }
        // return (range.max > range.min);


        //  My way using glm
        // auto t0 = (min - ray.pos) / ray.dir;
        // auto t1 = (max - ray.pos) / ray.dir;
        // auto t_min = glm::min(t0, t1);
        // auto t_max = glm::max(t0, t1);
        // auto range_min = compMax(f32x4{t_min, range.min});
        // auto range_max = compMin(f32x4{t_max, range.max});
        // return (range_max > range_min);


        // My way v2
        auto inv_dir = 1.0f / ray.dir;
        auto t0 = (min - ray.pos) * inv_dir;
        auto t1 = (max - ray.pos) * inv_dir;
        auto t_min = glm::min(t0, t1);
        auto t_max = glm::max(t0, t1);
        auto range_min = compMax(f32x4{t_min, ray.min});
        auto range_max = compMin(f32x4{t_max, ray.max});
        return (range_max >= range_min);
    }

    f32x3 center() const
    { return (min + max) * 0.5f; }

    static AABB merge(AABB const & a, AABB const & b)
    {
        return {
            .min = glm::min(a.min, b.min),
            .max = glm::max(a.max, b.max)
        };
    }

    static constexpr AABB identity()
    {
        return {
            .min = f32x3(+std::numeric_limits<f32>::max()),
            .max = f32x3(-std::numeric_limits<f32>::max()),
        };
    }
};

struct AABB_8wide
{
    AABB aabbs[8];

    std::bitset<8> hit(Ray const & ray) const
    {
        auto inv_dir = 1.f / ray.dir;

        // min and max values are reversed so the outcome is in correct order
        auto min_x = _mm256_set_ps(aabbs[7].min.x, aabbs[6].min.x, aabbs[5].min.x, aabbs[4].min.x, aabbs[3].min.x, aabbs[2].min.x, aabbs[1].min.x, aabbs[0].min.x);
        auto min_y = _mm256_set_ps(aabbs[7].min.y, aabbs[6].min.y, aabbs[5].min.y, aabbs[4].min.y, aabbs[3].min.y, aabbs[2].min.y, aabbs[1].min.y, aabbs[0].min.y);
        auto min_z = _mm256_set_ps(aabbs[7].min.z, aabbs[6].min.z, aabbs[5].min.z, aabbs[4].min.z, aabbs[3].min.z, aabbs[2].min.z, aabbs[1].min.z, aabbs[0].min.z);

        auto max_x = _mm256_set_ps(aabbs[7].max.x, aabbs[6].max.x, aabbs[5].max.x, aabbs[4].max.x, aabbs[3].max.x, aabbs[2].max.x, aabbs[1].max.x, aabbs[0].max.x);
        auto max_y = _mm256_set_ps(aabbs[7].max.y, aabbs[6].max.y, aabbs[5].max.y, aabbs[4].max.y, aabbs[3].max.y, aabbs[2].max.y, aabbs[1].max.y, aabbs[0].max.y);
        auto max_z = _mm256_set_ps(aabbs[7].max.z, aabbs[6].max.z, aabbs[5].max.z, aabbs[4].max.z, aabbs[3].max.z, aabbs[2].max.z, aabbs[1].max.z, aabbs[0].max.z);

        auto t0_x = _mm256_mul_ps(_mm256_sub_ps(min_x, _mm256_set1_ps(ray.pos.x)), _mm256_set1_ps(inv_dir.x));
        auto t0_y = _mm256_mul_ps(_mm256_sub_ps(min_y, _mm256_set1_ps(ray.pos.y)), _mm256_set1_ps(inv_dir.y));
        auto t0_z = _mm256_mul_ps(_mm256_sub_ps(min_z, _mm256_set1_ps(ray.pos.z)), _mm256_set1_ps(inv_dir.z));

        auto t1_x = _mm256_mul_ps(_mm256_sub_ps(max_x, _mm256_set1_ps(ray.pos.x)), _mm256_set1_ps(inv_dir.x));
        auto t1_y = _mm256_mul_ps(_mm256_sub_ps(max_y, _mm256_set1_ps(ray.pos.y)), _mm256_set1_ps(inv_dir.y));
        auto t1_z = _mm256_mul_ps(_mm256_sub_ps(max_z, _mm256_set1_ps(ray.pos.z)), _mm256_set1_ps(inv_dir.z));

        auto t_min_x = _mm256_min_ps(t0_x, t1_x);
        auto t_min_y = _mm256_min_ps(t0_y, t1_y);
        auto t_min_z = _mm256_min_ps(t0_z, t1_z);

        auto t_max_x = _mm256_max_ps(t0_x, t1_x);
        auto t_max_y = _mm256_max_ps(t0_y, t1_y);
        auto t_max_z = _mm256_max_ps(t0_z, t1_z);

        auto range_min = _mm256_max_ps(_mm256_max_ps(_mm256_max_ps(t_min_x, t_min_y), t_min_z), _mm256_set1_ps(ray.min));
        auto range_max = _mm256_min_ps(_mm256_min_ps(_mm256_min_ps(t_max_x, t_max_y), t_max_z), _mm256_set1_ps(ray.max));

        u8 results = _mm256_movemask_ps(_mm256_cmp_ps(range_max, range_min, _CMP_GE_OQ));
        return {results};
    }

    struct HitsSorted
    {
        u8 hittable_idx[8];
        u8 hit_count;
    };

    HitsSorted hit_sorted(Ray const & ray) const
    {
        auto inv_dir = 1.f / ray.dir;

        // min and max arder is reversed so the outcome is in correct order
        auto min_x = _mm256_set_ps(aabbs[7].min.x, aabbs[6].min.x, aabbs[5].min.x, aabbs[4].min.x, aabbs[3].min.x, aabbs[2].min.x, aabbs[1].min.x, aabbs[0].min.x);
        auto min_y = _mm256_set_ps(aabbs[7].min.y, aabbs[6].min.y, aabbs[5].min.y, aabbs[4].min.y, aabbs[3].min.y, aabbs[2].min.y, aabbs[1].min.y, aabbs[0].min.y);
        auto min_z = _mm256_set_ps(aabbs[7].min.z, aabbs[6].min.z, aabbs[5].min.z, aabbs[4].min.z, aabbs[3].min.z, aabbs[2].min.z, aabbs[1].min.z, aabbs[0].min.z);
 
        auto max_x = _mm256_set_ps(aabbs[7].max.x, aabbs[6].max.x, aabbs[5].max.x, aabbs[4].max.x, aabbs[3].max.x, aabbs[2].max.x, aabbs[1].max.x, aabbs[0].max.x);
        auto max_y = _mm256_set_ps(aabbs[7].max.y, aabbs[6].max.y, aabbs[5].max.y, aabbs[4].max.y, aabbs[3].max.y, aabbs[2].max.y, aabbs[1].max.y, aabbs[0].max.y);
        auto max_z = _mm256_set_ps(aabbs[7].max.z, aabbs[6].max.z, aabbs[5].max.z, aabbs[4].max.z, aabbs[3].max.z, aabbs[2].max.z, aabbs[1].max.z, aabbs[0].max.z);

        auto t0_x = _mm256_mul_ps(_mm256_sub_ps(min_x, _mm256_set1_ps(ray.pos.x)), _mm256_set1_ps(inv_dir.x));
        auto t0_y = _mm256_mul_ps(_mm256_sub_ps(min_y, _mm256_set1_ps(ray.pos.y)), _mm256_set1_ps(inv_dir.y));
        auto t0_z = _mm256_mul_ps(_mm256_sub_ps(min_z, _mm256_set1_ps(ray.pos.z)), _mm256_set1_ps(inv_dir.z));

        auto t1_x = _mm256_mul_ps(_mm256_sub_ps(max_x, _mm256_set1_ps(ray.pos.x)), _mm256_set1_ps(inv_dir.x));
        auto t1_y = _mm256_mul_ps(_mm256_sub_ps(max_y, _mm256_set1_ps(ray.pos.y)), _mm256_set1_ps(inv_dir.y));
        auto t1_z = _mm256_mul_ps(_mm256_sub_ps(max_z, _mm256_set1_ps(ray.pos.z)), _mm256_set1_ps(inv_dir.z));

        auto t_min_x = _mm256_min_ps(t0_x, t1_x);
        auto t_min_y = _mm256_min_ps(t0_y, t1_y);
        auto t_min_z = _mm256_min_ps(t0_z, t1_z);

        auto t_max_x = _mm256_max_ps(t0_x, t1_x);
        auto t_max_y = _mm256_max_ps(t0_y, t1_y);
        auto t_max_z = _mm256_max_ps(t0_z, t1_z);

        auto range_min = _mm256_max_ps(_mm256_max_ps(_mm256_max_ps(t_min_x, t_min_y), t_min_z), _mm256_set1_ps(ray.min));
        auto range_max = _mm256_min_ps(_mm256_min_ps(_mm256_min_ps(t_max_x, t_max_y), t_max_z), _mm256_set1_ps(ray.max));

        f32 range_mins[8];
        _mm256_storeu_ps(range_mins, range_min);

        u32 is_hits[8];
        _mm256_storeu_ps(reinterpret_cast<f32*>(&is_hits), _mm256_cmp_ps(range_max, range_min, _CMP_GE_OQ));

        HitsSorted hits{
            .hit_count = 0
        };

        //  only store idxs which is_hit[idx] == true
        for (auto i = 0; i < 8; ++i)
            if (is_hits[i])
                hits.hittable_idx[hits.hit_count++] = i;

        std::sort(
            hits.hittable_idx, hits.hittable_idx + hits.hit_count,
            [&range_mins](auto idx1, auto idx2){ return range_mins[idx1] < range_mins[idx2]; }
        );

        return hits;
    }
};

struct Hit
{
    bool is_hit;

    f32x3 pos;
    f32 dist;
    f32x3 normal;
    f32x3 barycentric;

    f32x3 attenuation;
};

struct Hittable
{
    virtual Hit hit(Ray &) const = 0;
    virtual AABB aabb() const  = 0;
};

struct AABBHeuristic
{
    Hittable * hittable;
    f32x3 value;
};

struct BoundingVolume final : Hittable, CTOR_Counter<BoundingVolume>
{
    array<Hittable *, 2> hittables;
    array<AABB, 2> aabbs;

    void create(span<AABBHeuristic> const & heuristics, u32 split_axis = 0)
    {
        if (heuristics.size() == 1)
        {
            assert(("please do not create a bounding box with 1 hittable", false));
        }
        else if (heuristics.size() == 2)
        {
            hittables[0] = heuristics[0].hittable;
            hittables[1] = heuristics[1].hittable;
        }
        else
        {
            std::ranges::sort(heuristics, {}, [split_axis](auto & ch){return ch.value[split_axis];});

            auto next_split_axis = (split_axis + 1) % 3;
            auto split = heuristics.size() / 2;
            span<AABBHeuristic> groups[2] = {
                heuristics.subspan(0, split),
                heuristics.subspan(split)
            };

            for (auto i = 0; i < 2; ++i)
            {
                if (groups[i].size() == 1)
                    hittables[i] = groups[i][0].hittable;
                else
                {
                    auto bv = new BoundingVolume();
                    bv->create(groups[i], next_split_axis);
                    hittables[i] = bv;
                }
            }
        }

        aabbs[0] = hittables[0]->aabb();
        aabbs[1] = hittables[1]->aabb();
    }

    Hit hit(Ray & ray) const final
    {
        Hit closest{.is_hit = false};

        for (auto i = 0; i < 2; ++i)
            if (aabbs[i].hit(ray))
            {
                auto hit = hittables[i]->hit(ray);
                if (hit.is_hit)
                {
                    closest = hit;
                    ray.max = hit.dist;
                }
            }

        return closest;
    }

    AABB aabb() const final
    { return AABB::merge(aabbs[0], aabbs[1]); }
};

struct BoundingVolume_8wide final : Hittable, CTOR_Counter<BoundingVolume_8wide>
{
    array<Hittable *, 8> hittables;
    AABB_8wide aabbs;

    void create(span<AABBHeuristic> const & heuristics)
    {
        if (heuristics.size() <= 8)
        {
            for (usize i = 0; i < 8; ++i)
                hittables[i] = heuristics[glm::min(i, heuristics.size() - 1)].hittable;
        }
        else
        {
            std::ranges::sort(heuristics, {}, [](auto & h) {return h.value.x; });
            auto split_x = heuristics.size() / 2;
            auto x0 = heuristics.subspan(0, split_x);
            auto x1 = heuristics.subspan(split_x);

            std::ranges::sort(x0, {}, [](auto & h) {return h.value.y; });
            auto split_x0_y = x0.size() / 2;
            auto x0_y0 = x0.subspan(0, split_x0_y);
            auto x0_y1 = x0.subspan(split_x0_y);
            std::ranges::sort(x1, {}, [](auto & h) {return h.value.y; });
            auto split_x1_y = x1.size() / 2;
            auto x1_y0 = x1.subspan(0, split_x1_y);
            auto x1_y1 = x1.subspan(split_x1_y);

            std::ranges::sort(x0_y0, {}, [](auto & h) {return h.value.z; });
            auto split_x0_y0_z = x0_y0.size() / 2;
            std::ranges::sort(x0_y1, {}, [](auto & h) {return h.value.z; });
            auto split_x0_y1_z = x0_y1.size() / 2;
            std::ranges::sort(x1_y0, {}, [](auto & h) {return h.value.z; });
            auto split_x1_y0_z = x1_y0.size() / 2;
            std::ranges::sort(x1_y1, {}, [](auto & h) {return h.value.z; });
            auto split_x1_y1_z = x1_y1.size() / 2;

            span<AABBHeuristic> groups[8] = {
                x0_y0.subspan(0, split_x0_y0_z),
                x0_y0.subspan(split_x0_y0_z),
                x0_y1.subspan(0, split_x0_y1_z),
                x0_y1.subspan(split_x0_y1_z),
                x1_y0.subspan(0, split_x1_y0_z),
                x1_y0.subspan(split_x1_y0_z),
                x1_y1.subspan(0, split_x1_y1_z),
                x1_y1.subspan(split_x1_y1_z),
            };

            for (auto i = 0; i < 8; ++i)
            {
                if (groups[i].size() == 1)
                {
                    hittables[i] = groups[i][0].hittable;
                }
                else if (groups[i].size() <= 6)
                {
                    auto bv = new BoundingVolume();
                    bv->create(groups[i]);
                    hittables[i] = bv;
                }
                else
                {
                    auto bv_8wide = new BoundingVolume_8wide();
                    bv_8wide->create(groups[i]);
                    hittables[i] = bv_8wide;
                }
            }
        }

        for (auto i = 0; i < 8; ++i)
            aabbs.aabbs[i] = hittables[i]->aabb();
    }

    Hit hit(Ray & ray) const final
    {
        Hit closest{.is_hit = false};

        if constexpr (true) // should use aabb hit distance order
        {
            auto aabb_hits = aabbs.hit_sorted(ray);
            for (auto i = 0; i < aabb_hits.hit_count; ++i)
            {
                auto hit = hittables[aabb_hits.hittable_idx[i]]->hit(ray);
                if (hit.is_hit)
                {
                    closest = hit;
                    ray.max = hit.dist;
                }
            }
        }
        else
        {
            auto aabb_hits = aabbs.hit(ray);
            for (auto i = 0; i < 8; ++i)
                if (aabb_hits[i])
                {
                    auto hit = hittables[i]->hit(ray);
                    if (hit.is_hit)
                    {
                        closest = hit;
                        ray.max = hit.dist;
                    }
                }
        }

        return closest;
    }

    AABB aabb() const final
    {
        auto aabb = aabbs.aabbs[0];
        for (auto i = 1; i < 8; ++i)
            aabb = AABB::merge(aabb, aabbs.aabbs[i]);
        return aabb;
    }
};

struct Triangle final : Hittable
{
    f32x3 vert[3];

    // https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/moller-trumbore-ray-triangle-intersection
    Hit hit(Ray & ray) const
    {
        auto edge_0to1 = vert[1] - vert[0];
        auto edge_0to2 = vert[2] - vert[0];

        auto _k = cross(ray.dir, edge_0to2);
        auto determinant = dot(edge_0to1, _k);

#if true // 2-sided
        if (glm::abs(determinant) < glm::epsilon<f32>())
            return {.is_hit = false};
#else
        if (determinant < glm::epsilon<f32>())
            return {.is_hit = false};
#endif

        auto _m = ray.pos - vert[0];
        auto u = dot(_m, _k) / determinant;
        if (u < 0 | u > 1)
            return {.is_hit = false};
        
        auto _n = cross(_m, edge_0to1);
        auto v = dot(ray.dir, _n) / determinant;
        if (v < 0 | u + v > 1)
            return {.is_hit = false};
        
        auto t = dot(edge_0to2, _n) / determinant;
        if ((t < ray.min) | (ray.max < t))
            return {.is_hit = false};

        auto pos = ray.at(t);
        auto normal = normalize(cross(edge_0to1, edge_0to2));
#if true // 2-sided
        normal *= glm::sign(determinant);
#endif
        auto barycentric = f32x3{1 - u - v, u, v};
        return {
            .is_hit = true,
            .pos = pos,
            .dist = t,
            .normal = normal,
            .barycentric = barycentric,
            // .attenuation = normal,
            .attenuation = barycentric,
            // .attenuation = {0.2, 0.4, 1},
        };
    }

    AABB aabb() const 
    {
        return {
            .min = min(min(vert[0], vert[1]), vert[2]),
            .max = max(max(vert[0], vert[1]), vert[2])
        };
    }
};

// JACCO BVH BEGIN
struct JACCO_BVH_Node
{
    AABB aabb;      // 24 bytes
    union           //  8 bytes            
    {
        // braanch node
        struct { u32 left_child; }; // right_child is always left_child + 1
        // leaf node
        struct { u32 tri_offset, tri_count; };
    };

    bool is_leaf() const
    { return tri_count != 0; }
};

struct JACCO_BVH final : Hittable
{
    vector<Triangle> tris;
    vector<JACCO_BVH_Node> bvh_pool;

    Hit hit(Ray & ray) const
    {
        Hit closest{.is_hit = false};

        vector<u32> possible_idxs;
        possible_idxs.reserve(log2(bvh_pool.size()));
        possible_idxs.push_back(0);

        while (not possible_idxs.empty())
        {
            auto node_idx = possible_idxs.back();
            possible_idxs.pop_back();

            auto & node = bvh_pool[node_idx];

            if (node.aabb.hit(ray))
            {
                if (node.is_leaf())
                {
                    for (auto i = 0; i < node.tri_count; i++)
                    {
                        auto hit = tris[node.tri_offset + i].hit(ray);
                        if (hit.is_hit)
                        {
                            closest = hit;
                            ray.max = hit.dist;
                        }
                    }
                }
                else
                {
                    possible_idxs.push_back(node.left_child);
                    possible_idxs.push_back(node.left_child + 1);
                }
            }
        }

        return closest;
    }

    AABB aabb() const
    { return bvh_pool[0].aabb; }
};

struct JACCO_BVH_Builder
{
    struct PrepData
    {
        AABB aabb;
        f32x3 centroid;
    };

    unique_array<u32> indices;
    unique_array<PrepData> prepdata;
    vector<JACCO_BVH_Node> bvh_pool;

    JACCO_BVH create(vector<Triangle> const & triangles)
    {
        fmt::print("Creating a JACCO_BVH with {} triangles\n", triangles.size());

        // preperation
        auto tri_count = triangles.size();

        indices = make_unique_array<u32>(tri_count);
        for (auto i = 0; i < tri_count; i++)
            indices[i] = i;

        prepdata = make_unique_array<PrepData>(tri_count);
        for (auto i = 0; i < tri_count; i++)
        {
            auto & t = triangles[i];
            auto & p = prepdata[i];
            p.centroid = (t.vert[0] + t.vert[1] + t.vert[2]) * (1.f / 3.f);
            p.aabb = {
                .min = min(min(t.vert[0], t.vert[1]), t.vert[2]),
                .max = max(max(t.vert[0], t.vert[1]), t.vert[2])
            };
        }

        //  bvh creation
        bvh_pool.reserve(2 * tri_count - 1);

        // add root node
        bvh_pool.push_back({
            .tri_offset = 0,
            .tri_count = u32(tri_count),
        });
        subdivide(0);

        // create bvh object
        JACCO_BVH bvh;

        bvh.tris.resize(tri_count);
        for (auto i = 0; i < tri_count; i++)
            bvh.tris[i] = triangles[indices[i]];

        bvh.bvh_pool = move(bvh_pool);

        return bvh;
    }

    void subdivide(u32 node_idx)
    {
        auto & node = bvh_pool[node_idx];

        // update bounds
        node.aabb = AABB::identity();
        for (auto i = 0; i < node.tri_count; i++)
            node.aabb = AABB::merge(node.aabb, prepdata[indices[node.tri_offset + i]].aabb);

        // subdivide
        if (node.tri_count <= 2)
            return;

        auto extents = node.aabb.max - node.aabb.min;
        auto split_axis = 0;
        if (extents.y > extents.x)
            split_axis = 1;
        if (extents.z > extents[split_axis])
            split_axis = 2;

        auto split_pos = node.aabb.min[split_axis] + 0.5f * extents[split_axis];

        u32 split_idx = std::partition(
            indices.get() + node.tri_offset,
            indices.get() + node.tri_offset + node.tri_count,
            [&](auto & i){ return prepdata[i].centroid[split_axis] < split_pos; }
        ) - indices.get();

        auto left_count = split_idx - node.tri_offset;

        if (left_count == 0 || left_count == node.tri_count)
            return;

        // fmt::print(
        //     "axis: {}, offset: {:_>4}, count: {:_>4}, split: {:_>4}, pos{{ min: {: 7.2f} | max: {: 7.2f} | split: {: 7.2f} }}\n",
        //     split_axis, node.tri_offset, node.tri_count, split_idx,
        //     node.aabb.min[split_axis], node.aabb.max[split_axis], split_pos
        // );

        auto left_child = bvh_pool.size();

        bvh_pool.push_back({
            .tri_offset = node.tri_offset,
            .tri_count = left_count,
        });
        bvh_pool.push_back({
            .tri_offset = split_idx,
            .tri_count = node.tri_count - left_count,
        });

        node.left_child = left_child;
        node.tri_count = 0;

        subdivide(node.left_child);
        subdivide(node.left_child + 1);
    }
};
// JACCO BVH END

struct Scene
{
    f32x3 background_color_up{0.63, 0.87, 0.99};
    f32x3 background_color_down{1, 1, 1};

    vector<Triangle> triangles;

    BoundingVolume bvh;
    BoundingVolume_8wide bvh_8wide;
    JACCO_BVH jacco_bvh;

    void create()
    {
        Timer timer;

        vector<AABBHeuristic> heuristics;
        heuristics.reserve(triangles.size());
        for (auto & triangle: triangles)
            heuristics.emplace_back((Hittable*)&triangle, triangle.aabb().center());

        timer.timeit(stderr, "BVH Heuristic calculated");

        BoundingVolume::ctor_count = BoundingVolume_8wide::ctor_count = 0;
        bvh.create(heuristics);
        timer.timeit(stderr, "BVH generated");
        fmt::print(
            "BoundingVolume: {} | BoundingVolume_8wide: {} | Mem: {}\n",
            BoundingVolume::ctor_count, BoundingVolume_8wide::ctor_count,
            BoundingVolume::ctor_count * sizeof(BoundingVolume) + BoundingVolume_8wide::ctor_count * sizeof(BoundingVolume_8wide) 
        );

        BoundingVolume::ctor_count = BoundingVolume_8wide::ctor_count = 0;
        bvh_8wide.create(heuristics);
        timer.timeit(stderr, "BVH_8wide generated");
        fmt::print(
            "BoundingVolume: {:L} | BoundingVolume_8wide: {:L} | Mem: {:L}\n",
            BoundingVolume::ctor_count, BoundingVolume_8wide::ctor_count,
            BoundingVolume::ctor_count * sizeof(BoundingVolume) + BoundingVolume_8wide::ctor_count * sizeof(BoundingVolume_8wide) 
        );

        jacco_bvh = JACCO_BVH_Builder().create(triangles);
        timer.timeit(stderr, "JACCO_BVH generated");
        fmt::print(
            "JACCO_BVH_Node: {:Ld} | Mem: {:Ld}\n",
            jacco_bvh.bvh_pool.size(), sizeof(JACCO_BVH_Node) * jacco_bvh.bvh_pool.size()
        );
    }

    f32x3 get_background_color(f32x3 const & dir)
    {
        auto t = 0.5f * (dir.z + 1); // map [-1, 1] -> [0, 1]
        return glm::lerp(background_color_down, background_color_up, t);
    }
};

}