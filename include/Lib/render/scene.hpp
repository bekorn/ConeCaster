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

    bool hit(Ray const & ray, Range<f32> const & range) const
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
        auto range_min = compMax(f32x4{t_min, range.min});
        auto range_max = compMin(f32x4{t_max, range.max});
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
};

struct AABB_8wide
{
    AABB aabbs[8];

    std::bitset<8> hit(Ray const & ray, Range<f32> const & range) const
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

        auto range_min = _mm256_max_ps(_mm256_max_ps(_mm256_max_ps(t_min_x, t_min_y), t_min_z), _mm256_set1_ps(range.min));
        auto range_max = _mm256_min_ps(_mm256_min_ps(_mm256_min_ps(t_max_x, t_max_y), t_max_z), _mm256_set1_ps(range.max));

        u8 results = _mm256_movemask_ps(_mm256_cmp_ps(range_max, range_min, _CMP_GE_OQ));
        return {results};
    }

    struct HitsSorted
    {
        u8 hittable_idx[8];
        u8 hit_count;
    };

    HitsSorted hit_sorted(Ray const & ray, Range<f32> const & range) const
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

        auto range_min = _mm256_max_ps(_mm256_max_ps(_mm256_max_ps(t_min_x, t_min_y), t_min_z), _mm256_set1_ps(range.min));
        auto range_max = _mm256_min_ps(_mm256_min_ps(_mm256_min_ps(t_max_x, t_max_y), t_max_z), _mm256_set1_ps(range.max));

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
    virtual Hit hit(Ray const &, Range<f32>) const = 0;
    virtual AABB aabb() const  = 0;
};

struct AABBHeuristic
{
    Hittable * hittable;
    f32x3 value;
};

struct BoundingVolume final : Hittable
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

    Hit hit(Ray const & ray, Range<f32> range) const final
    {
        Hit closest{.is_hit = false};

        for (auto i = 0; i < 2; ++i)
            if (aabbs[i].hit(ray, range))
            {
                auto hit = hittables[i]->hit(ray, range);
                if (hit.is_hit)
                {
                    closest = hit;
                    range.max = hit.dist;
                }
            }

        return closest;
    }

    AABB aabb() const final
    { return AABB::merge(aabbs[0], aabbs[1]); }
};

struct BoundingVolume_8wide final : Hittable
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
                else if (groups[i].size() < 6)
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

    Hit hit(Ray const & ray, Range<f32> range) const final
    {
        Hit closest{.is_hit = false};

        if constexpr (true) // should use aabb hit distance order
        {
            auto aabb_hits = aabbs.hit_sorted(ray, range);
            for (auto i = 0; i < aabb_hits.hit_count; ++i)
            {
                auto hit = hittables[aabb_hits.hittable_idx[i]]->hit(ray, range);
                if (hit.is_hit)
                {
                    closest = hit;
                    range.max = hit.dist;
                }
            }
        }
        else
        {
            auto aabb_hits = aabbs.hit(ray, range);
            for (auto i = 0; i < 8; ++i)
                if (aabb_hits[i])
                {
                    auto hit = hittables[i]->hit(ray, range);
                    if (hit.is_hit)
                    {
                        closest = hit;
                        range.max = hit.dist;
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

    Hit hit(Ray const & ray, Range<f32> range) const
    { return hit_mt(ray, range); }

    // https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/moller-trumbore-ray-triangle-intersection
    Hit hit_mt(Ray const & ray, Range<f32> range) const
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
        if (not range.contains(t))
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

    // https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/barycentric-coordinates
    Hit hit_barycentrics(Ray const & ray, Range<f32> range) const
    {
        // triangle's plane: dot(hit.pos, normal) + K = 0
        auto normal = cross(vert[2] - vert[0], vert[1] - vert[0]);
        auto area2 = length(normal); // length(cross product) yields area * 2
        normal = normalize(normal);
        auto K = -dot(vert[0], normal);

        // ray: ray.pos + ray.dir * t = hit.pos;
        // intersection: t = -(K + dot(ray.pos, normal)) / dot(ray.pos, normal)
        auto ray_dir_dot_normal = dot(ray.dir, normal);
        if (glm::abs(ray_dir_dot_normal) < glm::epsilon<f32>()) // ray and plane are parallel
            return {.is_hit = false};

        auto t = -(K + dot(ray.pos, normal)) / ray_dir_dot_normal;
        if (not range.contains(t))
            return {.is_hit = false};

        auto pos = ray.at(t);

        // calculate barycentric coordinates
        auto cross_0to1 = cross(pos - vert[0], vert[1] - vert[0]);
        auto cross_1to2 = cross(pos - vert[1], vert[2] - vert[1]);
        auto cross_2to0 = cross(pos - vert[2], vert[0] - vert[2]);
        auto barycentric = f32x3{
            length(cross_1to2) / area2,
            length(cross_2to0) / area2,
            0
        };
        barycentric.z = 1.f - barycentric.x - barycentric.y;

        // test if pos is inside the triangle
        auto on_left_of_edge01 = dot(normal, cross_0to1) >= 0.f;
        auto on_left_of_edge12 = dot(normal, cross_1to2) >= 0.f;
        auto on_left_of_edge20 = dot(normal, cross_2to0) >= 0.f;
        auto is_inside = on_left_of_edge01 & on_left_of_edge12 & on_left_of_edge20;

        if (not is_inside)
            return {.is_hit = false};

        return {
            .is_hit = true,
            .pos = pos,
            .dist = t,
            .normal = normal,
            .attenuation = barycentric,
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

struct HittableVector final : Hittable
{
    vector<Hittable *> hittables;

    Hit hit(Ray const & ray, Range<f32> range) const final
    {
        Hit closest{.is_hit = false};

        for (auto const & hittable: hittables)
        {
            auto hit = hittable->hit(ray, range);
            if (hit.is_hit)
            {
                closest = hit;
                range.max = hit.dist;
            }
        }

        return closest;
    }

    AABB aabb() const final
    {
        assert(not hittables.empty());

        auto aabb = hittables[0]->aabb();

        for (auto i = 1; i < hittables.size(); ++i)
            AABB::merge(aabb, hittables[i]->aabb());

        return aabb;
    }
};

struct Scene
{
    f32x3 background_color_up{0.63, 0.87, 0.99};
    f32x3 background_color_down{1, 1, 1};
    array<Triangle, 1000> tris;

    HittableVector hittables;
    BoundingVolume bvh;
    BoundingVolume_8wide bvh_8wide;

    void create()
    {
        for (auto & tri: tris)
        {
            auto center = Random::next<f32x3>(-200, 200);
            auto size = Random::next<f32>(30, 50);
            tri.vert[0] = center;
            tri.vert[1] = center + Random::next<f32x3>(-size, size);
            tri.vert[2] = center + Random::next<f32x3>(-size, size);
        }

        {
            f32x3 vertices[] = {
                {-120, 0, 0},
                {200, 0, 0},
                {0, 0, 100},
            };
            auto & tri = tris[0];
            for (auto i = 0; i < 3; ++i)
                tri.vert[i] = vertices[i];
        }

        hittables.hittables.reserve(tris.size());
        for (auto & tri: tris)
            hittables.hittables.push_back((Hittable*)&tri);

        Timer timer;

        vector<AABBHeuristic> heuristics;
        heuristics.reserve(hittables.hittables.size());
        for (auto & hittable: hittables.hittables)
            heuristics.emplace_back(hittable, hittable->aabb().center());

        timer.timeit(stderr, "BVH Heuristic calculated");

        bvh.create(heuristics);
        timer.timeit(stderr, "BVH generated");

        bvh_8wide.create(heuristics);
        timer.timeit(stderr, "BVH_8wide generated");
    }

    f32x3 get_background_color(f32x3 const & dir)
    {
        auto t = 0.5f * (dir.z + 1); // map [-1, 1] -> [0, 1]
        return glm::lerp(background_color_down, background_color_up, t);
    }
};

}