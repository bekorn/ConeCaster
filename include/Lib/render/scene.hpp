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


        // My way, optimized a tiny bit
        auto inv_dir = 1.0f / ray.dir;
        auto t0 = (min - ray.pos) * inv_dir;
        auto t1 = (max - ray.pos) * inv_dir;
        auto t_min = glm::min(t0, t1);
        auto t_max = glm::max(t0, t1);
        auto range_min = compMax(f32x4{t_min, range.min});
        auto range_max = compMin(f32x4{t_max, range.max});
        return (range_max > range_min);
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
    AABB volume;
    __m256 min_x, min_y, min_z;
    __m256 max_x, max_y, max_z;

    void create(array<AABB, 8> const & aabbs)
    {
        auto _min = aabbs[0].min, _max = aabbs[0].max;
        for (auto i = 1; i < 8; ++i)
            _min = glm::min(_min, aabbs[i].min), _max = glm::max(_max, aabbs[i].max);
        volume = AABB{.min = _min, .max = _max};

        min_x = _mm256_set_ps(aabbs[0].min.x, aabbs[1].min.x, aabbs[2].min.x, aabbs[3].min.x, aabbs[4].min.x, aabbs[5].min.x, aabbs[6].min.x, aabbs[7].min.x);
        min_y = _mm256_set_ps(aabbs[0].min.y, aabbs[1].min.y, aabbs[2].min.y, aabbs[3].min.y, aabbs[4].min.y, aabbs[5].min.y, aabbs[6].min.y, aabbs[7].min.y);
        min_z = _mm256_set_ps(aabbs[0].min.z, aabbs[1].min.z, aabbs[2].min.z, aabbs[3].min.z, aabbs[4].min.z, aabbs[5].min.z, aabbs[6].min.z, aabbs[7].min.z);

        max_x = _mm256_set_ps(aabbs[0].max.x, aabbs[1].max.x, aabbs[2].max.x, aabbs[3].max.x, aabbs[4].max.x, aabbs[5].max.x, aabbs[6].max.x, aabbs[7].max.x);
        max_y = _mm256_set_ps(aabbs[0].max.y, aabbs[1].max.y, aabbs[2].max.y, aabbs[3].max.y, aabbs[4].max.y, aabbs[5].max.y, aabbs[6].max.y, aabbs[7].max.y);
        max_z = _mm256_set_ps(aabbs[0].max.z, aabbs[1].max.z, aabbs[2].max.z, aabbs[3].max.z, aabbs[4].max.z, aabbs[5].max.z, aabbs[6].max.z, aabbs[7].max.z);
    }

    std::bitset<8> hit(Ray const & ray, Range<f32> const & range) const
    {
        auto inv_dir = 1.f / ray.dir;

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
        
        std::bitset<8> hits;
        auto result = _mm256_castps_si256(_mm256_cmp_ps(range_max, range_min, _CMP_GT_OQ));
        for (auto i = 0; i < 8; ++i)
            hits[i] = result.m256i_u32[7 - i]; // somehow the results are reversed
        return hits;
    }

    AABB aabb() const
    {
        return volume;
    }
};

struct Hit
{
    bool is_hit;

    f32x3 pos;
    f32 dist;
    f32x3 normal;

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
    Hittable * left, * right;
    AABB volume;

    void create(span<AABBHeuristic> const & heuristics, u32 split_axis = 0)
    {
        if (heuristics.size() == 1)
        {
            left = right = heuristics[0].hittable;
        }
        else if (heuristics.size() == 2)
        {
            left = heuristics[0].hittable;
            right = heuristics[1].hittable;
        }
        else
        {
            std::ranges::sort(heuristics, {}, [split_axis](auto & ch){return ch.value[split_axis];});

            auto next_split_axis = (split_axis + 1) % 3;
            auto split = heuristics.size() / 2;

            auto left_volume = new BoundingVolume();
            left_volume->create(heuristics.subspan(0, split), next_split_axis);
            left = left_volume;

            auto right_volume = new BoundingVolume();
            right_volume->create(heuristics.subspan(split), next_split_axis);
            right = right_volume;
        }

        volume = AABB::merge(left->aabb(), right->aabb());
    }

    Hit hit(Ray const & ray, Range<f32> range) const final
    {
        if (not volume.hit(ray, range))
            return {.is_hit = false};

        auto left_hit = left->hit(ray, range);
        if (left_hit.is_hit)
            range.max = left_hit.dist;
        
        auto right_hit = right->hit(ray, range);
        if (right_hit.is_hit) // right hit is closer
            return right_hit;
        else
            return left_hit;
    }

    AABB aabb() const final
    {
        return volume;
    }
};

struct BoundingVolume_8wide final : Hittable
{
    array<Hittable *,8> children;
    AABB_8wide volumes;

    void create(span<AABBHeuristic> const & heuristics)
    {
        if (heuristics.size() <= 8)
        {
            for (usize i = 0; i < 8; ++i)
                children[i] = heuristics[glm::min(i, heuristics.size() - 1)].hittable;
        }
        else
        {
            std::ranges::sort(heuristics, {}, [](auto & ch) {return ch.value.x; });
            auto split_x = heuristics.size() / 2;
            auto x0 = heuristics.subspan(0, split_x);
            auto x1 = heuristics.subspan(split_x);

            std::ranges::sort(x0, {}, [](auto & ch) {return ch.value.y; });
            auto split_x0_y = x0.size() / 2;
            auto x0_y0 = x0.subspan(0, split_x0_y);
            auto x0_y1 = x0.subspan(split_x0_y);
            std::ranges::sort(x1, {}, [](auto & ch) {return ch.value.y; });
            auto split_x1_y = x1.size() / 2;
            auto x1_y0 = x1.subspan(0, split_x1_y);
            auto x1_y1 = x1.subspan(split_x1_y);

            std::ranges::sort(x0_y0, {}, [](auto & ch) {return ch.value.z; });
            auto split_x0_y0_z = x0_y0.size() / 2;
            std::ranges::sort(x0_y1, {}, [](auto & ch) {return ch.value.z; });
            auto split_x0_y1_z = x0_y1.size() / 2;
            std::ranges::sort(x1_y0, {}, [](auto & ch) {return ch.value.z; });
            auto split_x1_y0_z = x1_y0.size() / 2;
            std::ranges::sort(x1_y1, {}, [](auto & ch) {return ch.value.z; });
            auto split_x1_y1_z = x1_y1.size() / 2;

            span<AABBHeuristic> children_groups[8] = {
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
                if (children_groups[i].size() == 1)
                {
                    children[i] = children_groups[i][0].hittable;
                }
                else if (children_groups[i].size() < 6)
                {
                    auto bv = new BoundingVolume();
                    bv->create(children_groups[i]);
                    children[i] = bv;
                }
                else
                {
                    auto bv_8wide = new BoundingVolume_8wide();
                    bv_8wide->create(children_groups[i]);
                    children[i] = bv_8wide;
                }
            }
        }

        array<AABB, 8> child_aabbs;
        for (auto i = 0; i < 8; ++i)
            child_aabbs[i] = children[i]->aabb();
        volumes.create(child_aabbs);
    }
 
    Hit hit(Ray const & ray, Range<f32> range) const final
    {
        auto volume_hits = volumes.hit(ray, range);
        if (volume_hits.none())
            return {.is_hit = false};

        Hit closest_hit{.is_hit = false}; // remove {.is_hit = false} to summon demons
        for (auto i = 0; i < 8; ++i)
            if (volume_hits[i])
            {
                auto hit = children[i]->hit(ray, range);
                if (hit.is_hit)
                {
                    closest_hit = hit;
                    range.max = hit.dist;
                }
            }

        return closest_hit;
    }

    AABB aabb() const final
    {
        return volumes.aabb();
    }
};

struct Sphere final : Hittable
{
    f32x3 pos;
    f32 r;
    f32x3 color;

    Hit hit(Ray const & ray, Range<f32> range) const final
    {
        // |(A + tD) - C| = r solutions for t are the intersections
        auto to_ray = ray.pos - pos;

        // a*t^2 + b*t + c = 0
        auto a = dot(ray.dir, ray.dir);
        auto b = 2 * dot(to_ray, ray.dir);
        auto c = dot(to_ray, to_ray) - r * r;
        auto discriminant = b * b - 4 * a * c;

        if (discriminant < 0) // No Solution
            return { .is_hit = false };

        auto sqrt_discr = glm::sqrt(discriminant);
        f32 distance;
        if (auto root = (-b - sqrt_discr) / (2 * a); range.contains(root))
            distance = root;
        else if (auto root = (-b + sqrt_discr) / (2 * a); range.contains(root))
            distance = root;
        else
            return { .is_hit = false };

        auto hit_pos = ray.at(distance);
        return Hit{
            .is_hit = true,
            .pos = hit_pos,
            .dist = distance,
            .normal = (hit_pos - pos) / r,
            .attenuation = color,
        };
    }

    AABB aabb() const final
    {
        return {
            .min = pos - r,
            .max = pos + r
        };
    }
};

struct HittableList final : Hittable
{
    vector<Hittable *> hittables;

    Hit hit(Ray const & ray, Range<f32> range) const final
    {
        Hit closest{
            .is_hit = false,
        };

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

struct Scene final : Hittable
{
    f32x3 background_color_up{0.63, 0.87, 0.99};
    f32x3 background_color_down{1, 1, 1};
    std::array<Sphere, 1000> spheres;

    HittableList hittable_list;
    BoundingVolume bvh;
    BoundingVolume_8wide bvh_8wide;

    void create()
    {
        for (auto & sphere: spheres)
        {
            sphere.r = Random::next<f32>(5, 10);
            sphere.pos = Random::next(f32x3(-200), f32x3(200));
            // sphere.color = sphere.pos / 512.f;
            sphere.color = Random::next(f32x3(0), f32x3(1));
        }

        spheres[0].r = 100;
        spheres[0].pos = {0, 0, 0};
        spheres[0].color = {0.8, 0.8, 0.1};

        hittable_list.hittables.reserve(spheres.size());
        for (auto & sphere: spheres)
            hittable_list.hittables.push_back(&sphere);

        vector<AABBHeuristic> heuristics;
        heuristics.reserve(hittable_list.hittables.size());
        for (auto & hittable: hittable_list.hittables)
            heuristics.emplace_back(hittable, hittable->aabb().center());

        bvh.create(heuristics);
        bvh_8wide.create(heuristics);
    }

    Hit hit(Ray const & ray, Range<f32> range) const final
    {
        return bvh.hit(ray, range);
    }

    AABB aabb() const final
    {
        return bvh.aabb();
    }

    f32x3 get_background_color(f32x3 const & dir)
    {
        auto t = 0.5f * (dir.z + 1); // map [-1, 1] -> [0, 1]
        return glm::lerp(background_color_down, background_color_up, t);
    }
};

}