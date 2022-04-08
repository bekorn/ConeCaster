#pragma once

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
        //     auto t0 = (min[a] - r.pos[a]) / r.dir[a];
        //     auto t1 = (max[a] - r.pos[a]) / r.dir[a];
        //     auto t_min = fmin(t0, t1);
        //     auto t_max = fmax(t0, t1);
        //     range.min = fmax(t_min, range.min);
        //     range.max = fmin(t_max, range.max);
        // }
        // return (range.max > range.min);

        // Book 2 Chapter 3.5
        // for (int a = 0; a < 3; a++) {
        //     auto invD = 1.0f / r.dir[a];
        //     auto t0 = (min[a] - r.pos[a]) * invD;
        //     auto t1 = (max[a] - r.pos[a]) * invD;
        //     if (invD < 0.0f)
        //         std::swap(t0, t1);
        //     range.min = fmax(t0, range.min);
        //     range.max = fmin(t1, range.max);
        // }
        // return (range.max > range.min);

        //  My way using glm
        auto t0 = (min - ray.pos) / ray.dir;
        auto t1 = (max - ray.pos) / ray.dir;
        auto t_min = glm::min(t0, t1);
        auto t_max = glm::max(t0, t1);

        auto range_min = compMax(f32x4{t_min, range.min});
        auto range_max = compMin(f32x4{t_max, range.max});
        return (range_max > range_min);
    }

    static AABB merge(AABB const & a, AABB const & b)
    {
        return {
            .min = glm::min(a.min, b.min),
            .max = glm::max(a.max, b.max)
        };
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

struct BoundingVolume final : Hittable
{
    Hittable * left, * right;
    AABB volume;

    void create(span<Hittable *> const & hittables, u32 split_axis = 0)
    {
        if (hittables.size() == 1)
        {
            left = right = hittables[0];
        }
        else if (hittables.size() == 2)
        {
            left = hittables[0];
            right = hittables[1];
        }
        else
        {
            struct Comparable
            {
                u32 idx;;
                f32 value;
            };
            vector<Comparable> comparables;
            comparables.reserve(hittables.size());
            for (auto idx = 0; auto const & hittable: hittables)
                comparables.emplace_back(idx++, hittable->aabb().min[split_axis]);
            std::ranges::sort(comparables, {}, &Comparable::value);

            vector<Hittable *> sorted_hittables;
            sorted_hittables.reserve(hittables.size());
            for (auto const & comparable: comparables)
                sorted_hittables.emplace_back(hittables[comparable.idx]);

            auto next_split_axis = split_axis + 1 % 3;

            auto left_size = hittables.size() / 2;
            auto left_volume = new BoundingVolume();
            left_volume->create({sorted_hittables.begin(), left_size}, next_split_axis);
            left = left_volume;

            auto right_size = hittables.size() - left_size;
            auto right_volume = new BoundingVolume();
            right_volume->create({sorted_hittables.begin() + left_size, right_size}, next_split_axis);
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

    AABB aabb() const  final
    {
        return volume;
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
    std::array<Sphere, 60> spheres;

    HittableList hittable_list;
    BoundingVolume bvh;

    void create()
    {
        for (auto & sphere: spheres)
        {
            sphere.r = Random::next<f32>(20, 40);
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

        bvh.create(hittable_list.hittables);
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