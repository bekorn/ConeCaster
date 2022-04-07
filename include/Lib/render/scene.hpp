#pragma once

#include "Lib/core/.hpp"

#include "core.hpp"

namespace Render
{

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
};

struct Scene final : Hittable
{
    f32x3 background_color_up{0.63, 0.87, 0.99};
    f32x3 background_color_down{1, 1, 1};
    std::array<Sphere, 60> spheres;

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
    }

    Hit hit(Ray const & ray, Range<f32> range) const final
    {
        Hit closest{
            .is_hit = false,
        };

        for (auto & sphere: spheres)
        {
            auto hit = sphere.hit(ray, range);
            if (hit.is_hit)
            {
                closest = hit;
                range.max = hit.dist;
            }
        }

        return closest;
    }

    f32x3 get_background_color(f32x3 const & dir)
    {
        auto t = 0.5f * (dir.y + 1); // map [-1, 1] -> [0, 1]
        return glm::lerp(background_color_down, background_color_up, t);
    }
};

}