#pragma once

#include "Lib/core/.hpp"
#include "Lib/opengl/framebuffer.hpp"

#include "core.hpp"
#include "scene.hpp"

namespace Render
{
struct Image
{
    u32x2 dimensions;
    u32 pixel_count;
    unique_array<f32x3> sample_sum; // color range [0-1]
    unique_array<f32> sample_count;
    unique_array<u8x4> result; // color range [0-255]

    void create(u32x2 dimensions)
    {
        this->dimensions = dimensions;
        pixel_count = dimensions.x * dimensions.y;

        sample_sum = make_unique_array<f32x3>(pixel_count);
        sample_count = make_unique_array<f32>(pixel_count);
        result = make_unique_array<u8x4>(pixel_count);
    }

    void add_sample(u32 x, u32 y, f32x3 const & color)
    {
        auto & sum = sample_sum[x + dimensions.x * y];
        auto & count = sample_count[x + dimensions.x * y];

        sum += color;
        count += 1;

        //const f32 threshhold = 4;
        //if (count >= threshhold)
        //{
        //    sum *= 1.f / threshhold;
        //    count = 1;
        //}
    }

    auto add_sample(u32x2 uv, f32x3 const & color)
    { return add_sample(uv.x, uv.y, color); }

    void calculate_result()
    {
        for (auto i = 0; i < pixel_count; ++i)
        {
            auto mean = sample_sum[i] / sample_count[i];
            mean *= 255.f;
            result[i] = u8x4(mean, 255);
        }
    }

    void calculate_result_with_gamma_correction()
    {
        for (auto i = 0; i < pixel_count; ++i)
        {
            auto mean = sample_sum[i] / sample_count[i];
            mean = glm::sqrt(mean); // gamma correction with g=2
            mean *= 255.f;
            result[i] = u8x4(mean, 255);
        }
    }
};

struct Camera
{
    f32x3 pos;
    f32x3 target;
    f32x2 screen_dimensions;
    f32x3 up, forward, right;
    f32x3 lower_left;

    struct Definition
    {
        f32x3 pos;
        f32x3 target;
        f32x3 up;
        f32 vfov;
    };
    Definition definition;

    Camera() = default;

    void create(Definition const & definition)
    {
        this->definition = definition;

        pos = definition.pos;
        target = definition.target;

        screen_dimensions.y = 2 * glm::tan(glm::radians(definition.vfov / 2));
        screen_dimensions.x = screen_dimensions.y;

        forward = normalize(target - pos);
        right = cross(forward, normalize(definition.up));
        up = cross(forward, -right);

        lower_left = pos
                   + forward
                   - right * (screen_dimensions.x * 0.5f)
                   - up    * (screen_dimensions.y * 0.5f);
    }

    Ray get_ray(f32x2 uv) const
    {
        auto out_pos = lower_left
                     + right * (screen_dimensions.x * uv.x)
                     + up * (screen_dimensions.y * uv.y);
        return Ray{
            .pos = pos,
            .dir = normalize(out_pos - pos),
            .bounce = 0
        };
    }
};

struct Renderer
{
    // settings
    u32 render_interval = 1;
    u32 sample_per_pixel = 16;
    u32 bounce_limit = 16;
    u32 scatter_count = 2;
    u32 rays_per_update = 500'000;
    u32 camera_rays_per_update = 0;

    bool activate_new_strategy = true;

    // resources
    GL::Texture2D gl_image;
    Image image;
    Scene scene;
    Camera camera;
    std::vector<Ray> rays;

    void create(u32x2 dimensions)
    {
        gl_image.create(GL::Texture2D::AttachmentDescription{
            .dimensions = dimensions,
            .internal_format = GL::GL_RGB8,
        });
        image.create(dimensions);
        scene.create();

        rays.reserve(dimensions.x * dimensions.y * 4);
        generate_camera_rays();
    }

    void reset_image()
    {
        image.create(image.dimensions);

        rays.clear();
        generate_camera_rays();
    }

    void generate_camera_rays()
    {
        // in the range [-0.5, 0.5]
        static const array sub_sample_points = []()
        {
            array<f32x2, 64> points;

            points[0] = {0, 0};

            for (auto i = 1; i < points.size(); ++i)
                points[i] = Random::next(f32x2(-0.5), f32x2(0.5));

            return points;
        }();

        auto const pos_normalizer = f32x2{1, 1} / f32x2(image.dimensions);

        for (auto x = 0; x < image.dimensions.x; ++x)
            for (auto y = 0; y < image.dimensions.y; ++y)
                for (auto i = 0; i < sample_per_pixel; ++i)
                {
                    auto ray = camera.get_ray((f32x2{x, y} + sub_sample_points[i]) * pos_normalizer);
                    ray.pixel = {x, y};
                    ray.color = {1, 1, 1};
                    rays.push_back(ray);
                }
    }

    void generate_random_camera_rays(u32 count)
    {
        auto const pos_normalizer = f32x2{1, 1} / f32x2(image.dimensions);

        for (auto _ = 0; _ < count; ++_)
        {
            auto pixel = Random::next(u32x2{0, 0}, image.dimensions - u32(1));
            auto ray = camera.get_ray(f32x2(pixel) * pos_normalizer);
            ray.color = {1, 1, 1};
            ray.pixel = pixel;
            rays.push_back(ray);
        }
    }

    bool is_done() const
    { return rays.empty(); }

    void update(FrameInfo const & frame_info)
    {
        static vector<Ray> generated_rays;

        generate_random_camera_rays(camera_rays_per_update);

        for (auto _ = glm::min(usize(rays_per_update), rays.size()); _ > 0; --_)
        {
            auto & ray = rays.back();

            Hit hit;
            
            if (activate_new_strategy)
                hit = scene.bvh.hit(ray, {0, std::numeric_limits<f32>::max()});
            else
                hit = scene.hittable_list.hit(ray, {0, std::numeric_limits<f32>::max()});

            if (hit.is_hit)
            {
                ray.color *= hit.attenuation;

                if (ray.bounce == bounce_limit)
                    image.add_sample(ray.pixel, ray.color);
                else
                    for (auto _ = ray.bounce < scatter_count ? 3 : 1; _ > 0; --_)
                    {
                        auto bounce_dir = normalize(hit.normal + Random::next_on_unit_sphere());
                        auto bounce_pos = hit.pos + bounce_dir * 0.001f;

                        generated_rays.push_back({
                            .pos = bounce_pos,
                            .dir = bounce_dir,
                            .color = ray.color,
                            .pixel = ray.pixel,
                            .bounce = ray.bounce + 1,
                        });
                    }
            }
            else
            {
                ray.color *= scene.get_background_color(ray.dir);
                image.add_sample(ray.pixel, ray.color);
            }

            rays.pop_back();
        }

        rays.insert(rays.end(), generated_rays.begin(), generated_rays.end());
        generated_rays.clear();
    }

    void render(FrameInfo const & frame_info)
    {
        using namespace GL;

        if (frame_info.idx % render_interval != 0)
            return;

        image.calculate_result_with_gamma_correction();

        glTextureSubImage2D(
            gl_image.id, 0,
            0, 0, image.dimensions.x, image.dimensions.y,
            GL_RGBA, GL_UNSIGNED_BYTE, image.result.get()
        );
    }
};
}
