#pragma once

#include "Lib/display/imgui.hpp"
#include "Lib/render/renderer.hpp"

struct Editor
{
    Render::Renderer & renderer;

    Editor(Render::Renderer & renderer) :
        renderer(renderer)
    {}

    void render(Render::FrameInfo const & frame_info)
    {
        using namespace ImGui;

        Begin("Render");
        {
            u64 id = renderer.gl_image.id;
            Image(
                (void*)id,
                {f32(renderer.image.dimensions.x), f32(renderer.image.dimensions.y)},
                {0, 1}, {1, 0}
            );
        }
        End();


        Begin("Metrics");
        Text("Resolution: %d x %d", renderer.image.dimensions.x, renderer.image.dimensions.y);
        Text("Rays Left: %d", renderer.rays.size());
        
        static f64 render_start = frame_info.seconds_since_start;
        Text("Render Start: %f", render_start);

        static f64 render_done;
        if (not renderer.is_done())
            render_done = frame_info.seconds_since_start;
        Text("Elapsed Time: %f", render_done - render_start);
        End();


        Begin("Settings");
        if (Button("Reset Image") or IsKeyPressed(ImGuiKey_R, false))
        {
            renderer.reset_image();
            render_start = frame_info.seconds_since_start;
        }
        Checkbox("New Strategy Active", &renderer.activate_new_strategy);
        {
            static array names = {
                "None", "BVH", "BVH_8WIDE"
            };
            if (BeginCombo("Acceleration", names[(u32)renderer.accelerator]))
            {
                for (auto i = 0; i < names.size(); ++i)
                    if (Selectable(names[i]))
                        renderer.accelerator = (Render::Renderer::Accelerator)i;
                EndCombo();
            }
        }
        {
            u32 min = 1, max = 64;
            SliderScalar("Samples per Pixel", ImGuiDataType_U32, &renderer.sample_per_pixel, &min, &max);
        }
        {
            u32 min = 0, max = 500, val = renderer.rays_per_update / 1'000;
            SliderScalar("Rays per Update", ImGuiDataType_U32, &val, &min, &max, "%d K");
            renderer.rays_per_update = val * 1'000;
        }
        {
            u32 min = 0, max = 100'000, val = renderer.camera_rays_per_update;
            SliderScalar("Camera Rays per Update", ImGuiDataType_U32, &val, &min, &max);
            renderer.camera_rays_per_update = val;
        }
        {
            u32 min = 1, max = 60;
            SliderScalar("Render Interval", ImGuiDataType_U32, &renderer.render_interval, &min, &max);
        }
        {
            u32 min = 0, max = 64;
            SliderScalar("Bounce Limit", ImGuiDataType_U32, &renderer.bounce_limit, &min, &max);
        }
        {
            u32 min = 0, max = 16;
            SliderScalar("Scatter First N", ImGuiDataType_U32, &renderer.scatter_count, &min, &max);
        }
        {
            // DragScalar("Focal Z", ImGuiDataType_Float, &renderer.focal_z);
        }
        End();


        Begin("Results");
        static array<GL::Texture2D, 10> results;

        static u32 result_idx = 0;
        static u32 const min = 0, max = results.size() - 1;
        SliderScalar("Result", ImGuiDataType_U32, &result_idx, &min, &max);

        // shortcut for keys [1-9] and [0]
        for (auto i = 0; i < results.size(); ++i)
            if (IsKeyPressed(ImGuiKey_0 + i, false))
                result_idx = i;

        auto & result = results[result_idx];

        if (Button("Save Render") or IsKeyPressed(ImGuiKey_S, false))
        {
            result = {};
            result.create(GL::Texture2D::AttachmentDescription{
                .dimensions = renderer.image.dimensions,
                .internal_format = GL::GL_RGB8,
            });
            GL::glCopyImageSubData(
                renderer.gl_image.id, GL::GL_TEXTURE_2D,
                0, 0, 0, 0,
                result.id, GL::GL_TEXTURE_2D,
                0, 0, 0, 0,
                renderer.image.dimensions.x, renderer.image.dimensions.y, 1
            );
        }

        if (result.id == 0)
            TextUnformatted("Empty");
        else
        {
            u64 id = result.id;
            Image(
                (void*)id,
                {f32(renderer.image.dimensions.x), f32(renderer.image.dimensions.y)},
                {0, 1}, {1, 0}
            );
        }

        End();

        
        Begin("Scene");
        ColorEdit3("Sky Up", begin(renderer.scene.background_color_up));
        ColorEdit3("Sky Down", begin(renderer.scene.background_color_down));

        {
            TextUnformatted("Camera");
            auto camera_definition = renderer.camera.definition;
            auto is_changed = false;
            PushID("Camera");
            is_changed |= DragFloat3("Position", begin(camera_definition.pos));
            is_changed |= DragFloat3("Target", begin(camera_definition.target));
            is_changed |= DragFloat3("Up", begin(camera_definition.up));
            is_changed |= DragFloat("VFOV", &camera_definition.vfov);
            PopID();
            if (is_changed)
                renderer.camera.create(camera_definition);
        }

        End();
    }
};