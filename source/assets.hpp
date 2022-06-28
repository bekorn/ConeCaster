#pragma once

#include "Lib/core/.hpp"

struct Assets
{
    std::filesystem::path assets_dir;

	Managed<GL::Texture2D> textures;
	Managed<unique_one<Render::IMaterial>> materials;
	Managed<Geometry::Primitive> primitives;
	Managed<Render::Mesh> meshes;

    Assets()
    {
        auto project_root = std::filesystem::current_path();
        while (project_root.filename() != "ConeCaster")
            project_root = project_root.parent_path();
        assets_dir = project_root / "assets";
    }

    void create()
    {
        Timer timer;

        // auto loaded_data = GLTF::Load({.name = "skull", .path = assets_dir / "skull/scene.gltf"});
        auto loaded_data = GLTF::Load({.name = "damaged helmet", .path = assets_dir / "damaged_helmet/DamagedHelmet.gltf"});
        timer.timeit(stderr, "Model Loaded");

        GLTF::Convert(loaded_data, textures, materials, primitives, meshes);
        timer.timeit(stderr, "Model Converted");
    }
};