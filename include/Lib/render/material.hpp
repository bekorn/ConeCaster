#pragma once

#include "Lib/opengl/core.hpp"

namespace Render
{
	struct IMaterial
	{
		// the functionality should be compatible with any shader,
		// but the data should be tightly coupled with individual shaders

		// TODO(bekorn): maybe generate the concrete Material structs with another program then compile this program.
		//  implication: no dynamic materials! for rendering purposes, ok. However, for an editor, this can be bad.
		//  maybe: the material class can be C-simple that it can be loaded dynamically on runtime (via a dll).
		//   so jit compile the generated Material into a dll and watch for changed dlls on this side.

		virtual GL::StorageBlock const & get_block() const = 0;
		virtual void write_to_buffer(byte * buffer) const = 0;
		virtual void read_from_buffer(const byte * buffer) = 0;
	};

	// for now, here is an example concrete Material
	struct Material_gltf_pbrMetallicRoughness final : IMaterial
	{
		static inline GL::StorageBlock block;

		GL::StorageBlock const & get_block() const final
		{ return block; }

		u64 base_color_texture_handle{0};
		f32x4 base_color_factor{0, 0, 0, 0};

		u64 metallic_roughness_texture_handle{0};
		f32x2 metallic_roughness_factor{0, 0};

		u64 emissive_texture_handle{0};
		f32x3 emissive_factor{0, 0, 0};

		u64 occlusion_texture_handle{0};
		u64 normal_texture_handle{0};

		void write_to_buffer(byte * buffer) const final
		{
			block.set(buffer, "Material.base_color_texture_handle", base_color_texture_handle);
			block.set(buffer, "Material.base_color_factor", base_color_factor);

			block.set(buffer, "Material.metallic_roughness_texture_handle", metallic_roughness_texture_handle);
			block.set(buffer, "Material.metallic_roughness_factor", metallic_roughness_factor);

			block.set(buffer, "Material.emissive_texture_handle", emissive_texture_handle);
			block.set(buffer, "Material.emissive_factor", emissive_factor);

			block.set(buffer, "Material.occlusion_texture_handle", occlusion_texture_handle);
			block.set(buffer, "Material.normal_texture_handle", normal_texture_handle);
		}

		void read_from_buffer(const byte * buffer) final
		{
			block.get(buffer, "Material.base_color_texture_handle", base_color_texture_handle);
			block.get(buffer, "Material.base_color_factor", base_color_factor);

			block.get(buffer, "Material.metallic_roughness_texture_handle", metallic_roughness_texture_handle);
			block.get(buffer, "Material.metallic_roughness_factor", metallic_roughness_factor);

			block.get(buffer, "Material.emissive_texture_handle", emissive_texture_handle);
			block.get(buffer, "Material.emissive_factor", emissive_factor);

			block.get(buffer, "Material.occlusion_texture_handle", occlusion_texture_handle);
			block.get(buffer, "Material.normal_texture_handle", normal_texture_handle);
		}
	};
}