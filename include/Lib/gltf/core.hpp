#pragma once

#include <string>

#include "Lib/core/core.hpp"

namespace GLTF
{
	// Spec: https://www.khronos.org/registry/glTF/specs/2.0/glTF-2.0.pdf

	using ::ByteBuffer;

	struct BufferView
	{
		u32 buffer_index;
		u32 offset;
		u32 length;
		optional<u32> stride;
	};

	struct Accessor
	{
		u32 buffer_view_index;
		u32 byte_offset;
		u32 vector_data_type;
		u32 vector_dimension;
		u32 count;
		bool normalized;
	};

	struct Image
	{
		ByteBuffer data;
		i32x2 dimensions;
		i32 channels;
	};

	struct Sampler
	{
		u32 min_filter;
		u32 mag_filter;
		u32 wrap_s;
		u32 wrap_t;
	};
	// TODO(bekorn): read about OpenGL Sampler Objects
	static Sampler constexpr SamplerDefault{
		.min_filter = 9729, // LINEAR
		.mag_filter = 9729, // LINEAR
		.wrap_s = 10497, // REPEAT
		.wrap_t = 10497 // REPEAT
	};

	struct Texture
	{
		std::string name;
		optional<u32> image_index;
		optional<u32> sampler_index;
	};

	struct Attribute
	{
		std::string name;
		u32 accessor_index;
	};

	// Equivalent of a draw call
	struct Primitive
	{
		std::string name;
		vector<Attribute> attributes;
		optional<u32> indices_accessor_index;
		optional<u32> material_index;
	};

	struct Mesh
	{
		std::string name;
		vector<Primitive> primitives;
	};

	struct Material
	{
		std::string name;

		struct TexInfo
		{
			u32 texture_index;
			u32 texcoord_index; // which TEXCOORD_n attribute to use
		};

		struct PbrMetallicRoughness
		{
			f32x4 base_color_factor;
			optional<TexInfo> base_color_texture;
			f32 metallic_factor;
			f32 roughness_factor;
			optional<TexInfo> metallic_roughness_texture;
		};
		optional<PbrMetallicRoughness> pbr_metallic_roughness;

		optional<TexInfo> normal_texture;
		f32 normal_texture_scale;

		optional<TexInfo> occlusion_texture;
		f32 occlusion_texture_strength;

		optional<TexInfo> emissive_texture;
		f32x3 emissive_factor;

		std::string alpha_mode;
		f32 alpha_cutoff;

		bool double_sided;
	};

	struct Node
	{
		std::string name;
		f32x3 translation;
		f32quat rotation;
		f32x3 scale;
		optional<u32> mesh_index;
		vector<u32> child_indices;
	};

	struct LoadedData
	{
		vector<ByteBuffer> buffers;
		vector<BufferView> buffer_views;
		vector<Accessor> accessors;

		vector<Image> images;
		vector<Sampler> samplers;
		vector<Texture> textures;

		vector<Mesh> meshes;
		vector<Material> materials;

		vector<Node> nodes;
	};

	auto const pbrMetallicRoughness_program_name = "gltf_pbrMetallicRoughness"_name;
}