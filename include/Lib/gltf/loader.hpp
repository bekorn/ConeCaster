#pragma once

#include <execution>

#include "Lib/core/file_io.hpp"

#include "core.hpp"

namespace GLTF
{
	struct Description
	{
		std::string name;
		std::filesystem::path path;
	};

	LoadedData Load(Description const & description)
	{
		using namespace rapidjson;
		using namespace File;
		using namespace File::JSON;

		LoadedData loaded;

		Document document;
		document.Parse(LoadAsString(description.path).c_str());

		auto const file_dir = description.path.parent_path();

		// Parse buffers
		for (auto const & item: document["buffers"].GetArray())
		{
			auto const & buffer = item.GetObject();

			auto file_size = buffer["byteLength"].GetUint64();
			auto file_name = buffer["uri"].GetString();
			// Limitation: only loads separate file binaries
			loaded.buffers.emplace_back(LoadAsBytes(file_dir / file_name, file_size));
		}

		// Parse buffer views
		for (auto const & item: document["bufferViews"].GetArray())
		{
			auto const & buffer_view = item.GetObject();

			loaded.buffer_views.push_back(
				{
					.buffer_index = GetU32(buffer_view, "buffer"),
					.offset = GetU32(buffer_view, "byteOffset", 0),
					.length = GetU32(buffer_view, "byteLength"),
					.stride = GetOptionalU32(buffer_view, "byteStride"),
				}
			);
		}

		// Parse images
		if (auto const member = document.FindMember("images"); member != document.MemberEnd())
		{
			auto const & items = member->value.GetArray();
			loaded.images.resize(items.Size());
			std::transform(
				std::execution::par_unseq,
				items.Begin(), items.End(),
				loaded.images.data(),
				[&file_dir](Document::Array::ValueType const & item) -> GLTF::Image
				{
					auto const & image = item.GetObject();

					auto const member = image.FindMember("uri");
					if (member == image.MemberEnd())
						throw std::runtime_error("images without a uri file path are not supported yet");

					auto uri = member->value.GetString();
					if (uri[5] == ':') // check for "data:" (base64 encoded data as a json string)
						throw std::runtime_error("images without a uri file path are not supported yet");

					auto const file_data = LoadAsBytes(file_dir / uri);
					i32x2 dimensions;
					i32 channels;
					void* raw_pixel_data = stbi_load_from_memory(
						file_data.data_as<const unsigned char>(), file_data.size,
						&dimensions.x, &dimensions.y,
						&channels, 0
					);

					auto data = ByteBuffer(
						move(raw_pixel_data),
						dimensions.x * dimensions.y * channels
					);

					return {
						.data = move(data),
						.dimensions = dimensions,
						.channels = channels,
					};
				}
			);
		}

		// Parse samplers
		if (auto const member = document.FindMember("samplers"); member != document.MemberEnd())
		{
			for (auto const & item: member->value.GetArray())
			{
				auto const & sampler = item.GetObject();

				// Min/Mag filters have no default values in the spec, I picked the values
				loaded.samplers.push_back(
					{
						.min_filter = GetU32(sampler, "minFilter", SamplerDefault.min_filter),
						.mag_filter = GetU32(sampler, "magFilter", SamplerDefault.mag_filter),
						.wrap_s = GetU32(sampler, "wrapS", SamplerDefault.wrap_s),
						.wrap_t = GetU32(sampler, "wrapT", SamplerDefault.wrap_t),
					}
				);
			}
		}

		// Parse textures
		NameGenerator texture_name_generator{.prefix = description.name + ":texture:"};
		if (auto const member = document.FindMember("textures"); member != document.MemberEnd())
		{
			for (auto const & item: member->value.GetArray())
			{
				auto const & texture = item.GetObject();

				loaded.textures.push_back(
					{
						.name = texture_name_generator.get(texture, "name"),
						.image_index = GetOptionalU32(texture, "source"),
						.sampler_index = GetOptionalU32(texture, "sampler"),
					}
				);
			}
		}

		// Parse accessors
		for (auto const & item: document["accessors"].GetArray())
		{
			auto const get_type_dimension = [](std::string const & type) -> u32
			{
				if (type == "SCALAR") return 1;
				if (type == "VEC2") return 2;
				if (type == "VEC3") return 3;
				if (type == "VEC4") return 4;
				if (type == "MAT2") return 4;
				if (type == "MAT3") return 9;
				/*if(type == "MAT4")*/ return 16;
			};

			auto const & accessor = item.GetObject();

			loaded.accessors.push_back(
				{
					.buffer_view_index = accessor["bufferView"].GetUint(),
					.byte_offset = GetU32(accessor, "byteOffset", 0),
					.vector_data_type = accessor["componentType"].GetUint(),
					.vector_dimension = get_type_dimension(accessor["type"].GetString()),
					.count = accessor["count"].GetUint(),
					.normalized = GetBool(accessor, "normalized", false),
				}
			);
		}

		// Parse meshes
		NameGenerator mesh_name_generator{.prefix = description.name + ":mesh:"};
		NameGenerator primitive_name_generator{.prefix = description.name + ":primitive:"};
		for (auto const & item: document["meshes"].GetArray())
		{
			auto const & mesh = item.GetObject();

			vector<Primitive> primitives;
			primitives.reserve(mesh["primitives"].Size());
			for (auto const & item: mesh["primitives"].GetArray())
			{
				auto const & primitive = item.GetObject();

				vector<Attribute> attributes;
				attributes.reserve(primitive["attributes"].MemberCount());
				for (auto const & attribute: primitive["attributes"].GetObject())
				{
					attributes.push_back(
						{
							.name = attribute.name.GetString(),
							.accessor_index = attribute.value.GetUint(),
						}
					);
				}

				primitives.push_back(
					{
						.name = primitive_name_generator.get(primitive, "name"),
						.attributes = attributes,
						.indices_accessor_index = GetOptionalU32(primitive, "indices"),
						.material_index = GetOptionalU32(primitive, "material"),
					}
				);
			}

			loaded.meshes.push_back(
				{
					.name = mesh_name_generator.get(mesh, "name"),
					.primitives = primitives,
				}
			);
		}

		// Parse materials
		NameGenerator material_name_generator{.prefix = description.name + ":material:"};
		for (auto const & item: document["materials"].GetArray())
		{
			auto const get_tex_info = [](JSONObj material, Key key) -> optional<Material::TexInfo>
			{
				auto member = material.FindMember(key.data());
				if (member != material.MemberEnd())
				{
					auto tex_info = member->value.GetObject();
					return Material::TexInfo{
						.texture_index = GetU32(tex_info, "index"),
						.texcoord_index = GetU32(tex_info, "texCoord", 0),
					};
				}
				else
					return {};
			};

			auto const & material = item.GetObject();
			Material mat{
				.name = material_name_generator.get(material, "name"),

				.normal_texture = get_tex_info(material, "normalTexture"),

				.occlusion_texture = get_tex_info(material, "occlusionTexture"),

				.emissive_texture = get_tex_info(material, "emissiveTexture"),
				.emissive_factor = GetF32x3(material, "emissiveFactor", {0, 0, 0}),

				.alpha_mode = GetString(material, "alphaMode", "OPAQUE"),
				.alpha_cutoff = GetF32(material, "alphaCutoff", 0.5),

				.double_sided = GetBool(material, "doubleSided", false),
			};

			if (auto member = material.FindMember("pbrMetallicRoughness"); member != material.MemberEnd())
			{
				auto const & pbrMetallicRoughness = member->value.GetObject();
				mat.pbr_metallic_roughness = {
					.base_color_factor = GetF32x4(pbrMetallicRoughness, "baseColorFactor", {1, 1, 1, 1}),
					.base_color_texture = get_tex_info(pbrMetallicRoughness, "baseColorTexture"),
					.metallic_factor = GetF32(pbrMetallicRoughness, "metallicFactor", 1),
					.roughness_factor = GetF32(pbrMetallicRoughness, "roughnessFactor", 1),
					.metallic_roughness_texture = get_tex_info(pbrMetallicRoughness, "metallicRoughnessTexture"),
				};
			}

			loaded.materials.push_back(mat);
		}

		// Parse nodes
		NameGenerator node_name_generator{.prefix = description.name + ":node:"};
		if (auto const member = document.FindMember("nodes"); member != document.MemberEnd())
		{
			for (auto const & item: member->value.GetArray())
			{
				auto const & gltf_node = item.GetObject();

				Node node{
					.name = node_name_generator.get(gltf_node, "name"),
					.mesh_index = GetOptionalU32(gltf_node, "mesh"),
				};

				if (auto member = gltf_node.FindMember("matrix"); member != gltf_node.MemberEnd())
				{
					f32x4x4 matrix;
					auto const & arr = member->value.GetArray();
					for (auto i = 0; i < 4; ++i)
						for (auto j = 0; j < 4; ++j)
							matrix[i][j] = arr[i * 4 + j].GetFloat();

					f32x3 skew;
					f32x4 perspective;
					glm::decompose(matrix, node.scale, node.rotation, node.translation, skew, perspective);
				}
				else
				{
					node.translation = GetF32x3(gltf_node, "translation", {0, 0, 0});
					auto rotation = GetF32x4(gltf_node, "rotation", {0, 0, 0, 1});			 // gltf order (x, y, z, w)
					node.rotation = f32quat(rotation.w, rotation.x, rotation.y, rotation.z); // glm  order (w, x, y, z)
					node.scale = GetF32x3(gltf_node, "scale", f32x3{1, 1, 1});
				}

				if (auto member = gltf_node.FindMember("children"); member != gltf_node.MemberEnd())
					for (auto & child_index : member->value.GetArray())
						node.child_indices.push_back(child_index.GetUint());

				loaded.nodes.push_back(node);
			}
		}

		return loaded;
	}
}
