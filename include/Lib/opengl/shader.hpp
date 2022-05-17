#pragma once

#include <string>
#include <regex>
#include <charconv>

#include "Lib/geometry/core.hpp"

#include "core.hpp"
#include "glue.hpp"

namespace GL
{
	// Pattern: String into Geometry::Attribute::Key
	Geometry::Attribute::Key IntoAttributeKey(std::string_view name)
	{
		using namespace Geometry::Attribute;
		
		if (name[0] == '_') // is custom
		{
			return KeyPool::get_or_create_key(name.substr(1));
		}
		else
		{
			if (name == "position") return Common::POSITION;
			if (name == "normal") return Common::NORMAL;
			if (name == "tangent") return Common::TANGENT;
			if (name == "texcoord") return Common::TEXCOORD;
			if (name == "color") return Common::COLOR;
			return KeyPool::get_or_create_key(name);
		}
	}

	struct ShaderStage : OpenGLObject
	{
		CTOR(ShaderStage, default)
		COPY(ShaderStage, delete)
		MOVE(ShaderStage, default)

		~ShaderStage()
		{
			glDeleteShader(id);
		}

		struct Description
		{
			GLenum stage;
			vector<char const*> const & sources;
		};

		void create(Description const & description)
		{
			id = glCreateShader(description.stage);

			glShaderSource(id, description.sources.size(), description.sources.data(), nullptr);
			glCompileShader(id);
		}

		bool is_compiled() const
		{
			GLsizei success;
			glGetShaderiv(id, GL_COMPILE_STATUS, &success);
			return success != 0;
		}

		std::string get_log() const
		{
			GLsizei log_size;
			glGetShaderiv(id, gl::GL_INFO_LOG_LENGTH, &log_size);

			std::string log(log_size, '\0');
			glGetShaderInfoLog(id, log_size, nullptr, log.data());

			return log;
		}

		std::string get_source() const
		{
			GLsizei source_size;
			glGetShaderiv(id, GL_SHADER_SOURCE_LENGTH, &source_size);

			std::string shader_source(source_size, '\0');
			glGetShaderSource(id, source_size, &source_size, shader_source.data());

			return shader_source;
		}
	};

	struct ShaderProgram : OpenGLObject
	{
		CTOR(ShaderProgram, default)
		COPY(ShaderProgram, delete)
		MOVE(ShaderProgram, default)

		~ShaderProgram()
		{
			glDeleteProgram(id);
		}

		struct Description
		{
			vector<ShaderStage const*> const & shader_stages;
		};

		void create(Description const & description)
		{
			id = glCreateProgram();

			for (auto & shader_stage: description.shader_stages)
				glAttachShader(id, shader_stage->id);

			glLinkProgram(id);

			for (auto & shader_stage: description.shader_stages)
				glDetachShader(id, shader_stage->id);
		}

		bool is_linked() const
		{
			i32 success;
			glGetProgramiv(id, GL_LINK_STATUS, &success);
			return success;
		}

		std::string get_log() const
		{
			i32 log_size;
			glGetProgramiv(id, gl::GL_INFO_LOG_LENGTH, &log_size);

			std::string log(log_size, '\0');
			glGetProgramInfoLog(id, log_size, nullptr, log.data());

			return log;
		}

		// https://www.khronos.org/opengl/wiki/Program_Introspection#Interface_query

		struct AttributeMapping
		{
			u32 location;
			GLenum glsl_type; // just for debug purposes
			bool per_patch;
			Geometry::Attribute::Key key;
		};
		vector<AttributeMapping> attribute_mappings;

		void update_attribute_mapping()
		{
			auto const interface = GL_PROGRAM_INPUT;

			array const query_props{
				GL_LOCATION,
				GL_TYPE,
				GL_IS_PER_PATCH,
			};
			array<i32, query_props.size()> query_results;

			attribute_mappings.clear();

			i32 max_name_size;
			glGetProgramInterfaceiv(id, interface, GL_MAX_NAME_LENGTH, &max_name_size);
			std::string name_buffer(max_name_size, 0);

			i32 resource_size;
			glGetProgramInterfaceiv(id, interface, GL_ACTIVE_RESOURCES, &resource_size);

			attribute_mappings.reserve(resource_size);
			for (auto i = 0; i < resource_size; ++i)
			{
				glGetProgramResourceiv(
					id, interface, i,
					query_props.size(), query_props.data(),
					query_props.size(), nullptr, query_results.data()
				);

				i32 name_size;
				glGetProgramResourceName(
					id, interface, i,
					name_buffer.size(), &name_size, name_buffer.data()
				);

				attribute_mappings.emplace_back(
					AttributeMapping{
						.location = static_cast<u32>(query_results[0]),
						.glsl_type = static_cast<GLenum>(query_results[1]),
						.per_patch = static_cast<bool>(query_results[2]),
						.key = IntoAttributeKey({name_buffer.data(), size_t(name_size)}),
					}
				);
			}
		};

		struct UniformMapping
		{
			u32 location;
			GLenum glsl_type; // just for debug purposes
			std::string key;
		};
		vector<UniformMapping> uniform_mappings;

		void update_uniform_mapping()
		{
			auto const interface = GL_UNIFORM;

			array const query_props{
				GL_LOCATION,
				GL_TYPE,
			};
			array<i32, query_props.size()> query_results;

			i32 max_name_size;
			glGetProgramInterfaceiv(id, interface, GL_MAX_NAME_LENGTH, &max_name_size);
			std::string name_buffer(max_name_size, 0);

			i32 resource_size;
			glGetProgramInterfaceiv(id, interface, GL_ACTIVE_RESOURCES, &resource_size);

			uniform_mappings.clear();

			uniform_mappings.reserve(resource_size);
			for (auto i = 0; i < resource_size; ++i)
			{
				glGetProgramResourceiv(
					id, interface, i,
					query_props.size(), query_props.data(),
					query_props.size(), nullptr, query_results.data()
				);

				// Skip block uniform
				if (query_results[0] == -1)
					continue;

				i32 name_size;
				glGetProgramResourceName(
					id, interface, i,
					name_buffer.size(), &name_size, name_buffer.data()
				);

				uniform_mappings.emplace_back(
					UniformMapping{
						.location = static_cast<u32>(query_results[0]),
						.glsl_type = static_cast<GLenum>(query_results[1]),
						.key = name_buffer.substr(0, name_size),
					}
				);
			}
		};

		struct UniformBlockMapping
		{
			u32 location;
			u32 data_size;
			std::string key;

			struct Variable
			{
				u32 offset;
				GLenum glsl_type; // just for debug purposes
				std::string key;
			};
			vector<Variable> variables;
		};
		vector<UniformBlockMapping> uniform_block_mappings;

		void update_uniform_block_mapping()
		{
			auto const interface = GL_UNIFORM_BLOCK;

			array const query_props{
				GL_BUFFER_BINDING,
				GL_BUFFER_DATA_SIZE,
				GL_NUM_ACTIVE_VARIABLES,
			};
			array<i32, query_props.size()> query_results;

			i32 max_name_size;
			glGetProgramInterfaceiv(id, interface, GL_MAX_NAME_LENGTH, &max_name_size);
			std::string name_buffer(max_name_size, 0);

			i32 resource_size;
			glGetProgramInterfaceiv(id, interface, GL_ACTIVE_RESOURCES, &resource_size);

			uniform_block_mappings.clear();
			uniform_block_mappings.reserve(resource_size);
			for (auto i = 0; i < resource_size; ++i)
			{
				glGetProgramResourceiv(
					id, interface, i,
					query_props.size(), query_props.data(),
					query_props.size(), nullptr, query_results.data()
				);

				i32 name_size;
				glGetProgramResourceName(
					id, interface, i,
					name_buffer.size(), &name_size, name_buffer.data()
				);

				uniform_block_mappings.emplace_back(
					UniformBlockMapping{
						.location = static_cast<u32>(query_results[0]),
						.data_size = static_cast<u32>(query_results[1]),
						.key = name_buffer.substr(0, name_size),
						.variables = query_uniform_block_variables(i, static_cast<u32>(query_results[2])),
					}
				);
			}
		};

		vector<UniformBlockMapping::Variable> query_uniform_block_variables(u32 block_idx, u32 variable_count)
		{
			vector<i32> variable_indices;
			variable_indices.resize(variable_count);

			glGetProgramResourceiv(
				id, GL_UNIFORM_BLOCK, block_idx,
				1, &GL_ACTIVE_VARIABLES,
				variable_count, nullptr, variable_indices.data()
			);

			auto const interface = GL_UNIFORM;

			array const query_props{
				GL_OFFSET,
				GL_TYPE,
			};
			array<i32, query_props.size()> query_results;

			i32 max_name_size;
			glGetProgramInterfaceiv(id, interface, GL_MAX_NAME_LENGTH, &max_name_size);
			std::string name_buffer(max_name_size, 0);

			vector<UniformBlockMapping::Variable> variables;
			variables.reserve(variable_count);

			for (auto i: variable_indices)
			{
				glGetProgramResourceiv(
					id, interface, i,
					query_props.size(), query_props.data(),
					query_props.size(), nullptr, query_results.data()
				);

				i32 name_size;
				glGetProgramResourceName(
					id, interface, i,
					name_buffer.size(), &name_size, name_buffer.data()
				);

				variables.emplace_back(UniformBlockMapping::Variable{
					.offset = static_cast<u32>(query_results[0]),
					.glsl_type = static_cast<GLenum>(query_results[1]),
					.key = name_buffer.substr(0, name_size),
				});
			}

			std::ranges::sort(variables, std::ranges::less{}, &UniformBlockMapping::Variable::offset);

			return variables;
		}

		struct StorageBlockMapping
		{
			u32 location;
			u32 data_size;
			std::string key;

			struct Variable
			{
				u32 offset;
				GLenum glsl_type; // just for debug purposes
				std::string key;
			};
			vector<Variable> variables;
		};
		vector<StorageBlockMapping> storage_block_mappings;

		void update_storage_block_mapping()
		{
			auto const interface = GL_SHADER_STORAGE_BLOCK;

			array const query_props{
				GL_BUFFER_BINDING,
				GL_BUFFER_DATA_SIZE,
				GL_NUM_ACTIVE_VARIABLES,
			};
			array<i32, query_props.size()> query_results;

			i32 max_name_size;
			glGetProgramInterfaceiv(id, interface, GL_MAX_NAME_LENGTH, &max_name_size);
			std::string name_buffer(max_name_size, 0);

			i32 resource_size;
			glGetProgramInterfaceiv(id, interface, GL_ACTIVE_RESOURCES, &resource_size);

			storage_block_mappings.clear();
			storage_block_mappings.reserve(resource_size);
			for (auto i = 0; i < resource_size; ++i)
			{
				glGetProgramResourceiv(
					id, interface, i,
					query_props.size(), query_props.data(),
					query_props.size(), nullptr, query_results.data()
				);

				i32 name_size;
				glGetProgramResourceName(
					id, interface, i,
					name_buffer.size(), &name_size, name_buffer.data()
				);

				storage_block_mappings.emplace_back(
					StorageBlockMapping{
						.location = static_cast<u32>(query_results[0]),
						.data_size = static_cast<u32>(query_results[1]),
						.key = name_buffer.substr(0, name_size),
						.variables = query_sotrage_block_variables(i, static_cast<u32>(query_results[2])),
					}
				);
			}
		};

		vector<StorageBlockMapping::Variable> query_sotrage_block_variables(u32 block_idx, u32 variable_count)
		{
			vector<i32> variable_indices;
			variable_indices.resize(variable_count);

			glGetProgramResourceiv(
				id, GL_SHADER_STORAGE_BLOCK, block_idx,
				1, &GL_ACTIVE_VARIABLES,
				variable_count, nullptr, variable_indices.data()
			);

			auto const interface = GL_BUFFER_VARIABLE;

			array const query_props{
				GL_OFFSET,
				GL_TYPE,
			};
			array<i32, query_props.size()> query_results;

			i32 max_name_size;
			glGetProgramInterfaceiv(id, interface, GL_MAX_NAME_LENGTH, &max_name_size);
			std::string name_buffer(max_name_size, 0);

			vector<StorageBlockMapping::Variable> variables;
			variables.reserve(variable_count);

			for (auto i: variable_indices)
			{
				glGetProgramResourceiv(
					id, interface, i,
					query_props.size(), query_props.data(),
					query_props.size(), nullptr, query_results.data()
				);

				i32 name_size;
				glGetProgramResourceName(
					id, interface, i,
					name_buffer.size(), &name_size, name_buffer.data()
				);

				variables.emplace_back(StorageBlockMapping::Variable{
					.offset = static_cast<u32>(query_results[0]),
					.glsl_type = static_cast<GLenum>(query_results[1]),
					.key = name_buffer.substr(0, name_size),
				});
			}

			std::ranges::sort(variables, std::ranges::less{}, &StorageBlockMapping::Variable::offset);

			return variables;
		}

		void update_interface_mapping()
		{
			update_attribute_mapping();
			update_uniform_mapping();
			update_uniform_block_mapping();
			update_storage_block_mapping();
		}
	};
}