#pragma once

#include "core.hpp"

namespace GL
{
	// TODO(bekorn): maybe merge this and UniformBlock into Block or InterfaceBlock
	//  because the only difference is the alignment size and that can be specified in the description
	struct StorageBlock
	{
		u32 binding;
		u32 data_size;
		usize aligned_size;
		std::string key;

		struct Variable
		{
			u32 offset;
			GLenum glsl_type;
		};
		// TODO(bekorn): Name can be utilized here to reduce string hashing
		std::unordered_map<std::string, Variable> variables;

		struct Description
		{
			GL::ShaderProgram::StorageBlockMapping const & layout;
		};

		void create(Description const & description)
		{
			i32 alignment;
			glGetIntegerv(GL_SHADER_STORAGE_BUFFER_OFFSET_ALIGNMENT, &alignment);

			binding = description.layout.location;
			data_size = description.layout.data_size;
			aligned_size = (data_size / alignment + 1) * alignment;
			key = description.layout.key;

			for (auto & variable: description.layout.variables)
				variables.try_emplace(variable.key, Variable{
					.offset = variable.offset,
					.glsl_type = variable.glsl_type,
				});
		}

		template<typename T>
		void set(byte * destination, std::string const & variable_key, T const & data) const
		{
			std::memcpy(
				destination + variables.at(variable_key).offset,
				&data, sizeof(T)
			);
		}

		template<typename T>
		void get(const byte * buffer, std::string const & variable_key, T & destination) const
		{
			std::memcpy(
				&destination,
				buffer + variables.at(variable_key).offset, sizeof(T)
			);
		}
	};
}