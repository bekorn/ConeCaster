#pragma once

#include "Lib/core/core.hpp"
#include "Lib/geometry/core.hpp"

#include "core.hpp"
#include "buffer.hpp"
#include "shader.hpp"
#include "glue.hpp"

namespace GL
{
	struct VertexArray : OpenGLObject
	{
		vector<Buffer> vertex_buffers;
		Buffer element_buffer;
		GLsizei element_count;

		CTOR(VertexArray, default);
		COPY(VertexArray, delete);
		MOVE(VertexArray, default);

		~VertexArray()
		{
			glDeleteVertexArrays(1, &id);
		}

		void create(Geometry::Primitive const & geometry, ShaderProgram const & shader)
		{
			// TODO(bekorn): This method creates a vbo for each attribute,
			//  it would be better to merge them in a single vbo

			glCreateVertexArrays(1, &id);

			for (auto const & attribute: shader.attribute_mappings)
			{
				auto const geo_attribute = geometry.attributes.find(attribute.key);
				if (geo_attribute == geometry.attributes.end())
					continue;

				auto const & [key, data] = *geo_attribute;

				Buffer vertex_buffer;
				vertex_buffer.create(
					{
						.data = data.buffer.span_as<byte>(),
					}
				);
				auto const buffer_bind_index = vertex_buffers.size();

				glVertexArrayVertexBuffer(
					id,
					buffer_bind_index,
					vertex_buffer.id, 0, data.type.size() * data.dimension
				);
				glVertexArrayAttribFormat(
					id,
					attribute.location,
					data.dimension, IntoGLenum(data.type),data.type.is_normalized(),
					0
				);
				glEnableVertexArrayAttrib(id, attribute.location);
				glVertexArrayAttribBinding(id, attribute.location, buffer_bind_index);

				vertex_buffers.emplace_back(move(vertex_buffer));
			}

			element_buffer.create(
				{
					.data = span((byte*)geometry.indices.data(), geometry.indices.size() * sizeof(u32)),
				}
			);

			glVertexArrayElementBuffer(id, element_buffer.id);

			element_count = geometry.indices.size();
		}
	};
}