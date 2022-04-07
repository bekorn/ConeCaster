#pragma once

#include "Lib/core/utils.hpp"

#include "core.hpp"
#include "uniform_block.hpp"
#include "storage_block.hpp"

namespace GL
{
	struct Buffer : OpenGLObject
	{
		CTOR(Buffer, default);
		COPY(Buffer, delete);
		MOVE(Buffer, default);

		~Buffer()
		{
			glDeleteBuffers(1, &id);
		}

		struct Description
		{
			GLenum usage = GL_STATIC_DRAW;
			span<byte> data = {};
		};

		void create(Description const & description)
		{
			glCreateBuffers(1, &id);
			glNamedBufferData(
				id,
				description.data.size(),
				description.data.data(),
				description.usage
			);
		}

		struct UniformBlockDescription
		{
			GLenum usage = GL_DYNAMIC_DRAW;
			GL::UniformBlock const & uniform_block;
			u32 array_size = 1;
		};

		void create(UniformBlockDescription const & description)
		{
			glCreateBuffers(1, &id);
			glNamedBufferData(
				id,
				description.uniform_block.aligned_size * description.array_size,
				nullptr,
				description.usage
			);
		}

		struct StorageBlockDescription
		{
			GLenum usage = GL_DYNAMIC_DRAW;
			GL::StorageBlock const & storage_block;
			usize array_size = 1;
		};

		void create(StorageBlockDescription const & description)
		{
			glCreateBuffers(1, &id);
			glNamedBufferData(
				id,
				description.storage_block.aligned_size * description.array_size,
				nullptr,
				description.usage
			);
		}
	};
}
