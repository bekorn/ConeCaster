#pragma once

#include "core.hpp"

namespace GL
{
	struct FrameBuffer : OpenGLObject
	{
		CTOR(FrameBuffer, default)
		COPY(FrameBuffer, delete)
		MOVE(FrameBuffer, default)

		~FrameBuffer()
		{
			glDeleteFramebuffers(1, &id);
		}

		struct Description
		{
			struct Attachment
			{
				GLenum type;
				Texture2D const & texture;
			};

			vector<Attachment> const & attachments;
		};

		void create(Description const & description)
		{
			glCreateFramebuffers(1, &id);

			for (auto const & attachment: description.attachments)
			{
				glNamedFramebufferTexture(
					id, attachment.type,
					attachment.texture.id, 0
				);
			}
		}
	};
}
