#pragma once

#include "Lib/core/core.hpp"
#include "Lib/opengl/vao.hpp"
#include "Lib/geometry/core.hpp"

namespace Render
{
    struct Drawable
	{
		Geometry::Primitive const & primitive;
		Named<unique_one<IMaterial>> const named_material;

		GL::VertexArray vertex_array;

		bool is_loaded() const
		{
			return vertex_array.id != 0;
		}

		// TODO(bekorn): ShaderProgram should be accessible from the material
		void load(GL::ShaderProgram const & program)
		{
			vertex_array.create(primitive, program);
		}

		void unload()
		{
			vertex_array = {};
		}
	};

	struct Mesh
	{
		vector<Drawable> drawables;

		CTOR(Mesh, default)
		COPY(Mesh, delete)
		MOVE(Mesh, default)
	};
}