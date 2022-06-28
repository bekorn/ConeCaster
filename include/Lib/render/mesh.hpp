#pragma once

#include "Lib/core/core.hpp"
#include "Lib/geometry/core.hpp"

namespace Render
{
    struct Drawable
	{
		Geometry::Primitive const & primitive;
		Named<unique_one<IMaterial>> const named_material;
	};

	struct Mesh
	{
		vector<Drawable> drawables;

		CTOR(Mesh, default)
		COPY(Mesh, delete)
		MOVE(Mesh, default)
	};
}