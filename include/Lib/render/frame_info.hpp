#pragma once

#include "Lib/core/core.hpp"

namespace Render
{
	struct FrameInfo
	{
		u64 idx;
		f64 seconds_since_start;
		f64 seconds_since_last_frame;
	};
}