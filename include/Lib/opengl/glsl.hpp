#pragma once

namespace GL
{
	auto const GLSL_VERSION_MACRO = "#version " + std::to_string(VERSION_MAJOR) + std::to_string(VERSION_MINOR) + "0\n";
}
