#pragma once

#include "Lib/core/core.hpp"

#include ".pch.hpp"

// To run the application with the external GPU, https://stackoverflow.com/a/39047129/2073225
extern "C"
{
__declspec(dllexport) i64 NvOptimusEnablement = 1;
__declspec(dllexport) i32 AmdPowerXpressRequestHighPerformance = 1;
}

namespace GL
{
	using namespace gl45core;
	using namespace gl45ext;

	auto constexpr VERSION_MAJOR = 4;
	auto constexpr VERSION_MINOR = 5;

	// following DSA API https://github.com/fendevel/Guide-to-Modern-OpenGL-Functions

	struct OpenGLObject
	{
		u32 id = 0;

		OpenGLObject() noexcept= default;

		OpenGLObject(OpenGLObject const &) noexcept = delete;
		OpenGLObject& operator=(OpenGLObject const &) noexcept = delete;

		OpenGLObject(OpenGLObject && other) noexcept
		{
			std::swap(id, other.id);
		}
		OpenGLObject& operator=(OpenGLObject && other) noexcept
		{
			std::swap(id, other.id);
			return *this;
		};
	};
}