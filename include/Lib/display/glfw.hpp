#pragma once

#include ".pch.hpp"

#include "Lib/core/expected.hpp"
#include "Lib/opengl/core.hpp"
#include "Lib/opengl/util.hpp"

namespace GLFW
{
	struct Context
	{
		static void error_callback(i32 error, const char* description)
		{
			fmt::print(stderr, "GLFW Error {}: \n", error, description);
		}

		CTOR(Context, default)
		COPY(Context, delete)
		MOVE(Context, default)

		optional<std::string> create()
		{
			glfwSetErrorCallback(error_callback);

			if (glfwInit())
				fmt::print(stderr, "GLFW initialized\n");
			else
				return "GLFW failed to initialize";

			return {};
		}

		~Context()
		{
			glfwTerminate();
		}
	};
	
	struct Window
	{
		GLFWwindow * glfw_window = nullptr;

		CTOR(Window, default)
		COPY(Window, delete)
		MOVE(Window, default)

		operator GLFWwindow *() const
		{
			return glfw_window;
		}

		static void resize_callback(GLFWwindow*, int width, int height)
		{
			GL::glViewport(0, 0, width, height);
		}

		struct Description
		{
			std::string_view title;
			i32x2 size = {720, 720};
			bool vsync = false;

			i32 const & gl_major;
			i32 const & gl_minor;
		};

		[[nodiscard]]
		optional<std::string> create(Description const & description)
		{
			glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, description.gl_major);
			glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, description.gl_minor);
			glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
			glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, true);
			if constexpr (GL_DEBUG)
				glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, true);

			glfw_window = glfwCreateWindow(
				description.size.x, description.size.y,
				description.title.data(), nullptr, nullptr
			);

			if (glfw_window != nullptr)
				fmt::print(stderr, "GLFW Window created\n");
			else
				return "GLFW Window failed to create";

			glfwMakeContextCurrent(glfw_window);

			glfwSwapInterval(description.vsync ? 1 : 0);

			// Load OpenGL API
			glbinding::initialize(glfwGetProcAddress);

			// Set callbacks
			if constexpr (GL_DEBUG)
			{
				using namespace GL;

				i32 flags = 0;
				glGetIntegerv(GL_CONTEXT_FLAGS, &flags);
				// TODO(bekorn): when run inside RenderDoc, can not create a debug context
				assert(flags & static_cast<i32>(GL_CONTEXT_FLAG_DEBUG_BIT));

				glEnable(GL_DEBUG_OUTPUT);
				glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
				glDebugMessageCallback(DebugCallback, nullptr);
				glDebugMessageControl(GL_DONT_CARE, GL_DONT_CARE, GL_DONT_CARE, 0, nullptr, GL_TRUE);
				
				fmt::print(stderr, "GL Debug initialized\n");
			}

			glfwSetWindowSizeCallback(glfw_window, resize_callback);

			return {};
		}

		~Window()
		{
			glfwDestroyWindow(glfw_window);
		}
	};
}
