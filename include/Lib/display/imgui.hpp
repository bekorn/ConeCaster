#pragma once

#include "Lib/opengl/glsl.hpp"

#include ".pch.hpp"

namespace Imgui
{
	struct Context
	{
		Context() = default;

		struct Description
		{
			GLFW::Window const & window;
		};

		void create(Description const & description)
		{
			IMGUI_CHECKVERSION();
			ImGui::CreateContext();

			// Setup Dear ImGui style
			ImGui::StyleColorsDark();

			// Setup Platform/Renderer backends
			ImGui_ImplGlfw_InitForOpenGL(description.window, true);

			ImGui_ImplOpenGL3_Init(GL::GLSL_VERSION_MACRO.data());
		}

		~Context()
		{
			ImGui_ImplOpenGL3_Shutdown();
			ImGui_ImplGlfw_Shutdown();
			ImGui::DestroyContext();
		}
	};
}

