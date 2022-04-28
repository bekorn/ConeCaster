#include <iostream>

#include "Lib/core/.hpp"
#include "Lib/display/.hpp"
#include "Lib/opengl/.hpp"
#include "Lib/render/.hpp"

#include "editor.hpp"

i32 main()
{
    GLFW::Context glfw_context;
    if (auto error = glfw_context.create())
    {
        fmt::print(stderr, "{}", error.value());
        return 1;
    }

    GLFW::Window window;
    if (auto error = window.create({
        .title = "Cone Caster",
        .size = {1600, 700},
        .vsync = true,
        .gl_major = GL::VERSION_MAJOR, .gl_minor = GL::VERSION_MINOR,
    }))
    {
        fmt::print(stderr, "{}", error.value());
        return 1;
    }

    fmt::print(stderr, "{}\n", GL::GetContextInfo());

    Imgui::Context imgui_context;
    imgui_context.create({.window = window});

	Render::FrameInfo frame_info;
	Render::FrameInfo previous_frame_info{
		.idx = 0,
		.seconds_since_start = glfwGetTime(),
		.seconds_since_last_frame = 0,
	};

	Render::Renderer renderer;
	renderer.camera.create({
		.pos = {0, -500, 0},
		.target = {0, 0, 0},
		.up = {0, 0, 1},
		.vfov = 60
	});
	renderer.create({512, 512});

	Editor editor(renderer);

	while (not glfwWindowShouldClose(window))
	{
		glfwPollEvents();

		frame_info.seconds_since_start = glfwGetTime();
		frame_info.seconds_since_last_frame = frame_info.seconds_since_start - previous_frame_info.seconds_since_start;
		frame_info.idx = previous_frame_info.idx + 1;

		{
			using namespace GL;

			// bind and clear defautl framebuffer
			glBindFramebuffer(GL_FRAMEBUFFER, 0);

			i32 w, h; glfwGetWindowSize(window, &w, &h);
			glViewport(0, 0, w, h);

			static f32x4 clear_color(0, 0, 0, 1);
			glClearNamedFramebufferfv(0, GL_COLOR, 0, begin(clear_color));
			static f32 clear_depth = 0;
			glClearNamedFramebufferfv(0, GL_DEPTH, 0, &clear_depth);
		}
		
		// Start ImGui frame
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();
		
		// Editor
		editor.render(frame_info);

		//	ConeCaster
		renderer.update(frame_info);
		renderer.render(frame_info);

		// Render Imgui frame
		ImGui::Render();
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		glfwSwapBuffers(window);

		previous_frame_info = frame_info;
	}

    return 0;
}