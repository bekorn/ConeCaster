#pragma once

#include <iostream>
#include <iomanip>
#include "fmt/core.h"
#include "fmt/chrono.h"

#include "core.hpp"

#include <chrono>

template<typename Duration = std::chrono::microseconds>
struct Timer
{
	// https://en.cppreference.com/w/cpp/chrono/high_resolution_clock see Notes paragraph 2
	using cpu_clock = std::chrono::steady_clock;
	using wall_clock = std::chrono::system_clock;

	cpu_clock::time_point cpu_last;
	wall_clock::time_point wall_last;

	Timer() :
		cpu_last(cpu_clock::now()),
		wall_last(wall_clock::now())
	{}

	struct TimeElapsed
	{
		Duration cpu;
		Duration wall;
	};

	TimeElapsed timeit()
	{
		auto cpu_now = cpu_clock::now();
		auto cpu_diff = std::chrono::duration_cast<Duration>(cpu_now - cpu_last);

		auto wall_now = wall_clock::now();
		auto wall_diff = std::chrono::duration_cast<Duration>(wall_now - wall_last);

		cpu_last = cpu_now;
		wall_last = wall_now;

		return {cpu_diff, wall_diff};
	};

	void timeit(std::ostream & out, std::string_view tag)
	{
		auto time_elapsed = timeit();

		out << std::right
			<< std::setw(32) << tag
			<< " | CPU:" << std::setw(16) << time_elapsed.cpu
			<< " | Wall:" << std::setw(16) << time_elapsed.wall
			<< '\n' << std::left;
	}

	void timeit(FILE * out, std::string_view tag)
	{
		auto time_elapsed = timeit();

		fmt::print(out,
			"{:>32} | CPU: {:>16} | Wall: {:>16}\n",
			tag, time_elapsed.cpu, time_elapsed.wall
		);
	}
};


// use this iterator to iterate an array of pointers as references
template<typename T>
struct PointerIterator
{
	T ** current;

	void operator++()
	{ current++; }

	T & operator*()
	{ return **current; }

	bool operator!=(PointerIterator const & other)
	{ return current != other.current; }
};


#include <random>

struct Random
{
	template<typename T>
	requires std::integral<T> || std::floating_point<T>
	static T next(T const & min, T const & max)
	{
		#if DEBUG
		static thread_local std::mt19937_64 generator(1337);
		#else
		static thread_local std::mt19937_64 generator(std::chrono::steady_clock::now().time_since_epoch().count());
		#endif
		
		if constexpr (std::is_integral_v<T>)
			return std::uniform_int_distribution<T>(min, max)(generator);
		else
			return std::uniform_real_distribution<T>(min, max)(generator);
	}

	template<glm::length_t L, typename T, glm::qualifier Q>
	static glm::vec<L, T, Q> next(glm::vec<L, T, Q> const & min, glm::vec<L, T, Q> const & max)
	{
		glm::vec<L, T, Q> result;
		for (auto i = 0; i < L; ++i)
			result[i] = next(min[i], max[i]);
		return result;
	}

	static f32x3 next_in_unit_cube()
	{
		return next(f32x3{-1, -1, -1}, f32x3{1, 1, 1});
	}

	// for the next 2 functions see https://datagenetics.com/blog/january32020/index.html
	static f32x3 next_in_unit_sphere()
	{
		// TODO(bekorn): compare with the discarding method 

		auto theta = next(f32(0), glm::two_pi<f32>());
		auto phi = acos(next<f32>(-1, 1));
		auto r = pow(next<f32>(0, 1), 1.f / 3.f);
		auto x = r * sin(phi) * cos(theta);
		auto y = r * sin(phi) * sin(theta);
		auto z = r * cos(phi);
		return {x, y, z};
	}

	static f32x3 next_on_unit_sphere()
	{
		// return normalize(next(f32x3(-1), f32x3(1)));
		return normalize(next_in_unit_sphere());
	}
};


#include <source_location>

// TODO(bekorn): this is demonstrative but still much better than depending on macros
void log(std::string_view message, std::source_location const location = std::source_location::current())
{
	std::cout
		<< "LOG "
		<< location.file_name() << "(" << location.line() << ") " << location.function_name() << ':'
		<< message
		<< '\n';
}
