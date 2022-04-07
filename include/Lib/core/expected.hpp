#pragma once

#include "core.hpp"

template<typename R, typename E>
struct Expected : variant<R, E>
{
	operator bool()
	{
		return std::holds_alternative<R>(*this);
	}

	R&& into_result()
	{
		return move(std::get<R>(*this));
	}

	E&& into_error()
	{
		return move(std::get<E>(*this));
	}
};
