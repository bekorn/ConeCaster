#pragma once

#include "meta.hpp"

[[noreturn]]
inline void constexpr unreachable() noexcept
{
#ifdef Compiler_MSVC
	// https://docs.microsoft.com/en-us/cpp/intrinsics/assume
	__assume(0);
#else
	// https://clang.llvm.org/docs/LanguageExtensions.html#builtin-unreachable
	__builtin_unreachable();
#endif
}