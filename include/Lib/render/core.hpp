#pragma once

#include "Lib/core/.hpp"

namespace Render
{

template<std::totally_ordered T>
struct Range
{
    T min, max;

    bool contains(T const & val) const
    { return min < val & val < max; }
};

struct Ray
{
    f32x3 pos;
    f32x3 dir;

    f32x3 color;
    u32x2 pixel;

    u32 bounce;

    f32x3 at(f32 t) const
    { return pos + dir * t; }
};

}