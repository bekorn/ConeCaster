#pragma once

#include "Lib/core/.hpp"

namespace Render
{

struct Ray
{
    f32x3 pos, dir;
    f32 min, max;

    f32x3 color;
    u32x2 pixel;

    u32 bounce;

    f32x3 at(f32 t) const
    { return pos + dir * t; }
};

}