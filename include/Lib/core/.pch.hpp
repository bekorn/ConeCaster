#pragma once

// {fmt}
#include "fmt/core.h"
#include "fmt/format.h"
#include "fmt/chrono.h"

// Configure GLM
#define GLM_FORCE_NO_CTOR_INIT
#define GLM_FORCE_EXPLICIT_CTOR

// forward declarations to reduce compile times
#include <glm/fwd.hpp>

// source file
#include <glm/glm.hpp>

// extensions (http://glm.g-truc.net/glm.pdf)
#include <glm/gtc/quaternion.hpp>
#include <glm/gtc/type_precision.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/gtx/euler_angles.hpp>
#include <glm/gtx/matrix_decompose.hpp>
#include <glm/gtx/range.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtx/string_cast.hpp>
#include <glm/gtc/random.hpp>
#include <glm/gtx/norm.hpp>

// File IO
#include <fstream>
#include <filesystem>
#include <rapidjson/document.h>

// Handy std classes
#include <array>
#include <memory>
#include <span>
#include <optional>
#include <variant>
#include <vector>
#include <unordered_map>
