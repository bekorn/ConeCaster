#pragma once

#include "Lib/geometry/core.hpp"

namespace GL
{
	std::string_view GLSLTypeToString(GLenum type)
	{
		using namespace std::string_view_literals;

		switch (type)
		{
		// https://www.khronos.org/registry/OpenGL-Refpages/gl4/html/glGetActiveUniform.xhtml
		case (GL_FLOAT): return "float"sv;
		case (GL_FLOAT_VEC2): return "vec2"sv;
		case (GL_FLOAT_VEC3): return "vec3"sv;
		case (GL_FLOAT_VEC4): return "vec4"sv;
		case (GL_DOUBLE): return "double"sv;
		case (GL_DOUBLE_VEC2): return "dvec2"sv;
		case (GL_DOUBLE_VEC3): return "dvec3"sv;
		case (GL_DOUBLE_VEC4): return "dvec4"sv;
		case (GL_INT): return "int"sv;
		case (GL_INT_VEC2): return "ivec2"sv;
		case (GL_INT_VEC3): return "ivec3"sv;
		case (GL_INT_VEC4): return "ivec4"sv;
		case (GL_UNSIGNED_INT): return "unsigned int"sv;
		case (GL_UNSIGNED_INT_VEC2): return "uvec2"sv;
		case (GL_UNSIGNED_INT_VEC3): return "uvec3"sv;
		case (GL_UNSIGNED_INT_VEC4): return "uvec4"sv;
		case (GL_BOOL): return "bool"sv;
		case (GL_BOOL_VEC2): return "bvec2"sv;
		case (GL_BOOL_VEC3): return "bvec3"sv;
		case (GL_BOOL_VEC4): return "bvec4"sv;
		case (GL_FLOAT_MAT2): return "mat2"sv;
		case (GL_FLOAT_MAT3): return "mat3"sv;
		case (GL_FLOAT_MAT4): return "mat4"sv;
		case (GL_FLOAT_MAT2x3): return "mat2x3"sv;
		case (GL_FLOAT_MAT2x4): return "mat2x4"sv;
		case (GL_FLOAT_MAT3x2): return "mat3x2"sv;
		case (GL_FLOAT_MAT3x4): return "mat3x4"sv;
		case (GL_FLOAT_MAT4x2): return "mat4x2"sv;
		case (GL_FLOAT_MAT4x3): return "mat4x3"sv;
		case (GL_DOUBLE_MAT2): return "dmat2"sv;
		case (GL_DOUBLE_MAT3): return "dmat3"sv;
		case (GL_DOUBLE_MAT4): return "dmat4"sv;
		case (GL_DOUBLE_MAT2x3): return "dmat2x3"sv;
		case (GL_DOUBLE_MAT2x4): return "dmat2x4"sv;
		case (GL_DOUBLE_MAT3x2): return "dmat3x2"sv;
		case (GL_DOUBLE_MAT3x4): return "dmat3x4"sv;
		case (GL_DOUBLE_MAT4x2): return "dmat4x2"sv;
		case (GL_DOUBLE_MAT4x3): return "dmat4x3"sv;
		case (GL_SAMPLER_1D): return "sampler1D"sv;
		case (GL_SAMPLER_2D): return "sampler2D"sv;
		case (GL_SAMPLER_3D): return "sampler3D"sv;
		case (GL_SAMPLER_CUBE): return "samplerCube"sv;
		case (GL_SAMPLER_1D_SHADOW): return "sampler1DShadow"sv;
		case (GL_SAMPLER_2D_SHADOW): return "sampler2DShadow"sv;
		case (GL_SAMPLER_1D_ARRAY): return "sampler1DArray"sv;
		case (GL_SAMPLER_2D_ARRAY): return "sampler2DArray"sv;
		case (GL_SAMPLER_1D_ARRAY_SHADOW): return "sampler1DArrayShadow"sv;
		case (GL_SAMPLER_2D_ARRAY_SHADOW): return "sampler2DArrayShadow"sv;
		case (GL_SAMPLER_2D_MULTISAMPLE): return "sampler2DMS"sv;
		case (GL_SAMPLER_2D_MULTISAMPLE_ARRAY): return "sampler2DMSArray"sv;
		case (GL_SAMPLER_CUBE_SHADOW): return "samplerCubeShadow"sv;
		case (GL_SAMPLER_BUFFER): return "samplerBuffer"sv;
		case (GL_SAMPLER_2D_RECT): return "sampler2DRect"sv;
		case (GL_SAMPLER_2D_RECT_SHADOW): return "sampler2DRectShadow"sv;
		case (GL_INT_SAMPLER_1D): return "isampler1D"sv;
		case (GL_INT_SAMPLER_2D): return "isampler2D"sv;
		case (GL_INT_SAMPLER_3D): return "isampler3D"sv;
		case (GL_INT_SAMPLER_CUBE): return "isamplerCube"sv;
		case (GL_INT_SAMPLER_1D_ARRAY): return "isampler1DArray"sv;
		case (GL_INT_SAMPLER_2D_ARRAY): return "isampler2DArray"sv;
		case (GL_INT_SAMPLER_2D_MULTISAMPLE): return "isampler2DMS"sv;
		case (GL_INT_SAMPLER_2D_MULTISAMPLE_ARRAY): return "isampler2DMSArray"sv;
		case (GL_INT_SAMPLER_BUFFER): return "isamplerBuffer"sv;
		case (GL_INT_SAMPLER_2D_RECT): return "isampler2DRect"sv;
		case (GL_UNSIGNED_INT_SAMPLER_1D): return "usampler1D"sv;
		case (GL_UNSIGNED_INT_SAMPLER_2D): return "usampler2D"sv;
		case (GL_UNSIGNED_INT_SAMPLER_3D): return "usampler3D"sv;
		case (GL_UNSIGNED_INT_SAMPLER_CUBE): return "usamplerCube"sv;
		case (GL_UNSIGNED_INT_SAMPLER_1D_ARRAY): return "usampler2DArray"sv;
		case (GL_UNSIGNED_INT_SAMPLER_2D_ARRAY): return "usampler2DArray"sv;
		case (GL_UNSIGNED_INT_SAMPLER_2D_MULTISAMPLE): return "usampler2DMS"sv;
		case (GL_UNSIGNED_INT_SAMPLER_2D_MULTISAMPLE_ARRAY): return "usampler2DMSArray"sv;
		case (GL_UNSIGNED_INT_SAMPLER_BUFFER): return "usamplerBuffer"sv;
		case (GL_UNSIGNED_INT_SAMPLER_2D_RECT): return "usampler2DRect"sv;
		case (GL_IMAGE_1D): return "image1D"sv;
		case (GL_IMAGE_2D): return "image2D"sv;
		case (GL_IMAGE_3D): return "image3D"sv;
		case (GL_IMAGE_2D_RECT): return "image2DRect"sv;
		case (GL_IMAGE_CUBE): return "imageCube"sv;
		case (GL_IMAGE_BUFFER): return "imageBuffer"sv;
		case (GL_IMAGE_1D_ARRAY): return "image1DArray"sv;
		case (GL_IMAGE_2D_ARRAY): return "image2DArray"sv;
		case (GL_IMAGE_2D_MULTISAMPLE): return "image2DMS"sv;
		case (GL_IMAGE_2D_MULTISAMPLE_ARRAY): return "image2DMSArray"sv;
		case (GL_INT_IMAGE_1D): return "iimage1D"sv;
		case (GL_INT_IMAGE_2D): return "iimage2D"sv;
		case (GL_INT_IMAGE_3D): return "iimage3D"sv;
		case (GL_INT_IMAGE_2D_RECT): return "iimage2DRect"sv;
		case (GL_INT_IMAGE_CUBE): return "iimageCube"sv;
		case (GL_INT_IMAGE_BUFFER): return "iimageBuffer"sv;
		case (GL_INT_IMAGE_1D_ARRAY): return "iimage1DArray"sv;
		case (GL_INT_IMAGE_2D_ARRAY): return "iimage2DArray"sv;
		case (GL_INT_IMAGE_2D_MULTISAMPLE): return "iimage2DMS"sv;
		case (GL_INT_IMAGE_2D_MULTISAMPLE_ARRAY): return "iimage2DMSArray"sv;
		case (GL_UNSIGNED_INT_IMAGE_1D): return "uimage1D"sv;
		case (GL_UNSIGNED_INT_IMAGE_2D): return "uimage2D"sv;
		case (GL_UNSIGNED_INT_IMAGE_3D): return "uimage3D"sv;
		case (GL_UNSIGNED_INT_IMAGE_2D_RECT): return "uimage2DRect"sv;
		case (GL_UNSIGNED_INT_IMAGE_CUBE): return "uimageCube"sv;
		case (GL_UNSIGNED_INT_IMAGE_BUFFER): return "uimageBuffer"sv;
		case (GL_UNSIGNED_INT_IMAGE_1D_ARRAY): return "uimage1DArray"sv;
		case (GL_UNSIGNED_INT_IMAGE_2D_ARRAY): return "uimage2DArray"sv;
		case (GL_UNSIGNED_INT_IMAGE_2D_MULTISAMPLE): return "uimage2DMS"sv;
		case (GL_UNSIGNED_INT_IMAGE_2D_MULTISAMPLE_ARRAY): return "uimage2DMSArray"sv;
		case (GL_UNSIGNED_INT_ATOMIC_COUNTER): return "atomic_uint"sv;

		// https://www.khronos.org/registry/OpenGL/extensions/ARB/ARB_gpu_shader_int64.txt
		case (GL_INT64_ARB): return "int64"sv;
		case (GL_INT64_VEC2_ARB): return "i64vec2"sv;
		case (GL_INT64_VEC3_ARB): return "i64vec3"sv;
		case (GL_INT64_VEC4_ARB): return "i64vec4"sv;
		case (GL_UNSIGNED_INT64_ARB): return "uint64"sv;
		case (GL_UNSIGNED_INT64_VEC2_ARB): return "u64vec2"sv;
		case (GL_UNSIGNED_INT64_VEC3_ARB): return "u64vec3"sv;
		case (GL_UNSIGNED_INT64_VEC4_ARB): return "u64vec4"sv;

		default: return "unknown type :("sv;
		}
	}

	u32 ComponentTypeSize(GLenum type)
	{
		// https://www.khronos.org/opengl/wiki/Vertex_Specification#Component_type
		switch (type)
		{
		case GL_BYTE: return sizeof(GLbyte);
		case GL_SHORT: return sizeof(GLshort);
		case GL_INT: return sizeof(GLint);

		case GL_UNSIGNED_BYTE: return sizeof(GLubyte);
		case GL_UNSIGNED_SHORT: return sizeof(GLushort);
		case GL_UNSIGNED_INT: return sizeof(GLuint);

		case GL_FIXED: return sizeof(GLfixed);
		case GL_HALF_FLOAT: return sizeof(GLhalf);
		case GL_FLOAT: return sizeof(GLfloat);
		case GL_DOUBLE: return sizeof(GLdouble);

		default: throw std::runtime_error("Unknown component Type :(");
		}
	}

	GLenum IntoGLenum(Geometry::Attribute::Type::Value type)
	{
		using enum Geometry::Attribute::Type::Value;
		switch (type)
		{
		case F32: return GL_FLOAT;
		case I8:
		case I8NORM: return GL_BYTE;
		case I16:
		case I16NORM: return GL_SHORT;
		case I32:
		case I32NORM: return GL_INT;
		case U8:
		case U8NORM: return GL_UNSIGNED_BYTE;
		case U16:
		case U16NORM: return GL_UNSIGNED_SHORT;
		case U32:
		case U32NORM: return GL_UNSIGNED_INT;
		}
		unreachable();
	}

	Geometry::Attribute::Type IntoAttributeType(GLenum type, bool is_normalized)
	{
		using enum Geometry::Attribute::Type::Value;

		if (is_normalized)
			switch (type)
			{
			case GL_BYTE: return I8NORM;
			case GL_SHORT: return I16NORM;
			case GL_INT: return I32NORM;
			case GL_UNSIGNED_BYTE: return U8NORM;
			case GL_UNSIGNED_SHORT: return U16NORM;
			case GL_UNSIGNED_INT: return U32NORM;
			}
		else
			switch (type)
			{
			case GL_FLOAT: return F32;
			case GL_BYTE: return I8;
			case GL_SHORT: return I16;
			case GL_INT: return I32;
			case GL_UNSIGNED_BYTE: return U8;
			case GL_UNSIGNED_SHORT: return U16;
			case GL_UNSIGNED_INT: return U32;
			}

		throw std::runtime_error("Unknown type or combination :(");
	}
}