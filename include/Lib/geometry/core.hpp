#pragma once

#include <unordered_map>

#include "Lib/core/core.hpp"
#include "Lib/core/named.hpp"
#include "Lib/core/intrinsics.hpp"

namespace Geometry
{
	namespace Attribute
	{
		struct Key
		{
			enum class Common : u8
			{
				POSITION,
				NORMAL,
				TANGENT,
				TEXCOORD,
				COLOR,
				JOINTS,
				WEIGHTS,
			};

			variant<Common, std::string> name;
			u8 layer;

			bool operator==(Key const & other) const
			{
				return layer == other.layer and name == other.name;
			}

			struct Hasher
			{
				usize operator()(Key const & key) const
				{
					// TODO(bekorn): is this a good hash???
					return std::hash<decltype(key.layer)>{}(key.layer) ^ std::hash<decltype(key.name)>{}(key.name);
				}
			};
		};

		struct Type
		{
			enum class Value : u16
			{
				F32,
				I8, I16, I32, I8NORM, I16NORM, I32NORM,
				U8, U16, U32, U8NORM, U16NORM, U32NORM,
			};
			using enum Value;

			Value value;

			Type(Value value) :
				value(value)
			{}

			operator Value() const
			{ return value; }

			u8 size() const
			{
				switch (value)
				{
				case I8:
				case U8:
				case I8NORM:
				case U8NORM: return 1;
				case I16:
				case U16:
				case I16NORM:
				case U16NORM: return 2;
				case F32:
				case I32:
				case U32:
				case I32NORM:
				case U32NORM: return 4;
				}
				unreachable();
			}

			bool is_normalized() const
			{
				switch (value)
				{
				case I8NORM:
				case U8NORM:
				case I16NORM:
				case U16NORM:
				case I32NORM:
				case U32NORM: return true;
				case I8:
				case U8:
				case I16:
				case U16:
				case F32:
				case I32:
				case U32: return false;
				}
				unreachable();
			}
		};

		struct Data
		{
			Type type;
			u8 dimension;

			ByteBuffer buffer;
		};
	}

	struct Primitive
	{
		std::unordered_map<Attribute::Key, Attribute::Data, Attribute::Key::Hasher> attributes;
		vector<u32> indices;

		CTOR(Primitive, default);
		COPY(Primitive, delete);
		MOVE(Primitive, default);
	};
}

template<>
struct fmt::formatter<Geometry::Attribute::Key::Common> : formatter<std::string_view>
{
	template<typename FormatContext>
	auto format(Geometry::Attribute::Key::Common const & semantic, FormatContext & ctx)
	{
		using enum Geometry::Attribute::Key::Common;
		switch (semantic)
		{
		case POSITION: return formatter<std::string_view>::format("POSITION", ctx);
		case NORMAL: return formatter<std::string_view>::format("NORMAL", ctx);
		case TANGENT: return formatter<std::string_view>::format("TANGENT", ctx);
		case TEXCOORD: return formatter<std::string_view>::format("TEXCOORD", ctx);
		case COLOR: return formatter<std::string_view>::format("COLOR", ctx);
		case JOINTS: return formatter<std::string_view>::format("JOINTS", ctx);
		case WEIGHTS: return formatter<std::string_view>::format("WEIGHTS", ctx);
		}
		unreachable();
	}
};

template<>
struct fmt::formatter<Geometry::Attribute::Key>
{
	template <typename ParseContext>
	constexpr auto parse(ParseContext & ctx)
	{ return ctx.end(); }

	template <typename FormatContext>
	auto format(Geometry::Attribute::Key const & key, FormatContext & ctx)
	{
		if (std::holds_alternative<std::string>(key.name))
			return fmt::format_to(ctx.out(), "{}:{}", std::get<std::string>(key.name), key.layer);
		else
			return fmt::format_to(ctx.out(), "{}:{}", std::get<Geometry::Attribute::Key::Common>(key.name), key.layer);

		// auto out = ctx.out();
		// std::visit(
		// 	[&](auto const & name) { out = fmt::format_to(ctx.out(), "{}:{}", name, key.layer); },
		// 	key.name
		// );
		// return out;
	}
};

template<>
struct fmt::formatter<Geometry::Attribute::Type::Value> : formatter<std::string_view>
{
	template<typename FormatContext>
	auto format(Geometry::Attribute::Type::Value const & type, FormatContext & ctx)
	{
		using enum Geometry::Attribute::Type::Value;
		switch (type)
		{
		case F32: return formatter<std::string_view>::format("F32", ctx);
		case I8: return formatter<std::string_view>::format("I8", ctx);
		case I16: return formatter<std::string_view>::format("I16", ctx);
		case I32: return formatter<std::string_view>::format("I32", ctx);
		case I8NORM: return formatter<std::string_view>::format("I8NORM", ctx);
		case I16NORM: return formatter<std::string_view>::format("I16NORM", ctx);
		case I32NORM: return formatter<std::string_view>::format("I32NORM", ctx);
		case U8: return formatter<std::string_view>::format("U8", ctx);
		case U16: return formatter<std::string_view>::format("U16", ctx);
		case U32: return formatter<std::string_view>::format("U32", ctx);
		case U8NORM: return formatter<std::string_view>::format("U8NORM", ctx);
		case U16NORM: return formatter<std::string_view>::format("U16NORM", ctx);
		case U32NORM: return formatter<std::string_view>::format("U32NORM", ctx);
		}
	}
};
