#pragma once

#include <unordered_map>
#include <ostream>

#include "Lib/core/core.hpp"
#include "Lib/core/intrinsics.hpp"
#include "Lib/core/named.hpp"

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

		std::ostream & operator<<(std::ostream & os, Key::Common semantic)
		{
			using enum Key::Common;
			switch (semantic)
			{
			case POSITION: return os << "POSITION";
			case NORMAL: return os << "NORMAL";
			case TANGENT: return os << "TANGENT";
			case TEXCOORD: return os << "TEXCOORD";
			case COLOR: return os << "COLOR";
			case JOINTS: return os << "JOINTS";
			case WEIGHTS: return os << "WEIGHTS";
			}
			unreachable();
		}

		std::ostream & operator<<(std::ostream & os, Key const & key)
		{
			std::visit([&os](auto const & s) { os << s; }, key.name);
			return os << ':' << static_cast<u32>(key.layer);
		}

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
			}
		};

		std::ostream & operator<<(std::ostream & os, Type::Value type)
		{
			using enum Type::Value;
			switch (type)
			{
			case F32: return os << "F32";
			case I8: return os << "I8";
			case I16: return os << "I16";
			case I32: return os << "I32";
			case I8NORM: return os << "I8NORM";
			case I16NORM: return os << "I16NORM";
			case I32NORM: return os << "I32NORM";
			case U8: return os << "U8";
			case U16: return os << "U16";
			case U32: return os << "U32";
			case U8NORM: return os << "U8NORM";
			case U16NORM: return os << "U16NORM";
			case U32NORM: return os << "U32NORM";
			}
			unreachable();
		}

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

