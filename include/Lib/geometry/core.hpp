#pragma once

#include <unordered_map>
#include <map>

#include "Lib/core/core.hpp"
#include "Lib/core/named.hpp"
#include "Lib/core/intrinsics.hpp"

namespace Geometry
{
	namespace Attribute
	{
		struct Key
		{
			u32 begin_idx;
			u32 size;

			bool operator==(Key const & other) const
			{ return begin_idx == other.begin_idx; }

			struct Hasher
			{
				usize operator()(Key const & key) const
				{ return std::hash<u32>{}(key.begin_idx); }
			};
		};

		struct KeyPool
		{
			inline static vector<char> pool;
			inline static std::map<u64, Key> hash2key;

			static Key get_or_create_key(std::string_view const & name)
			{
				auto hash = std::hash<std::string_view>{}(name);

				auto it = hash2key.find(hash);
				if (it != hash2key.end())
					return it->second;

				Key key{
					.begin_idx = u32(pool.size()),
					.size = u32(name.size())
				};

				pool.reserve(pool.size() + name.size());
				pool.insert(pool.end(), name.begin(), name.end());

				hash2key.emplace(hash, key);
				return key;
			}

			static std::string_view get_name(Key key)
			{ return {pool.data() + key.begin_idx, key.size}; }
		};

		struct Common
		{
			inline static Key POSITION;
			inline static Key NORMAL;
			inline static Key TANGENT;
			inline static Key TEXCOORD;
			inline static Key COLOR;
			inline static Key JOINTS;
			inline static Key WEIGHTS;

			static void create()
			{
				Common::POSITION = KeyPool::get_or_create_key("POSITION");
				Common::NORMAL = KeyPool::get_or_create_key("NORMAL");
				Common::TANGENT = KeyPool::get_or_create_key("TANGENT");
				Common::TEXCOORD = KeyPool::get_or_create_key("TEXCOORD");
				Common::COLOR = KeyPool::get_or_create_key("COLOR");
				Common::JOINTS = KeyPool::get_or_create_key("JOINTS");
				Common::WEIGHTS = KeyPool::get_or_create_key("WEIGHTS");
			}
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
struct fmt::formatter<Geometry::Attribute::Key> : fmt::formatter<string_view>
{
	template <typename FormatContext>
	auto format(Geometry::Attribute::Key const & key, FormatContext & ctx)
	{
		auto pool_begin = Geometry::Attribute::KeyPool::pool.begin();
		string_view key_name(pool_begin + key.begin_idx, key.size);
		return formatter<string_view>::format(key_name, ctx);
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
