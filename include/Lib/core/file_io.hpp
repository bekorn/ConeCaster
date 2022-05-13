#pragma once

#include "core.hpp"

// This does not play well with precompiled headers, so it will be included in a regular header
#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

namespace File
{
	ByteBuffer LoadAsBytes(std::filesystem::path const & path)
	{
		assert(std::filesystem::exists(path));
		std::basic_ifstream<byte> file(path, std::ios::in | std::ios::binary | std::ios::ate);

		usize file_size = file.tellg();
		ByteBuffer buffer(file_size);

		file.seekg(0);
		file.read(buffer.data.get(), file_size);

		return buffer;
	}

	ByteBuffer LoadAsBytes(std::filesystem::path const & path, usize file_size)
	{
		assert(std::filesystem::exists(path));
		std::basic_ifstream<byte> file(path, std::ios::in | std::ios::binary);

		ByteBuffer buffer(file_size);
		file.read(buffer.data.get(), file_size);

		return buffer;
	}

	std::string LoadAsString(std::filesystem::path const & path)
	{
		assert(std::filesystem::exists(path));
		std::basic_ifstream<char> file(path, std::ios::in | std::ios::binary | std::ios::ate);

		usize file_size = file.tellg();
		std::string buffer(file_size, '\0');

		file.seekg(0);
		file.read(buffer.data(), file_size);

		return buffer;
	}

    namespace JSON
	{
		using JSONObj = rapidjson::Document::ConstObject const &;
		using Key = std::string_view const &;

		u32 GetU32(JSONObj obj, Key key)
		{
			return obj[key.data()].GetUint();
		}

		u32 GetU32(JSONObj obj, Key key, u32 def_value)
		{
			auto member = obj.FindMember(key.data());
			if (member != obj.MemberEnd())
				return member->value.GetUint();
			else
				return def_value;
		}

		optional<u32> GetOptionalU32(JSONObj obj, Key key)
		{
			auto member = obj.FindMember(key.data());
			if (member != obj.MemberEnd())
				return member->value.GetUint();
			else
				return nullopt;
		}

		std::string GetString(JSONObj obj, Key key, std::string const & def_value)
		{
			auto member = obj.FindMember(key.data());
			if (member != obj.MemberEnd())
				return member->value.GetString();
			else
				return def_value;
		}

		bool GetBool(JSONObj obj, Key key, bool def_value)
		{
			auto member = obj.FindMember(key.data());
			if (member != obj.MemberEnd())
				return member->value.GetBool();
			else
				return def_value;
		}

		f32 GetF32(JSONObj obj, Key key, f32 def_value)
		{
			auto member = obj.FindMember(key.data());
			if (member != obj.MemberEnd())
				return member->value.GetFloat();
			else
				return def_value;
		}

		f32x3 GetF32x3(JSONObj obj, Key key, f32x3 def_value)
		{
			auto member = obj.FindMember(key.data());
			if (member != obj.MemberEnd())
			{
				auto const & arr = member->value.GetArray();
				f32x3 val;
				for (auto i = 0; i < 3; ++i)
					val[i] = arr[i].GetFloat();
				return val;
			}
			else
				return def_value;
		}

		f32x4 GetF32x4(JSONObj obj, Key key, f32x4 def_value)
		{
			auto member = obj.FindMember(key.data());
			if (member != obj.MemberEnd())
			{
				auto const & arr = member->value.GetArray();
				f32x4 val;
				for (auto i = 0; i < 4; ++i)
					val[i] = arr[i].GetFloat();
				return val;
			}
			else
				return def_value;
		}

		struct NameGenerator
		{
			std::string const prefix;
			u64 counter = 0;

			std::string get(JSONObj obj, Key key)
			{
				auto member = obj.FindMember(key.data());
				if (member != obj.MemberEnd())
					return prefix + std::to_string(counter++) + ':' + member->value.GetString();
				else
					return prefix + std::to_string(counter++);
			}
		};
	}
}
