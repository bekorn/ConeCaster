#pragma once

#include <unordered_map>
#include <iomanip>

#include "core.hpp"

// inspired by https://github.com/skypjack/entt/blob/master/src/entt/core/hashed_string.hpp
struct Name
{
	usize hash;
	std::string string;

	CTOR(Name, default)
	COPY(Name, default)
	MOVE(Name, default)

	Name(std::string_view const & sv) :
		hash(std::hash<std::string_view>{}(sv)), string(sv)
	{}

	Name(std::string const & str) :
		hash(std::hash<std::string>{}(str)), string(str)
	{}

	bool operator==(Name const & other) const
	{ return hash == other.hash; }

	struct Hasher
	{
		usize operator()(Name const & name) const
		{ return name.hash; }
	};

	friend std::ostream & operator<<(std::ostream & out, Name const & name)
	{ return out << std::right << std::setw(20) << name.hash << std::left << '|' << name.string; }
};

Name operator ""_name(char const * literal, usize size)
{
	return std::string_view(literal, size);
}

template<typename T>
struct Named
{
	Name const & name;
	T & data;

	bool operator==(Named const & other) const
	{ return name == other.name; }
};

template<typename T>
struct Managed
{
	std::unordered_map<Name, T, Name::Hasher> resources;

	template<typename... Args>
	Named<T> generate(Name const & name, Args && ... args)
	{
		auto[it, is_emplaced] = resources.try_emplace(name, std::forward<Args>(args)...);
		if (not is_emplaced)
			std::cerr << "!! Resource is not emplaced: " << name << '\n';
		return {it->first, it->second};
	}

	void erase(Name const & name)
	{ resources.erase(name); }

	T & get(Name const & name)
	{ return resources.at(name); }

	T const & get(Name const & name) const
	{ return resources.at(name); }

	Named<T> get_named(Name const & name)
	{
		auto it = resources.find(name);
		return {it->first, it->second};
	}

	T & get_or_generate(Name const & name)
	{ return resources[name]; }

	auto contains(Name const & name) const
	{ return resources.contains(name); }

	auto find(Name const & name) const
	{ return resources.find(name); }

	auto begin()
	{ return resources.begin(); }

	auto end()
	{ return resources.end(); }
};