//
// Created by vfernandez on 19/03/19.
//

#ifndef TRIMAL_HASHER_H
#define TRIMAL_HASHER_H
#include <cstdint>

// https://gist.github.com/Lee-R/3839813

namespace detail
{
    // FNV-1a 32bit hashing algorithm.
    constexpr std::uint32_t fnv1a_32(char const* s, std::size_t count)
    {
        return ((count ? fnv1a_32(s, count - 1) : 2166136261u) ^ s[count]) * 16777619u;
    }
}    // namespace detail

constexpr std::uint32_t operator"" _hash(char const* s, std::size_t count)
{
    return detail::fnv1a_32(s, count);
}



#endif //TRIMAL_HASHER_H
