#pragma once

namespace Gx {
	using nullptr_t = decltype(nullptr);
	using size_t = decltype(sizeof 0);
	using int8_t = signed char;
	using uint8_t = unsigned char;
	using int16_t = signed short;
	using uint16_t = unsigned short;
	using int32_t = signed int;
	using uint32_t = unsigned int;
	using int64_t = signed long long;
	using uint64_t = unsigned long long;

	template<size_t size> struct signed_size_t;
	template<size_t size> struct unsigned_size_t;
	template<> struct signed_size_t<1> { using type = int8_t; };
	template<> struct signed_size_t<2> { using type = int16_t; };
	template<> struct signed_size_t<4> { using type = int32_t; };
	template<> struct signed_size_t<8> { using type = int64_t; };

	template<> struct unsigned_size_t<1> { using type = uint8_t; };
	template<> struct unsigned_size_t<2> { using type = uint16_t; };
	template<> struct unsigned_size_t<4> { using type = uint32_t; };
	template<> struct unsigned_size_t<8> { using type = uint64_t; };

	using ssize_t = signed_size_t<sizeof(size_t)>::type;

	using intptr_t = signed_size_t<sizeof(void*)>::type;
	using uintptr_t = unsigned_size_t<sizeof(void*)>::type;

	constexpr size_t BitsPerByte = 8u;

	static_assert(sizeof(int8_t) == 1, "");
	static_assert(sizeof(uint8_t) == 1, "");

	static_assert(sizeof(int16_t) == 2, "");
	static_assert(sizeof(uint16_t) == 2, "");

	static_assert(sizeof(int32_t) == 4, "");
	static_assert(sizeof(uint32_t) == 4, "");

	static_assert(sizeof(int64_t) == 8, "");
	static_assert(sizeof(uint64_t) == 8, "");

	template<typename T, typename U>
	constexpr auto offsetOf(U const T::*member) {
		return reinterpret_cast<size_t >(&(reinterpret_cast<T *>(0u)->*member));
	}

	template<typename T> constexpr
	T const& max(T const& a, T const& b) {
		return a > b ? a : b;
	}

	template<typename T> constexpr
	T const& min(T const& a, T const& b) {
		return a < b ? a : b;
	}
}
