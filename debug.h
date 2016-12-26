#pragma once

#include "basetype.h"

namespace Gx {
     inline uint64_t rdtscp() {
        union {
            Gx::uint64_t as64;
            struct {
                Gx::uint32_t low;
                Gx::uint32_t high;
            };
        } ret{.as64 = 0ull};
        __asm__ __volatile__ ("rdtscp;": "=a" (ret.low), "=d" (ret.high) : : );
        return ret.as64;
    } 
}