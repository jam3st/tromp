#pragma once

#include "basetype.h"
#include "debug.h"
#include <immintrin.h>
#include <memory.h>

namespace Gx {
       namespace Blake2b {
              static constexpr u_int8_t alignmentBits = 64u;
              static constexpr auto BlakeHashSize64Bits = 8u; 
              static constexpr auto BlakeBlockSize64Bits = 16u;
              static constexpr auto HeaderLen64Bits = 18u;

              template<uint32_t size64Bits>
              class Chunk {
                      public:
                            uint64_t const& operator[](size_t const index) const {
                                   return v64[index];
                            }                 

                            uint64_t& operator[](size_t const index) {
                                   return v64[index];
                            }

                            uint32_t get32(size_t const index) const {
                                   return v32[index];
                            }

                            void set32(size_t const index, uint32_t const val)  {
                                   v32[index] = val;
                            }

                            uint32_t get8(size_t const index) const {
                                   return v8[index];
                            }
                            
                     private:
                            union __attribute__((__aligned__(alignmentBits))) {
                                   uint64_t v64[size64Bits];
                                   uint32_t v32[size64Bits * 2u];
                                   uint16_t v16[size64Bits * 4u];
                                   uint8_t v8[size64Bits * 8u];
                            };
              };

              using Hash = Chunk<BlakeHashSize64Bits>;
              using Block = Chunk<BlakeBlockSize64Bits>;
              using Header = Chunk<HeaderLen64Bits>;

              static constexpr auto HashLengthBlake = 0x01010000ull | 50ull;
              static constexpr auto ZchashPow = 0x576f50687361635aull;
              static constexpr auto ZchashParam = 0x9000000c8ull;
              static constexpr auto MidCount = 128ull;
              static constexpr auto FullCount = 144ull;

              static constexpr uint8_t Rounds = 12u;
              static constexpr uint8_t Sigma[Rounds][BlakeBlockSize64Bits] = {
                     {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15 },
                     { 14, 10,  4,  8,  9, 15, 13,  6,  1, 12,  0,  2, 11,  7,  5,  3 },
                     { 11,  8, 12,  0,  5,  2, 15, 13, 10, 14,  3,  6,  7,  1,  9,  4 },
                     {  7,  9,  3,  1, 13, 12, 11, 14,  2,  6,  5, 10,  4,  0, 15,  8 },
                     {  9,  0,  5,  7,  2,  4, 10, 15, 14,  1, 11, 12,  6,  8,  3, 13 },
                     {  2, 12,  6, 10,  0, 11,  8,  3,  4, 13,  7,  5, 15, 14,  1,  9 },
                     { 12,  5,  1, 15, 14, 13,  4, 10,  0,  7,  6,  3,  9,  2,  8, 11 },
                     { 13, 11,  7, 14, 12,  1,  3,  9,  5,  0, 15,  4,  8,  6,  2, 10 },
                     {  6, 15, 14,  9, 11,  3,  0,  8, 12,  2, 13,  7,  1,  4, 10,  5 },
                     { 10,  2,  8,  4,  7,  6,  1,  5, 15, 11,  9, 14,  3, 12, 13,  0 },
                     {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15 },
                     { 14, 10,  4,  8,  9, 15, 13,  6,  1, 12,  0,  2, 11,  7,  5,  3 }
              };

              static constexpr uint64_t IV[] = { 
                     0x6a09e667f3bcc908ull, 0xbb67ae8584caa73bull, 0x3c6ef372fe94f82bull, 0xa54ff53a5f1d36f1ull,
                     0x510e527fade682d1ull, 0x9b05688c2b3e6c1full, 0x1f83d9abfb41bd6bull, 0x5be0cd19137e2179ull
              };

              static inline uint64_t rotr64(uint64_t a, uint8_t bits) {
                     return (a >> bits) | (a << (64u - bits));
              }

              static inline void mix(uint64_t *va, uint64_t *vb, uint64_t *vc, uint64_t *vd, uint64_t const x, uint64_t const y) {
                     *va = (*va + *vb + x);
                     *vd = rotr64(*vd ^ *va, 32);
                     *vc = (*vc + *vd);
                     *vb = rotr64(*vb ^ *vc, 24);
                     *va = (*va + *vb + y);
                     *vd = rotr64(*vd ^ *va, 16);
                     *vc = (*vc + *vd);
                     *vb = rotr64(*vb ^ *vc, 63);
              }

              Hash initMidstate(Block const& m) {
                     Hash tmp;
                     tmp[0] = IV[0] ^ HashLengthBlake;
                     tmp[1] = IV[1];
                     tmp[2] = IV[2];
                     tmp[3] = IV[3];
                     tmp[4] = IV[4];
                     tmp[5] = IV[5];
                     tmp[6] = IV[6] ^ ZchashPow;
                     tmp[7] = IV[7] ^ ZchashParam;
    
                     uint64_t            v[16];
                     memcpy(&v[0], &tmp[0], sizeof(IV));
                     memcpy(&v[8], &IV[0], sizeof(IV));
                     v[12] ^= MidCount;
                     for (auto i = 0u; i < Rounds; ++i) {
                            uint8_t const* const s = Sigma[i];
                            mix(v + 0, v + 4, v + 8,  v + 12, m[s[0]],  m[s[1]]);
                            mix(v + 1, v + 5, v + 9,  v + 13, m[s[2]],  m[s[3]]);
                            mix(v + 2, v + 6, v + 10, v + 14, m[s[4]],  m[s[5]]);
                            mix(v + 3, v + 7, v + 11, v + 15, m[s[6]],  m[s[7]]);
                            mix(v + 0, v + 5, v + 10, v + 15, m[s[8]],  m[s[9]]);
                            mix(v + 1, v + 6, v + 11, v + 12, m[s[10]], m[s[11]]);
                            mix(v + 2, v + 7, v + 8,  v + 13, m[s[12]], m[s[13]]);
                            mix(v + 3, v + 4, v + 9,  v + 14, m[s[14]], m[s[15]]);
                     }
                     for (auto i = 0u; i < 8u; i++) {
                            tmp[i] ^= v[i] ^ v[i + 8u];
                     }
                     return tmp;
              }

              void finalise(Hash& hash, Block const& m) {
                     uint64_t            v[16];
                     memcpy(&v[0], &hash[0], sizeof(IV));
                     memcpy(&v[8], &IV[0], sizeof(IV));
                     v[12] ^= FullCount;
                     v[14] ^= 0xffffffffffffffffull;
                     for (auto i = 0u; i < Rounds; ++i) {
                            uint8_t const* const s = Sigma[i];
                            mix(v + 0, v + 4, v + 8,  v + 12, m[s[0]],  m[s[1]]);
                            mix(v + 1, v + 5, v + 9,  v + 13, m[s[2]],  m[s[3]]);
                            mix(v + 2, v + 6, v + 10, v + 14, m[s[4]],  m[s[5]]);
                            mix(v + 3, v + 7, v + 11, v + 15, m[s[6]],  m[s[7]]);
                            mix(v + 0, v + 5, v + 10, v + 15, m[s[8]],  m[s[9]]);
                            mix(v + 1, v + 6, v + 11, v + 12, m[s[10]], m[s[11]]);
                            mix(v + 2, v + 7, v + 8,  v + 13, m[s[12]], m[s[13]]);
                            mix(v + 3, v + 4, v + 9,  v + 14, m[s[14]], m[s[15]]);
                     }
                     for (auto i = 0u; i < 8u; i++) {
                            hash[i] ^= v[i] ^ v[i + 8];
                     }
              }
       };
}