#pragma once

#include "blake2b.h"
#include <iostream>
#include <iomanip>

namespace Gx {
       namespace Equihash {
              static constexpr auto MaxNumSolutions = 8u; // mean ~ 2, stdev ~ 2 so 2 + 2 x 3
              static constexpr auto NumHashPerSolution64Bits = 256u;
              static constexpr auto NumberOfIndices64Bits = 32768u;
              static constexpr auto HeaderNonceLen8Bits = 140u;
              static constexpr auto HashStringLen64Bits = 5u;
              using HashString = Blake2b::Chunk<HashStringLen64Bits>;
              using Solution = Blake2b::Chunk<NumHashPerSolution64Bits>;
              using Solutions = struct {
                     Solution sol0;
                     Solution sol1;
                     Solution sol2;
                     Solution sol3;
                     uint64_t numSol = 0ull;
              };
              using Duplicates = Blake2b::Chunk<NumberOfIndices64Bits>;

              inline Blake2b::Block getMs(Blake2b::Header const& header) {
                     Blake2b::Block tmp;
                     for(auto i = 0u; i < Blake2b::BlakeBlockSize64Bits; ++i) {
                            tmp[i] = header[i];
                     }
                     return tmp;
              }

              inline Blake2b::Block getFinal(Blake2b::Header const& header) {
                     Blake2b::Block tmp;
                     static constexpr auto MaxHeaderLen64Bits = (HeaderNonceLen8Bits + 7u) / 8u - Blake2b::BlakeBlockSize64Bits;
                     for(auto i = Blake2b::BlakeBlockSize64Bits; i < Blake2b::BlakeBlockSize64Bits + MaxHeaderLen64Bits; ++i) {
                            tmp[i - Blake2b::BlakeBlockSize64Bits] = header[i];
                     }
                     for(auto i = MaxHeaderLen64Bits; i < Blake2b::BlakeBlockSize64Bits; ++i) {
                            tmp[i] = 0ull;
                     }
                     return tmp;
              }

              
              inline bool verifySolution(Solution const& solution, Blake2b::Header const& header) {
                     static HashString outs[512u];
                     static HashString ins[512u]; 
                     static Duplicates dups;

                     auto msBlk = getMs(header);
                     auto const ms = Gx::Blake2b::initMidstate(msBlk);
                     auto finalBlk = getFinal(header);
                     
                     for(auto i = 0u; i < NumberOfIndices64Bits; ++i) {
                            dups[i] = 0ull;
                     }
                     auto prev = solution.get32(0u);;
                     for(auto i = 0u; i < 512u; ++i) {
                            auto tmp = finalBlk;
                            auto const index = solution.get32(i);
                            auto const indexDupPos = index >> 6u;
                            auto const indexBitPos = 1ull <<  (index & ( (1 << 6u) - 1u));
                            if((dups[indexDupPos] & indexBitPos) != 0u) {
                                   return false;
                            } 
                            if((i & 0x1) == 0x1u) {
                                   if(prev >=  solution.get32(i)) {
                                          return false;
                                   }
                            } else {
                                   prev = solution.get32(i);
                            }
                            dups[indexDupPos] ^= indexBitPos;
                            tmp.set32(3u, htole32(index >> 1u) );
                            auto copyMs = ms;
                            finalise(copyMs, tmp);
                            auto const startPos = ((solution.get32(i) & 0x1u) == 0x1u) ? 25u : 0u;
                            for(auto j = 0u; j < 5u; ++j) {
                                   const auto currPos = startPos + 5u * j;
                                   outs[i].set32(j << 1u,         copyMs.get8(currPos)              | (copyMs.get8(currPos + 1u) << 8u) | (copyMs.get8(currPos + 2u) >> 4u) << 16u);
                                   outs[i].set32((j << 1u) + 1u,  (copyMs.get8(currPos + 2u) & 0xf) | (copyMs.get8(currPos + 3u) << 4u) | (copyMs.get8(currPos + 4u) << 12u));
                            }
                     }
                     
                     for(auto i = 9u; i != 0u ; --i) {
                            HashString* const src = ((i & 0x1u) == 0x1u) ? outs : ins;
                            HashString* const dst = ((i & 0x1u) == 0x1u) ? ins : outs;
                            for(auto j = 0u; j < (1u << i); j += 2) {
                                   for(auto k = 0u; k < 10u; ++k) {
                                          dst[j >> 1u].set32(k, src[j].get32(k) ^ src[j + 1u].get32(k));
                                   }
                            }
                     }
                     HashString* const res = ins;
                     auto finalValue = res[0u].get32(0u);
                     for(auto i = 1u; i < 10u ; ++i) {
                            finalValue ^= res[0].get32(i);
                     }
                     return finalValue == 0x0u;
              }

              inline Solutions solve(Blake2b::Header const& header) {
                     auto msBlk = getMs(header);
                     auto const ms = Gx::Blake2b::initMidstate(msBlk);
                     auto finalBlk = getFinal(header);
                     

                     initializeLists(header);
                     
                     Solutions s;
                     return s;
              }
       }
}