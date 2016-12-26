#include "equihash.h"

#include "equi_miner.h"
#include "ctype.h"

static int hextobyte(const char *x) {
        u32 b = 0;
        for (int i = 0; i < 2; i++) {
                uchar c = tolower(x[i]);
                assert(isxdigit(c));
                b = (b << 4) | (c - (c >= '0' && c <= '9' ? '0' : ('a' - 10)));
        }
        return b;
}

int main(int argc, char **argv) {
        Gx::Blake2b::Header hdr;
        for(auto i = 0u; i < Gx::Blake2b::HeaderLen64Bits; ++i) {
                hdr[i] = 0ull;
        }
        // for(auto i = 0u; i < 20u; ++i) {
        //         hdr[i] = reinterpret_cast<Gx::uint64_t const* const>(&headernonce[1])[i];
        // }
        auto sols = Gx::Equihash::solve(hdr);        
        Gx::Equihash::verifySolution(sols.sol0, hdr);
        return 0;
}
