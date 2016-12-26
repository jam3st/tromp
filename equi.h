#include "blake/blake2.h"
#include <endian.h>
#include <stdint.h> // for types uint32_t,uint64_t
#include <stdlib.h> // for function qsort
#include <string.h> // for functions memset
#include <iostream>
#include <iomanip>

typedef uint32_t u32;
typedef unsigned char uchar;

// algorithm parameters, prefixed with W (for Wagner) to reduce include file
// conflicts

static const u32 WN = 25u;

static const u32 WK = 9u;

static const u32 HEADERNONCELEN = 140u;

static const u32 NDIGITS = (WK + 1); // 10
static const u32 DIGITBITS = (WN * 8u / (NDIGITS)); // 20

const static u32 PROOFSIZE = 1 << WK;
const static u32 BASE = 1 << DIGITBITS; // 1024
const static u32 NHASHES = 2 * BASE; // 2048
const static u32 HASHESPERBLAKE = 512 / (WN * 8u);
const static u32 HASHOUT = HASHESPERBLAKE * WN;

typedef u32 proof[PROOFSIZE];

void setheader(blake2b_state *ctx, char const headernonce[]) {
        uint32_t le_N = htole32(WN * 8u);
        uint32_t le_K = htole32(WK);
        uchar personal[] = "ZcashPoW01230123";
        memcpy(personal + 8, &le_N, 4);
        memcpy(personal + 12, &le_K, 4);
        blake2b_param P[1];
        P->digest_length = HASHOUT; // 50
        P->key_length = 0;
        P->fanout = 1;
        P->depth = 1;
        P->leaf_length = 0;
        P->node_offset = 0;
        P->node_depth = 0;
        P->inner_length = 0;
        memset(P->reserved, 0, sizeof(P->reserved));
        memset(P->salt, 0, sizeof(P->salt));
        memcpy(P->personal, (const uint8_t *)personal, 16);
        blake2b_init_param(ctx, P);
        blake2b_update(ctx, (const uchar *)headernonce, HEADERNONCELEN);
}

enum verify_code {      
        POW_OK,
        POW_HEADER_LENGTH,
        POW_DUPLICATE,
        POW_OUT_OF_ORDER,
        POW_NONZERO_XOR
};
const char *errstr[] = {"OK", "wrong header length", "duplicate index",
                        "indices out of order", "nonzero xor"};

void genhash(blake2b_state const* const ctx, u32 idx, uchar hash[]) {
        blake2b_state state = *ctx;
        u32 leb = htole32(idx / HASHESPERBLAKE);
        blake2b_update(&state, (uchar *)&leb, sizeof(u32));
        uchar blakehash[HASHOUT];
        blake2b_final(&state, blakehash, HASHOUT);
        memcpy(hash, blakehash + (idx % HASHESPERBLAKE) * WN, WN);
}

int verifyrec(blake2b_state const* const ctx, u32 const indices[], uchar hash[WN], u32 r) {
        if (r == 0) {
                std::cerr << "GH " <<  std::hex << std::setfill('0') << std::setw(8) << indices[0] << ":";
                genhash(ctx, *indices, hash);
                for(int i = 0; i < 25; ++i) {
                        std::cerr << " "  << std::hex << std::setfill('0') << std::setw(2) << (u32)hash[i]      ;
                }
                std::cerr << std::endl;
                return POW_OK;
        }
        u32 const* const indices1 = &indices[0] + (1 << (r - 1u));
        if (*indices >= *indices1) {
                return POW_OUT_OF_ORDER;
        }
        uchar hash0[WN], hash1[WN];
        int vrf0 = verifyrec(ctx, indices, hash0, r - 1);
        if (vrf0 != POW_OK) {
                return vrf0;
        }
        int vrf1 = verifyrec(ctx, indices1, hash1, r - 1);
        if (vrf1 != POW_OK) {
                return vrf1;
        }
        for (u32 i = 0u; i < WN; i++)
                hash[i] = hash0[i] ^ hash1[i];
        u32 b = r < WK ? r * DIGITBITS : WN * 8u;
        std::cerr << "B: " << std::dec   << b/8 << std::endl;   
        u32 i = 0;
        for (; i < b / 8; i++) {
                if (hash[i]) {
                        return POW_NONZERO_XOR;
                }
        }
        if ((b % 8) && hash[i] >> (8 - (b % 8))) {
                return POW_NONZERO_XOR;
        }
        return POW_OK;
}

int compu32(const void *pa, const void *pb) {
        u32 a = *(u32 *)pa, b = *(u32 *)pb;
        return a < b ? -1 : a == b ? 0 : +1;
}

bool duped(proof const prf) {
        proof sortprf;
        memcpy(sortprf, prf, sizeof(proof));
        qsort(sortprf, PROOFSIZE, sizeof(u32), &compu32);
        for (u32 i = 1; i < PROOFSIZE; i++)
                if (sortprf[i] <= sortprf[i - 1])
                        return true;
        return false;
}

// verify Wagner conditions
int verify(u32 const indices[PROOFSIZE], char const headernonce[], const u32 headerlen) {
        if (headerlen != HEADERNONCELEN)
                return POW_HEADER_LENGTH;
        if (duped(indices))
                return POW_DUPLICATE;
        blake2b_state ctx;
        setheader(&ctx, headernonce);
        uchar hash[WN];
        return verifyrec(&ctx, indices, hash, WK);
}
