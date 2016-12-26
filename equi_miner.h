// Equihash solver
// Copyright (c) 2016 John Tromp

// Equihash presents the following problem
//
// Fix N, K, such that N is a multiple of K+1
// Let integer n = N/(K+1), and view N-bit words
// as having K+1 "digits" of n bits each
// Fix M = 2^{n+1} N-bit hashes H_0, ... , H_{M-1}
// as outputs of a hash function applied to an (n+1)-bit index
//
// Problem: find a binary tree on 2^K distinct indices,
// for which the exclusive-or of leaf hashes is all 0s
// Additionally, it should satisfy the Wagner conditions:
// 1) for each height i subtree, the exclusive-or
// of its 2^i leaf hashes starts with i*n 0 bits,
// 2) the leftmost leaf of any left subtree is less
// than the leftmost leaf of the corresponding right subtree
//
// The algorithm below solves this by storing trees
// as a directed acyclic graph of K layers
// The n digit bits are split into
// n-RESTBITS bucket bits and RESTBITS leftover bits
// Each layer i, consisting of height i subtrees
// whose xor starts with i*n 0s, is partitioned into
// 2^{n-RESTBITS} buckets according to the next n-RESTBITS
// in the xor
// Within each bucket, trees whose xor match in the
// next RESTBITS bits are combined to produce trees
// in the next layer
// To eliminate trees with duplicated indices,
// we simply test if the last 32 bits of the xor are 0,
// and if so, assume that this is due to index duplication
// In practice this works very well to avoid bucket overflow
// and produces negligible false positives

#include <iostream>
#include "equi.h"
#include <assert.h>
#include <stdio.h>

#include "blake2-avx2/blake2bip.h"
typedef blake2b_state blake_state;

// u32 already defined in equi.h
typedef uint16_t u16;
typedef uint64_t u64;

// required for avoiding multio-threading race conflicts
typedef u32 au32;

static const u32 RESTBITS = 8u;

// 2_log of number of buckets
static const u32 BUCKBITS = (DIGITBITS - RESTBITS); // 12


static const u32 NBUCKETS = 1 << BUCKBITS;      // number of buckets
static const u32 BUCKMASK = NBUCKETS - 1;       // corresponding bucket mask
static const u32 SLOTBITS = RESTBITS + 1 + 1;   // 2_log of number of slots per bucket
static const u32 SLOTRANGE = 1 << SLOTBITS;     // default bucket capacity
static const u32 SLOTMASK = SLOTRANGE - 1;      // corresponding SLOTBITS mask
static const u32 SLOTMSB = 1 << (SLOTBITS - 1); // most significat bit in SLOTMASK
static const u32 NSLOTS = SLOTRANGE;  // number of slots per bucket
static const u32 NRESTS = 1 << RESTBITS;        // number of possible values of RESTBITS bits
static const u32 MAXSOLS = 8;                   // more than 8 solutions are rare

// tree node identifying its children as two different slots in
// a bucket on previous layer with matching rest bits (x-tra hash)
struct tree {
        // formerly i had these bitfields
        // unsigned bucketid : BUCKBITS;
        // unsigned slotid0  : SLOTBITS;
        // unsigned slotid1 : SLOTBITS;
        // but these were poorly optimized by the compiler
        // so now we do things "manually"
        u32 bid_s0_s1;

        static_assert(BUCKBITS + 2 * SLOTBITS <= 32, "cantor throws a fit");

        // constructor for height 0 trees stores index instead
        tree(const u32 idx) { bid_s0_s1 = idx; }
        static u32 cantor(u32 s0, u32 s1) { return s1 * (s1 + 1) / 2 + s0; }
        tree(const u32 bid, const u32 s0, const u32 s1) {
                bid_s0_s1 = (((bid << SLOTBITS) | s0) << SLOTBITS) | s1;
        }
        // retrieve hash index from tree(const u32 idx) constructor
        u32 getindex() const { return bid_s0_s1; }
        // retrieve bucket index
        u32 bucketid() const {
                return bid_s0_s1 >> (2 * SLOTBITS);
        }
// retrieve first slot index
        u32 slotid0() const { return (bid_s0_s1 >> SLOTBITS) & SLOTMASK; }
        // retrieve second slot index
        u32 slotid1() const {
                return bid_s0_s1 & SLOTMASK;
        }
        bool prob_disjoint(const tree other) const {
                tree xort(bid_s0_s1 ^ other.bid_s0_s1);
                return xort.bucketid() || (xort.slotid0() && xort.slotid1());
// next two tests catch much fewer cases and are therefore skipped
// && slotid0() != other.slotid1() && slotid1() != other.slotid0()
        }
};

// each bucket slot occupies a variable number of hash/tree units,
// all but the last of which hold the xor over all leaf hashes,
// or what's left of it after stripping the initial i*n 0s
// the last unit holds the tree node itself
// the hash is sometimes accessed 32 bits at a time (word)
// and sometimes 8 bits at a time (bytes)
union htunit {
        tree tag;
        u32 word;
        uchar bytes[sizeof(u32)];
};

constexpr u32 WORDS(u32 const bits) {
         return (bits + 31u) / 32u;
}

u32 const HASHWORDS0 = WORDS(WN * 8u - DIGITBITS + RESTBITS);
u32 const HASHWORDS1 = WORDS(WN * 8u - 2 * DIGITBITS + RESTBITS);

// A slot is up to HASHWORDS0 hash units followed by a tag
typedef htunit slot0[HASHWORDS0 + 1];
typedef htunit slot1[HASHWORDS1 + 1];
// a bucket is NSLOTS treenodes
typedef slot0 bucket0[NSLOTS];
typedef slot1 bucket1[NSLOTS];
// the N-bit hash consists of K+1 n-bit "digits"
// each of which corresponds to a layer of NBUCKETS buckets
typedef bucket0 digit0[NBUCKETS];
typedef bucket1 digit1[NBUCKETS];
typedef au32 bsizes[NBUCKETS];

// The algorithm proceeds in K+1 rounds, one for each digit
// All data is stored in two heaps,
// heap0 of type digit0, and heap1 of type digit1
// The following table shows the layout of these heaps
// in each round, which is an optimized version
// of xenoncat's fixed memory layout, avoiding any waste
// Each line shows only a single slot, which is actually
// replicated NSLOTS * NBUCKETS times
//
//             heap0         heap1
// round  hashes   tree   hashes tree
// 0      A A A A A A 0   . . . . . .
// 1      A A A A A A 0   B B B B B 1
// 2      C C C C C 2 0   B B B B B 1
// 3      C C C C C 2 0   D D D D 3 1
// 4      E E E E 4 2 0   D D D D 3 1
// 5      E E E E 4 2 0   F F F 5 3 1
// 6      G G 6 . 4 2 0   F F F 5 3 1
// 7      G G 6 . 4 2 0   H H 7 5 3 1
// 8      I 8 6 . 4 2 0   H H 7 5 3 1
//
// Round 0 generates hashes and stores them in the buckets
// of heap0 according to the initial n-RESTBITS bits
// These hashes are denoted A above and followed by the
// tree tag denoted 0
// In round 1 we combine each pair of slots in the same bucket
// with matching RESTBITS of digit 0 and store the resulting
// 1-tree in heap1 with its xor hash denoted B
// Upon finishing round 1, the A space is no longer needed,
// and is re-used in round 2 to store both the shorter C hashes,
// and their tree tags denoted 2
// Continuing in this manner, each round reads buckets from one
// heap, and writes buckets in the other heap.
// In the final round K, all pairs leading to 0 xors are identified
// and their leafs recovered through the DAG of tree nodes

// convenience function
u32 min(const u32 a, const u32 b) { return a < b ? a : b; }

// size (in bytes) of hash in round 0 <= r < WK
u32 hashsize(const u32 r) {
        const u32 hashbits = (WN * 8) - (r + 1) * DIGITBITS + RESTBITS;
        return (hashbits + 7) / 8;
}

// convert bytes into words,rounding up
u32 hashwords(u32 bytes) { return (bytes + 3) / 4; }

// manages hash and tree data
struct htalloc {
        bucket0 *heap0;
        bucket1 *heap1;
        u32 alloced;
        htalloc() { alloced = 0; }
        void alloctrees() {
                static_assert(DIGITBITS >= 16, "needed to ensure hashes shorten by 1 unit every 2 digits");
                heap0 = (bucket0 *)alloc(NBUCKETS, sizeof(bucket0));
                heap1 = (bucket1 *)alloc(NBUCKETS, sizeof(bucket1));
        }
        void dealloctrees() {
                free(heap0);
                free(heap1);
        }
        void *alloc(const u32 n, const u32 sz) {
                void *mem = calloc(n, sz);
                assert(mem);
                alloced += n * sz;
                return mem;
        }
};

// main solver object, shared between all threads
struct equi {
        blake_state blake_ctx; // holds blake2b midstate after call to setheadernounce
        htalloc hta;           // holds allocated heaps
        bsizes *nslots;        // counts number of slots used in buckets
        proof *sols;           // store found solutions here (only first MAXSOLS)
        au32 nsols;            // number of solutions found
        u32 bfull;               // count number of times bucket can't fit new item
        u32 hfull;               // count number of xor-ed hash with last 32 bits zero
        equi() {
                static_assert(sizeof(htunit) == 4, "");
                static_assert(WK & 1, "K assumed odd in candidate() calling indices1()");
                hta.alloctrees();
                nslots = (bsizes *)hta.alloc(2 * NBUCKETS, sizeof(au32));
                sols = (proof *)hta.alloc(MAXSOLS, sizeof(proof));
        }
        ~equi() {
                hta.dealloctrees();
                free(nslots);
                free(sols);
        }
        // prepare blake2b midstate for new run and initialize counters
        void setheadernonce(const char *headernonce, const u32 len) {
                setheader(&blake_ctx, headernonce);
                nsols = bfull = hfull = 0;
        }
        // get heap0 bucket size in threadsafe manner
        u32 getslot0(const u32 bucketi) {
                return nslots[0][bucketi]++;
        }
        // get heap1 bucket size in threadsafe manner
        u32 getslot1(const u32 bucketi) {
                return nslots[1][bucketi]++;
        }
        // get old heap0 bucket size and clear it for next round
        u32 getnslots0(const u32 bid) {
                au32 &nslot = nslots[0][bid];
                const u32 n = min(nslot, NSLOTS);
                nslot = 0;
                return n;
        }
        // get old heap1 bucket size and clear it for next round
        u32 getnslots1(const u32 bid) {
                au32 &nslot = nslots[1][bid];
                const u32 n = min(nslot, NSLOTS);
                nslot = 0;
                return n;
        }
        // recognize most (but not all) remaining dupes while Wagner-ordering
        // the indices
        bool orderindices(u32 *indices, u32 size) {
                if (indices[0] > indices[size]) {
                        for (u32 i = 0; i < size; i++) {
                                const u32 tmp = indices[i];
                                indices[i] = indices[size + i];
                                indices[size + i] = tmp;
                        }
                }
                return false;
        }
        // listindices combines index tree reconstruction with probably dupe
        // test
        bool listindices0(u32 r, const tree t, u32 *indices) {
                if (r == 0) {
                        *indices = t.getindex();
                        return false;
                }
                const slot1 *buck = hta.heap1[t.bucketid()];
                const u32 size = 1 << --r;
                u32 tagi = hashwords(hashsize(r));
                u32 s1 = t.slotid1(), s0 = t.slotid0();
                tree t0 = buck[s0][tagi].tag, t1 = buck[s1][tagi].tag;
                return !t0.prob_disjoint(t1) || listindices1(r, t0, indices) || listindices1(r, t1, indices + size) || orderindices(indices, size) ||
                       indices[0] == indices[size];
        }
        // need separate instance for accessing (differently typed) heap1
        bool listindices1(u32 r, const tree t, u32 *indices) {
                const slot0 *buck = hta.heap0[t.bucketid()];
                const u32 size = 1 << --r;
                u32 tagi = hashwords(hashsize(r));
                u32 s1 = t.slotid1(), s0 = t.slotid0();
                tree t0 = buck[s0][tagi].tag, t1 = buck[s1][tagi].tag;
                return listindices0(r, t0, indices) || listindices0(r, t1, indices + size) || orderindices(indices, size) || indices[0] == indices[size];
        }
        // check a candidate that resulted in 0 xor
        // add as solution, with proper subtree ordering, if it has unique
        // indices
        void candidate(const tree t) {
                proof prf;
                // listindices combines index tree reconstruction with probably
                // dupe test
                if (listindices1(WK, t, prf) || duped(prf))
                        return; // assume WK odd
                                // and now we have ourselves a genuine solution
                u32 soli = nsols++;
                // copy solution into final place
                if (soli < MAXSOLS)
                        memcpy(sols[soli], prf, sizeof(proof));
        }
        // thread-local object that precomputes various slot metrics for each
        // round
        // facilitating access to various bits in the variable size slots
        struct htlayout {
                htalloc hta;
                u32 prevhtunits;
                u32 nexthtunits;
                u32 dunits;
                u32 prevbo;

                htlayout(equi *eq, u32 r) : hta(eq->hta), prevhtunits(0), dunits(0) {
                        u32 nexthashbytes = hashsize(r);        // number of bytes occupied by round r hash
                        nexthtunits = hashwords(nexthashbytes); // number of 32bit words
                                                                // taken up by those bytes
                        prevbo = 0;                             // byte offset for accessing hash form
                                                                // previous round
                        if (r) {                                // similar measure for previous round
                                u32 prevhashbytes = hashsize(r - 1);
                                prevhtunits = hashwords(prevhashbytes);
                                prevbo = prevhtunits * sizeof(htunit) - prevhashbytes; // 0-3
                                dunits = prevhtunits - nexthtunits;                    // number of words by
                                                                                       // which hash shrinks
                        }
                }
                // extract remaining bits in digit slots in same bucket still
                // need to collide on
                u32 getxhash0(const htunit *slot) const {
                        return (slot->bytes[prevbo] & 0xf) << 4 | slot->bytes[prevbo + 1] >> 4;
                }
                // similar but accounting for possible change in hashsize modulo
                // 4 bits
                u32 getxhash1(const htunit *slot) const {
                        return slot->bytes[prevbo];
                }
                // test whether two hashes match in last 32 bits
                bool equal(const htunit *hash0, const htunit *hash1) const { return hash0[prevhtunits - 1].word == hash1[prevhtunits - 1].word; }
        };

        // this thread-local object performs in-bucket collisions
        // by linking together slots that have identical rest bits
        // (which is in essense a 2nd stage bucket sort)
        struct collisiondata {
// This maintains NRESTS = 2^RESTBITS lists whose starting slot
// are in xhashslots[] and where subsequent (next-lower-numbered)
// slots in each list are found through nextxhashslot[]
// since 0 is already a valid slot number, use ~0 as nil value
                typedef u16 xslot;
                static const xslot xnil = ~0;
                xslot xhashslots[NRESTS];
                xslot nextxhashslot[NSLOTS];
                xslot nextslot;
                u32 s0;

                void clear() {
                        memset(xhashslots, xnil, NRESTS * sizeof(xslot));
                        memset(nextxhashslot, xnil, NSLOTS * sizeof(xslot));
                }
                void addslot(u32 s1, u32 xh) {
                        nextslot = xhashslots[xh];
                        nextxhashslot[s1] = nextslot;
                        xhashslots[xh] = s1;
                }
                bool nextcollision() const {
                        return nextslot != xnil;
                }
                u32 slot() {
                        nextslot = nextxhashslot[s0 = nextslot];
                        return s0;
                }
        };


        // number of hashes extracted from NBLAKES blake2b outputs
        static const u32 HASHESPERBLOCK = HASHESPERBLAKE;
        // number of blocks of parallel blake2b calls
        static const u32 NBLOCKS = (NHASHES + HASHESPERBLOCK - 1) / HASHESPERBLOCK;

        void digit0() {
                htlayout htl(this, 0);
                const u32 hashbytes = hashsize(0);
                uchar hashes[64];
                blake_state state0 = blake_ctx; // local copy on stack can be copied faster
                std::cerr << "Num blocks " << NBLOCKS << std::endl;
                for (u32 block = 0; block < NBLOCKS; ++block) {
                        blake_state state = state0; // make another copy since
                                                    // blake2b_final modifies it
                        u32 leb = htole32(block);
                        blake2b_update(&state, (uchar *)&leb, sizeof(u32));
                        blake2b_final(&state, hashes, HASHOUT);
                        for (u32 j = 0; j < HASHESPERBLAKE; j++) {
                                const uchar *ph = hashes + j * WN;
                                const u32 bucketid = ((u32)ph[0] << 4) | ph[1] >> 4;
                                // grab next available slot in that
                                // bucket
                                const u32 slot = getslot0(bucketid);    
                                if (slot >= NSLOTS) {
                                        bfull++; // this actually never
                                                        // seems to happen in
                                                        // round 0 due to
                                                        // uniformity
                                        continue;
                                }
                                // location for slot's tag
                                htunit *s = hta.heap0[bucketid][slot] + htl.nexthtunits;
                                // hash should end right before tag
                                memcpy(s->bytes - hashbytes, ph + WN - hashbytes, hashbytes);
                                // round 0 tags store hash-generating
                                // index
                                s->tag = tree(block * HASHESPERBLAKE + j);
                        }
                }
        }

        void digitodd(u32 const r) {
                htlayout htl(this, r);
                collisiondata cd;
                // threads process buckets in round-robin fashion
                for (u32 bucketid = 0; bucketid < NBUCKETS; ++bucketid) {
                        cd.clear();                            // could have made this the constructor, and
                                                               // declare here
                        slot0 *buck = htl.hta.heap0[bucketid]; // point to first
                                                               // slot of this
                                                               // bucket
                        u32 bsize = getnslots0(bucketid);      // grab and reset bucket size
                        for (u32 s1 = 0; s1 < bsize; s1++) {   // loop over slots
                                const htunit *slot1 = buck[s1];
                                cd.addslot(s1, htl.getxhash0(slot1)); // identify list
                                                                      // of previous
                                                                      // colliding
                                                                      // slots
                                for (; cd.nextcollision();) {
                                        const u32 s0 = cd.slot();
                                        const htunit *slot0 = buck[s0];
                                        if (htl.equal(slot0, slot1)) { // expect
                                                                       // difference in
                                                                       // last 32 bits
                                                                       // unless duped
                                                hfull++;               // record discarding
                                                continue;
                                        }
                                        u32 xorbucketid; // determine bucket for
                                                         // s0 xor s1
                                        const uchar *bytes0 = slot0->bytes, *bytes1 = slot1->bytes;
                                        xorbucketid = (((u32)(bytes0[htl.prevbo + 1] ^ bytes1[htl.prevbo + 1]) & 0xf) << 8) |
                                                      (bytes0[htl.prevbo + 2] ^ bytes1[htl.prevbo + 2]);
                                        // grab next available slot in that
                                        // bucket
                                        const u32 xorslot = getslot1(xorbucketid);
                                        if (xorslot >= NSLOTS) {
                                                bfull++; // SAVEMEM determines
                                                         // how often this
                                                         // happens
                                                continue;
                                        }
                                        // start of slot for s0 ^ s1
                                        htunit *xs = htl.hta.heap1[xorbucketid][xorslot];
                                        // store xor of hashes possibly minus
                                        // initial 0 word due to collision
                                        for (u32 i = htl.dunits; i < htl.prevhtunits; i++)
                                                xs++->word = slot0[i].word ^ slot1[i].word;
                                        // store tree node right after hash
                                        xs->tag = tree(bucketid, s0, s1);
                                }
                        }
                }
        }

        void digiteven(u32 const r) {
                htlayout htl(this, r);
                collisiondata cd;
                for (u32 bucketid = 0; bucketid < NBUCKETS; ++bucketid) {
                        cd.clear();
                        slot1 *buck = htl.hta.heap1[bucketid];
                        u32 bsize = getnslots1(bucketid);
                        for (u32 s1 = 0; s1 < bsize; s1++) {
                                const htunit *slot1 = buck[s1];
                                cd.addslot(s1, htl.getxhash1(slot1));
                                for (; cd.nextcollision();) {
                                        const u32 s0 = cd.slot();
                                        const htunit *slot0 = buck[s0];
                                        if (htl.equal(slot0, slot1)) {
                                                hfull++;
                                                continue;
                                        }
                                        u32 xorbucketid;
                                        const uchar *bytes0 = slot0->bytes, *bytes1 = slot1->bytes;
                                        xorbucketid = ((u32)(bytes0[htl.prevbo + 1] ^ bytes1[htl.prevbo + 1]) << 4) |
                                                      (bytes0[htl.prevbo + 2] ^ bytes1[htl.prevbo + 2]) >> 4;
                                        const u32 xorslot = getslot0(xorbucketid);
                                        if (xorslot >= NSLOTS) {
                                                bfull++;
                                                continue;
                                        }
                                        htunit *xs = htl.hta.heap0[xorbucketid][xorslot];
                                        for (u32 i = htl.dunits; i < htl.prevhtunits; i++)
                                                xs++->word = slot0[i].word ^ slot1[i].word;
                                        xs->tag = tree(bucketid, s0, s1);
                                }
                        }
                }
        }

        // final round looks simpler
        void digitK() {
                collisiondata cd;
                htlayout htl(this, WK);
                u32 nc = 0;
                for (u32 bucketid = 0; bucketid < NBUCKETS; ++bucketid) {
                        cd.clear();
                        slot0 *buck = htl.hta.heap0[bucketid]; // assume WK odd
                        u32 bsize = getnslots0(bucketid);      // assume WK odd
                        for (u32 s1 = 0; s1 < bsize; s1++) {
                                const htunit *slot1 = buck[s1];
                                cd.addslot(s1, htl.getxhash0(slot1)); // assume WK odd
                                for (; cd.nextcollision();) {
                                        const u32 s0 = cd.slot();
                                        const htunit *slot0 = buck[s0];
                                        // there is only 1 word of hash left
                                        if (htl.equal(slot0, slot1) && slot0[1].tag.prob_disjoint(slot1[1].tag)) {
                                                candidate(tree(bucketid, s0, s1)); // so a match gives a
                                                                                   // solution candidate
                                                nc++;
                                        }
                                }
                        }
                }
                // printf(" %d candidates ", nc);  // this gets uncommented a
                // lot for debugging
        }
        
        void solve() {
                printf("Digit 0");
                digit0();
                for (u32 r = 1; r < WK; ++r) {
                        r & 1u ? digitodd(r) : digiteven(r);
                }       
                digitK();
        }
};
