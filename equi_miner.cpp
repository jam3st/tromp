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
        int nonce = 0;
        bool showsol = true;
        const char *header = "";
        const char *hex = "";
        equi eq;

        printf("Looking for wagner-tree on (\"%s\",%d", hex ? "0x..." : header, nonce);
        printf(") with %d %d-bit digits\n", NDIGITS, DIGITBITS);
        printf("Using %dMB of memory\n", 1 + eq.hta.alloced / 0x100000);
        u32 sumnsols = 0;
        char headernonce[HEADERNONCELEN];
        u32 hdrlen = strlen(header);
        if (*hex) {
                assert(strlen(hex) == 2 * HEADERNONCELEN);
                for (u32 i = 0; i < HEADERNONCELEN; i++)
                        headernonce[i] = hextobyte(&hex[2 * i]);
        } else {
                memcpy(headernonce, header, hdrlen);
                memset(headernonce + hdrlen, 0, sizeof(headernonce) - hdrlen);
        }
        ((u32 *)headernonce)[32] = htole32(nonce);
        eq.setheadernonce(headernonce, sizeof(headernonce));
        eq.solve();
        u32 nsols, maxsols = min(MAXSOLS, eq.nsols);
        for (nsols = 0; nsols < maxsols; nsols++) {
                if (showsol) {
                        printf("\nSolution");
                        for (u32 i = 0; i < PROOFSIZE; i++)
                                printf(" %jx",
                                        (uintmax_t)eq.sols[nsols][i]);
                }
                int pow_rc = verify(eq.sols[nsols], headernonce, sizeof(headernonce));
                if (pow_rc == POW_OK)
                        printf("\n is Verified\n");
                else
                        printf("FAILED due to %s\n", errstr[pow_rc]);

                printf("\n%d solutions\n", nsols);
                sumnsols += nsols;
        }
        printf("%d total solutions\n", sumnsols);
        return 0;
}
