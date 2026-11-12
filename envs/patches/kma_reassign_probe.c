/* Diagnostic harness for the unmodified upstream matching API. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "assembly.h"
#include "compdna.h"
#include "reassign.h"

static CompDNA *encode(const char *text) {
    const char *alphabet = "ACGT";
    size_t n = strlen(text);
    unsigned char *bases = malloc(n + 1);
    CompDNA *value = setComp(n + 64);
    if(!bases || !value) exit(2);
    for(size_t i = 0; i < n; ++i) {
        const char *base = strchr(alphabet, text[i]);
        bases[i] = base ? (unsigned char)(base - alphabet) : 4;
    }
    compDNA(value, bases, n);
    free(bases);
    return value;
}

int main(int argc, char **argv) {
    if(argc != 3) return 2;
    CompDNA *consensus = encode(argv[1]);
    CompDNA *candidate = encode(argv[2]);
    printf("%d\n", reassign_matchseqs(consensus, candidate));
    destroyComp(consensus);
    destroyComp(candidate);
    return 0;
}
