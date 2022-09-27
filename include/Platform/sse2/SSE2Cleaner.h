#ifndef TRIMAL_SSE2CLEANER_H
#define TRIMAL_SSE2CLEANER_H

#ifdef HAVE_SSE2

#include "Cleaner.h"

class SSECleaner: public Cleaner {
private:
    // SIMD index for which residues can be skipped in `skipResidues`
    unsigned char* skipResidues;
    unsigned char* skipResidues_unaligned;
    // temporary counters for `calculateSpuriousVector`
    uint8_t* hits_u8;
    uint8_t* hits_u8_unaligned;
    uint32_t* hits;
    uint32_t* hits_unaligned;
public:
    SSECleaner(Alignment* parent);
    ~SSECleaner();
    void calculateSeqIdentity() override;
    bool calculateSpuriousVector(float overlap, float *spuriousVector);
};

#endif
#endif
