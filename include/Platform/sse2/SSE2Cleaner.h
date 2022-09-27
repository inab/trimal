#ifndef TRIMAL_SSE2CLEANER_H
#define TRIMAL_SSE2CLEANER_H

#ifdef HAVE_SSE2

#include "Cleaner.h"

class SSE2Cleaner: public Cleaner {
public:
    SSE2Cleaner(Alignment* parent);
    SSE2Cleaner(Alignment* parent, Cleaner* mold);
    void calculateSeqIdentity() override;
    bool calculateSpuriousVector(float overlap, float *spuriousVector) override;
};

#endif
#endif
