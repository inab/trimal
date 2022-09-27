#ifdef HAVE_SSE2

#include <immintrin.h>
#include <climits>
#include <cstdint>

#include "Alignment/Alignment.h"
#include "Statistics/Manager.h"
#include "Statistics/Similarity.h"
#include "InternalBenchmarker.h"
#include "reportsystem.h"
#include "defines.h"
#include "utils.h"

#include "Platform/sse2/SSE2Cleaner.h"

#define NLANES_8   sizeof(__m128i) / sizeof(uint8_t)   // number of 8-bit lanes in __m128i
#define NLANES_32  sizeof(__m128i) / sizeof(uint32_t)  // number of 32-bit lanes in __m128i

#define ALLOC_MASK           (sizeof(__m128i) - 1)
#define ALLOC_SIZE(N, T)     ((N * sizeof(T) + ALLOC_MASK) & (~ALLOC_MASK))
#define ALIGNED_ALLOC(N, T)  (static_cast<T*>(aligned_alloc(sizeof(__m128i), ALLOC_SIZE(N, T))))

static inline uint32_t _mm_hsum_epi8(__m128i a) {
    __m128i vsum = _mm_sad_epu8(a, _mm_setzero_si128());
    return _mm_extract_epi16(vsum, 0) + _mm_extract_epi16(vsum, 4);
}

SSE2Cleaner::SSE2Cleaner(Alignment* parent): Cleaner(parent) {}

SSE2Cleaner::SSE2Cleaner(Alignment* parent, Cleaner* mold): Cleaner(parent, mold) {}

void SSE2Cleaner::calculateSeqIdentity() {
  // Create a timer that will report times upon its destruction
  //	which means the end of the current scope.
  StartTiming("void SSE2Cleaner::calculateSeqIdentity(void) ");

  // allocate aligned memory for faster SIMD loads
  uint8_t* skipResidues = ALIGNED_ALLOC(alig->originalNumberOfResidues, uint8_t);

  // create a bitmask for residues to skip
  for(int i = 0; i < alig->originalNumberOfResidues; i++) {
      skipResidues[i] = alig->saveResidues[i] == -1 ? 0xFF : 0;
  }

  // declare indices
  int i, j, k, l;

  // Depending on alignment type, indetermination symbol will be one or other
  char indet = (alig->getAlignmentType() & SequenceTypes::AA) ? 'X' : 'N';

  // prepare constant SIMD vectors
  const __m128i allindet = _mm_set1_epi8(indet);
  const __m128i allgap   = _mm_set1_epi8('-');
  const __m128i ONES     = _mm_set1_epi8(1);

  // Create identities matrix to store identities scores
  alig->identities = new float*[alig->originalNumberOfSequences];
  for(i = 0; i < alig->originalNumberOfSequences; i++) {
      if (alig->saveSequences[i] == -1) continue;
      alig->identities[i] = new float[alig->originalNumberOfSequences];
      alig->identities[i][i] = 0;
  }

  // For each seq, compute its identity score against the others in the MSA
  for (i = 0; i < alig->originalNumberOfSequences; i++) {
      if (alig->saveSequences[i] == -1) continue;

      // Compute identity scores for the current sequence against the rest
      for (j = i + 1; j < alig->originalNumberOfSequences; j++) {
          if (alig->saveSequences[j] == -1) continue;

          const char* datai = alig->sequences[i].data();
          const char* dataj = alig->sequences[j].data();

          __m128i seqi;
          __m128i seqj;
          __m128i skip;
          __m128i eq;
          __m128i gapsi;
          __m128i gapsj;
          __m128i mask;
          __m128i dst_acc = _mm_setzero_si128();
          __m128i hit_acc = _mm_setzero_si128();

          int hit = 0;
          int dst = 0;

          // run with unrolled loops of UINT8_MAX iterations first
          for (k = 0; ((int) (k + NLANES_8*UINT8_MAX)) < alig->originalNumberOfResidues;) {
              for (l = 0; l < UINT8_MAX; l++, k += NLANES_8) {
                  // load data for the sequences
                  seqi = _mm_loadu_si128( (const __m128i*) (&datai[k]));
                  seqj = _mm_loadu_si128( (const __m128i*) (&dataj[k]));
                  skip = _mm_load_si128(  (const __m128i*) (&skipResidues[k]));
                  eq = _mm_cmpeq_epi8(seqi, seqj);
                  // find which sequence characters are gap or indet
                  gapsi = _mm_or_si128(_mm_cmpeq_epi8(seqi, allgap), _mm_cmpeq_epi8(seqi, allindet));
                  gapsj = _mm_or_si128(_mm_cmpeq_epi8(seqj, allgap), _mm_cmpeq_epi8(seqj, allindet));
                  // find position where not both characters are gap
                  mask = _mm_andnot_si128(skip, _mm_andnot_si128(_mm_and_si128(gapsi, gapsj), ONES));
                  // update counters
                  dst_acc = _mm_add_epi8(dst_acc,                   mask);
                  hit_acc = _mm_add_epi8(hit_acc, _mm_and_si128(eq, mask));
              }
              // merge accumulators
              dst += _mm_hsum_epi8(dst_acc);
              hit += _mm_hsum_epi8(hit_acc);
              dst_acc = _mm_setzero_si128();
              hit_acc = _mm_setzero_si128();
          }

          // run remaining iterations in SIMD while possible
          for (; ((int) (k + NLANES_8)) < alig->originalNumberOfResidues; k += NLANES_8) {
              // load data for the sequences; load is unaligned because
              // string data is not guaranteed to be aligned.
              seqi = _mm_loadu_si128( (const __m128i*) (&datai[k]));
              seqj = _mm_loadu_si128( (const __m128i*) (&dataj[k]));
              skip = _mm_load_si128(  (const __m128i*) (&skipResidues[k]));
              eq = _mm_cmpeq_epi8(seqi, seqj);
              // find which sequence characters are either a gap or
              // an indeterminate character
              gapsi = _mm_or_si128(_mm_cmpeq_epi8(seqi, allgap), _mm_cmpeq_epi8(seqi, allindet));
              gapsj = _mm_or_si128(_mm_cmpeq_epi8(seqj, allgap), _mm_cmpeq_epi8(seqj, allindet));
              // find position where not both characters are gap
              mask = _mm_andnot_si128(skip, _mm_andnot_si128(_mm_and_si128(gapsi, gapsj), ONES));
              // update counters: update dst if any of the two sequence
              // characters is non-gap/indeterminate; update length if
              // any of the characters is non-gappy/indeterminate.
              dst_acc = _mm_add_epi8(dst_acc,                   mask);
              hit_acc = _mm_add_epi8(hit_acc, _mm_and_si128(eq, mask));
          }

          // update counters after last loop
          hit += _mm_hsum_epi8(hit_acc);
          dst += _mm_hsum_epi8(dst_acc);

          // process the tail elements when there remain less than
          // can be fitted in a SIMD vector
          for (; k < alig->originalNumberOfResidues; k++) {
              int gapi = (datai[k] == indet) || (datai[k] == '-');
              int gapj = (dataj[k] == indet) || (dataj[k] == '-');
              dst += (!(gapi && gapj)) && (!skipResidues[k]);
              hit += (!(gapi && gapj)) && (!skipResidues[k]) && (datai[k] == dataj[k]);
          }

          if (dst == 0) {
              debug.report(ErrorCode::NoResidueSequences,
                  new std::string[2]
                      {
                          alig->seqsName[i],
                          alig->seqsName[j]
                      });
              alig->identities[i][j] = 0;
          } else {
              // Identity score between two sequences is the ratio of identical residues
              // by the total length (common and no-common residues) among them
              alig->identities[i][j] = (float) hit / dst;
          }

          alig->identities[j][i] = alig->identities[i][j];
      }
  }
  
  // free SIMD buffers
  free(skipResidues);
}

bool SSE2Cleaner::calculateSpuriousVector(float overlap, float *spuriousVector) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("bool SSE2Cleaner::calculateSpuriousVector(float overlap, float *spuriousVector) ");

    // abort if there is not output vector to write to
    if (spuriousVector == nullptr)
        return false;

    // allocate aligned memory for faster SIMD loads of partial column sums
    uint32_t* hits    = ALIGNED_ALLOC(alig->originalNumberOfResidues, uint32_t);
    uint8_t*  hits_u8 = ALIGNED_ALLOC(alig->originalNumberOfResidues, uint8_t);

    // compute number of sequences from overlap threshold
    uint32_t  ovrlap  = uint32_t(ceil(overlap * float(alig->originalNumberOfSequences - 1)));

    // Depending on alignment type, indetermination symbol will be one or other
    char indet = (alig->getAlignmentType() & SequenceTypes::AA) ? 'X' : 'N';

    // prepare constant SIMD vectors
    const __m128i allindet  = _mm_set1_epi8(indet);
    const __m128i allgap    = _mm_set1_epi8('-');
    const __m128i ONES      = _mm_set1_epi8(1);

    // for each sequence in the alignment, computes its overlap
    for (int i = 0; i < alig->originalNumberOfSequences; i++) {

        // reset hits count
        memset(&hits[0],    0, alig->originalNumberOfResidues*sizeof(uint32_t));
        memset(&hits_u8[0], 0, alig->originalNumberOfResidues*sizeof(uint8_t));

        // compare sequence to other sequences for every position
        for (int j = 0; j < alig->originalNumberOfSequences; j++) {

            // don't compare sequence to itself
            if (j == i)
                continue;

            const char* datai = alig->sequences[i].data();
            const char* dataj = alig->sequences[j].data();

            __m128i seqi;
            __m128i seqj;
            __m128i eq;
            __m128i gapsi;
            __m128i gapsj;
            __m128i gaps;
            __m128i n;
            __m128i hit;

            int k = 0;

            // run iterations in SIMD while possible
            for (; ((int) (k + NLANES_8)) <= alig->originalNumberOfResidues; k += NLANES_8) {
                // load data for the sequences
                seqi = _mm_loadu_si128((const __m128i*) (&datai[k]));
                seqj = _mm_loadu_si128((const __m128i*) (&dataj[k]));
                // find which sequence characters are gap or indet
                gapsi = _mm_or_si128(_mm_cmpeq_epi8(seqi, allgap), _mm_cmpeq_epi8(seqi, allindet));
                gapsj = _mm_or_si128(_mm_cmpeq_epi8(seqj, allgap), _mm_cmpeq_epi8(seqj, allindet));
                gaps  = _mm_andnot_si128( _mm_or_si128(gapsi, gapsj), _mm_set1_epi8(0xFF));
                // find which sequence characters match
                eq = _mm_cmpeq_epi8(seqi, seqj);
                // find position where either not both characters are gap, or they are equal
                n = _mm_and_si128(_mm_or_si128(eq, gaps), ONES);
                // update counters
                hit = _mm_load_si128((const __m128i*) (&hits_u8[k]));
                hit = _mm_add_epi8(hit, n);
                _mm_store_si128((__m128i*) (&hits_u8[k]), hit);
            }

            // process the tail elements when there remain less than
            // can be fitted in a SIMD vector
            for (; k < alig->originalNumberOfResidues; k++) {
                int nongapi = (datai[k] != indet) && (datai[k] != '-');
                int nongapj = (dataj[k] != indet) && (dataj[k] != '-');
                hits_u8[k] += ((nongapi && nongapj) || (datai[k] == dataj[k]));
            }

            // we can process up to UINT8_MAX sequences, otherwise hits_u8[k]
            // may overflow, so every UINT8_MAX iterations we transfer the
            // partial hit counts from `hits_u8` to `hits`
            if ((j % UINT8_MAX) == 0) {
                for (k = 0; k < alig->originalNumberOfResidues; k++) hits[k] += hits_u8[k];
                memset(hits_u8, 0, alig->originalNumberOfResidues*sizeof(uint8_t));
            }
        }

        // update counters after last loop
        for (int k = 0; k < alig->originalNumberOfResidues; k++) hits[k] += hits_u8[k];

        // compute number of good positions in for sequence i
        uint32_t seqValue = 0;
        for (int k = 0; k < alig->originalNumberOfResidues; k++)
            if (hits[k] >= ovrlap)
                seqValue++;

        // compute overlap of current sequence as the fraction of columns
        // above overlap threshold
        spuriousVector[i] = ((float) seqValue / alig->originalNumberOfResidues);
    }

    // free temporary memory
    free(hits);
    free(hits_u8);

    // If there is not problem in the method, return true
    return true;
}

#endif
