/* *****************************************************************************

    trimAl v2.0: a tool for automated alignment trimming in large-scale
                 phylogenetics analyses.

    readAl v2.0: a tool for automated alignment conversion among different
                 formats.

    2022
        Larralde M.             (martin.larralde@embl.de)

    2009-2019
        Fernandez-Rodriguez V.  (victor.fernandez@bsc.es)
        Capella-Gutierrez S.    (salvador.capella@bsc.es)
        Gabaldon, T.            (tgabaldon@crg.es)

    This file is part of trimAl/readAl.

    trimAl/readAl are free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, the last available version.

    trimAl/readAl are distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with trimAl/readAl. If not, see <http://www.gnu.org/licenses/>.

***************************************************************************** */

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

#include "Platform/sse2/SSE2Similarity.h"

#define NLANES_8   sizeof(__m128i) / sizeof(uint8_t)   // number of 8-bit lanes in __m128i
#define NLANES_32  sizeof(__m128i) / sizeof(uint32_t)  // number of 32-bit lanes in __m128i

#define ALLOC_MASK           (sizeof(__m128i) - 1)
#define ALLOC_SIZE(N, T)     ((N * sizeof(T) + ALLOC_MASK) & (~ALLOC_MASK))
#define ALIGNED_ALLOC(N, T)  (static_cast<T*>(aligned_alloc(sizeof(__m128i), ALLOC_SIZE(N, T))))

static inline uint32_t _mm_hsum_epi8(__m128i a) {
    __m128i vsum = _mm_sad_epu8(a, _mm_setzero_si128());
    return _mm_extract_epi16(vsum, 0) + _mm_extract_epi16(vsum, 4);
}

using namespace statistics;

SSE2Similarity::SSE2Similarity(Alignment* parentAlignment): Similarity(parentAlignment) {}

void SSE2Similarity::calculateMatrixIdentity() {
    // Create a timerLevel that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("void SSE2Similarity::calculateMatrixIdentity() ");

    // We don't want to calculate the matrix identity
    // if it has been previously calculated
    if (matrixIdentity != nullptr)
        return;

    // declare indices
    int i, j, k, l;

    // Allocate memory for the matrix identity
    matrixIdentity = new float *[alig->originalNumberOfSequences];
    for (i = 0; i < alig->originalNumberOfSequences; i++) {
        matrixIdentity[i] = new float[alig->originalNumberOfSequences];
    }

    // Depending on alignment type, indetermination symbol will be one or other
    const char indet = alig->getAlignmentType() & SequenceTypes::AA ? 'X' : 'N';

    // prepare constant SIMD vectors
    const __m128i allindet = _mm_set1_epi8(indet);
    const __m128i allgap   = _mm_set1_epi8('-');
    const __m128i ONES     = _mm_set1_epi8(1);

    // For each sequences' pair, compare identity
    // #pragma omp parallel for num_threads(NUMTHREADS) private(j, k, l) if(alig->originalNumberOfSequences>MINPARALLELSIZE)
    for (i = 0; i < alig->originalNumberOfSequences; i++) {
        for (j = i + 1; j < alig->originalNumberOfSequences; j++) {

            const char* datai = alig->sequences[i].data();
            const char* dataj = alig->sequences[j].data();

            __m128i seqi;
            __m128i seqj;
            __m128i eq;
            __m128i gapsi;
            __m128i gapsj;
            __m128i len_acc = _mm_setzero_si128();
            __m128i sum_acc = _mm_setzero_si128();

            int sum = 0;
            int length = 0;

            // run with unrolled loops of UCHAR_MAX iterations first
            for (k = 0; ((int) (k + NLANES_8*UCHAR_MAX)) < alig->originalNumberOfResidues;) {
                // unroll the internal loop
                for (l = 0; l < UCHAR_MAX; l++, k += NLANES_8) {
                    // load data for the sequences
                    seqi = _mm_loadu_si128( (const __m128i*) (&datai[k]));
                    seqj = _mm_loadu_si128( (const __m128i*) (&dataj[k]));
                    // find which sequence characters are gap or indet
                    gapsi = _mm_or_si128(_mm_cmpeq_epi8(seqi, allgap), _mm_cmpeq_epi8(seqi, allindet));
                    gapsj = _mm_or_si128(_mm_cmpeq_epi8(seqj, allgap), _mm_cmpeq_epi8(seqj, allindet));
                    // find which sequence characters are equal
                    eq    = _mm_cmpeq_epi8(seqi, seqj);
                    // update counters
                    sum_acc = _mm_add_epi8(sum_acc, _mm_and_si128(eq, _mm_andnot_si128( _mm_or_si128(gapsi, gapsj), ONES)));
                    len_acc = _mm_add_epi8(len_acc,                   _mm_andnot_si128(_mm_and_si128(gapsi, gapsj), ONES));
                }
                // merge and reset accumulators
                sum    += _mm_hsum_epi8(sum_acc);
                length += _mm_hsum_epi8(len_acc);
                sum_acc = _mm_setzero_si128();
                len_acc = _mm_setzero_si128();
            }

            // run remaining iterations in SIMD while possible
            for (; ((int) (k + NLANES_8)) < alig->originalNumberOfResidues; k += NLANES_8) {
                // load data for the sequences; load is unaligned because
                // string data is not guaranteed to be aligned.
                seqi = _mm_loadu_si128( (const __m128i*) (&datai[k]));
                seqj = _mm_loadu_si128( (const __m128i*) (&dataj[k]));
                // find which sequence characters are either a gap or
                // an indeterminate character
                gapsi = _mm_or_si128(_mm_cmpeq_epi8(seqi, allgap), _mm_cmpeq_epi8(seqi, allindet));
                gapsj = _mm_or_si128(_mm_cmpeq_epi8(seqj, allgap), _mm_cmpeq_epi8(seqj, allindet));
                // find which sequence characters are equal
                eq    = _mm_cmpeq_epi8(seqi, seqj);
                // update counters: update sum if both sequence characters
                // are non-gap/indeterminate and equal; update length if
                // any of the characters is non-gap/indeterminate.
                sum_acc = _mm_add_epi8(sum_acc, _mm_and_si128(eq, _mm_andnot_si128( _mm_or_si128(gapsi, gapsj), ONES)));
                len_acc = _mm_add_epi8(len_acc,                   _mm_andnot_si128(_mm_and_si128(gapsi, gapsj), ONES));
            }

            // merge accumulators
            sum    += _mm_hsum_epi8(sum_acc);
            length += _mm_hsum_epi8(len_acc);

            // process the tail elements when there remain less than
            // can be fitted in a SIMD vector
            for (; k < alig->originalNumberOfResidues; k++) {
                int gapi = (datai[k] == '-') || (datai[k] == indet);
                int gapj = (dataj[k] == '-') || (dataj[k] == indet);
                sum    += (!gapi) && (!gapj) && (datai[k] == dataj[k]);
                length += (!gapi) || (!gapj);
            }

            // Calculate the value of matrix idn for columns j and i
            matrixIdentity[i][j] = matrixIdentity[j][i] = (1.0F - ((float)sum / length));
        }
    }
}

#endif