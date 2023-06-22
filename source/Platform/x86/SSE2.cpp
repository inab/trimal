/* *****************************************************************************

    trimAl v2.0: a tool for automated alignment trimming in large-scale
                 phylogenetics analyses.

    readAl v2.0: a tool for automated alignment conversion among different
                 formats.

    2022-2023
        Larralde, M.            (martin.larralde@embl.de)

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

#include <climits>
#include <cstdint>
#include <emmintrin.h>

#include "Alignment/Alignment.h"
#include "InternalBenchmarker.h"
#include "Statistics/Gaps.h"
#include "Statistics/Manager.h"
#include "Statistics/Similarity.h"
#include "Platform/x86/SSE2.h"
#include "Platform/template.h"
#include "defines.h"
#include "reportsystem.h"
#include "utils.h"


class SSE2Vector {
private:
  __m128i vector;
  inline SSE2Vector(__m128i vec) : vector(vec) {}

public:
  const static size_t LANES = 16;
  const static size_t SIZE = sizeof(__m128i);

  inline SSE2Vector() : vector(_mm_setzero_si128()) {}

  inline static SSE2Vector duplicate(const uint8_t value) {
    return SSE2Vector(_mm_set1_epi8(value));
  }

  inline static SSE2Vector load(const uint8_t *data) {
    return SSE2Vector(_mm_load_si128((const __m128i *)data));
  }

  inline static SSE2Vector loadu(const uint8_t *data) {
    return SSE2Vector(_mm_loadu_si128((const __m128i *)data));
  }

  inline void store(uint8_t *data) const {
    _mm_store_si128((__m128i *)data, vector);
  }

  inline void storeu(uint8_t *data) const {
    _mm_storeu_si128((__m128i *)data, vector);
  }

  inline SSE2Vector &operator+=(const SSE2Vector &rhs) {
    vector = _mm_add_epi8(vector, rhs.vector);
    return *this;
  }

  inline SSE2Vector operator==(const SSE2Vector &rhs) const {
    return SSE2Vector(_mm_cmpeq_epi8(vector, rhs.vector));
  }

  inline SSE2Vector operator&(const SSE2Vector &rhs) const {
    return SSE2Vector(_mm_and_si128(vector, rhs.vector));
  }

  inline SSE2Vector operator|(const SSE2Vector &rhs) const {
    return SSE2Vector(_mm_or_si128(vector, rhs.vector));
  }

  inline SSE2Vector operator!() const {
    return SSE2Vector(_mm_andnot_si128(vector, _mm_set1_epi8(0xFF)));
  }

  inline SSE2Vector andnot(const SSE2Vector &rhs) const {
    return SSE2Vector(_mm_andnot_si128(rhs.vector, vector));
  }

  inline uint16_t sum() const {
    __m128i vsum = _mm_sad_epu8(vector, _mm_setzero_si128());
    return _mm_extract_epi16(vsum, 0) + _mm_extract_epi16(vsum, 4);
  }

  inline void clear() { vector = _mm_setzero_si128(); }
};

namespace statistics {

bool SSE2Similarity::calculateVectors(bool cutByGap) {
  StartTiming("bool SSE2Similarity::calculateVectors(bool cutByGap) ");
  return simd::calculateSimilarityVectors<SSE2Vector>(*this, cutByGap);
}

void SSE2Gaps::CalculateVectors() {
  StartTiming("bool SSE2Gaps::CalculateVectors() ");
  simd::calculateGapVectors<SSE2Vector>(*this);
}

void SSE2Identity::calculateSeqIdentity() {
  StartTiming("void SSE2Identity::calculateSeqIdentity() ");
  simd::calculateSeqIdentity<SSE2Vector>(*this);
}

bool SSE2Overlap::calculateSpuriousVector(float overlap, float *spuriousVector) {
  StartTiming("bool SSE2Overlap::calculateSpuriousVector(float overlap, float "
              "*spuriousVector) ");
  return simd::calculateSpuriousVector<SSE2Vector>(*this, overlap,
                                                  spuriousVector);
}

} // namespace statistics
