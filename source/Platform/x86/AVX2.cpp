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
#include <immintrin.h>

#include "Alignment/Alignment.h"
#include "InternalBenchmarker.h"
#include "Statistics/Gaps.h"
#include "Statistics/Manager.h"
#include "Statistics/Similarity.h"
#include "Platform/x86/AVX2.h"
#include "Platform/template.h"
#include "defines.h"
#include "reportsystem.h"
#include "utils.h"
#include "cpuinfo_x86.h"

using namespace cpu_features;


class AVX2Vector {
private:
  __m256i vector;
  inline AVX2Vector(__m256i vec) : vector(vec) {}

public:
  const static size_t LANES = 32;
  const static size_t SIZE = sizeof(__m256i);

  inline AVX2Vector() : vector(_mm256_setzero_si256()) {}

  inline static AVX2Vector duplicate(const uint8_t value) {
    return AVX2Vector(_mm256_set1_epi8(value));
  }

  inline static AVX2Vector load(const uint8_t *data) {
    return AVX2Vector(_mm256_load_si256((const __m256i *)data));
  }

  inline static AVX2Vector loadu(const uint8_t *data) {
    return AVX2Vector(_mm256_loadu_si256((const __m256i *)data));
  }

  inline void store(uint8_t *data) const {
    _mm256_store_si256((__m256i *)data, vector);
  }

  inline void storeu(uint8_t *data) const {
    _mm256_storeu_si256((__m256i *)data, vector);
  }

  inline AVX2Vector &operator+=(const AVX2Vector &rhs) {
    vector = _mm256_add_epi8(vector, rhs.vector);
    return *this;
  }

  inline AVX2Vector operator==(const AVX2Vector &rhs) const {
    return AVX2Vector(_mm256_cmpeq_epi8(vector, rhs.vector));
  }

  inline AVX2Vector operator&(const AVX2Vector &rhs) const {
    return AVX2Vector(_mm256_and_si256(vector, rhs.vector));
  }

  inline AVX2Vector operator|(const AVX2Vector &rhs) const {
    return AVX2Vector(_mm256_or_si256(vector, rhs.vector));
  }

  inline AVX2Vector operator!() const {
    return AVX2Vector(_mm256_andnot_si256(vector, _mm256_set1_epi8(0xFF)));
  }

  inline AVX2Vector andnot(const AVX2Vector &rhs) const {
    return AVX2Vector(_mm256_andnot_si256(rhs.vector, vector));
  }

#ifdef CPU_FEATURES_ARCH_X86_64
  inline uint16_t sum() const {
    __m256i sum256 = _mm256_sad_epu8(vector, _mm256_setzero_si256());
    __m128i sum128 = _mm_add_epi64(_mm256_extractf128_si256(sum256, 1), _mm256_castsi256_si128(sum256));
    return _mm_extract_epi64(sum128, 0) + _mm_extract_epi64(sum128, 1);
  }
#else
  inline uint16_t sum() const {
    __m256i sum256 = _mm256_sad_epu8(vector, _mm256_setzero_si256());
    __m128i sum128 = _mm_add_epi64(_mm256_extractf128_si256(sum256, 1), _mm256_castsi256_si128(sum256));
    return _mm_extract_epi32(sum128, 0) + _mm_extract_epi32(sum128, 2);
  }
#endif


  inline void clear() { vector = _mm256_setzero_si256(); }
};

namespace statistics {

bool AVX2Similarity::calculateVectors(bool cutByGap) {
  StartTiming("bool AVX2Similarity::calculateVectors(bool cutByGap) ");
  return simd::calculateSimilarityVectors<AVX2Vector>(*this, cutByGap);
}

void AVX2Gaps::CalculateVectors() {
  StartTiming("bool AVX2Gaps::CalculateVectors() ");
  simd::calculateGapVectors<AVX2Vector>(*this);
}

void AVX2Identity::calculateSeqIdentity() {
  StartTiming("void AVX2Identity::calculateSeqIdentity() ");
  simd::calculateSeqIdentity<AVX2Vector>(*this);
}

bool AVX2Overlap::calculateSpuriousVector(float overlap, float *spuriousVector) {
  StartTiming("bool AVX2Overlap::calculateSpuriousVector(float overlap, float "
              "*spuriousVector) ");
  return simd::calculateSpuriousVector<AVX2Vector>(*this, overlap,
                                                  spuriousVector);
}

} // namespace statistics

