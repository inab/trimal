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

#include <arm_neon.h>
#include <climits>
#include <cstdint>

#include "Alignment/Alignment.h"
#include "InternalBenchmarker.h"
#include "Statistics/Manager.h"
#include "Statistics/Similarity.h"
#include "Platform/Arm/NEON.h"
#include "Platform/template.h"
#include "defines.h"
#include "reportsystem.h"
#include "utils.h"


class NEONVector {
private:
  uint8x16_t vector;
  inline NEONVector(uint8x16_t vec) : vector(vec) {}

public:
  const static size_t LANES = 16;
  const static size_t SIZE = sizeof(uint8x16_t);

  inline NEONVector() : vector(vdupq_n_u8(0)) {}

  inline static NEONVector duplicate(const uint8_t value) {
    return NEONVector(vdupq_n_u8(value));
  }

  inline static NEONVector load(const uint8_t *data) {
    return NEONVector(vld1q_u8(data));
  }

  inline static NEONVector loadu(const uint8_t *data) {
    return NEONVector(vld1q_u8(data));
  }

  inline void store(uint8_t *data) const { vst1q_u8(data, vector); }

  inline void storeu(uint8_t *data) const { vst1q_u8(data, vector); }

  inline NEONVector &operator+=(const NEONVector &rhs) {
    vector = vaddq_u8(vector, rhs.vector);
    return *this;
  }

  inline NEONVector operator==(const NEONVector &rhs) const {
    return NEONVector(vceqq_u8(vector, rhs.vector));
  }

  inline NEONVector operator&(const NEONVector &rhs) const {
    return NEONVector(vandq_u8(vector, rhs.vector));
  }

  inline NEONVector operator|(const NEONVector &rhs) const {
    return NEONVector(vorrq_u8(vector, rhs.vector));
  }

  inline NEONVector operator!() const { return NEONVector(vmvnq_u8(vector)); }

  inline NEONVector andnot(const NEONVector &rhs) const {
    return NEONVector(vbicq_u8(vector, rhs.vector));
  }

  inline uint16_t sum() const {
#ifdef __aarch64__
    // Don't use `vaddvq_u8` because it can overflow, first add pairwise
    // into 16-bit lanes to avoid overflows
    return vaddvq_u16(vpaddlq_u8(vector));
#else
    uint64x2_t paired = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(vector)));
    return vgetq_lane_u64(paired, 0) + vgetq_lane_u64(paired, 1);
#endif
  }

  inline void clear() { vector = vdupq_n_u8(0); }
};

namespace statistics {

bool NEONSimilarity::calculateVectors(bool cutByGap) {
  StartTiming("bool NEONSimilarity::calculateVectors(bool cutByGap) ");
  return simd::calculateSimilarityVectors<NEONVector>(*this, cutByGap);
}

void NEONGaps::CalculateVectors() {
  StartTiming("bool NEONGaps::CalculateVectors() ");
  simd::calculateGapVectors<NEONVector>(*this);
}

void NEONIdentity::calculateSeqIdentity() {
  StartTiming("void NEONIdentity::calculateSeqIdentity() ");
  simd::calculateSeqIdentity<NEONVector>(*this);
}

bool NEONOverlap::calculateSpuriousVector(float overlap, float *spuriousVector) {
  StartTiming("bool NEONOverlap::calculateSpuriousVector(float overlap, float "
              "*spuriousVector) ");
  return simd::calculateSpuriousVector<NEONVector>(*this, overlap,
                                                  spuriousVector);
}

} // namespace statistics