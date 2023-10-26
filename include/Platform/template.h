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

#ifndef TRIMAL_PLATFORM_TEMPLATE_H
#define TRIMAL_PLATFORM_TEMPLATE_H

#include <cstdlib>
#include <climits>

#include "Alignment/Alignment.h"
#include "InternalBenchmarker.h"
#include "Statistics/Manager.h"
#include "Statistics/Similarity.h"
#include "Statistics/similarityMatrix.h"
#include "Statistics/Identity.h"
#include "defines.h"
#include "reportsystem.h"
#include "utils.h"

namespace simd {

template <class T, class Vector> 
T *aligned_array(size_t n) {
  T* ptr = nullptr;
  if (Vector::SIZE > 1) {
    const size_t mask = Vector::SIZE - 1;
    const size_t size = (n * sizeof(T) + mask) & (~mask);
    void* memptr = nullptr;
    if (posix_memalign(&memptr, Vector::SIZE, size) == 0)
      ptr = static_cast<T *>(memptr);
  } else {
    const size_t size = (n * sizeof(T));
    ptr = static_cast<T *>(malloc(size));
  }
  if (ptr == nullptr)
    throw std::bad_alloc();
  return ptr;
}

template <class Vector>
inline bool calculateSimilarityVectors(statistics::Similarity &s,
                                       bool cutByGap) {
  // A similarity matrix must be defined. If not, return false
  if (s.simMatrix == nullptr)
    return false;

  // Get the similarity matrix raw storage
  const float** distMat = s.simMatrix->getDistanceMatrix();

  s.alig->Statistics->calculateSeqIdentity();
  float *identities = s.alig->Statistics->identity->identities;

  // Create the variable gaps, in case we want to cut by gaps
  int *gaps = nullptr;

  // Retrieve the gaps values in case we want to set to 0 the similarity value
  // in case the gaps value for that column is bigger or equal to 0.8F
  if (cutByGap) {
    if (s.alig->Statistics->gaps == nullptr)
      s.alig->Statistics->calculateGapStats();
    gaps = s.alig->Statistics->gaps->getGapsWindow();
  }

  // Initialize the variables used
  int i, j, k;
  float num, den;
  size_t arrayIdentityPosition;

  // Depending on alignment type, indetermination symbol will be one or other
  char indet = s.alig->getAlignmentType() & SequenceTypes::AA ? 'X' : 'N';

  // Q temporal value
  float Q;
  // Temporal chars that will contain the residues to compare by pair.
  int numA, numB;

  // Calculate the maximum number of gaps a column can have to calculate it's
  //      similarity
  float gapThreshold = 0.8F * s.alig->numberOfResidues;

  // Cache pointers to matrix rows to avoid dereferencing in inner loops
  const float *distRow;

  // Create buffers to store column data
  std::vector<char> colnum =
      std::vector<char>(s.alig->originalNumberOfSequences);
  std::vector<char> colgap =
      std::vector<char>(s.alig->originalNumberOfSequences);

  // For each column calculate the Q value and the MD value using an equation
  for (i = 0; i < s.alig->originalNumberOfResidues; i++) {
    // Set MDK for columns with gaps values bigger or equal to threshold
    if ((gaps != nullptr) && gaps[i] >= gapThreshold) {
      s.MDK[i] = 0.F;
      continue;
    }

    // Fill the column buffer with the current column and check
    // characters are well-defined with respect to the similarity matrix
    for (j = 0; j < s.alig->originalNumberOfSequences; j++) {
      char letter = utils::toUpper(s.alig->sequences[j][i]);
      if ((letter == indet) || (letter == '-')) {
        colgap[j] = 1;
      } else {
        colgap[j] = 0;
        if ((letter < 'A') || (letter > 'Z')) {
          debug.report(ErrorCode::IncorrectSymbol,
                       new std::string[1]{std::string(1, letter)});
          return false;
        } else {
            int index = s.simMatrix->getLetterIndex(letter);
            if (index == -1) {
                debug.report(ErrorCode::UndefinedSymbol,
                            new std::string[1]{std::string(1, letter)});
                return false;
            } else {
                colnum[j] = index;
            }
        }
      }
    }

    // For each AAs/Nucleotides' pair in the column we compute its distance
    arrayIdentityPosition = 0;
    for (j = 0, num = 0, den = 0; j < s.alig->originalNumberOfSequences; j++) {
      // We don't compute the distance if the first element is
      // a indeterminate (XN) or a gap (-) element.
      if (colgap[j]) {
        arrayIdentityPosition += s.alig->originalNumberOfSequences - j - 1;
        continue;
      }

      // Get the index of the first residue
      // and cache pointers to matrix rows
      numA = colnum[j];
      distRow = distMat[numA];

      for (k = j + 1; k < s.alig->originalNumberOfSequences; k++) {
        // We don't compute the distance if the second element is
        //      a indeterminate (XN) or a gap (-) element
        if (colgap[k]) {
          arrayIdentityPosition++;
          continue;
        }

        // Get the index of the second residue and compute
        // fraction with identity value for the two pairs and
        // its distance based on similarity matrix's value.
        numB = colnum[k];
        num += (1.0F - identities[arrayIdentityPosition]) * distRow[numB];
        den += (1.0F - identities[arrayIdentityPosition]);
        arrayIdentityPosition++;
      }
    }

    // If we are processing a column with only one AA/nucleotide, MDK = 0
    if (den == 0) {
      s.MDK[i] = 0;
    } else {
      Q = num / den;
      // If the MDK value is more than 1, we normalized this value to 1.
      //      Only numbers higher than 0 yield exponents higher than 1
      //      Using this we can test if the result is going to be higher than 1.
      //      And thus, prevent calculating the exp.
      // Take in mind that the Q is negative, so we must test if Q is LESSER
      //      than one, not bigger.
      if (Q < 0)
        s.MDK[i] = 1.F;
      else
        s.MDK[i] = exp(-Q);
    }
  }

  return true;
}

template <class Vector>
inline bool calculateSpuriousVector(statistics::Overlap &o, const float overlap,
                                    float *spuriousVector) {
  // abort if there is not output vector to write to
  if (spuriousVector == nullptr)
    return false;

  // Get the alignment from the Cleaner object
  const Alignment* alig = o.alig;

  // compute number of sequences from overlap threshold
  uint32_t ovrlap =
      uint32_t(ceil(overlap * float(alig->originalNumberOfSequences - 1)));

  // Depending on alignment type, indetermination symbol will be one or other
  char indet = (alig->getAlignmentType() & SequenceTypes::AA) ? 'X' : 'N';

  // prepare constant SIMD vectors
  const Vector ALLINDET = Vector::duplicate(indet);
  const Vector ALLGAP = Vector::duplicate('-');
  const Vector ONES = Vector::duplicate(1);

  // allocate aligned memory for faster SIMD loads
  uint32_t *hits =
      aligned_array<uint32_t, Vector>(alig->originalNumberOfResidues);
  uint8_t *hits_u8 =
      aligned_array<uint8_t, Vector>(alig->originalNumberOfResidues);

  // for each sequence in the alignment, computes its overlap
  for (int i = 0; i < alig->originalNumberOfSequences; i++) {

    // reset hits count
    memset(&hits[0], 0, alig->originalNumberOfResidues * sizeof(uint32_t));
    memset(&hits_u8[0], 0, alig->originalNumberOfResidues * sizeof(uint8_t));

    unsigned int processedSequences = 0;
    const uint8_t *datai =
        reinterpret_cast<const uint8_t *>(alig->sequences[i].data());

    // compare sequence to other sequences for every position
    for (int j = 0; j < alig->originalNumberOfSequences; j++) {

      // don't compare sequence to itself
      if (j == i)
        continue;

      const uint8_t *dataj =
          reinterpret_cast<const uint8_t *>(alig->sequences[j].data());

      int k = 0;

      // run iterations in SIMD while possible
      for (; ((int)(k + Vector::LANES)) <= alig->originalNumberOfResidues;
           k += Vector::LANES) {
        // load data for the sequences
        const Vector seqi = Vector::loadu(&datai[k]);
        const Vector seqj = Vector::loadu(&dataj[k]);
        // find which sequence characters are gap or indet
        const Vector gapsi = (seqi == ALLGAP) | (seqi == ALLINDET);
        const Vector gapsj = (seqj == ALLGAP) | (seqj == ALLINDET);
        const Vector gaps = !(gapsi | gapsj);
        // find which sequence characters match
        const Vector eq = (seqi == seqj);
        // find position where either not both characters are gap, or they are
        // equal
        const Vector n = (eq | gaps) & ONES;
        // update counters
        Vector hit = Vector::load(&hits_u8[k]);
        hit += n;
        hit.store(&hits_u8[k]);
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
      if ((processedSequences % UINT8_MAX) == 0) {
        for (k = 0; k < alig->originalNumberOfResidues; k++)
          hits[k] += hits_u8[k];
        memset(hits_u8, 0, alig->originalNumberOfResidues * sizeof(uint8_t));
      }
      processedSequences++;
    }

    // update counters after last loop
    for (int k = 0; k < alig->originalNumberOfResidues; k++)
      hits[k] += hits_u8[k];

    // compute number of good positions in for sequence i
    uint32_t seqValue = 0;
    for (int k = 0; k < alig->originalNumberOfResidues; k++)
      if (hits[k] >= ovrlap)
        seqValue++;

    // compute overlap of current sequence as the fraction of columns
    // above overlap threshold
    spuriousVector[i] = ((float)seqValue / alig->originalNumberOfResidues);
  }

  // free allocated memory
  free(hits);
  free(hits_u8);

  // If there is not problem in the method, return true
  return true;
}

template <class Vector> 
inline void calculateSeqIdentity(statistics::Identity &id) {
  int i, j, k, l;
  size_t arrayIdentitySize, arrayIdentityPosition;
  const Alignment* alig = id.alig;

  // create identities matrix to store identities scores
  arrayIdentitySize = ((float) alig->originalNumberOfSequences * alig->originalNumberOfSequences + alig->originalNumberOfSequences) / 2;
  id.identities = new float[arrayIdentitySize];

  // Depending on alignment type, indetermination symbol will be one or other
  char indet = (alig->getAlignmentType() & SequenceTypes::AA) ? 'X' : 'N';

  // prepare constant SIMD vectors
  const Vector ALLINDET = Vector::duplicate(indet);
  const Vector ALLGAP = Vector::duplicate('-');
  const Vector ONES = Vector::duplicate(1);

  // create an index of residues to skip
  uint8_t *skipResidues = aligned_array<uint8_t, Vector>(
      alig->originalNumberOfResidues); // ALIGNED_ALLOC(alig->originalNumberOfResidues,
                                         // uint8_t);
  for (int k = 0; k < alig->originalNumberOfResidues; k++) {
    skipResidues[k] = alig->saveResidues[k] == -1 ? 0xFF : 0;
  }

  // For each seq, compute its identity score against the others in the MSA
  arrayIdentityPosition = 0;
  for (i = 0; i < alig->originalNumberOfSequences; i++) {
    if (alig->saveSequences[i] == -1)
      continue;

    const uint8_t *datai =
        reinterpret_cast<const uint8_t *>(alig->sequences[i].data());

    // Compute identity scores for the current sequence against the rest
    for (j = i + 1; j < alig->originalNumberOfSequences; j++) {
      if (alig->saveSequences[j] == -1)
        continue;

      const uint8_t *dataj =
          reinterpret_cast<const uint8_t *>(alig->sequences[j].data());

      Vector dst_acc = Vector();
      Vector hit_acc = Vector();

      int hit = 0;
      int dst = 0;

      // run with unrolled loops of UINT8_MAX iterations first
      for (k = 0; ((int)(k + Vector::LANES * UINT8_MAX)) <
                  alig->originalNumberOfResidues;) {
        for (l = 0; l < UINT8_MAX; l++, k += Vector::LANES) {
          // load data for the sequences
          Vector seqi = Vector::loadu(&datai[k]);
          Vector seqj = Vector::loadu(&dataj[k]);
          Vector skip = Vector::load(&skipResidues[k]);
          Vector eq = (seqi == seqj);
          // find which sequence characters are gap or indet
          Vector gapsi = ((seqi == ALLGAP) | (seqi == ALLINDET));
          Vector gapsj = ((seqj == ALLGAP) | (seqj == ALLINDET));
          // find position where not both characters are gap
          Vector mask = ONES.andnot(gapsi & gapsj).andnot(skip);
          // update counters
          dst_acc += mask;
          hit_acc += (eq & mask);
        }
        // merge accumulators
        dst += dst_acc.sum();
        hit += hit_acc.sum();
        dst_acc.clear();
        hit_acc.clear();
      }

      // run remaining iterations in SIMD while possible
      for (; ((int)(k + Vector::LANES)) < alig->originalNumberOfResidues;
           k += Vector::LANES) {
        // load data for the sequences; load is unaligned because
        // string data is not guaranteed to be aligned.
        Vector seqi = Vector::loadu(&datai[k]);
        Vector seqj = Vector::loadu(&dataj[k]);
        Vector skip = Vector::load(&skipResidues[k]);
        Vector eq = (seqi == seqj);
        // find which sequence characters are gap or indet
        Vector gapsi = ((seqi == ALLGAP) | (seqi == ALLINDET));
        Vector gapsj = ((seqj == ALLGAP) | (seqj == ALLINDET));
        // find position where not both characters are gap
        Vector mask = ONES.andnot(gapsi & gapsj).andnot(skip);
        // update counters
        dst_acc += mask;
        hit_acc += (eq & mask);
      }

      // update counters after last loop
      hit += hit_acc.sum();
      dst += dst_acc.sum();

      // process the tail elements when there remain less than
      // can be fitted in a SIMD vector
      for (; k < alig->originalNumberOfResidues; k++) {
        int gapi = (datai[k] == indet) || (datai[k] == '-');
        int gapj = (dataj[k] == indet) || (dataj[k] == '-');
        dst += (!(gapi && gapj)) && (!skipResidues[k]);
        hit +=
            (!(gapi && gapj)) && (!skipResidues[k]) && (datai[k] == dataj[k]);
      }

      if (dst == 0) {
        id.identities[arrayIdentityPosition] = 0;
      } else {
        // Identity score between two sequences is the ratio of identical
        // residues by the total length (common and no-common residues) among
        // them
        id.identities[arrayIdentityPosition] = (float)hit / dst;
      }

      arrayIdentityPosition++;
    }
  }

  // free allocated memory
  free(skipResidues);
}

template <class Vector> 
inline void calculateGapVectors(statistics::Gaps &g) {
  int i, j;
  const Vector ALLGAP = Vector::duplicate('-');
  const Vector ONES = Vector::duplicate(1);
  const Alignment* alig = g.alig;


  // use temporary buffer for storing 8-bit partial sums
  uint8_t *gapsInColumn_u8 =
      aligned_array<uint8_t, Vector>(alig->originalNumberOfResidues);
  memset(g.gapsInColumn, 0, sizeof(int) * alig->originalNumberOfResidues);
  memset(gapsInColumn_u8, 0,
         sizeof(uint8_t) * alig->originalNumberOfResidues);

  // count gaps per column
  for (j = 0; j < alig->originalNumberOfSequences; j++) {
    // skip sequences not retained in alignment
    if (alig->saveSequences[j] == -1)
      continue;
    // process the whole sequence, LANES lanes at a time
    const uint8_t *data =
        reinterpret_cast<const uint8_t *>(alig->sequences[j].data());
    for (i = 0; ((int)(i + Vector::LANES)) < alig->originalNumberOfResidues;
         i += Vector::LANES) {
      Vector letters = Vector::loadu(&data[i]);
      Vector counts = Vector::load(&gapsInColumn_u8[i]);
      counts += (ONES & (letters == ALLGAP));
      counts.store(&gapsInColumn_u8[i]);
    }
    // count the remaining gap elements without SIMD
    for (; i < alig->originalNumberOfResidues; i++)
      if (data[i] == '-')
        gapsInColumn_u8[i]++;
    // every UINT8_MAX iterations the accumulator may overflow, so the
    // temporary counts are moved into the final counter array, and the
    // accumulator is reset
    if (j % UINT8_MAX == 0) {
      for (i = 0; i < alig->originalNumberOfResidues; i++)
        g.gapsInColumn[i] += gapsInColumn_u8[i];
      memset(gapsInColumn_u8, 0,
             sizeof(uint8_t) * alig->originalNumberOfResidues);
    }
  }
  // collect the remaining partial counts into the final counter array
  for (i = 0; i < alig->originalNumberOfResidues; i++)
    g.gapsInColumn[i] += gapsInColumn_u8[i];

  // free temporary buffer
  free(gapsInColumn_u8);

  // build histogram and find largest number of gaps
  for (int i = 0; i < alig->originalNumberOfResidues; i++) {
    // g.totalGaps += g.gapsInColumn[i];
    g.numColumnsWithGaps[g.gapsInColumn[i]]++;
    if (g.gapsInColumn[i] > g.maxGaps)
      g.maxGaps = g.gapsInColumn[i];
  }
}
} // namespace simd

#endif