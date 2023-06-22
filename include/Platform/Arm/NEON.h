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

#ifndef TRIMAL_PLATFORM_X86_NEON_H
#define TRIMAL_PLATFORM_X86_NEON_H

#include "Statistics/Identity.h"
#include "Statistics/Gaps.h"
#include "Statistics/Overlap.h"
#include "Statistics/Similarity.h"

namespace statistics {

class NEONSimilarity : public Similarity {
public:
  NEONSimilarity(Alignment *parentAlignment) : Similarity(parentAlignment) {}
  NEONSimilarity(Alignment *parentAlignment, Similarity *parentSimilarity) :
    Similarity(parentAlignment, parentSimilarity) {}
  bool calculateVectors(bool cutByGap) override;
};
class NEONGaps : public Gaps {
public:
  NEONGaps(Alignment *parentAlignment) : Gaps(parentAlignment) {}
  NEONGaps(Alignment *parentAlignment, Gaps *parentGaps):
    Gaps(parentAlignment, parentGaps) {}
  void CalculateVectors() override;
};

class NEONOverlap : public Overlap {
public:
  NEONOverlap(Alignment *parent) : Overlap(parent) {}
  NEONOverlap(Alignment *parent, Overlap *parentOverlap) : Overlap(parent, parentOverlap) {}
  bool calculateSpuriousVector(float overlap, float *spuriousVector) override;
};

class NEONIdentity : public Identity {
public:
  NEONIdentity(Alignment *parent) : Identity(parent) {}
  NEONIdentity(Alignment *parent, Identity *parentIdentity) : Identity(parent, parentIdentity) {}
  void calculateSeqIdentity() override;
};

} // namespace statistics

#endif