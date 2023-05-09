/* *****************************************************************************

    trimAl v2.0: a tool for automated alignment trimming in large-scale
                 phylogenetics analyses.

    readAl v2.0: a tool for automated alignment conversion among different
                 formats.

    2022-2023
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

#include "Statistics/Overlap.h"
#include "InternalBenchmarker.h"
#include "Statistics/Manager.h"
#include "reportsystem.h"
#include "Alignment/Alignment.h"
#include "defines.h"
#include "utils.h"

namespace statistics {

    Overlap::Overlap(Alignment *parentAlignment) {
        StartTiming("Overlap::Overlap(Alignment *parentAlignment) ");
       
        alig = parentAlignment;
        refCounter = new int(1);
    }

    Overlap::Overlap(Alignment *parentAlignment, Overlap *mold) {
        StartTiming("Overlap::Overlap(Alignment *parentAlignment, Overlap *mold) ");

        alig = parentAlignment;

        overlaps = mold->overlaps;
        refCounter = mold->refCounter;
        (*refCounter)++;
    }

    Overlap::~Overlap() {
        if (refCounter == nullptr || --(*refCounter) < 1) {
            if (overlaps != nullptr) {
                for (int i = 0; i < alig->numberOfSequences; i++)
                    delete[] overlaps[i];
                delete[] overlaps;
            }
            delete refCounter;
            refCounter = nullptr;
        }
    }

    void Overlap::calculateSeqOverlap() {
        // Create a timer that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("void Overlap::calculateSeqOverlap(void) ");

        // Compute the overlap between sequences taken each of them as the reference
        // to compute such scores. It will lead to a non-symmetric matrix.
        int i, j, k, shared, referenceLength;
        char indet;

        // Depending on alignment type, indetermination symbol will be one or other
        indet = alig->getAlignmentType() & SequenceTypes::AA ? 'X' : 'N';

        // Create overlap matrix to store overlap scores
        overlaps = new float *[alig->numberOfSequences];

        // For each seq, compute its overlap score against the others in the MSA
        for (i = 0; i < alig->numberOfSequences; i++) {
            overlaps[i] = new float[alig->numberOfSequences];

            for (j = 0; j < alig->numberOfSequences; j++) {
                for (k = 0, shared = 0, referenceLength = 0; k < alig->numberOfResidues; k++) {
                    // If there a valid residue for the reference sequence, then see if
                    // there is a valid residue for the other sequence.
                    if ((alig->sequences[i][k] != indet) && (alig->sequences[i][k] != '-')) {
                        referenceLength++;
                        if ((alig->sequences[j][k] != indet) && (alig->sequences[j][k] != '-'))
                            shared++;
                    }
                }
                // Overlap score between two sequences is the ratio of shared valid
                // residues divided by the sequence length taken as reference. The
                // overlaps matrix, therefore, will be not symmetric.
                overlaps[i][j] = (float) shared / referenceLength;
            }
        }
    }

    bool Overlap::calculateSpuriousVector(float overlap, float *spuriousVector) {
        // Create a timer that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("bool Overlap::calculateSpuriousVector(float overlap, float *spuriousVector) ");

        int i, j, k, seqValue, ovrlap, hit;
        char indet;

        float floatOverlap = overlap * float(alig->originalNumberOfSequences - 1);
        ovrlap = int(overlap * (alig->originalNumberOfSequences - 1));
        if (floatOverlap > float(ovrlap))
            ovrlap++;

        if (spuriousVector == nullptr)
            return false;
        // Depending on the kind of Alignment, we have
        // different indetermination symbol
        if (alig->getAlignmentType() & SequenceTypes::AA)
            indet = 'X';
        else
            indet = 'N';
        // For each Alignment's sequence, computes its overlap
        for (i = 0, seqValue = 0; i < alig->originalNumberOfSequences; i++, seqValue = 0) {
            // For each Alignment's column, computes the overlap
            // between the selected sequence and the other ones
            for (j = 0; j < alig->originalNumberOfResidues; j++) {
                for (k = 0, hit = 0; k < alig->originalNumberOfSequences; k++) {
                    // Do not compare a sequence to itself
                    if (k == i) continue;
                    // If the element of sequence selected is the same
                    // that the element of sequence considered, computes
                    // a hit
                    if (alig->sequences[i][j] == alig->sequences[k][j])
                        hit++;
                        // If the element of sequence selected isn't a 'X' nor
                        // 'N' (indetermination) or a '-' (gap) and the element
                        // of sequence considered isn't a  a 'X' nor 'N'
                        // (indetermination) or a '-' (gap), computes a hit
                    else if ((alig->sequences[i][j] != indet) && (alig->sequences[i][j] != '-')
                            && (alig->sequences[k][j] != indet) && (alig->sequences[k][j] != '-'))
                        hit++;
                }
                // Finally, if the hit's number divided by number of
                // sequences minus one is greater or equal than
                // overlap's value, compute a sequence hit
                if (hit >= ovrlap)
                    seqValue++;
            }

            // For each Alignment's sequence, computes its spurious's
            // or overlap's value as the column's hits -for that
            // sequence- divided by column's number.
            spuriousVector[i] = ((float) seqValue / alig->originalNumberOfResidues);
        }

        // If there is not problem in the method, return true
        return true;
    }
}