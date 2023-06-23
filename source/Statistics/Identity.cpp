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

#include "Statistics/Identity.h"
#include "InternalBenchmarker.h"
#include "Statistics/Manager.h"
#include "reportsystem.h"
#include "Alignment/Alignment.h"
#include "defines.h"
#include "utils.h"

namespace statistics {

    Identity::Identity(Alignment *parentAlignment) {
        StartTiming("Identity::Identity(Alignment *parentAlignment) ");
       
        alig = parentAlignment;
        refCounter = new int(1);
    }

    Identity::Identity(Alignment *parentAlignment, Identity *mold) {
        StartTiming("Identity::Identity(Alignment *parentAlignment, Identity *mold) ");

        alig = parentAlignment;

        identities = mold->identities;
        refCounter = mold->refCounter;
        (*refCounter)++;
    }

    void Identity::calculateSeqIdentity() {
        // Create a timer that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("void Identity::calculateSeqIdentity(void) ");

        int i, j, k, hit, dst;
        char indet;

        // Depending on alignment type, indetermination symbol will be one or other
        indet = (alig->getAlignmentType() & SequenceTypes::AA) ? 'X' : 'N';

        // Create identities matrix to store identities scores
        identities = new float *[alig->originalNumberOfSequences];

        // For each seq, compute its identity score against the others in the MSA
        for (i = 0; i < alig->originalNumberOfSequences; i++) {
            if (alig->saveSequences[i] == -1) continue;
            identities[i] = new float[alig->originalNumberOfSequences];

            // It's a symmetric matrix, copy values that have been already computed
            for (j = 0; j < i; j++) {
                if (alig->saveSequences[j] == -1) continue;
                identities[i][j] = identities[j][i];
            }
            identities[i][i] = 0;

            // Compute identity scores for the current sequence against the rest
            for (j = i + 1; j < alig->originalNumberOfSequences; j++) {
                if (alig->saveSequences[j] == -1) continue;
                for (k = 0, hit = 0, dst = 0; k < alig->numberOfResidues; k++) {
                    if (alig->saveResidues[k] == -1) continue;
                    // If one of the two positions is a valid residue,
                    // count it for the common length
                    if (((alig->sequences[i][k] != indet) && (alig->sequences[i][k] != '-')) ||
                        ((alig->sequences[j][k] != indet) && (alig->sequences[j][k] != '-'))) {
                        dst++;
                        // If both positions are the same, count a hit
                        if (alig->sequences[i][k] == alig->sequences[j][k])
                            hit++;
                    }
                }

                // Identity score between two sequences is the ratio of identical residues
                // by the total length (common and no-common residues) among them
                identities[i][j] = (float) hit / dst;
            }
        }
    }

    void Identity::calculateRelaxedSeqIdentity() {
        // Create a timer that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("void Identity::calculateRelaxedSeqIdentity(void) ");
        // Raw approximation of sequence identity computation designed for reducing
        // comparisons for huge alignemnts

        int i, j, k, hit;

        float inverseResidNumber = 1.F / alig->originalNumberOfResidues;

        // Create identities matrix to store identities scores
        identities = new float *[alig->originalNumberOfSequences];

        // For each seq, compute its identity score against the others in the MSA
        for (i = 0; i < alig->originalNumberOfSequences; i++) {
            if (alig->saveSequences[i] == -1) continue;
            identities[i] = new float[alig->originalNumberOfSequences];

            // It's a symmetric matrix, copy values that have been already computed
            for (j = 0; j < i; j++) {
                if (alig->saveSequences[j] == -1) continue;
                identities[i][j] = identities[j][i];
            }
            identities[i][i] = 0;

            // Compute identity score between the selected sequence and the others
            for (j = i + 1; j < alig->originalNumberOfSequences; j++) {
                if (alig->saveSequences[j] == -1) continue;
                for (k = 0, hit = 0; k < alig->originalNumberOfResidues; k++) {
                    if (alig->saveResidues[k] == -1) continue;
                    // If both positions are the same, count a hit
                    if (alig->sequences[i][k] == alig->sequences[j][k])
                        hit++;
                }
                // Raw identity score is computed as the ratio of identical residues between
                // alignment length 
                identities[i][j] = hit * inverseResidNumber;
            }
        }
    }

    Identity::~Identity() {
        if (refCounter == nullptr || --(*refCounter) < 1) {
            if (identities != nullptr) {
                for (int i = 0; i < alig->numberOfSequences; i++)
                    delete[] identities[i];
                delete[] identities;
            }
            delete refCounter;
            refCounter = nullptr;
        }
    }
}