/* *****************************************************************************

    trimAl v2.0: a tool for automated alignment trimming in large-scale
                 phylogenetics analyses.

    readAl v2.0: a tool for automated alignment conversion among different
                 formats.

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

#include "Statistics/Similarity.h"
#include "Statistics/Consistency.h"
#include "Statistics/Identity.h"
#include "Statistics/Manager.h"
#include "Statistics/Overlap.h"
#include "InternalBenchmarker.h"

#if defined(HAVE_AVX2)
#include "Platform/x86/AVX2.h"
#endif
#if defined(HAVE_NEON)
#include "Platform/Arm/NEON.h"
#endif
#if defined(HAVE_SSE2)
#include "Platform/x86/SSE2.h"
#endif

namespace statistics {
    bool Manager::calculateConservationStats() {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("bool Manager::calculateConservationStats(void) ");


        // It the gaps statistics object has not been created
        // we create it
        if (!calculateGapStats())
            return false;

        // It the similarity statistics object has not been
        // created we create it
        if (similarity == nullptr) {
            similarity = new Similarity(alig);
            similarity->setSimilarityMatrix(_similarityMatrix);
            similarity->applyWindow(shWindow);
        }

        // Ask for the similarity matrix
        if (!similarity->isSimMatrixDef())
            return false;

        // Compute the similarity statistics from the input
        // Alignment
        if (!similarity->calculateVectors(true))
            return false;

        // Ask to know if it is necessary to apply any window
        // method. If it's necessary, we apply it
        if (alig->Statistics->similarity->isDefinedWindow())
            return true;
        else
            return alig->Statistics->similarity->applyWindow(shWindow);
    }


    void Manager::printStatisticsConservationColumns() {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("void Manager::printStatisticsConservationColumns(void) ");

        // Check if the similarity statistics object has been
        // created */
        if (calculateConservationStats())
            /* then prints the information */
            alig->Statistics->similarity->printConservationColumns();

    }

    void Manager::printStatisticsConservationTotal() {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("void Manager::printStatisticsConservationTotal(void) ");


        // Check if the similarity statistics object has been
        // created
        if (calculateConservationStats())
            // then prints the information
            alig->Statistics->similarity->printConservationAcl();

    }

    bool Manager::setSimilarityMatrix(similarityMatrix *sm) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("bool Manager::setSimilarityMatrix(similarityMatrix *sm) ");

        _similarityMatrix = sm;

        // If scons object is not created, we create them
        if (alig->Statistics->similarity == nullptr)
#if defined(HAVE_AVX2)
            alig->Statistics->similarity = new AVX2Similarity(alig);
#elif defined(HAVE_SSE2)
            alig->Statistics->similarity = new SSE2Similarity(alig);
#elif defined(HAVE_NEON)
            alig->Statistics->similarity = new NEONSimilarity(alig);
#else
            alig->Statistics->similarity = new Similarity(alig);
#endif

        // Associate the matrix to the similarity statistics object
        // If it's OK, we return true
        return alig->Statistics->similarity->setSimilarityMatrix(sm);


    }

    void Manager::printStatisticsGapsColumns() {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("void Manager::printStatisticsGapsColumns(void) ");


        // Check if there is computed the gaps statistics
        if (calculateGapStats())
            // then prints the information
            alig->Statistics->gaps->printGapsColumns();

    }

    void Manager::printStatisticsGapsTotal() {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("void Manager::printStatisticsGapsTotal(void) ");


        // Check it there is computed the gaps statistics
        if (calculateGapStats())
            // then prints the information
            alig->Statistics->gaps->printGapsAcl();

    }

    void Manager::printCorrespondence() {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("void Manager::printCorrespondence(void) ");
        int i;

        std::cout << "#ColumnsMap\t";

        // Print the first residue
        for (i = 0; i < alig->originalNumberOfResidues - 1; i++) {
            if (alig->saveResidues[i] == -1) continue;
            std::cout << alig->saveResidues[i];
            break;
        }
        // Print the rest of the residues, using the separator ", "
        for (i++; i < alig->originalNumberOfResidues; i++) {
            if (alig->saveResidues[i] == -1) continue;
            std::cout << ", " << alig->saveResidues[i];
        }

    }

    bool Manager::calculateGapStats() {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("bool Manager::calculateGapStats(void) ");

        // If Alignment matrix is not created, return false
        if (alig->sequences == nullptr)
            return false;

        // If sgaps object is not created, we create them
        // and calculate the statistics
        if (gaps == nullptr) {
#if defined(HAVE_AVX2)
            gaps = new AVX2Gaps(alig);
#elif defined(HAVE_SSE2)
            gaps = new SSE2Gaps(alig);
#elif defined(HAVE_NEON)
            gaps = new NEONGaps(alig);
#else
            gaps = new Gaps(alig);
#endif
            gaps->CalculateVectors();
        }
        return gaps->applyWindow(ghWindow);
    }

    bool Manager::calculateSeqIdentity() {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("bool Manager::calculateSeqIdentity(void) ");

        // If Alignment matrix is not created, return false
        if (alig->sequences == nullptr)
            return false;

        // If sgaps object is not created, we create them
        // and calculate the statistics
        if (identity == nullptr) {
            identity = new Identity(alig);
            identity->calculateSeqIdentity();
        }
        return true;
    }

    bool Manager::calculateSeqOverlap() {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("bool Manager::calculateSeqOverlap(void) ");

        // If Alignment matrix is not created, return false
        if (alig->sequences == nullptr)
            return false;

        // If sgaps object is not created, we create them
        // and calculate the statistics
        if (overlap == nullptr) {
#if defined(HAVE_AVX2)
            overlap = new AVX2Overlap(alig);
#elif defined(HAVE_SSE2)
            overlap = new SSE2Overlap(alig);
#elif defined(HAVE_NEON)
            overlap = new NEONOverlap(alig);
#else
            overlap = new Overlap(alig);
#endif
            overlap->calculateSeqOverlap();
        }
        return true;
    }

    bool Manager::calculateSpuriousVector(float overlapColumn, float *spuriousVector) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("bool Manager::calculateSpuriousVector(float, float *) ");

        // If Alignment matrix is not created, return false
        if (alig->sequences == nullptr)
            return false;

        // If sgaps object is not created, we create them
        // and calculate the statistics
        if (overlap == nullptr) {
#if defined(HAVE_AVX2)
            overlap = new AVX2Overlap(alig);
#elif defined(HAVE_SSE2)
            overlap = new SSE2Overlap(alig);
#elif defined(HAVE_NEON)
            overlap = new NEONOverlap(alig);
#else
            overlap = new Overlap(alig);
#endif
        }

        return overlap->calculateSpuriousVector(overlapColumn, spuriousVector);
    }

    Manager::Manager(Alignment *parent) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("Manager::Manager(Alignment *parent) ");
        alig = parent;

        ghWindow = 0;
        shWindow = 0;
    }

    Manager::Manager(Alignment *parent, Manager *mold) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("Manager::Manager(Alignment *parent, Manager *mold) ");
        alig = parent;

        _similarityMatrix = mold->_similarityMatrix;

        ghWindow = mold->ghWindow;
        shWindow = mold->shWindow;

        if (mold->similarity)
#if defined(HAVE_AVX2)
            similarity = new AVX2Similarity(parent, mold->similarity);
#elif defined(HAVE_SSE2)
            similarity = new SSE2Similarity(parent, mold->similarity);
#elif defined(HAVE_NEON)
            similarity = new NEONSimilarity(parent, mold->similarity);
#else
            similarity = new Similarity(parent, mold->similarity);
#endif

        if (mold->consistency)
            consistency = new Consistency(parent, mold->consistency);

        if (mold->gaps)
#if defined(HAVE_AVX2)
            gaps = new AVX2Gaps(parent, mold->gaps);
#elif defined(HAVE_SSE2)
            gaps = new SSE2Gaps(parent, mold->gaps);
#elif defined(HAVE_NEON)
            gaps = new NEONGaps(parent, mold->gaps);
#else
            gaps = new Gaps(parent, mold->gaps);
#endif

        if (mold->identity)
#if defined(HAVE_AVX2)
            identity = new AVX2Identity(parent, mold->identity);
#elif defined(HAVE_SSE2)
            identity = new SSE2Identity(parent, mold->identity);
#elif defined(HAVE_NEON)
            identity = new NEONIdentity(parent, mold->identity);
#else
            identity = new Identity(parent, mold->identity);
#endif

        if (mold->overlap)
#if defined(HAVE_AVX2)
        overlap = new AVX2Overlap(parent, mold->overlap);
#elif defined(HAVE_SSE2)
        overlap = new SSE2Overlap(parent, mold->overlap);
#elif defined(HAVE_NEON)
        overlap = new NEONOverlap(parent, mold->overlap);
#else
        overlap = new Overlap(parent, mold->overlap);
#endif

    }

    Manager::~Manager() {
        delete gaps;
        gaps = nullptr;

        delete similarity;
        similarity = nullptr;

        delete consistency;
        consistency = nullptr;

        delete identity;
        identity = nullptr;

        delete overlap;
        overlap = nullptr;
    }
}
