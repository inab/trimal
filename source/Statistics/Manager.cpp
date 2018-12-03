//
// Created by bioinfo on 5/06/17.
//

#include "Statistics/Similarity.h"
#include "Statistics/Consistency.h"
#include "InternalBenchmarker.h"
#include "Statistics/Manager.h"
#include "Statistics/Gaps.h"
#include "Alignment/Alignment.h"

namespace statistics {
    bool Manager::calculateConservationStats(void) {
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


    void Manager::printStatisticsConservationColumns(void) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("void Manager::printStatisticsConservationColumns(void) ");

        // Check if the similarity statistics object has been
        // created */
        if (calculateConservationStats())
            /* then prints the information */
            alig->Statistics->similarity->printConservationColumns();

    }

    void Manager::printStatisticsConservationTotal(void) {
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
            alig->Statistics->similarity = new Similarity(alig);

        // Associate the matrix to the similarity statistics object
        // If it's OK, we return true
        return alig->Statistics->similarity->setSimilarityMatrix(sm);


    }

    void Manager::printStatisticsGapsColumns(void) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("void Manager::printStatisticsGapsColumns(void) ");


        // Check if there is computed the gaps statistics
        if (calculateGapStats())
            // then prints the information
            alig->Statistics->gaps->printGapsColumns();

    }

    void Manager::printStatisticsGapsTotal(void) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("void Manager::printStatisticsGapsTotal(void) ");


        // Check it there is computed the gaps statistics
        if (calculateGapStats())
            // then prints the information
            alig->Statistics->gaps->printGapsAcl();

    }

    void Manager::printCorrespondence(void) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("void Manager::printCorrespondence(void) ");
        int i;

        std::cout << "#ColumnsMap\t";
        // Print the saveResidues relathionship
        for (i = 0; i < alig->originalNumberOfResidues - 1; i++)
            std::cout << alig->saveResidues[i] << ", ";
        std::cout << alig->saveResidues[i] << std::endl;

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
            gaps = new Gaps(alig);
        }
        gaps->CalculateVectors();
        return gaps->applyWindow(ghWindow);
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
            similarity = new Similarity(parent, mold->similarity);

        if (mold->consistency)
            consistency = new Consistency(parent, mold->consistency);

        if (mold->gaps)
            gaps = new Gaps(parent, mold->gaps);
    }

    Manager::~Manager() {
        delete gaps;
        gaps = nullptr;

        delete similarity;
        similarity = nullptr;

        delete consistency;
        consistency = nullptr;
    }
}