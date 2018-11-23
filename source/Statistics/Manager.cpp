//
// Created by bioinfo on 5/06/17.
//

#include "Statistics/Conservation.h"
#include "Statistics/Consistency.h"
#include "Statistics/Manager.h"
#include "Statistics/Gaps.h"
#include "Alignment.h"
#include "InternalBenchmarker.h"

namespace Statistics {
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
        if (conservation == nullptr) {
            conservation = new Conservation(_alignment);
            conservation->setSimilarityMatrix(_similarityMatrix);
            conservation->applyWindow(shWindow);
        }

        // Ask for the similarity matrix
        if (!conservation->isSimMatrixDef())
            return false;

        // Compute the similarity statistics from the input
        // Alignment
        if (!conservation->calculateVectors(true))
            return false;

        // Ask to know if it is necessary to apply any window
        // method. If it's necessary, we apply it
        if (_alignment->Statistics->conservation->isDefinedWindow())
            return true;
        else
            return _alignment->Statistics->conservation->applyWindow(shWindow);
    }


    void Manager::printStatisticsConservationColumns(void) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("void Manager::printStatisticsConservationColumns(void) ");

        // Check if the similarity statistics object has been
        // created */
        if (calculateConservationStats())
            /* then prints the information */
            _alignment->Statistics->conservation->printConservationColumns();

    }

    void Manager::printStatisticsConservationTotal(void) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("void Manager::printStatisticsConservationTotal(void) ");


        // Check if the similarity statistics object has been
        // created
        if (calculateConservationStats())
            // then prints the information
            _alignment->Statistics->conservation->printConservationAcl();

    }

    bool Manager::setSimilarityMatrix(similarityMatrix *sm) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("bool Manager::setSimilarityMatrix(similarityMatrix *sm) ");

        _similarityMatrix = sm;

        // If scons object is not created, we create them
        if (_alignment->Statistics->conservation == nullptr)
            _alignment->Statistics->conservation = new Conservation(_alignment);

        // Associate the matrix to the similarity statistics object
        // If it's OK, we return true
        return _alignment->Statistics->conservation->setSimilarityMatrix(sm);


    }

    void Manager::printStatisticsGapsColumns(void) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("void Manager::printStatisticsGapsColumns(void) ");


        // Check if there is computed the gaps statistics
        if (calculateGapStats())
            // then prints the information
            _alignment->Statistics->gaps->printGapsColumns();

    }

    void Manager::printStatisticsGapsTotal(void) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("void Manager::printStatisticsGapsTotal(void) ");


        // Check it there is computed the gaps statistics
        if (calculateGapStats())
            // then prints the information
            _alignment->Statistics->gaps->printGapsAcl();

    }

    void Manager::printCorrespondence(void) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("void Manager::printCorrespondence(void) ");
        int i;

        std::cout << "#ColumnsMap\t";
        // Print the saveResidues relathionship
        for (i = 0; i < _alignment->numberOfResidues - 1; i++)
            std::cout << _alignment->saveResidues[i] << ", ";
        std::cout << _alignment->saveResidues[i] << std::endl;

    }

    bool Manager::calculateGapStats() {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("bool Manager::calculateGapStats(void) ");

        // If Alignment matrix is not created, return false
        if (_alignment->sequences == nullptr)
            return false;

        // If sgaps object is not created, we create them
        // and calculate the statistics
        if (gaps == nullptr) {
            gaps = new Gaps(_alignment);
        }
        gaps->CalculateVectors();
        return gaps->applyWindow(ghWindow);
    }

    Manager::Manager(Alignment *parent) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("Manager::Manager(Alignment *parent) ");
        _alignment = parent;

        ghWindow = 0;
        shWindow = 0;
    }

    Manager::Manager(Alignment *parent, Manager *mold) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("Manager::Manager(Alignment *parent, Manager *mold) ");
        _alignment = parent;

        _similarityMatrix = mold->_similarityMatrix;

        ghWindow = mold->ghWindow;
        shWindow = mold->shWindow;

        if (mold->conservation)
            conservation = new Conservation(parent, mold->conservation);

        if (mold->consistency)
            consistency = new Consistency(parent, mold->consistency);

        if (mold->gaps)
            gaps = new Gaps(parent, mold->gaps);
    }

    Manager::~Manager() {
        delete gaps;
        gaps = nullptr;

        delete conservation;
        conservation = nullptr;

        delete consistency;
        consistency = nullptr;
    }
}