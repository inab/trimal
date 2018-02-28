#include <TimerFactory.h>
//
// Created by bioinfo on 5/06/17.
//

#include "../include/Statistics/StatisticsManager.h"
#include "../include/newAlignment.h"
#include "../include/Statistics/statisticsGaps.h"
#include "../include/Statistics/statisticsConservation.h"
#include "../include/Statistics/statisticsConsistency.h"


bool StatisticsManager::calculateConservationStats(void) {
	 // Create a timer that will report times upon its destruction
	 //	which means the end of the current scope.
	StartTiming("bool StatisticsManager::calculateConservationStats(void) ");


    // It the gaps statistics object has not been created
    // we create it 
    if (calculateGapStats() != true)
        return false;

    // It the similarity statistics object has not been
    // created we create it
    if (_alignment->Statistics->conservation == NULL)
        return false;

    // Ask for the similarity matrix
    if (_alignment->Statistics->conservation->isSimMatrixDef() != true)
        return false;

    // Compute the similarity statistics from the input
    // newAlignment
    if (!_alignment->Statistics->conservation->calculateVectors(_alignment->Statistics->gaps->getGapsWindow()))
        return false;

    // Ask to know if it is necessary to apply any window
    // method. If it's necessary, we apply it
    if (_alignment->Statistics->conservation->isDefinedWindow())
        return true;
    else
        return _alignment->Statistics->conservation->applyWindow(shWindow);
}


void StatisticsManager::printStatisticsConservationColumns(void) {
	 // Create a timer that will report times upon its destruction
	 //	which means the end of the current scope.
	StartTiming("void StatisticsManager::printStatisticsConservationColumns(void) ");

    // Check if the similarity statistics object has been
    // created */
    if (calculateConservationStats())
        /* then prints the information */
        _alignment->Statistics->conservation->printConservationColumns();

}

void StatisticsManager::printStatisticsConservationTotal(void) {
	 // Create a timer that will report times upon its destruction
	 //	which means the end of the current scope.
	StartTiming("void StatisticsManager::printStatisticsConservationTotal(void) ");


    // Check if the similarity statistics object has been
    // created 
    if (calculateConservationStats())
        // then prints the information
        _alignment->Statistics->conservation->printConservationAcl();

}

bool StatisticsManager::setSimilarityMatrix(similarityMatrix *sm) {
	 // Create a timer that will report times upon its destruction
	 //	which means the end of the current scope.
	StartTiming("bool StatisticsManager::setSimilarityMatrix(similarityMatrix *sm) ");

    // If scons object is not created, we create them
    if (_alignment->Statistics->conservation == NULL)
        _alignment->Statistics->conservation = new statisticsConservation2(_alignment);

    // Associate the matrix to the similarity statistics
    // object
    if (!_alignment->Statistics->conservation->setSimilarityMatrix(sm))
        return false;

    // If it's OK, we return true
    return true;

}

void StatisticsManager::printStatisticsGapsColumns(void) {
	 // Create a timer that will report times upon its destruction
	 //	which means the end of the current scope.
	StartTiming("void StatisticsManager::printStatisticsGapsColumns(void) ");


    // Check if there is computed the gaps statistics
    if (calculateGapStats())
        // then prints the information
        _alignment->Statistics->gaps->printGapsColumns();

}

void StatisticsManager::printStatisticsGapsTotal(void) {
	 // Create a timer that will report times upon its destruction
	 //	which means the end of the current scope.
	StartTiming("void StatisticsManager::printStatisticsGapsTotal(void) ");


    // Check it there is computed the gaps statistics 
    if (calculateGapStats())
        // then prints the information
        _alignment->Statistics->gaps->printGapsAcl();

}

void StatisticsManager::printCorrespondence(void) {
	 // Create a timer that will report times upon its destruction
	 //	which means the end of the current scope.
	StartTiming("void StatisticsManager::printCorrespondence(void) ");
    int i;

    cout << "#ColumnsMap\t";
    // Print the saveResidues relathionship
    for (i = 0; i < _alignment->residNumber - 1; i++)
        cout << _alignment->saveResidues[i] << ", ";
    cout << _alignment->saveResidues[i] << endl;

}

bool StatisticsManager::calculateGapStats(void) {
	 // Create a timer that will report times upon its destruction
	 //	which means the end of the current scope.
	StartTiming("bool StatisticsManager::calculateGapStats(void) ");

    // If newAlignment matrix is not created, return false
    if (_alignment->sequences == NULL)
        return false;

    // If sgaps object is not created, we create them
    // and calculate the statistics
    if (_alignment->Statistics->gaps == NULL) {
        _alignment->Statistics->gaps = new statisticsGaps(_alignment);
        _alignment->Statistics->gaps->applyWindow(ghWindow);
    }

    return true;
}

StatisticsManager::StatisticsManager(newAlignment *parent) {
	 // Create a timer that will report times upon its destruction
	 //	which means the end of the current scope.
	StartTiming("StatisticsManager::StatisticsManager(newAlignment *parent) ");
    _alignment = parent;

    ghWindow = 0;
    shWindow = 0;
}

StatisticsManager::StatisticsManager(newAlignment *parent, StatisticsManager *mold) {
	 // Create a timer that will report times upon its destruction
	 //	which means the end of the current scope.
	StartTiming("StatisticsManager::StatisticsManager(newAlignment *parent, StatisticsManager *mold) ");
    _alignment = parent;

    ghWindow = mold->ghWindow;
    shWindow = mold->shWindow;
}

StatisticsManager::~StatisticsManager() {
    delete gaps;
    delete conservation;
}
