//
// Created by bioinfo on 5/06/17.
//

#include <StatisticsManager.h>
#include <newAlignment.h>
#include <Cleaner.h>

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This method compute the similarity values from the input newAlignment if it
 * has not been computed before */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
bool StatisticsManager::calculateConservationStats(void) {

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* It the gaps statistics object has not been created
     * we create it */
    if(calculateGapStats() != true)
        return false;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* It the similarity statistics object has not been
     * created we create it */
    if(_alignment->scons == NULL)
        return false;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Ask for the similarity matrix */
    if(_alignment->scons -> isSimMatrixDef() != true)
        return false;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Compute the similarity statistics from the input
     * newAlignment */
    if(!_alignment->scons -> calculateVectors(_alignment->sequences, _alignment->sgaps->getGapsWindow()))
        return false;

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Ask to know if it is necessary to apply any window
     * method. If it's necessary, we apply it */
    if(_alignment->scons->isDefinedWindow())
        return true;
    else
        return _alignment->scons->applyWindow(_alignment->shWindow);
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Print the similarity value for each column from the newAlignment */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void StatisticsManager::printStatisticsConservationColumns(void) {

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Check if the similarity statistics object has been
     * created */
    if(calculateConservationStats())
        /* then prints the information */
        _alignment->scons -> printConservationColumns();
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Print the accumulative similarity distribution values from the newAlignment */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void StatisticsManager::printStatisticsConservationTotal(void) {

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Check if the similarity statistics object has been
     * created */
    if(calculateConservationStats())
        /* then prints the information */
        _alignment->scons -> printConservationAcl();
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Set the similarity matrix. This matrix is necessary for some methods in
 * the program */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
bool StatisticsManager::setSimilarityMatrix(similarityMatrix *sm) {

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* If scons object is not created, we create them */
    if(_alignment->scons == NULL)
        _alignment->scons = new statisticsConservation(_alignment->sequences,
                                                       _alignment->sequenNumber,
                                                       _alignment->residNumber,
                                                       _alignment->dataType);
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Associate the matrix to the similarity statistics
     * object */
    if(!_alignment->scons -> setSimilarityMatrix(sm))
        return false;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* If it's OK, we return true */
    return true;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Print the gaps value for each column in the newAlignment */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void StatisticsManager::printStatisticsGapsColumns(void) {

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Check if there is computed the gaps statistics */
    if(calculateGapStats())
        /* then prints the information */
        _alignment->sgaps -> printGapsColumns();
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Print the acumulative gaps distribution value from the input newAlignment */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void StatisticsManager::printStatisticsGapsTotal(void) {

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Check it there is computed the gaps statistics */
    if(calculateGapStats())
        /* then prints the information */
        _alignment->sgaps -> printGapsAcl();
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This method prints the correspondece between the columns in the original
 * and in the trimmed alignment */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void StatisticsManager::printCorrespondence(void) {
    int i;

    cout << "#ColumnsMap\t";
    /* Print the saveResidues relathionship */
    for(i = 0; i < _alignment->residNumber - 1; i++)
        cout << _alignment->saveResidues[i] << ", ";
    cout << _alignment->saveResidues[i] << endl;

}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This function computes the gaps statistics for the input newAlignment. */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
bool StatisticsManager::calculateGapStats(void) {

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* If newAlignment matrix is not created, return false */
    if(_alignment->sequences == NULL)
        return false;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* If sgaps object is not created, we create them
       and calculate the statistics */
    if(_alignment->sgaps == NULL) {
        _alignment->sgaps = new statisticsGaps(_alignment->sequences,
                                               _alignment->sequenNumber,
                                               _alignment->residNumber,
                                               _alignment->dataType);
        _alignment->sgaps -> applyWindow(_alignment->ghWindow);
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    return true;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

StatisticsManager::StatisticsManager(newAlignment *parent) {
    _alignment = parent;
}

