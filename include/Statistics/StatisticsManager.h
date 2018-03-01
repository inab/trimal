//
// Created by bioinfo on 5/06/17.
//

#ifndef TRIMAL_STATISTICSMANAGER_H
#define TRIMAL_STATISTICSMANAGER_H
#include "similarityMatrix.h"
//#include "statisticsGaps.h"
//#include "statisticsConservation.h"
//#include "statisticsConsistency.h"

class statisticsGaps;
class statisticsConservation2;
class statisticsConsistency;

class newAlignment;
/// \brief Class to handle the interaction with statistics and statistics objects.
class StatisticsManager {
public:

    statisticsGaps * gaps                   = nullptr;
    statisticsConservation2 * conservation  = nullptr;
    statisticsConsistency * consistency     = nullptr;

    int ghWindow;
    int shWindow;

    bool calculateConservationStats();

    bool setSimilarityMatrix(similarityMatrix *sm);

    bool calculateGapStats();

    void printStatisticsGapsColumns();

    void printStatisticsGapsTotal();

    void printStatisticsConservationColumns();

    void printStatisticsConservationTotal();

    void printCorrespondence();

//    void saveStatistics(similarityMatrix *sm);
//
//    void saveStatistics(similarityMatrix *, int, int);

private:
    friend class newAlignment;

    newAlignment* _alignment;

    StatisticsManager(newAlignment* parent);

    StatisticsManager(newAlignment* parent, StatisticsManager* mold);

    ~StatisticsManager();
};


#endif //TRIMAL_STATISTICSMANAGER_H
