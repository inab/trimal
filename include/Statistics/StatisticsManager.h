//
// Created by bioinfo on 5/06/17.
//

#ifndef TRIMAL_STATISTICSMANAGER_H
#define TRIMAL_STATISTICSMANAGER_H
#include "../../include/similarityMatrix.h"

// Forward declarations
class statisticsGaps;
class statisticsConservation;
class statisticsConsistency;

class newAlignment;

/// \brief Class to handle the interaction with statistics and statistics objects.
class StatisticsManager {
public:

    statisticsGaps * gaps                   = nullptr;
    statisticsConservation * conservation   = nullptr;
    statisticsConsistency * consistency     = nullptr;

    similarityMatrix * _similarityMatrix    = nullptr;

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

private:
    friend class newAlignment;

    newAlignment* _alignment;

    explicit StatisticsManager(newAlignment* parent);

    StatisticsManager(newAlignment* parent, StatisticsManager* mold);

    ~StatisticsManager();
};


#endif //TRIMAL_STATISTICSMANAGER_H
