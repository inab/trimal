//
// Created by bioinfo on 5/06/17.
//

#ifndef TRIMAL_STATISTICSMANAGER_H
#define TRIMAL_STATISTICSMANAGER_H
#include "../include/similarityMatrix.h"
#include "../include/Statistics/statisticsGaps.h"
#include "../include/Statistics/statisticsConservation2.h"

class newAlignment;
/// \brief Class to handle the interaction with statistics and statistics objects.
class StatisticsManager {
public:

    statisticsGaps * gaps = NULL;
    statisticsConservation2 * conservation = NULL;

    int ghWindow;
    int shWindow;

    bool calculateConservationStats(void);

    bool setSimilarityMatrix(similarityMatrix *sm);

    bool calculateGapStats(void);

    void printStatisticsGapsColumns(void);

    void printStatisticsGapsTotal(void);

    void printStatisticsConservationColumns(void);

    void printStatisticsConservationTotal(void);

    void printCorrespondence(void);

    void saveStatistics(similarityMatrix *sm);

    void saveStatistics(similarityMatrix *, int, int);

private:
    friend class newAlignment;

    newAlignment* _alignment;

    StatisticsManager(newAlignment* parent);

    StatisticsManager(newAlignment* parent, StatisticsManager* mold);

    ~StatisticsManager();
};


#endif //TRIMAL_STATISTICSMANAGER_H
