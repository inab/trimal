//
// Created by bioinfo on 5/06/17.
//

#ifndef TRIMAL_STATISTICSMANAGER_H
#define TRIMAL_STATISTICSMANAGER_H
#include <similarityMatrix.h>

class newAlignment;
class StatisticsManager {
public:
    /* Statistics calculation */

    /** \brief Basic conservation statistics calculation.
     * \return \e true if all is ok, \e false if there were errors (i.e. there is no similarity matrix defined
     * in conservation statistics).
     *
     * This method calculates conservation statistics with the previously defined similarity matrix.
     */
    bool calculateConservationStats(void);

    /** \brief Conservation statistics calculation.
     * \param sm similarity matrix used for the statistics calculation.
     * \return \e true if all is ok, \e false if there were errors.
     *
     * This method calculates conservation statistics using the \b sm similarity matrix.
     */
    bool setSimilarityMatrix(similarityMatrix *sm);

    /** \brief Gap statistics calculation.
     * \return \e true if all is ok, \e false if there were errors.
     *
     * This method calculates gap statistics without window calculation (half window value = 0).
     */
    bool calculateGapStats(void);

    /* Output Statistics */

    /** \brief Printing normal gap statistics method.
     *
     * This method prints gap statistics for each column.
     */
    void printStatisticsGapsColumns(void);

    /** \brief Printing accumulated gap statistics method.
     *
     * This method prints accumulated gap statistics.
     */
    void printStatisticsGapsTotal(void);

    /** \brief Printing conservation values method.
     *
     * This method prints conservation value for each column of the alignment.
     */
    void printStatisticsConservationColumns(void);

    /** \brief Printing accumulated conservation statistics method.
     * This method prints the accumulated number of columns for each conservation
     * value from the the alignment.
     */
    void printStatisticsConservationTotal(void);

    void printCorrespondence(void);

private:
    friend class newAlignment;

    newAlignment* _alignment;

    StatisticsManager(newAlignment* parent);
};


#endif //TRIMAL_STATISTICSMANAGER_H
