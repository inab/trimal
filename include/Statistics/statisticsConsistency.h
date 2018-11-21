/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
   ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

    trimAl v1.4: a tool for automated alignment trimming in large-scale
                 phylogenetics analyses.

    2009-2015 Capella-Gutierrez S. and Gabaldon, T.

              [scapella, tgabaldon]@crg.es

    This file is part of trimAl

    trimAl is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, the last available version.

    trimAl is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with trimAl. If not, see <http://www.gnu.org/licenses/>.

***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

#ifndef COMPAREFILES_H
#define COMPAREFILES_H


#include "../../include/FormatHandling/FormatManager.h"
#include "../../include/trimalManager.h"

#include <iostream>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iomanip>


// Forward declaration
class newAlignment;

/// \brief Helper class to handle the comparison between alignments. It also handles Consistency staticstics.
class statisticsConsistency {

  public:

    void perform(char *comparesetFilePath, FormatManager &FormatManager, trimAlManager &manager, char *forceFile);

/**
 \brief Print the consistency value for each column from the selected alignment.
 \param _alignment Alignment used to obtain the accumulated consistency value
 \param compareVect Vector containing the consistency value for each column.
 */
    static void printStatisticsFileColumns(newAlignment& _alignment, float * compareVect);
/**
 \brief Print the accumulated consistency value from the selected alignment.
 \param _alignment Alignment used to obtain the accumulated consistency value
 \param compareVect Vector containing the consistency value for each column.
 */
    static void printStatisticsFileAcl(newAlignment& _alignment, float * compareVect);
/**
 \brief Applies a new window to the alignment.
 \param _halfWindow Half size of window to apply.
 \return \b True
 */
    bool applyWindow(int _halfWindow);
/**
 \brief Method to compare a set of alignments to select the most consistent one respect the others.\n 
 To compute the consistency values we use the proportion of residue pairs per columns in the alignments to compare.
 \param vectAlignments Alignment vector to compare and select the most consistent.
 \param fileNames Vector containing all the filenames. Useful only if verbosity==True.
 \param[out] columnsValue Consistency values of selected alignment.
 \param numAlignments Number of alignments to compare.
 \param verbosity Wether or not report by printing some results.
 \return \b -1 if there was any error.\n <b> Alignment index </b> of the selected algorithm otherwise.
 */
    static int compareAndChoose(newAlignment **vectAlignments, char **fileNames, float *columnsValue, int numAlignments, bool verbosity);
/**
 \brief Method to obtain the consistency values vector for a given alignment against a set of alignments with the same sequences.
 \param vectAlignments Alignment vector to compare against the selected alignment
 \param numAlignments Number of alignments to compare
 \param selected Alignment to compare against the set of alignments.
 \param[out] columnsValue Vector to fill with the consistency values.
 \return Wether or not the method went ok.
 */
    static bool forceComparison(newAlignment **vectAlignments, int numAlignments, newAlignment *selected, float *columnsValue);

    statisticsConsistency(newAlignment *pAlignment, statisticsConsistency *pConsistency);
    statisticsConsistency();

    ~statisticsConsistency();

    float * getValues();

private:

    newAlignment * _alignment               = nullptr;
    newAlignment ** compareAlignmentsArray  = nullptr;

    std::ifstream compare;

    float *values          = nullptr;
    float *values_windowed = nullptr;

    int
        numFiles            = 0,
        i                   = 0,
        maxAminos           = 0,
        prevType            = -1,
        referFile           = -1,
        halfWindow          = -1,
        residues            = -1;

    int * refCounter;

    char
        c,
        * line              = nullptr,
        ** filesToCompare   = nullptr;

    std::string nline;


    bool appearErrors       = false;

    void delete_variables();

    bool isDefinedWindow();
};
#endif
