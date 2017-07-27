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

#include <stdlib.h>

#include <string>
#include <iostream>

#include "alignment.h"
#include "newAlignment.h"

class compareFiles {

  public:
/**
 \brief Print the consistency value for each column from the selected alignment.
 \param numAminos Number of aminos in the alignment
 \param compareVect Vector containing the consistency value for each column.
 */
    static void printStatisticsFileColumns(int numAminos, float * compareVect);
/**
 \brief Print the accumulated consistency value from the selected alignment.
 \param numAminos Number of aminos in the alignment
 \param compareVect Vector containing the consistency value for each column.
 */
    static void printStatisticsFileAcl(int numAminos, float * compareVect);
/**
 \brief Applies a new window to the alignment.
 \param columns Number of columns present in the alignment.
 \param halfWindow Half size of window to apply.
 \param[in,out] columnsValue Vector containing the values that should be averaged and, at the same time, vector that will contain the new averaged values.
 \return \b True
 */
    static bool applyWindow(int columns, int halfWindow, float * columnsValue);
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
    static int algorithm(alignment **vectAlignments, char **fileNames, float *columnsValue, int numAlignments, bool verbosity);
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
    static int algorithm(newAlignment **vectAlignments, char **fileNames, float *columnsValue, int numAlignments, bool verbosity);
/**
 \brief Method to obtain the consistency values vector for a given alignment against a set of alignments with the same sequences.
 \param vectAlignments Alignment vector to compare against the selected alignment
 \param numAlignments Number of alignments to compare
 \param selected Alignment to compare against the set of alignments.
 \param[out] columnsValue Vector to fill with the consistency values.
 \return Wether or not the method went ok.
 */
    static bool forceComparison(alignment **vectAlignments, int numAlignments, alignment *selected, float *columnsValue);
/**
 \brief Method to obtain the consistency values vector for a given alignment against a set of alignments with the same sequences.
 \param vectAlignments Alignment vector to compare against the selected alignment
 \param numAlignments Number of alignments to compare
 \param selected Alignment to compare against the set of alignments.
 \param[out] columnsValue Vector to fill with the consistency values.
 \return Wether or not the method went ok.
 */
    static bool forceComparison(newAlignment **vectAlignments, int numAlignments, newAlignment *selected, float *columnsValue);

};
#endif
