/* *****************************************************************************

    trimAl v2.0: a tool for automated alignment trimming in large-scale
                 phylogenetics analyses.

    readAl v2.0: a tool for automated alignment conversion among different
                 formats.

    2009-2019
        Fernandez-Rodriguez V.  (victor.fernandez@bsc.es)
        Capella-Gutierrez S.    (salvador.capella@bsc.es)
        Gabaldon, T.            (tgabaldon@crg.es)

    This file is part of trimAl/readAl.

    trimAl/readAl are free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, the last available version.

    trimAl/readAl are distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with trimAl/readAl. If not, see <http://www.gnu.org/licenses/>.

***************************************************************************** */

#ifndef COMPAREFILES_H
#define COMPAREFILES_H


#include "FormatHandling/FormatManager.h"
#include "trimalManager.h"

#include <iostream>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iomanip>

// Forward declaration
class Alignment;

namespace statistics {


    /**
     * \brief Class to calculate the consistency between several MSA containing the same sequences, differently aligned.\n
     * Using this statistics, the class is able to select the most consistent alignment between all alignments provided.\n
     * It is possible to forcefully select an alignment, but to calculate the statistics for latter use.\n
     * After selecting an alignment (most consistent or manually selected),
     * it is possible to use this statistic to trim the alignment, removing
     *      columns that are not consistent enough with the other alignments.
    */
    class Consistency {

    public:

        /**
         * \brief Method to compare a set of MSA,
         *      all containing the same sequences and residues.\n
         * The number of residues must be the same, but gaps are not taken into account.\n
         * This is due to the same sequence being aligned in different ways,
         *      which changes the gap patterns.
         * \param comparesetFilePath Path to the file containing paths for each alignment to compare.\n One per line
         * \param formatManager Format manager, to load and save the alignments.
         * \param manager trimAl manager, to store the choosen alignment in trimAlManager::origAlig
         * \param forceFile path to file to forcefully select. If nullptr, the most consistent alignment will be selected.
         *
         * \note The method does not return anything, as the alignment is stored on trimAlManager::origAlig
         *      It also stores the current instance of statisticsConsistency into the consistent alignment.
         *      This allows us to use this information to trim, or represent it on the HTML/SVG reports.
        */
        bool perform(char *comparesetFilePath, FormatHandling::FormatManager &formatManager, trimAlManager &manager, char *forceFile);

        /**
         * \brief Print the consistency value for each column from the selected alignment.
         * \param alig Alignment used to obtain the accumulated consistency value
         * \param compareVect Vector containing the consistency value for each column.
        */
        static void printStatisticsFileColumns(Alignment &alig, float *compareVect);

        /**
         * \brief Print the accumulated consistency value from the selected alignment.
         * \param alig Alignment used to obtain the accumulated consistency value
         * \param compareVect Vector containing the consistency value for each column.
        */
        static void printStatisticsFileAcl(Alignment &alig, float *compareVect);

        /**
         * \brief Applies a new window to the alignment.
         * \param halfW Half size of window to apply.
         * \return \b True if correct \b Exits if false
        */
        bool applyWindow(int halfW);

        /**
         * \brief Copy constructor
        */
        Consistency(Alignment *pAlignment, Consistency *pConsistency);

        /**
         * \brief Default Construtor
        */
        Consistency();

        /**
         * \brief Default Destructor
        */
        ~Consistency();

        /**
         * \brief Stat Getter \n
         * \return
         * If a window has been set and applied, it will return the windowed values
         * If a window has been set, but not applied, it wil apply it and return the windowed values.
         * If no window has been set, it will return the raw values.
        */
        float *getValues();

    private:

        /** \brief Original alignment for which the stat was calculated */
        Alignment *alig = nullptr;

        /** \brief Array of alignments to compare */
        Alignment **compareAlignmentsArray = nullptr;

        /** \brief Raw consistency values */
        float *values = nullptr;

        /** \brief Windowed consistency values */
        float *values_windowed = nullptr;

        int
        /** \brief Number of files to compare */
                numFiles = 0,
        /** \brief Temporary variable used on loops. */
                i = 0,
        /** \brief Maximum number of residues on the whole dataset */
                maxResidues = 0,
        /** \brief Variable to store the type of the alignment from the last alignment. */
                halfWindow = -1,
        /** \brief Number of residues of the selected alignment */
                residues = -1;

        /** \brief Counter of how many statisticsConsistency share the same MDK values. */
        int *refCounter;

        /** \brief Intermediate variable to keep track of the progress status */
        bool appearErrors = false;

        /**
         * \brief Method to check wether or not a window has been applied.
         * \returns \b True halfWindow is different than \b -1
        */
        bool isWindowDefined();

        /**
         * \brief Method to compare a set of alignments to select the most consistent one respect the others.\n
         * To compute the consistency values we use the proportion of residue pairs per columns in the alignments to compare.
         * \param vectAlignments Alignment vector to compare and select the most consistent.
         * \param fileNames Vector containing all the filenames. Useful only if verbosity==True.
         * \param[out] columnsValue Consistency values of selected alignment.
         * \param numAlignments Number of alignments to compare.
         * \param verbosity Wether or not report by printing some results.
         * \return \b -1 if there was any error.\n <b> Alignment index </b> of the selected algorithm otherwise.
        */
        static int compareAndChoose(Alignment **vectAlignments, char **fileNames, float *columnsValue, int numAlignments, bool verbosity);

        /**
         * \brief Method to obtain the consistency values vector for a given alignment against a set of alignments with the same sequences.
         * \param vectAlignments Alignment vector to compare against the selected alignment
         * \param numAlignments Number of alignments to compare
         * \param selected Alignment to compare against the set of alignments.
         * \param[out] columnsValue Vector to fill with the consistency values.
         * \return Wether or not the method went ok.
        */
        static bool forceComparison(Alignment **vectAlignments, int numAlignments, Alignment *selected, float *columnsValue);
    };
}
#endif
