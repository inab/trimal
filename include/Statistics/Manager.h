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

#ifndef TRIMAL_STATISTICSMANAGER_H
#define TRIMAL_STATISTICSMANAGER_H
#include "similarityMatrix.h"

// Forward declarations
class Alignment;



/// Namespace containing all classes related to statistics handling
namespace statistics {
// Forward declarations
    class Gaps;

    class Similarity;

    class Consistency;

    /**
     * \brief Class to handle the interaction with Alignment and statistics objects.\n
     * It serves as a wrapper or intermediate between the alignment and each specific stat.\n
     * It also encapsulates the similarityMatrix.
     */
    class Manager {
    public:

        /**
         * \brief Gaps submodule
         * */
        Gaps *gaps = nullptr;
        /**
         * \brief Similarity submodule
         * */
        Similarity *similarity = nullptr;
        /**
         * \brief Consistency submodule
         * */
        Consistency *consistency = nullptr;

        /**
         * \brief SimilarityMatrix used by Similarity
         * */
        similarityMatrix *_similarityMatrix = nullptr;

        /**
         * \brief Gap window
         * */
        int ghWindow;
        /**
         * \brief Similarity window
         * */
        int shWindow;

        /**
         * \brief Method to set a similarity matrix
         * */
        bool setSimilarityMatrix(similarityMatrix *sm);

        /**
         * \brief Method to handle gap stat calculation\n
         * It checks if the #gaps submodule has been created, otherwise, creates it
         * */
        bool calculateGapStats();

        /**
         * \brief Wrapper to Statistics::Gaps::printGapsColumns()\n
         * It calls to calculateGapStats() to make sure the information is available before reporting the requested values
         * */
        void printStatisticsGapsColumns();

        /**
         * \brief Wrapper to Statistics::Gaps::printGapsAcl()\n
         * It calls to calculateGapStats() to make sure the information is available before reporting the requested values
         * */
        void printStatisticsGapsTotal();

        /**
         * \brief Method to handle similarity stat calculation\n
         * It checks if the #similarity submodule has been created, otherwise, creates it
         * */
        bool calculateConservationStats();

        /**
         * \brief Wrapper to Statistics::Similarity::printConservationAcl()\n
         * It calls to calculateConservationStats() to make sure the information is available before reporting the requested values
         * */
        void printStatisticsConservationColumns();

        /**
         * \brief Wrapper to Statistics::Similarity::printGapsTotal()\n
         * It calls to calculateConservationStats() to make sure the information is available before reporting the requested values
         * */
        void printStatisticsConservationTotal();

        /**
         * \brief Method to print the vector containing the keep/reject (Alignment::saveResidues) values of the associated Alignment\n
         * If \b -1 , the residue would be rejected, otherwise, it prints the position on the alignment of the kept residue
         * */
        void printCorrespondence();

    private:
        /// Making Alignment a friend class allows us to have private constructor and destructors.\n
        /// This gives us more control over the statistics classes, hinting that they shouldn't be created outside the scope of a newAlignmet.
        friend class ::Alignment;

        Alignment *alig;

        explicit Manager(Alignment *parent);

        Manager(Alignment *parent, Manager *mold);

        ~Manager();
    };

}
#endif //TRIMAL_STATISTICSMANAGER_H
