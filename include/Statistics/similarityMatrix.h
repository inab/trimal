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

#ifndef SIMILARITYMATRIX_H
#define SIMILARITYMATRIX_H

#include <iostream>
#include <cstring>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cctype>
#include <cmath>


namespace statistics {
/** 
 * \brief Class that contains information of similarity matrices.\n
 * These are used to calculate the similarity between residues on the same column.
 * <br>
 * Default matrices for AA, NT and DEG NT are provided, along with method for loading custom matrices.
 */
    class similarityMatrix {

        int *vhash;
        float **simMat;
        float **distMat;
        int numPositions;

    private:
        /**
         * \brief Method to allocate memory for the similiarity matrix
         * \param nPos Number of different possible residues in the alignment.\n
         * This are, on the default matrices:
         *     -# AA:      20 residues
         *     -# NT:       5 residues
         *     -# DEG NT:  15 residues
        */
        void memoryAllocation(int nPos);

        /**
         * \brief Method to deallocate memory allocated on the similarityMatrix::memoryAllocation method.\n
         * It makes use of the #numPositions to effectively remove the memory.
        */
        void memoryDeletion();

    public:
        /** \brief Constructor **/
        similarityMatrix();

        /** \brief Destructor **/
        ~similarityMatrix();

        /**
         * \brief Method to load a custom matrix
         * \param filename Path to file containing the matrix to load
         * \return \b True if loaded \n \b False if an error ocurred**/
        bool loadSimMatrix(char *filename);

        /**
         * \brief Method to load the default AA similarity matrix
        */
        void defaultAASimMatrix();

        /**
         * \brief Method to load the default NT similarity matrix
        */
        void defaultNTSimMatrix();

        /**
         * \brief Method to load the default DEG NT similarity matrix
        */
        void defaultNTDegeneratedSimMatrix();

        /**
         * \brief Method to load alternative similarity matrices also included on the suite.
         * Currently, only one type of alternative matrix is available: \n
         * \b matrix_code: 1 \b datatype SequenceTypes::AA
         * \param matrix_code ID of the matrix
         * \param datatype Numberical representation of the data type.
         *  See #SequenceTypes
        */
        void alternativeSimilarityMatrices(int matrix_code, int datatype);

        /**
         * \brief Method to get the similarity distance between two residues, A and B\n
         * Characters provided must be both uppercase, please, refer to utils::toUpper
         * \param a First residue to compare
         * \param b Second residue to compare
         * \return Distance between A and B based on the similarity matrix loaded.
        */
        float getDistance(char a, char b);

        /**
         * \brief Method to print the loaded matrix.
        */
        void printMatrix();
    };
}
#endif
