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

#ifndef STATISTICS_CONSERVATION2_H
#define STATISTICS_CONSERVATION2_H

#include "Gaps.h"
#include "similarityMatrix.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <math.h>

// Forward declaration
class Alignment;

namespace statistics {


    /**
     * \brief Class to handle the calculation relative to similarity.\n
     * This class is narrowly connected to similarityMatrix, as the latter contains the information
     *      to calculate the similarity of the alignment.
    */
    class Similarity {
    public:

        Similarity(Alignment *parentAlignment, Similarity *mold);

        Alignment * alig;

        /** \brief Half Window used on the calculation of the similarity */
        int halfWindow          = -1;

        /* Similarity vectors */
        /** \brief Raw similarity values */
        float *MDK                  = nullptr;
        /** \brief Windowed convervation values */
        float *MDK_Window           = nullptr;

        /** \brief Identity weight matrix between alignment rows */
        float **matrixIdentity      = nullptr;

        /** \brief Similarity matrix used to similarity calculations */
        similarityMatrix *simMatrix = nullptr;

        /** \brief Counter of how many other instances share the same information */
        int * refCounter;

    public:
        /** \brief Computes the matrix identity between alignment's columns. */
        void calculateMatrixIdentity();

        /** \brief Constructor without any parameters */
        explicit Similarity(Alignment * parentAlignment);

        /** \brief Destructor */
        ~Similarity();

        /**
            \brief Method to calculate the similarity values of a alignment matrix.
            \param cutByGap Wheter to cut by gap or not
        */
        bool calculateVectors(bool cutByGap = true);

        /**
            \brief Allows us compute the conservationWindow's values.
            \param halfW Half window size to apply.
            \return \b False if there is a previously computed vector for this window size or half window size is greater than 1/4 of the alignment length.
        */
        bool applyWindow(int halfW);

        /**
            \brief Returns if a windows size value has been defined or not.
            \return \b True if a windows size has been defined.\b False otherwise.
        */
        bool isDefinedWindow();

        /**
            \brief This methods returns a pointer to conservationWindow's vector
            \return Similarity window vector.
        */
        float *getMdkWindowedVector();

        /**
            \brief Stores a valid similarity matrix point to use.
            \param sm Similarity matrix pointer to associate.
            \return \b True if sm is valid, \b False if it's null
        */
        bool setSimilarityMatrix(similarityMatrix * sm);

        /**
            \brief Returns if a similarity matrix is being used or not.
            \return \b True if there is a similarity matrix set, \b False otherwise.
        */
        bool isSimMatrixDef();

        /**
            \brief Computes and selects the cut point values based on similarity's values.
            \todo This description seems a little vague.
            \param baseLine Percentage of columns desired.
            \param conservationPct Percentage of similarity desired.
        */
        double calcCutPoint(float baseLine, float conservationPct);

        /** \brief Prints the similarity's value for each alignment's column. */
        void printConservationColumns();

        /** \brief Computes and prints the accumulative statistics associated to the alignment. */
        void printConservationAcl();

    };

}

#endif
