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

#ifndef TRIMAL_ENTROPY_H
#define TRIMAL_ENTROPY_H

#include "Gaps.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <math.h>

// Forward declaration
class Alignment;

namespace statistics {


    /**
     * \brief Class to handle the calculation relative to Entropy.\n
     * This class is narrowly connected to EntropyMatrix, as the latter contains the information
     *      to calculate the Entropy of the alignment.
    */
    class Entropy {
    public:

        Entropy(Alignment *parentAlignment, Entropy *mold);

        Alignment * alig;

        /** \brief Half Window used on the calculation of the Entropy */
        int halfWindow          = -1;

        /* Entropy vectors */
        /** \brief Raw Entropy values */
        float *EntropyValues                  = nullptr;
        /** \brief Windowed entropy values */
        float *EntropyValues_Window           = nullptr;

        /** \brief Counter of how many other instances share the same information */
        int * refCounter;

    public:

        /** \brief Constructor without any parameters */
        explicit Entropy(Alignment * parentAlignment);

        /** \brief Destructor */
        ~Entropy();

        /**
            \brief Method to calculate the Entropy values of a alignment matrix.
        */
        bool calculateVectors(bool cutByGap);

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
            \return Entropy window vector.
        */
        float *getValues();

        /**
            \brief Computes and selects the cut point values based on Entropy's values.
            \todo This description seems a little vague.
            \param baseLine Percentage of columns desired.
            \param conservationPct Percentage of Entropy desired.
        */
        double calcCutPoint(float baseLine, float conservationPct);

        /** \brief Prints the Entropy's value for each alignment's column. */
        void printConservationColumns();

        /** \brief Computes and prints the accumulative statistics associated to the alignment. */
        void printConservationAcl();

    };

}

#endif //TRIMAL_ENTROPY_H