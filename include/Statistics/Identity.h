/* *****************************************************************************

    trimAl v2.0: a tool for automated alignment trimming in large-scale
                 phylogenetics analyses.

    readAl v2.0: a tool for automated alignment conversion among different
                 formats.

    2022-2023
        Larralde M.             (martin.larralde@embl.de)

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

#ifndef STATISTICS_IDENTITY_H
#define STATISTICS_IDENTITY_H

class Alignment;

namespace statistics {

    /**
     * \brief Class to handle the calculation relative to identity.\n
     * Computation was previously done by Cleaner directly on Alignment 
     * objects, but was moved to a dedicated class for consistency with
     * other statistics.
    */
    class Identity {
    public:
        Identity(Alignment *parentAlignment, Identity *mold);

        Alignment *alig;

        /** \brief Identity between sequences */
        float *identities = nullptr;

        /** \brief Counter of how many other instances share the same information */
        int * refCounter;

    public:
        /** \brief Constructor without any parameters */
        explicit Identity(Alignment* parentAlignment);

        /** \brief Destructor */
        virtual ~Identity();

        /**
        \brief Method to calculate identities between the sequences from the alignment.\n
        */
        virtual void calculateSeqIdentity();
    };
}

#endif