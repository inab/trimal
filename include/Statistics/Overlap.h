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

#ifndef STATISTICS_OVERLAP_H
#define STATISTICS_OVERLAP_H

class Alignment;

namespace statistics {

    /**
     * \brief Class to handle the calculation relative to overlap.\n
     * Computation was previously done by Cleaner directly on Alignment 
     * objects, but was moved to a dedicated class for consistency with
     * other statistics.
    */
    class Overlap {
    public:
        Overlap(Alignment *parentAlignment, Overlap *mold);

        Alignment *alig;

        /** \brief Overlap between sequences */
        float **overlaps = nullptr;

        /** \brief Counter of how many other instances share the same information */
        int * refCounter;

    public:
        /** \brief Constructor without any parameters */
        explicit Overlap(Alignment* parentAlignment);

        /** \brief Destructor */
        virtual ~Overlap();

        /**
        \brief Method to calculate overlap between the sequences from the alignment.\n
        */
        virtual void calculateSeqOverlap();

        /**
        \brief Method to compute the overlap values.\n
        \param overlap
        Overlap threshold.
        \param[out] spuriousVector
        Pointer to the spuriousVector to fill.
        \return   <b> True </b> if the calculation went ok.\n
                    <b> False </b> otherwise.\n
        This should happen only if you pass a null pointer instead of a spuriousVector.
        */
        virtual bool calculateSpuriousVector(float overlapColumn, float *spuriousVector);
    };
}

#endif