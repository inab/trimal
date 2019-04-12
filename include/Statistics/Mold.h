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

#ifndef STATISTICS_MOLD_H
#define STATISTICS_MOLD_H

#include <sstream>
#include <iomanip>
#include <string>

// Forward declaration
class Alignment;

namespace statistics {

    /** \brief <b> Currently not in use </b>\n Class to serve as a mold to multiple stats.\n
     * This would help to maintain coherence and to ease the understanding of each of the stats. \n
     * As the number of stats grows, this standarization will become more valuable.\n
    */
    class Mold {
    public:

        Alignment *alig;

        int
                residues = -1,
                sequences = -1,
                halfWindowApplied = -1;

        /* Similarity vectors */
        float *values = nullptr;
        float *valuesWindow = nullptr;

        int *refCounter = nullptr;

        const static std::string statName;


    public:
        Mold(Alignment *parentAlignment, Mold *mold);

        explicit Mold(Alignment *parentAlignment);

        ~Mold();

        virtual bool calculate() { return false; };

        virtual bool applyWindow(int halfW);

        virtual bool isDefinedWindow();

        virtual float *getVector();

        virtual void printByColumn(bool calculateRelative = false); // NOLINT ->
        // Default values for virtual methods is 'prohibited'

        virtual void printAccumulative(bool calculateRelative = false); // NOLINT ->
        // Default values for virtual methods is 'prohibited'

    };
}
#endif
