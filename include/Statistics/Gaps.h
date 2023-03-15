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
#ifndef STATISTICSGAPS_H
#define STATISTICSGAPS_H

#include <iostream>
#include <iomanip>
#include <sstream>

// Forward declaration
class Alignment;

namespace statistics {

    /// \brief Class to handle gaps statistics
    class Gaps {
    public:

        Gaps(Alignment *pAlignment, Gaps *pGaps);

        Alignment *alig;

        int maxGaps = -1;
        int halfWindow = -1;
        int totalGaps = -1;

        int *gapsInColumn = nullptr;
        int *numColumnsWithGaps = nullptr;
        int *gapsWindow = nullptr;
        int *refCounter;

    public:

        // Class constructor without parameters.
        explicit Gaps(Alignment *parent);

        // Class destroyer.
        ~Gaps();

        // Methods allows us compute the gapWindows' values.
        bool applyWindow(int);

        // This methods returns a gaps' vector reference.
        int *getGapsWindow();

        // Allows compute and select the cut point value.
        double calcCutPoint(float, float);

        // Automatic method to find a cut point value
        // using the first and the second slopes.
        int calcCutPointMixSlope();

        // Automatic method to compute a cut point value
        // using the second slope approach.
        int calcCutPoint2ndSlope();

        // This methods print the gaps' percentage
        // of each column in the alignment.
        void printGapsColumns() const;

        // This methods prints the statistics
        // for the alignment relates to gaps.
        void printGapsAcl() const;

        bool isDefinedWindow();

        void CalculateVectors();
    };
}
#endif
