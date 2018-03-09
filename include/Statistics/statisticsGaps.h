/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
   ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

    trimAl v1.4: a tool for automated alignment trimming in large-scale
                 phylogenetics analyses.

    2009-2015 Capella-Gutierrez S. and Gabaldon, T.
              [scapella, tgabaldon]@crg.es

    This file is part of trimAl.

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
#ifndef STATISTICSGAPS_H
#define STATISTICSGAPS_H

#include <iostream>
#include <iomanip>

using namespace std;

class newAlignment;

/// \brief Class to handle Gaps statistics
class statisticsGaps {
public:

    statisticsGaps(newAlignment *pAlignment, statisticsGaps *pGaps);

    newAlignment *_alignment;

    int residNumber = -1;
    int sequenNumber = -1;
    int maxGaps = -1;
    int halfWindow = -1;

    int *gapsInColumn = nullptr;
    int *numColumnsWithGaps = nullptr;
    int *aminosXInColumn = nullptr;
    int *gapsWindow = nullptr;
    int *refCounter;

public:

    // Class constructor without parameters.
    explicit statisticsGaps(newAlignment *parent);

    // Class destroyer.
    ~statisticsGaps();

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
    void printGapsColumns();

    // This methods prints the statistics
    // for the alignment relates to gaps.
    void printGapsAcl();

    bool isDefinedWindow();
};

#endif
