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

#include "Statistics/statisticsMold.h"
#include "newAlignment.h"
#include "reportsystem.h"
#include "TimerFactory.h"

statisticsMold::statisticsMold(newAlignment *parentAlignment, statisticsMold *mold) {

}

statisticsMold::statisticsMold(newAlignment *parentAlignment) {

}

statisticsMold::~statisticsMold() {
    if (--(*refCounter) == 0)
    {
        delete [] values;
    }
    delete [] valuesWindow;
    _alignment = nullptr;
}

bool statisticsMold::applyWindow(int _halfWindow) {
    // Calculate the values array if it has not been calculated previously
    if (values == nullptr)
        calculateVectors();

    // Check is the half window value passed is in the valid range
    if (_halfWindow > residues / 4) {
        debug.report(ErrorCode::SimilarityWindowTooBig);
        return false;
    }

    // If the current half window is the same as the last one, don't do anything
    if (halfWindowApplied == _halfWindow) return true;

    // Save the requested half window. This is useful when making a copy of the
    // alignment, as the window values are not valid anymore but don't want to
    // calculate them if not needed anymore
    halfWindowRequested = _halfWindow;

    // If the half window requested is 0 or a negative number
    // we simply delete the window values.
    if (_halfWindow < 1) {
        if (halfWindowApplied > 0)
            delete[] valuesWindow;

        valuesWindow = nullptr;
        return true;
    }

    // Initialize the values used in the calculation
    int i, j, window;

    // Initialize the values window array if it's null
    if (valuesWindow == nullptr)
        valuesWindow = new float[residues + 1];


    halfWindowApplied = _halfWindow;
    window = 2 * halfWindowApplied + 1;

    // Do the average window calculations 
    for (i = 0; i < residues; i++) {
        valuesWindow[i] = 0.F;
        for (j = i - halfWindowApplied; j <= i + halfWindowApplied; j++) {
            if (j < 0) valuesWindow[i] += values[-j];
            else if (j >= residues) valuesWindow[i] += values[((2 * residues - j) - 2)];
            else valuesWindow[i] += values[j];
        }

        // Calculate the average value, by dividing the values
        valuesWindow[i] = valuesWindow[i] / (float) window;
    }

    return true;
}

bool statisticsMold::isDefinedWindow() {
    return (halfWindowRequested != -1);
}

float *statisticsMold::getVector() {
    if (isDefinedWindow()) {
        // Check if the window has been applied
        if (halfWindowRequested != halfWindowApplied)
            applyWindow(halfWindowRequested);
        // Return the windowed value
        return valuesWindow;
    }
        // Return the original values
    else return values;
}

void statisticsMold::printByColumn() {

}

void statisticsMold::printAccumulative() {

}
