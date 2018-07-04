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

#include <sstream>
#include "Statistics/statisticsMold.h"
#include "newAlignment.h"
#include "reportsystem.h"
#include "TimerFactory.h"

const std::string statisticsMold::statName = "Gaps"; // NOLINT

statisticsMold::statisticsMold(newAlignment *parentAlignment, statisticsMold *mold) {
    _alignment = parentAlignment;
    refCounter = mold->refCounter;
    refCounter++;
}

statisticsMold::statisticsMold(newAlignment *parentAlignment) {
    _alignment = parentAlignment;
    refCounter = new int(1);
}

statisticsMold::~statisticsMold() {
    if (--(*refCounter) == 0) {
        delete[] values;
        delete refCounter;
    }
    delete[] valuesWindow;
    _alignment = nullptr;
}

bool statisticsMold::applyWindow(int _halfWindow) {
    // Calculate the values array if it has not been calculated previously
    if (values == nullptr)
        calculate();

    // Check is the half window value passed is in the valid range
    if (_halfWindow > residues / 4) {
        debug.report(ErrorCode::WindowTooBig);
        return false;
    }

    // If the current half window is the same as the last one, don't do anything
    if (halfWindowApplied == _halfWindow) return true;

    // If the half window requested is 0 or a negative number
    // we simply delete the window values.
    if (_halfWindow < 1) {
        if (halfWindowApplied > 0)
            delete[] valuesWindow;

        valuesWindow = nullptr;
        return true;
    }

    // Initialize the values used in the calculation
    int i, j;

    // Initialize the values window array if it's null
    if (valuesWindow == nullptr)
        valuesWindow = new float[residues];

    // Special value. We don't compute the average
    // if there is no neighbour on one side
    valuesWindow[0] = values[0];
    int currentWindowSize = 1;

    // First loop: While index is lower than the window value,
    // We apply a window as big as the index
    for (i = 1; i < _halfWindow + 1; i++) {

        // Increase the window size by 2, one residue on each side.
        currentWindowSize += 2;

        valuesWindow[i] = 0.F;

        for (j = 0; j < currentWindowSize; j++) {
            valuesWindow[i] += values[j];
        }

        valuesWindow[i] = valuesWindow[i] / currentWindowSize;
    }

    // Second loop: Apply requested window size until
    // i + windowsize > array size
    {
        float window = 2 * _halfWindow + 1;

        for (; i != residues - _halfWindow; i++) {

            valuesWindow[i] = values[i];

            for (j = 1; j < _halfWindow + 1; j++) {
                valuesWindow[i] += values[i + j];
                valuesWindow[i] += values[i - j];
            }

            valuesWindow[i] = valuesWindow[i] / window;
        }
    }

    // Third loop: As index is bigger than the number of residues - halfwindow
    // Apply a symetric window as big as possible
    for (; i < residues - 1; i++) {

        // Decrease the window size by 2, one residue on each side.
        currentWindowSize -= 2;

        valuesWindow[i] = values[i];

        for (j = 1; (i + j) < residues; j++) {
            valuesWindow[i] += values[i + j];
            valuesWindow[i] += values[i - j];
        }

        // Calculate the average value, by dividing the current value (sumatory)
        // by the current window size
        valuesWindow[i] = valuesWindow[i] / currentWindowSize;

    }

    // Special value. We don't compute the average
    // if there is no neighbour on one side
    valuesWindow[residues - 1] = values[residues - 1];

    // New window has been applied
    halfWindowApplied = _halfWindow;

    return true;
}

bool statisticsMold::isDefinedWindow() {
    return valuesWindow != nullptr;
}

float *statisticsMold::getVector() {
    // Return the window values if it is allocated.
    // Else, return the raw values
    return valuesWindow == nullptr ? values : valuesWindow;
}

void statisticsMold::printByColumn(bool calculateRelative) {

    int i, size = 20;

    std::string &filename = _alignment->filename;

    cout << std::fixed << std::setw(filename.length())
         << std::setfill(' ') << std::left << "" << endl

         << "#\33[0;31m File:\33[0;1m "
         << filename << "\33[0m"

         << std::setw(filename.length())
         << std::setfill(' ')
         << std::left << "" << endl;

    cout << "#\33[0;36m BlockSize: \33[0;1m"
         << size << "\33[0m" << endl

         << "#\33[0;32m Statistic:\33[0;1m "
         << statName << " per column\33[0m" << endl;

    int headerSize = std::max(8 + filename.length(), 24 + statName.length());

    cout << std::setw(headerSize) << std::setfill('-')
         << std::left << "" << std::setfill(' ') << endl;

    cout << "\33[0;33;1m"
         << std::setw(size) << std::left << " Residue"
         << std::left << " " << statName << endl

         << std::setw(size) << std::left << " Number"
         << std::left << " Value" << endl
         << std::setfill('-') << "\33[0;m"
         << std::setw(size) << std::right << "  "
         << std::setw(size) << std::right << "  " << endl
         << std::setfill(' ');

    cout.precision(5);

    float *reportValues;

    // If vector is defined,
    // we use it to print the conservation's values.
    reportValues = getVector();

    for (i = 0; i < _alignment->originalResidNumber; i++)
        cout << setw(size) << std::left << i
             << setw(2) << std::right << reportValues[i] << endl;
}

void statisticsMold::printAccumulative(bool calculateRelative) {
    float refer, *vectAux;
    int i, num, acm;
    int size = std::max(20, (int)statName.length() + 4);

    // Allocate memory
    vectAux = new float[_alignment->originalResidNumber];

    // Copy the current vector (windowed or not) to a new array
    utils::copyVect(getVector(), vectAux, _alignment->originalResidNumber);

    // Sort the copied vector
    utils::quicksort(vectAux, 0, _alignment->originalResidNumber - 1);

    std::string &filename = _alignment->filename;

    cout << std::fixed << std::setw(filename.length())
         << std::setfill(' ') << std::left << "" << endl

         << "#\33[0;31m File:\33[0;1m "
         << filename << "\33[0m"

         << std::setw(filename.length())
         << std::setfill(' ')
         << std::left << "" << endl;

    cout << "#\33[0;36m BlockSize: \33[0;1m"
         << size << "\33[0m" << endl

         << "#\33[0;32m Statistic:\33[0;1m "
         << statName << " total\33[0m" << endl;

    int headerSize = std::max(8 + filename.length(), 19 + statName.length());

    cout << std::setw(headerSize) << std::setfill('-')
         << std::left << "" << std::setfill(' ') << endl;

    std::stringstream firstLine;
    std::stringstream secondLine;
    std::stringstream thirdLine;

    firstLine << std::setw(size) << std::left << "";
    secondLine << std::setw(size) << std::left << " Number of";
    thirdLine << std::setw(size) << std::left << " Residues";

    firstLine << std::setw(size) << std::left << "";
    secondLine << std::setw(size) << std::left << " Percentage";
    thirdLine << std::setw(size) << std::left << " of Alignment";

    firstLine << std::setw(size) << std::left << " Accumulative";
    secondLine << std::setw(size) << std::left << " Number of";
    thirdLine << std::setw(size) << std::left << " Residues";

    firstLine << std::setw(size) << std::left << " Accumulative";
    secondLine << std::setw(size) << std::left << " percentage";
    thirdLine << std::setw(size) << std::left << " of alignment";

    firstLine << std::setw(size) << std::left << " ";
    secondLine << std::setw(size) << std::left << " " + statName;
    thirdLine << std::setw(size) << std::left << " Value";

    if (calculateRelative)
    {
        firstLine << std::setw(size) << std::left << " Percentage of";
        secondLine << std::setw(size) << std::left << " " + statName;
        thirdLine << std::setw(size) << std::left << " per column";

        firstLine << std::setw(size) << std::left << " " + statName;
        secondLine << std::setw(size) << std::left << " score";
        thirdLine << std::setw(size) << std::left << " per column";
    }

    cout << "\33[0;33;1m"
         << firstLine.rdbuf() << endl
         << secondLine.rdbuf() << endl
         << thirdLine.rdbuf() << endl
         << "\33[0;m"
         << std::setfill('-');

    for (i = 0; i < (calculateRelative ? 7 : 5); i++)
        cout << setw(size) << std::right << "   ";

    cout << endl << setfill(' ');
    cout.precision(10);


    // Initialize temporal values
    refer = vectAux[_alignment->originalResidNumber - 1];
    acm = 0;
    num = 1;

    // Count the columns with the same value and compute this information
    // to show the accumulative stat in the alignment.
    for (i = _alignment->originalResidNumber - 2; i > -1; i--) {
        acm++;

        if (refer != vectAux[i]) {

            cout
                << setw(size) << std::left << num

                << setw(size) << std::left
                << setw(size - 6) << std::right << ((float) num / _alignment->originalResidNumber * 100.0F)
                << setw(6) << std::right << " "

                << setw(size) << std::left << acm

                << setw(size) << std::left
                << setw(size - 6) << std::right << ((float) acm / _alignment->originalResidNumber * 100.0F)
                << setw(6) << std::right << " "

                << setw(size) << std::left << refer;

            if (calculateRelative)
            {
                cout
                    << setw(size) << std::left
                    << setw(size - 6) << std::right << (vectAux[i] * 100.0F) / _alignment->originalResidNumber
                    << setw(6) << std::right << " "

                    << setw(size) << std::left
                    << setw(size - 6) << std::right << 1.0F - (vectAux[i] / _alignment->originalResidNumber)
                    << setw(6) << std::right << " ";
            }

            cout << endl;
            refer = vectAux[i];
            num = 1;
        } else num++;
    }
    acm++;

    // Repeat the body of the last for, to show the last result
    cout
        << setw(size) << std::left << num

        << setw(size) << std::left
        << setw(size - 6) << std::right << ((float) num / _alignment->originalResidNumber * 100.0F)
        << setw(6) << std::right << " "

        << setw(size) << std::left << acm

        << setw(size) << std::left
        << setw(size - 6) << std::right << ((float) acm / _alignment->originalResidNumber * 100.0F)
        << setw(6) << std::right << " "

        << setw(size) << std::left << refer;

    if (calculateRelative)
    {
        cout
                << setw(size) << std::left
                << setw(size - 6) << std::right << (vectAux[i] * 100.0F) / _alignment->originalResidNumber
                << setw(6) << std::right << " "

                << setw(size) << std::left
                << setw(size - 6) << std::right << 1.0F - (vectAux[i] / _alignment->originalResidNumber)
                << setw(6) << std::right << " ";
    }

    cout << endl;

    // Deallocate the reserved memory.
    delete[] vectAux;
}
