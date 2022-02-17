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
#include "InternalBenchmarker.h"
#include "Statistics/Gaps.h"
#include "reportsystem.h"
#include "Alignment/Alignment.h"
#include "defines.h"
#include "utils.h"

namespace statistics {

    Gaps::Gaps(Alignment *parent) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("Gaps::Gaps(Alignment *parent) ");

        alig = parent;

        maxGaps = 0;
        halfWindow = 0;
        totalGaps = 0;

        // Memory allocation for the vectors and its initialization
        gapsInColumn = new int[alig->originalNumberOfResidues];
        utils::initlVect(gapsInColumn, alig->originalNumberOfResidues, 0);

        numColumnsWithGaps = new int[alig->originalNumberOfSequences + 1];
        utils::initlVect(numColumnsWithGaps, alig->originalNumberOfSequences + 1, 0);

        refCounter = new int(1);
    }

    Gaps::Gaps(Alignment *pAlignment,
                                   Gaps *pGaps) {
        alig = pAlignment;

        maxGaps = 0;
        totalGaps = 0;

        // Pointer initialization
        gapsInColumn = pGaps->gapsInColumn;

        numColumnsWithGaps = pGaps->numColumnsWithGaps;

        gapsWindow = pGaps->gapsWindow;

        // Count the gaps and indeterminations of each columns
        maxGaps = pGaps->maxGaps;

        refCounter = pGaps->refCounter;
        (*refCounter)++;
    }


    Gaps::~Gaps() {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("Gaps::Gaps");

        // Only free memory if there is previous memory allocation
        if (--(*refCounter) == 0) {
            delete[] gapsInColumn;
            delete[] numColumnsWithGaps;
            delete[] gapsWindow;
            delete refCounter;
        }
    }

    bool Gaps::applyWindow(int halfW) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("bool Gaps::applyWindow(int halfW) ");

        if (halfW > alig->originalNumberOfResidues / 4) {
            debug.report(ErrorCode::GapWindowTooBig);
            return false;
        }

        // Save the requested half window. This is useful when making a copy of the
        // alignment, as the window values are not valid anymore but don't want to
        // calculate them if not needed anymore
        halfWindow = halfW;

        // If the half window requested is 0 or a negative number
        // we simply delete the window values.
        if (halfW < 1) {

            delete[] gapsWindow;
            gapsWindow = nullptr;
            return true;
        }

        int i, j, window;

        if (gapsWindow == nullptr)
            gapsWindow = new int[alig->originalNumberOfResidues];

        // Initialize to 0 the vector that will store the number of gaps of each column
        // and the vector that will store window processing results
        utils::initlVect(numColumnsWithGaps, alig->originalNumberOfSequences + 1, 0);

        // Initialize maximum gaps number per column value and store the mediumWinSize value in the object.
        maxGaps = 0;
        window = (2 * halfWindow + 1);

        // We calculate some statistics for every column in the alignment,and the maximum gaps' number value
        for (i = 0; i < alig->originalNumberOfResidues; i++) {
            gapsWindow[i] = 0;
            // Sum the total number of gaps for the considered window
            for (j = i - halfWindow, gapsWindow[i] = 0; j <= i + halfWindow; j++) {
                if (j < 0)
                    gapsWindow[i] += gapsInColumn[-j];
                else if (j >= alig->originalNumberOfResidues)
                    gapsWindow[i] += gapsInColumn[((2 * alig->originalNumberOfResidues - j) - 2)];
                else
                    gapsWindow[i] += gapsInColumn[j];
            }

            // Calculate, and round to the nearest integer, the number of gaps for the i column
            gapsWindow[i] = utils::roundInt(((double) gapsWindow[i] / window));
            // Increase in 1 the number of colums with the same number of gaps than column i
            numColumnsWithGaps[gapsWindow[i]]++;

            // Update the max number of gaps in the alignment, if neccesary
            if (gapsWindow[i] > maxGaps)
                maxGaps = gapsWindow[i];
        }
        return true;
    }

    bool Gaps::isDefinedWindow() {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("bool Similarity::isDefinedWindow(void) ");

        return (halfWindow > 0);
    }

    int *Gaps::getGapsWindow() {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("int *Gaps::getGapsWindow(void) ");

        // If a window is defined
        if (isDefinedWindow()) {
            // Check if the window has been applied
            if (gapsWindow == nullptr)
                // Apply if not
                applyWindow(halfWindow);
            // Return the windowed value
            return gapsWindow;
        }
            // Return the original values
        else return gapsInColumn;
    }

    double Gaps::calcCutPoint(float minInputAlignment, float gapThreshold) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("double Gaps::calcCutPoint(float minInputAlignment, float gapThreshold) ");

        // Method to select the cutting point based on gaps values from the input
        // alignment. The cutting point is selected as the maximum gaps number allowed
        // in the output alignment given a minimum percentage of the input alignment
        // to be kept and a maximum gaps number. In case of both values set different
        // cutting points, the minimum percentage of the input alignment prevails.

        double cuttingPoint_MinimumConserv, cuttingPoint_gapThreshold;
        int i, acum;

        // Calculate the gap number represented by the gaps threshold. This gap number
        // represents the maximum gap number in any column in the output alignment
        cuttingPoint_gapThreshold = (double) alig->numberOfSequences * gapThreshold;

        // Compute the minimum columns number to be kept from the input alignment
        cuttingPoint_MinimumConserv = utils::roundInt(((double) (alig->originalNumberOfResidues * minInputAlignment) / 100.0));
        if (cuttingPoint_MinimumConserv > alig->originalNumberOfResidues)
            cuttingPoint_MinimumConserv = alig->originalNumberOfResidues;

        // We look at the number of gaps which allows us to keep the minimum columns
        // number from the input alignment
        for (i = 0, acum = 0; i < alig->originalNumberOfSequences; i++) {
            acum += numColumnsWithGaps[i];
            if (acum >= cuttingPoint_MinimumConserv)
                break;
        }

        // If there is not an exact number for the gaps cutting point, compute such
        // value as the inmediate superior gap number minus the proportion of columns
        // necessary to respect the minimum percentage from the input alignment to be
        // kept
        if (numColumnsWithGaps[i])
            cuttingPoint_MinimumConserv =
                    (double) (i - ((float) (acum - cuttingPoint_MinimumConserv) / numColumnsWithGaps[i]));
        else
            cuttingPoint_MinimumConserv = 0;

        // Return the maximum gap number of the two computed cutting points.
        return (utils::max(cuttingPoint_MinimumConserv, cuttingPoint_gapThreshold));
    }

    int Gaps::calcCutPointMixSlope() {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("int Gaps::calcCutPointMixSlope(void) ");

        float delta = 0, maxSlope = -1, *firstSlopeVector, *secondSlopeVector;
        int prev, pprev, maxIter, row = 1, act = 0, max = 0;

        // We build two slope vectors, one vector for the first slope and another one for the second.
        firstSlopeVector = new float[maxGaps + 1];
        secondSlopeVector = new float[maxGaps + 1];

        // We initialize them with -1.0 value and fix the maximum iteractions' number as maximun gaps' number plus 1
        utils::initlVect(firstSlopeVector, maxGaps, -1.0F);
        utils::initlVect(secondSlopeVector, maxGaps, -1.0F);
        maxIter = maxGaps + 1;

        // Until to achieve the maximum iteractions' number.
        while (act < maxIter) {

            // We look for a first point to second slope.
            while ((numColumnsWithGaps[act]) == 0) act++;
            pprev = act;
            if ((act + 1) >= maxIter) break;

            // We look for a first point to first slope.
            do {
                act++;
            } while ((numColumnsWithGaps[act]) == 0);
            prev = act;
            if ((act + 1) >= maxIter) break;

            // We look for a second point to first and second slope.
            do {
                act++;
            } while ((numColumnsWithGaps[act]) == 0);
            if (act >= maxIter) break;

            // Calculate the first slope between the earlier previus and previus points.
            firstSlopeVector[prev] = ((float) (prev - pprev) / alig->originalNumberOfSequences);
            firstSlopeVector[prev] /= ((float) numColumnsWithGaps[prev] / alig->originalNumberOfResidues);

            // Calculate the first slope between the previus and current points.
            firstSlopeVector[act] = ((float) (act - prev) / alig->originalNumberOfSequences);
            firstSlopeVector[act] /= ((float) numColumnsWithGaps[act] / alig->originalNumberOfResidues);

            // Calculate the second slope between the earlier previus and current points.
            secondSlopeVector[act] = ((float) (act - pprev) / alig->originalNumberOfSequences);
            secondSlopeVector[act] /= ((float) (numColumnsWithGaps[act] + numColumnsWithGaps[prev]) / alig->originalNumberOfResidues);

            // If the ratio between ...
            if ((secondSlopeVector[pprev] != -1.0) || (firstSlopeVector[pprev] != -1.0)) {

                // first slope previus and first slope earlier previus points.
                if (firstSlopeVector[pprev] != -1.0) {
                    delta = firstSlopeVector[prev] / firstSlopeVector[pprev];
                    row = pprev;
                }

                // first slope current and first slope previus points.
                if (delta < (firstSlopeVector[act] / firstSlopeVector[prev])) {
                    delta = firstSlopeVector[act] / firstSlopeVector[prev];
                    row = prev;
                }

                // second slope current and second slope earlier previus points.
                if (secondSlopeVector[pprev] != -1.0) {
                    if (delta < (secondSlopeVector[act] / secondSlopeVector[pprev])) {
                        delta = secondSlopeVector[act] / secondSlopeVector[pprev];
                        row = pprev;
                    }
                }

                // ... is better that current maxSlope then we modify the maxSlope with the best ratio.
                if (delta > maxSlope) {
                    maxSlope = delta;
                    max = row;
                }
            }
            act = prev;
        }

        // We delete the local memory.
        delete[] firstSlopeVector;
        delete[] secondSlopeVector;

        // and, finally, we return the cut point.
        return max;
    }

    int Gaps::calcCutPoint2ndSlope() {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("int Gaps::calcCutPoint2ndSlope(void) ");

        float maxSlope = -1, *secondSlopeVector;
        int prev, pprev, maxIter, act = 0, max = 0;

        // We build one slope vector and fix the maximum iterations' number
        // as the gaps'number plus 1.
        secondSlopeVector = new float[maxGaps + 1];
        utils::initlVect(secondSlopeVector, maxGaps + 1, -1.0F);
        maxIter = maxGaps + 1;

        // Find the lowest number of gaps into the input alignment. If there are few
        // points, it is possible that lowest number of gaps is returned as the threshold.
        // It could happen input alignment does not have columns with no-gaps
        for (act = 0, max = 0; numColumnsWithGaps[act] == 0; act++)
            max = act + 1;

        act = 0;
        while (act < maxIter) {

            // We look for a first point to second slope.
            while ((numColumnsWithGaps[act]) == 0)
                act++;
            pprev = act;
            if ((act + 1) >= maxIter)
                break;

            // We look for a first point to first slope.
            do {
                act++;
            } while ((numColumnsWithGaps[act]) == 0);
            prev = act;
            if ((act + 1) >= maxIter)
                break;

            // We look for a second point to first and second slope.
            do {
                act++;
            } while ((numColumnsWithGaps[act]) == 0);
            if (act >= maxIter)
                break;

            // Calculate the second slope between the earlier previous and current points.
            secondSlopeVector[act] = ((float) (act - pprev) / alig->numberOfSequences);
            secondSlopeVector[act] /= ((float) (numColumnsWithGaps[act] + numColumnsWithGaps[prev]) / alig->originalNumberOfResidues);

            // If the ratio between second slope current and second slope earlier previous points.
            if (secondSlopeVector[pprev] != -1.0) {
                if ((secondSlopeVector[act] / secondSlopeVector[pprev]) > maxSlope) {
                    maxSlope = (secondSlopeVector[act] / secondSlopeVector[pprev]);
                    max = pprev;
                }
            } else if (secondSlopeVector[prev] != -1.0) {
                if ((secondSlopeVector[act] / secondSlopeVector[prev]) > maxSlope) {
                    maxSlope = (secondSlopeVector[act] / secondSlopeVector[prev]);
                    max = pprev;
                }
            }
            act = prev;
        }

        // We deallocate local memory.
        delete[] secondSlopeVector;

        // Finally, we return the selected cut point.
        return max;
    }

    void Gaps::printGapsColumns() const {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("void Gaps::printGapsColumns(void) ");

        // Colors can be removed with
        // ' sed -r "s:\x1B\[[0-9;]*[mK]::g" '

        int *vectAux;

        // We allocate a local vector to recovery information on it
        vectAux = new int[alig->originalNumberOfResidues];

        // We decide about the information's source then we get the information.
        if (halfWindow == 0)
            utils::copyVect(gapsInColumn, vectAux, alig->originalNumberOfResidues);
        else
            utils::copyVect(gapsWindow, vectAux, alig->originalNumberOfResidues);

        int size = 20;

        std::string fname = alig->filename;

        std::cout
                << std::setw(fname.length() + 7)
                << std::setfill(' ')
                << std::left << "" << "\n";

        std::cout << "#\33[0;31m File :\33[0;1m" << fname << "\33[0m";

        fname = std::to_string(size);

        std::cout
                << std::setw(fname.length() + 7)
                << std::setfill(' ')
                << std::left << "" << "\n";

        std::cout << "#\33[0;36m BlockSize : \33[0;1m" << fname << "\33[0m" << "\n";

        fname = " Gaps per Column";

        std::cout << "#\33[0;32m Statistic :\33[0;1m" << fname << "\33[0m" << "\n";

        std::cout << std::setw(alig->filename.size())
                  << std::setfill('-')
                  << std::left << ""
                  << std::setfill(' ')
                  << "\n";

        std::cout << std::setfill(' ') << "\33[0;33;1m"
                  << std::setw(size) << std::left << " Residue"
                  << std::setw(size) << std::left << " % Gaps"
                  << std::setw(size) << std::left << " Gap Score"
                  << "\33[0m" << "\n";

        std::cout << std::setfill('-');

        std::cout << std::setw(size) << std::right << "  "
                  << std::setw(size) << std::right << "  "
                  << std::setw(size) << std::right << "  " << "\n";

        std::cout << std::setfill(' ') << std::fixed;
        std::cout.precision(10);

        // Show the information that have been requered
        for (int i = 0; i < alig->originalNumberOfResidues; i++)
            std::cout << std::setw(size) << std::setfill(' ') << std::left << i
                      << std::setw(size) << std::setfill(' ') << std::left
                      << std::setw(size - 6) << std::setfill(' ') << std::right
                      << (vectAux[i] * 100.0) / alig->originalNumberOfSequences
                      << std::setw(size) << std::setfill(' ') << std::right
                      << 1.F - (float(vectAux[i]) / alig->originalNumberOfSequences) << "\n";

        // Finally, we deallocate the local memory
        delete[] vectAux;
    }

    void Gaps::printGapsAcl() const {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("void Gaps::printGapsAcl(void) ");

        // Colors can be removed with
        // ' sed -r "s:\x1B\[[0-9;]*[mK]::g" '

        int acm, i;
        int size = 20;

        std::string fname = alig->filename;

        std::cout
                << std::setw(fname.length() + 7)
                << std::setfill(' ')
                << std::left << "" << "\n";

        std::cout << "#\33[0;31m File :\33[0;1m" << fname << "\33[0m";

        fname = std::to_string(size);

        std::cout
                << std::setw(fname.length() + 7)
                << std::setfill(' ')
                << std::left << "" << "\n";

        std::cout << "#\33[0;36m BlockSize : \33[0;1m" << fname << "\33[0m" << "\n";

        fname = " Gaps Total";

        std::cout << "#\33[0;32m Statistic :\33[0;1m" << fname << "\33[0m" << "\n";

        std::cout << std::setw(alig->filename.size())
                  << std::setfill('-')
                  << std::left << ""
                  << std::setfill(' ')
                  << "\n";


        std::stringstream firstLine;
        std::stringstream secondLine;
        std::stringstream thirdLine;

        firstLine << std::setw(size) << std::left << "";
        secondLine << std::setw(size) << std::left << " Number of";
        thirdLine << std::setw(size) << std::left << " residues";

        firstLine << std::setw(size) << std::left << "";
        secondLine << std::setw(size) << std::left << " Percentage";
        thirdLine << std::setw(size) << std::left << " of alignment";

        firstLine << std::setw(size) << std::left << "";
        secondLine << std::setw(size) << std::left << " Accumulative";
        thirdLine << std::setw(size) << std::left << " residues";

        firstLine << std::setw(size) << std::left << " Accumulative";
        secondLine << std::setw(size) << std::left << " percent of";
        thirdLine << std::setw(size) << std::left << " alignment";

        firstLine << std::setw(size) << std::left << " Number";
        secondLine << std::setw(size) << std::left << " of gaps";
        thirdLine << std::setw(size) << std::left << " per column";

        firstLine << std::setw(size) << std::left << " Percentage";
        secondLine << std::setw(size) << std::left << " of gaps";
        thirdLine << std::setw(size) << std::left << " per column";

        firstLine << std::setw(size) << std::left << "";
        secondLine << std::setw(size) << std::left << " Gaps score";
        thirdLine << std::setw(size) << std::left << " per column";

        std::cout << "\33[0;33;1m"
                  << firstLine.rdbuf() << "\n"
                  << secondLine.rdbuf() << "\n"
                  << thirdLine.rdbuf() << "\n"
                  << "\33[0;m"
                  << std::setfill('-');

        for (i = 0; i < 7; i++) {
            std::cout << std::setw(size) << std::right << "   ";
        }

        std::cout << "\n"
                  << std::fixed
                  << std::setfill(' ');

        std::cout.precision(10);

        // Count for each gaps' number the columns' number with that gaps' number.
        for (i = 0, acm = 0; i <= maxGaps; i++) {

            // If the columns' number with this gaps' number is not equal to zero, we will count the columns.
            if (numColumnsWithGaps[i] != 0) {

                // Compute and prints the accumulative values for the gaps in the alignment.
                acm += numColumnsWithGaps[i];

                // Number of residues
                std::cout << std::setw(size) << std::left << numColumnsWithGaps[i];

                // Percentage of alignment
                std::cout << std::setw(size) << std::left
                          << std::setw(size - 6) << std::right << (numColumnsWithGaps[i] * 100.0F) / alig->originalNumberOfResidues
                          << std::setw(6) << std::left << " ";

                // Accumulative residues
                std::cout << std::setw(size) << std::left << acm;

                // Accumulative percent of alignment
                std::cout << std::setw(size) << std::left
                          << std::setw(size - 6) << std::right << (acm * 100.0F) / alig->originalNumberOfResidues
                          << std::setw(6) << std::left << " ";

                // Number of gaps per column
                std::cout << std::setw(size) << std::left << i;

                // Percentage of gaps per column
                std::cout << std::setw(size) << std::left << (i * 1.0F) / alig->originalNumberOfSequences;

                // Gaps score per column
                std::cout << std::setw(size) << std::left << 1.F - (((float) i) / alig->originalNumberOfSequences);

                // End line
                std::cout << "\n";
            }
        }

        std::cout << std::setw(alig->filename.size())
                  << std::setfill('-')
                  << std::left << ""
                  << std::setfill(' ')
                  << "\n";

        std::cout << "\n"
                  << std::fixed
                  << std::setfill(' ');

        std::cout << "## AverageGaps\t" << float(totalGaps)/(alig->originalNumberOfResidues * alig->originalNumberOfSequences) << "\n";
        std::cout << "#> AverageGaps\tAverage gaps of the alignment" << "\n";
    }

    void Gaps::CalculateVectors() {
        int i, j;
        // Count the gaps and indeterminations of each columns
        for (i = 0; i < alig->originalNumberOfResidues; i++) {
            gapsInColumn[i] = 0;
            for (j = 0; j < alig->originalNumberOfSequences; j++) {
                if (alig->saveSequences[j] == -1) continue;
                if (alig->sequences[j][i] == '-') {
                    gapsInColumn[i]++;
                    totalGaps++;
                }
            }

            // Increase the number of colums with the number of gaps of the last processed column
            numColumnsWithGaps[gapsInColumn[i]]++;

            if (gapsInColumn[i] > maxGaps)
                maxGaps = gapsInColumn[i];
        }
    }
}