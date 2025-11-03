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

#include "Statistics/Entropy.h"
#include "InternalBenchmarker.h"
#include "Statistics/Manager.h"
#include "reportsystem.h"
#include "Alignment/Alignment.h"
#include "defines.h"
#include "utils.h"
#include "math.h"

namespace statistics {

    Entropy::Entropy(Alignment *parentAlignment) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("Entropy::Entropy(Alignment *parentAlignment) ");

        alig = parentAlignment;

        EntropyValues = new float[alig->originalNumberOfResidues];
        utils::initlVect(EntropyValues, alig->originalNumberOfResidues, 0);

        refCounter = new int(1);
    }

    Entropy::Entropy(Alignment *parentAlignment,
                           Entropy *mold) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("Entropy::Entropy(Alignment *parentAlignment) ");

        alig = parentAlignment;

        halfWindow = 0;

        EntropyValues = mold->EntropyValues;
        EntropyValues_Window = mold->EntropyValues_Window;

        refCounter = mold->refCounter;
        (*refCounter)++;
    }

    Entropy::~Entropy() {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("Entropy::Entropy");
        // Deallocate common variables only in case there is no module that has a
        // reference to them
        if (refCounter == nullptr || --(*refCounter) < 1) {
            delete[] EntropyValues;
            EntropyValues = nullptr;
            delete[] EntropyValues_Window;
            EntropyValues_Window = nullptr;
            delete refCounter;
            refCounter = nullptr;
        }
    }


    bool Entropy::calculateVectors([[maybe_unused]] bool cutByGap = false) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("bool Entropy::calculateVectors(bool cutByGap) ");

        std::map<char, float> residueCounter {};

        // Initialize the variables used
        int residueIndex, sequenceIndex;

        // Depending on alignment type, indetermination symbol will be one or other
        char indet = alig->getAlignmentType() & SequenceTypes::AA ? 'X' : 'N';

        // Calculate the max entropy value for the current dataType, to normalize data
        float maxEntropyDivisor;
        int N = 0;
        switch (alig->getAlignmentType())
        {
            case DNA: [[fallthrough]];
            case RNA:
                N = 4;  break; // 4 Elements
            case (DNA | DEG): [[fallthrough]];
            case (RNA | DEG):
                N = 14; break; // 4 Elements + 6 Pairs + 4 Triplets
            case AA:
                N = 22; break; // 22 Elements
            case (AA | DEG):
                N = 24; break; // 24 Elements
            default:
                debug.log(VerboseLevel::ERROR)
                    << "Contact with developer: Entropy::calculateVectors\n";
                break;
        }
        // Max entropy is obtained when all entries have the same proportion: 1/N
        // Max Shannon entropy = 1.0F/N * -log2(1.0F/N) * N -> log2(1.0F/N)
        // We make 1F/log2(1.0F/N) so we can use it as multiplier and not divisor.
        //      This makes the calculation faster
        maxEntropyDivisor = 1.0F / log2(1.0F / N);

        float sequenNumberDivisor = 1.0F / alig->originalNumberOfSequences;
        char currentResidue;
        // this->alig->Statistics->calculateGapStats(); To be deleted if not considering gap stats

        for (residueIndex = 0; residueIndex < alig->originalNumberOfResidues; residueIndex++)
        {
            residueCounter.clear();
            for (sequenceIndex = 0; sequenceIndex < alig->originalNumberOfSequences; sequenceIndex++)
            {
                currentResidue = alig->sequences[sequenceIndex][residueIndex];
                currentResidue = utils::toUpper(currentResidue);

                if (residueCounter.count(currentResidue) == 0) {
                    residueCounter[currentResidue] = 1;
                } else {
                    residueCounter[currentResidue]++;
                }
            }

            // residueCounter[indet] = 1; ??
            // residueCounter['-'] = 1; ??

            EntropyValues[residueIndex] = 0;
            float residueFrequency;
            for (auto & residueCount : residueCounter) {
                // Divide the raw count by the number of sequences to obtain the freq
                residueFrequency = residueCount.second * sequenNumberDivisor;
                // Add the entropy for each element
                EntropyValues[residueIndex] += residueFrequency * log2(residueFrequency);
            }
            // Although Shannon entropy is the negative value of the sum of (x*log2(x))
            // We remove the negative values
            //  by multiplying with maxEntropyDivisor, which is also negative.
            EntropyValues[residueIndex] *= maxEntropyDivisor;
        }

        return true;
    }

    bool Entropy::applyWindow(int halfW) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("bool Entropy::applyWindow(int halfW) ");

        if (EntropyValues == nullptr)
            calculateVectors(false);

        if (halfW > alig->originalNumberOfResidues / 4) {
            debug.report(ErrorCode::EntropyWindowTooBig);
            return false;
        }

        // If the current half window is the same as the last one, don't do anything
        if (halfWindow == halfW) return true;

        // Save the requested half window. This is useful when making a copy of the
        // alignment, as the window values are not valid anymore but don't want to
        // calculate them if not needed anymore
        halfWindow = halfW;

        // If the half window requested is 0 or a negative number
        // we simply delete the window values.
        if (halfW < 1) {
            delete[] EntropyValues_Window;

            EntropyValues_Window = nullptr;
            return true;
        }

        // Initialize the values used in the calculation
        int i, j, window;

        // Initialize the MDK window array if it's null
        if (EntropyValues_Window == nullptr)
            EntropyValues_Window = new float[alig->originalNumberOfResidues + 1];

        window = 2 * halfWindow + 1;

        // Do the average window calculations
        for (i = 0; i < alig->originalNumberOfResidues; i++) {
            EntropyValues_Window[i] = 0.F;
            for (j = i - halfWindow; j <= i + halfWindow; j++) {
                if (j < 0)
                    EntropyValues_Window[i] += EntropyValues[-j];
                else if (j >= alig->originalNumberOfResidues)
                    EntropyValues_Window[i] += EntropyValues[((2 * alig->originalNumberOfResidues - j) - 2)];
                else
                    EntropyValues_Window[i] += EntropyValues[j];
            }

            // Calculate the average value, by dividing the values
            EntropyValues_Window[i] = EntropyValues_Window[i] / (float) window;
        }
        return true;
    }

    bool Entropy::isDefinedWindow() {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("bool Entropy::isDefinedWindow(void) ");

        return (halfWindow > 0);
    }

    float *Entropy::getValues() {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("float *Entropy::getMdkWindowedVector(void) ");

        // If a window is defined
        if (isDefinedWindow()) {
            // Check if the window has been applied
            if (EntropyValues_Window == nullptr)
                applyWindow(halfWindow);
            // Return the windowed value
            return EntropyValues_Window;
        }
            // Return the original values
        else return EntropyValues;
    }

    double Entropy::calcCutPoint(float baseLine, float conservationPct) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("double Entropy::calcCutPoint(float baseLine, float conservationPct) ");
        // It computes the cutting point based on alignment's Entropy values -
        // the so-called 'Entropy'. It also takes into account the minimum percentage
        // from the input alignment to be kept. Depending on those two values, the
        // method will select a different cutting-point.

        double cuttingPoint_MinimumConserv, cuttingPoint_SimilThreshold;
        int i, highestPos;
        float *vectAux;

        vectAux = new float[alig->originalNumberOfResidues];

        // Sort a copy of the vector containing the Entropy values after applying
        // any windows methods. Take the columns value that it lower than the minimum
        // Entropy threshold set by the user
        utils::copyVect(getValues(), vectAux, alig->originalNumberOfResidues);
        utils::quicksort(vectAux, 0, alig->originalNumberOfResidues - 1);

        for (i = alig->originalNumberOfResidues - 1; i >= 0; i--)
            if (vectAux[i] < conservationPct)
                break;
        cuttingPoint_SimilThreshold = vectAux[i];

        // It is possible that due to number casting, we get a number out of the
        // vector containing the Entropy values - it is not reporting an overflow
        // situation but giving back a 0 when it should be a number equal (or closer)
        // to 1.
        highestPos = (int) ((double) (alig->originalNumberOfResidues - 1) * (100.0 - baseLine) / 100.0);
        highestPos = highestPos < (alig->originalNumberOfResidues - 1) ? highestPos : alig->originalNumberOfResidues - 1;
        cuttingPoint_MinimumConserv = vectAux[highestPos];

        delete[] vectAux;

        // Return the minimum cutting point between the one set by the threshold and
        // the one set by the minimum percentage of the input alignment to be kept
        return (cuttingPoint_MinimumConserv < cuttingPoint_SimilThreshold ?
                cuttingPoint_MinimumConserv : cuttingPoint_SimilThreshold);
    }

    void Entropy::printConservationColumns() {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("void Entropy::printConservationColumns(void) ");

        int i, size = 20;

        std::string fname = alig->filename;


        std::cout << std::fixed
                  << std::setw(fname.length() + 7)
                  << std::setfill(' ')
                  << std::left << "" << std::endl;

        std::cout << "#\33[0;31m File :\33[0;1m" << fname << "\33[0m";

        fname = std::to_string(size);

        std::cout
                << std::setw(fname.length() + 7)
                << std::setfill(' ')
                << std::left << "" << std::endl;

        std::cout << "#\33[0;36m BlockSize : \33[0;1m" << fname << "\33[0m" << std::endl;

        fname = " Entropy per Column";

        std::cout << "#\33[0;32m Statistic :\33[0;1m" << fname << "\33[0m" << std::endl;

        std::cout << std::setw(alig->filename.substr(6, alig->filename.size() - 7).length() + 7)
                  << std::setfill('-')
                  << std::left << ""
                  << std::setfill(' ')
                  << std::endl;

        std::cout << "\33[0;33;1m"
                  << std::setw(size) << std::left << " Residue" << std::left << " Entropy" << std::endl
                  << std::setw(size) << std::left << " Number" << std::left << " Value" << std::endl
                  << std::setfill('-')
                  << "\33[0;m"
                  << std::setw(size) << std::right << "  "
                  << std::setw(size) << std::right << "  " << std::endl
                  << std::setfill(' ');

        std::cout.precision(10);

        float *values;

        // If MDK_Window vector is defined, we use it to print the Entropy's values.
        if (EntropyValues_Window != nullptr)
            values = EntropyValues_Window;
            // In others cases, we uses the MDK vector to print the Entropy's vlaues.
        else
            values = EntropyValues;

        for (i = 0; i < alig->originalNumberOfResidues; i++)
            std::cout << std::setw(size) << std::left << i << values[i] << std::endl;
    }

    void Entropy::printConservationAcl() {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("void Entropy::printConservationAcl(void) ");

        float refer, *vectAux;
        int i, num, acm;
        int size = 20;
        // Allocate memory
        vectAux = new float[alig->originalNumberOfResidues];

        // Select the Entropy's value source and copy that vector in a auxiliar vector
        if (EntropyValues_Window != nullptr) utils::copyVect(EntropyValues_Window, vectAux, alig->originalNumberOfResidues);
        else utils::copyVect(EntropyValues, vectAux, alig->originalNumberOfResidues);

        // Sort the auxiliar vector.
        utils::quicksort(vectAux, 0, alig->originalNumberOfResidues - 1);

        // Print filename
        std::string fname = alig->filename;

        std::cout << std::fixed
                  << std::setw(fname.length() + 7)
                  << std::setfill(' ')
                  << std::left << "" << std::endl;

        std::cout << "#\33[0;31m File :\33[0;1m" << fname << "\33[0m";

        fname = std::to_string(size);

        std::cout
                << std::setw(fname.length() + 7)
                << std::setfill(' ')
                << std::left << "" << std::endl;

        std::cout << "#\33[0;36m BlockSize : \33[0;1m" << fname << "\33[0m" << std::endl;

        fname = " Entropy Total";

        std::cout << "#\33[0;32m Statistic :\33[0;1m" << fname << "\33[0m" << std::endl;

        std::cout << std::setw(alig->filename.substr(6, alig->filename.size() - 7).length() + 7)
                  << std::setfill('-')
                  << std::left << ""
                  << std::setfill(' ')
                  << std::endl;


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
        secondLine << std::setw(size) << std::left << " Entropy";
        thirdLine << std::setw(size) << std::left << " Value";

        std::cout << "\33[0;33;1m"
                  << firstLine.rdbuf() << std::endl
                  << secondLine.rdbuf() << std::endl
                  << thirdLine.rdbuf() << std::endl
                  << "\33[0;m"
                  << std::setfill('-');

        for (i = 0; i < 5; i++)
            std::cout << std::setw(size) << std::right << "   ";

        std::cout << std::endl << std::setfill(' ');
        std::cout.precision(10);


        // Initializate some values
        refer = vectAux[alig->originalNumberOfResidues - 1];
        acm = 0;
        num = 1;

        // Count the columns with the same Entropy's value and compute this information to shows the accumulative
        // statistics in the alignment.
        for (i = alig->originalNumberOfResidues - 2; i >= 0; i--) {
            acm++;

            if (refer != vectAux[i]) {

                std::cout
                        << std::setw(size) << std::left << num

                        << std::setw(size) << std::left
                        << std::setw(size - 6) << std::right << ((float) num / alig->originalNumberOfResidues * 100.0F)
                        << std::setw(6) << std::right << " "

                        << std::setw(size) << std::left << acm

                        << std::setw(size) << std::left
                        << std::setw(size - 6) << std::right << ((float) acm / alig->originalNumberOfResidues * 100.0F)
                        << std::setw(6) << std::right << " "

                        << std::setw(size) << std::left << refer

                        << std::endl;
                refer = vectAux[i];
                num = 1;
            } else num++;
        }
        acm++;

        std::cout
                << std::setw(size) << std::left << num

                << std::setw(size) << std::left
                << std::setw(size - 6) << std::right << ((float) num / alig->originalNumberOfResidues * 100.0F)
                << std::setw(6) << std::right << " "

                << std::setw(size) << std::left << acm

                << std::setw(size) << std::left
                << std::setw(size - 6) << std::right << ((float) acm / alig->originalNumberOfResidues * 100.0F)
                << std::setw(6) << std::right << " "

                << std::setw(size) << std::left << refer

                << std::endl;

        // Deallocate the reserved memory.
        delete[] vectAux;
    }

}