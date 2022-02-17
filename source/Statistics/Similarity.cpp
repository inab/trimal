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

#include "Statistics/Similarity.h"
#include "InternalBenchmarker.h"
#include "Statistics/Manager.h"
#include "reportsystem.h"
#include "Alignment/Alignment.h"
#include "defines.h"
#include "utils.h"
#include "omp.h"

namespace statistics {

    Similarity::Similarity(Alignment *parentAlignment) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("Similarity::Similarity(Alignment *parentAlignment) ");

        alig = parentAlignment;

        MDK = new float[alig->originalNumberOfResidues];
        utils::initlVect(MDK, alig->originalNumberOfResidues, 0);

        // Initialize the similarity matrix to nullptr.
        simMatrix = nullptr;

        refCounter = new int(1);
    }

    Similarity::Similarity(Alignment *parentAlignment,
                                                   Similarity *mold) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("Similarity::Similarity(Alignment *parentAlignment) ");

        alig = parentAlignment;

        halfWindow = 0;

        MDK = mold->MDK;
        MDK_Window = mold->MDK_Window;

        simMatrix = mold->simMatrix;

        refCounter = mold->refCounter;
        (*refCounter)++;
    }

    Similarity::~Similarity() {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("Similarity::Similarity");
        // Deallocate common variables only in case there is no module that has a
        // reference to them
        if (refCounter == nullptr || --(*refCounter) < 1) {
            delete[] MDK;
            MDK = nullptr;
            delete[] MDK_Window;
            MDK_Window = nullptr;

            if (matrixIdentity != nullptr)
                for (int i = 0; i < alig->numberOfSequences; i++)
                    delete[] matrixIdentity[i];
            delete[] matrixIdentity;
            matrixIdentity = nullptr;
            delete refCounter;
            refCounter = nullptr;
        }
    }

    void Similarity::calculateMatrixIdentity() {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("void Similarity::calculateMatrixIdentity() ");

        // We don't want to calculate the matrix identity
        // if it has been previously calculated
        if (matrixIdentity != nullptr)
            return;

        // Allocate temporal variables
        char indet;
        int i, j, k;
        int sum, length;

        // Allocate memory for the matrix identity
        matrixIdentity = new float *[alig->originalNumberOfSequences];
        #pragma omp parallel for num_threads(NUMTHREADS) if(alig->originalNumberOfSequences>MINPARALLELSIZE)
        for (i = 0; i < alig->originalNumberOfSequences; i++) {
            matrixIdentity[i] = new float[alig->originalNumberOfSequences];
        }

        // Depending on alignment type, indetermination symbol will be one or other
        indet = alig->getAlignmentType() & SequenceTypes::AA ? 'X' : 'N';

        // For each sequences' pair
        #pragma omp parallel for num_threads(NUMTHREADS) private(j, k, sum, length) if(alig->originalNumberOfSequences>MINPARALLELSIZE)
        for (i = 0; i < alig->originalNumberOfSequences; i++) {
            for (j = i + 1; j < alig->originalNumberOfSequences; j++) {
                // For each position in the alignment of that pair than we are processing
                for (k = 0, sum = 0, length = 0; k < alig->originalNumberOfResidues; k++) {
                    // If we find a element that is not a gap or an X aminoacid in the first sequence of the pair
                    if ((alig->sequences[i][k] != '-') &&
                        (alig->sequences[i][k] != indet)) {

                        // If we also find a valid element in the second sequence
                        if ((alig->sequences[j][k] != '-') &&
                            (alig->sequences[j][k] != indet))
                        {
                            // If the two valid elements are the same increase the sum
                            if (alig->sequences[j][k] == alig->sequences[i][k])
                                sum++;
                        }
                        // Increase the length of the sequence free of gaps and X elements
                        length++;
                    }
                    // If the first processed element is invalid and in the second we find a valid element increase the length of
                    // the sequence free of gaps and X elements
                    else if ((alig->sequences[j][k] != '-') &&
                             (alig->sequences[j][k] != indet))
                    {
                        length++;
                    }
                }

                // Calculate the value of matrix idn for columns j and i
                matrixIdentity[j][i] = (1.0F - ((float)sum / length) );
                matrixIdentity[i][j] = matrixIdentity[j][i];
            }
        }
    }

    bool Similarity::calculateVectors(bool cutByGap) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("bool Similarity::calculateVectors(int *gaps) ");

        // A similarity matrix must be defined. If not, return false
        if (simMatrix == nullptr)
            return false;

        // Calculate the matrix identity in case it's not done before
        if (matrixIdentity == nullptr)
            calculateMatrixIdentity();

        // Create the variable gaps, in case we want to cut by gaps
        int *gaps = nullptr;

        // Retrieve the gaps values in case we want to set to 0 the similarity value
        // in case the gaps value for that column is bigger or equal to 0.8F
        if (cutByGap)
        {
            if (alig->Statistics->gaps == nullptr)
                alig->Statistics->calculateGapStats();
            gaps = alig->Statistics->gaps->getGapsWindow();
        }

        // Initialize the variables used
        int i, j, k;
        float num, den;

        // Depending on alignment type, indetermination symbol will be one or other
        char indet = alig->getAlignmentType() & SequenceTypes::AA ? 'X' : 'N';

        // Q temporal value
        float Q;
        // Temporal chars that will contain the residues to compare by pair.
        char chA, chB;

        // Calculate the maximum number of gaps a column can have to calculate it's
        //      similarity
        float gapThreshold = 0.8F * alig->numberOfResidues;

        // For each column calculate the Q value and the MD value using an equation
        for (i = 0; i < alig->originalNumberOfResidues; i++) {
            // Set MDK for columns with gaps values bigger or equal to 0.8F
            if (cutByGap && gaps[i] >= gapThreshold) {
                MDK[i] = 0.F;
                continue;
            }
            // For each AAs/Nucleotides' pair in the column we compute its distance
            for (j = 0, num = 0, den = 0; j < alig->originalNumberOfSequences; j++) {


                // Calculate the upper value of the residue,
                //      to use in simMatrix->getDistance
                // This is faster than calculating the upper on that method
                //      as this is done before entering the loop
                // Doing this before checking if the element is indeterminate or gap
                //      allows to check if the indetermination is not capitalized
                chA = utils::toUpper(alig->sequences[j][i]);

                // We don't compute the distance if the first element is
                // a indeterminate (XN) or a gap (-) element.
                if ((chA == '-') || (chA == indet))
                    continue;

                for (k = j + 1; k < alig->originalNumberOfSequences; k++) {
                    // We calculate the upper value of the residue,
                    //      to use in simMatrix->getDistance
                    // This is equally faster as if it was done inside the method
                    //      but to prevent errors, the method doesn't 'upper'
                    //      the given chars.
                    // Doing this before checking if the element is indeterminate or gap
                    //      allows to check if the indetermination is not capitalized
                    chB = utils::toUpper(alig->sequences[k][i]);

                    // We don't compute the distance if the second element is
                    //      a indeterminate (XN) or a gap (-) element
                    if ((chB == '-') || (chB == indet))
                        continue;

                    // We use the identity value for the two pairs and
                    //      its distance based on similarity matrix's value.
                    num += matrixIdentity[j][k] * simMatrix->getDistance(chA, chB);
                    den += matrixIdentity[j][k];
                }
            }

            // If we are processing a column with only one AA/nucleotide, MDK = 0
            if (den == 0)
                MDK[i] = 0;
            else
            {
                Q = num / den;
                // If the MDK value is more than 1, we normalized this value to 1.
                //      Only numbers higher than 0 yield exponents higher than 1
                //      Using this we can test if the result is going to be higher than 1.
                //      And thus, prevent calculating the exp.
                // Take in mind that the Q is negative, so we must test if Q is LESSER
                //      than one, not bigger.
                if (Q < 0)
                    MDK[i] = 1.F;
                else
                    MDK[i] = exp(-Q);
            }
        }

        for (i = 0; i < alig->originalNumberOfSequences; i++)
            delete[] matrixIdentity[i];
        delete[] matrixIdentity;
        matrixIdentity = nullptr;

        return true;
    }

    bool Similarity::applyWindow(int halfW) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("bool Similarity::applyWindow(int halfW) ");

        // Calculate the MDK array if it has not been calculated previously
        if (MDK == nullptr)
            calculateVectors(true);

        // Check is the half window value passed is in the valid range
        if (halfW > alig->originalNumberOfResidues / 4) {
            debug.report(ErrorCode::SimilarityWindowTooBig);
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
            delete[] MDK_Window;

            MDK_Window = nullptr;
            return true;
        }

        // Initialize the values used in the calculation
        int i, j, window;

        // Initialize the MDK window array if it's null
        if (MDK_Window == nullptr)
            MDK_Window = new float[alig->originalNumberOfResidues + 1];

        window = 2 * halfWindow + 1;

        // Do the average window calculations
        for (i = 0; i < alig->originalNumberOfResidues; i++) {
            MDK_Window[i] = 0.F;
            for (j = i - halfWindow; j <= i + halfWindow; j++) {
                if (j < 0)
                    MDK_Window[i] += MDK[-j];
                else if (j >= alig->originalNumberOfResidues)
                    MDK_Window[i] += MDK[((2 * alig->originalNumberOfResidues - j) - 2)];
                else
                    MDK_Window[i] += MDK[j];
            }

            // Calculate the average value, by dividing the values
            MDK_Window[i] = MDK_Window[i] / (float) window;
        }
        return true;
    }

    bool Similarity::isDefinedWindow() {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("bool Similarity::isDefinedWindow(void) ");

        return (halfWindow > 0);
    }

    float *Similarity::getMdkWindowedVector() {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("float *Similarity::getMdkWindowedVector(void) ");

        // If a window is defined
        if (isDefinedWindow()) {
            // Check if the window has been applied
            if (MDK_Window == nullptr)
                applyWindow(halfWindow);
            // Return the windowed value
            return MDK_Window;
        }
            // Return the original values
        else return MDK;
    }

    bool Similarity::setSimilarityMatrix(similarityMatrix *sm) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("bool Similarity::setSimilarityMatrix(similarityMatrix *sm) ");

        // Checks if a similarity matrix is already being used.
        if (sm == nullptr)
            return false;

        if (simMatrix == sm)
            return true;

        delete simMatrix;

        simMatrix = sm;
        return true;
    }

    bool Similarity::isSimMatrixDef() {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("bool Similarity::isSimMatrixDef(void) ");

        return simMatrix != nullptr;
    }

    double Similarity::calcCutPoint(float baseLine, float conservationPct) {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("double Similarity::calcCutPoint(float baseLine, float conservationPct) ");
        // It computes the cutting point based on alignment's similarity values -
        // the so-called 'similarity'. It also takes into account the minimum percentage
        // from the input alignment to be kept. Depending on those two values, the
        // method will select a different cutting-point.

        double cuttingPoint_MinimumConserv, cuttingPoint_SimilThreshold;
        int i, highestPos;
        float *vectAux;

        vectAux = new float[alig->originalNumberOfResidues];

        // Sort a copy of the vector containing the similarity values after applying
        // any windows methods. Take the columns value that it lower than the minimum
        // similarity threshold set by the user
        utils::copyVect(getMdkWindowedVector(), vectAux, alig->originalNumberOfResidues);
        utils::quicksort(vectAux, 0, alig->originalNumberOfResidues - 1);

        for (i = alig->originalNumberOfResidues - 1; i >= 0; i--)
            if (vectAux[i] < conservationPct)
                break;
        cuttingPoint_SimilThreshold = vectAux[i];

        // It is possible that due to number casting, we get a number out of the
        // vector containing the similarity values - it is not reporting an overflow
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

    void Similarity::printConservationColumns() {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("void Similarity::printConservationColumns(void) ");

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

        fname = " Similarity per Column";

        std::cout << "#\33[0;32m Statistic :\33[0;1m" << fname << "\33[0m" << std::endl;

        std::cout << std::setw(alig->filename.size())
             << std::setfill('-')
             << std::left << ""
             << std::setfill(' ')
             << std::endl;

        std::cout << "\33[0;33;1m"
             << std::setw(size) << std::left << " Residue" << std::left << " Similarity" << std::endl
             << std::setw(size) << std::left << " Number" << std::left << " Value" << std::endl
             << std::setfill('-')
             << "\33[0;m"
             << std::setw(size) << std::right << "  "
             << std::setw(size) << std::right << "  " << std::endl
             << std::setfill(' ');

        std::cout.precision(10);

        float *values;

        // If MDK_Window vector is defined, we use it to print the similarity's values.
        if (MDK_Window != nullptr)
            values = MDK_Window;
            // In others cases, we uses the MDK vector to print the similarity's vlaues.
        else
            values = MDK;

        for (i = 0; i < alig->originalNumberOfResidues; i++)
            std::cout << std::setw(size) << std::left << i << values[i] << std::endl;
    }

    void Similarity::printConservationAcl() {
        // Create a timerLevel that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("void Similarity::printConservationAcl(void) ");

        float refer, *vectAux;
        int i, num, acm;
        int size = 20;
        // Allocate memory
        vectAux = new float[alig->originalNumberOfResidues];

        // Select the similarity's value source and copy that vector in a auxiliar vector
        if (MDK_Window != nullptr) utils::copyVect(MDK_Window, vectAux, alig->originalNumberOfResidues);
        else utils::copyVect(MDK, vectAux, alig->originalNumberOfResidues);

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

        fname = " Similarity Total";

        std::cout << "#\33[0;32m Statistic :\33[0;1m" << fname << "\33[0m" << std::endl;

        std::cout << std::setw(alig->filename.size())
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
        secondLine << std::setw(size) << std::left << " Similarity";
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

        // Count the columns with the same similarity's value and compute this information to shows the accumulative
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

