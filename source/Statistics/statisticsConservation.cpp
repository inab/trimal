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

#include "Statistics/statisticsConservation.h"
#include "defines.h"
#include "newAlignment.h"
#include "reportsystem.h"
#include <sstream>
#include <TimerFactory.h>

statisticsConservation::statisticsConservation(newAlignment *parentAlignment) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("statisticsConservation::statisticsConservation(newAlignment *parentAlignment) ");

    _alignment = parentAlignment;

    residues = _alignment->originalResidNumber;

    Q = new float[residues];
    utils::initlVect(Q, residues, 0);

    MDK = new float[residues];
    utils::initlVect(MDK, residues, 0);

    sequences = _alignment->sequenNumber;

    // Initialize the similarity matrix to nullptr.
    simMatrix = nullptr;

    refCounter = new int(1);
}

statisticsConservation::statisticsConservation(newAlignment *parentAlignment,
                                               statisticsConservation *mold) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("statisticsConservation::statisticsConservation(newAlignment *parentAlignment) ");

    _alignment = parentAlignment;

    residues = _alignment->originalResidNumber;

    sequences = _alignment->sequenNumber;

    halfWindowApplied = -1;

    halfWindowRequested = mold->halfWindowRequested;

    Q = mold->Q;
    MDK = mold->MDK;

    simMatrix = mold->simMatrix;

    refCounter = mold->refCounter;
    (*refCounter)++;
}

statisticsConservation::~statisticsConservation() {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("statisticsConservation::~statisticsConservation() ");

    // Deallocate common variables only in case there is no module that has a
    // reference to them
    if (--(*refCounter) == 0) {
        delete[] Q;
        delete[] MDK;

        if (matrixIdentity != nullptr)
            for (int i = 0; i < _alignment->sequenNumber; i++)
                delete[] matrixIdentity[i];

        delete[] matrixIdentity;
        delete refCounter;
    }

    // We always want to delete the windows values as they are related to the
    // specific alignment and not the whole set of derived alignments
    delete[] MDK_Window;
}

void statisticsConservation::calculateMatrixIdentity() {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("void statisticsConservation::calculateMatrixIdentity() ");

    // We don't want to calculate the matrix identity
    // if it has been previously calculated
    if (matrixIdentity != nullptr)
        return;

    // Allocate temporal variables
    char indet;
    int i, ii, j, jj, k, sum, length;

    // Allocate memory for the matrix identity
    matrixIdentity = new float *[sequences + 1];
    for (i = 0; i < sequences; i++) {
        matrixIdentity[i] = new float[sequences];
        utils::initlVect(matrixIdentity[i], sequences, 0);
    }

    // Depending on alignment type, indetermination symbol will be one or other 
    indet = (_alignment->getAlignmentType() & SequenceTypes::AA) ? 'X' : 'N';

    // For each sequences' pair 
    for (i = 0, ii = -1; i < _alignment->originalSequenNumber; i++) {
        if (_alignment->saveSequences[i] == -1) continue;
        ii++;
        for (j = i + 1, jj = ii; j < _alignment->originalSequenNumber; j++) {
            if (_alignment->saveSequences[j] == -1) continue;
            jj++;
            // For each position in the alignment of that pair than we are processing 
            for (k = 0, sum = 0, length = 0; k < _alignment->originalResidNumber; k++) {
                if (_alignment->saveResidues[k] == -1) continue;

                // If we find a element that is not a gap or an X aminoacid in the first sequence of the pair 
                if ((_alignment->sequences[i][k] != '-') && (_alignment->sequences[i][k] != indet)) {

                    // If we also find a valid element in the second sequence  
                    if ((_alignment->sequences[j][k] != '-') && (_alignment->sequences[j][k] != indet))

                        // If the two valid elements are the same increase the sum 
                        if (_alignment->sequences[j][k] == _alignment->sequences[i][k])
                            sum++;

                    // Increase the length of the sequence free of gaps and X elements 
                    length++;
                }

                    // If the first processed element is invalid and in the second we find a valid element increase the length of
                    // the sequence free of gaps and X elements
                else if ((_alignment->sequences[j][k] != '-') && (_alignment->sequences[j][k] != indet))
                    length++;
            }

            // Calculate the value of matrix idn for columns j and i
            matrixIdentity[jj][ii] = (100.0F - ((float) sum / length) * 100.0F);
            matrixIdentity[ii][jj] = matrixIdentity[jj][ii];

        }
    }
}

bool statisticsConservation::calculateVectors(bool cutByGap) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("bool statisticsConservation::calculateVectors(int *gaps) ");

    // A conservation matrix must be defined. If not, return false
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
        gaps = _alignment->Statistics->gaps->getGapsWindow();

    // Initialize the variables used
    char indet;
    int i, j, jj, k, kk;
    float num, den;

    // Depending on alignment type, indetermination symbol will be one or other 
    indet = (_alignment->getAlignmentType() & SequenceTypes::AA) ? 'X' : 'N';

    // For each column calculate the Q value and the MD value using an equation 
    for (i = 0; i < _alignment->originalResidNumber; i++) {
        if (_alignment->saveResidues[i] == -1) continue;

        // Set MDK for columns with gaps values bigger or equal to 0.8F
        if (cutByGap) {
            if ((float) gaps[i] / sequences >= 0.8F) {
                MDK[i] = 0.F;
                continue;
            }
        }

        // For each AAs/Nucleotides' pair in the column we compute its distance
        for (j = 0, jj = -1, num = 0, den = 0; j < _alignment->originalSequenNumber; j++) {
            if (_alignment->saveSequences[j] == -1) continue;
            jj++;
            // We don't compute the distance if the first element is a indeterminate (X) or a gap (-) element.
            if ((_alignment->sequences[j][i] != '-') && (_alignment->sequences[j][i] != indet))
                for (k = j + 1, kk = jj; k < _alignment->originalSequenNumber; k++) {
                    if (_alignment->saveSequences[k] == -1) continue;
                    kk++;
                    // We don't compute the distance between the pair if the second element is a indeterminate or a gap element
                    if ((_alignment->sequences[k][i] != '-') && (_alignment->sequences[k][i] != indet)) {
                        // We use the identity value for the two pairs and its distance based on similarity matrix's value. 
                        num += matrixIdentity[jj][kk] *
                               simMatrix->getDistance(
                                       _alignment->sequences[j][i],
                                       _alignment->sequences[k][i]
                               );
                        den += matrixIdentity[jj][kk];
                    }
                }
        }

        // If we are processing a column with only one AA/nucleotide,
        // the denominator is 0 and we don't execute the division
        // and we set the Q[i] value to 0. 
        Q[i] = (den == 0) ? 0 : num / den;
        MDK[i] = exp(-Q[i]);

        // If the MDK value is more than 1, we normalized this value to 1. 
        if (MDK[i] > 1) MDK[i] = 1;

    }

    for (i = 0; i < _alignment->sequenNumber; i++)
        delete[] matrixIdentity[i];
    delete[] matrixIdentity;
    matrixIdentity = nullptr;

    return true;
}

bool statisticsConservation::applyWindow(int _halfWindow) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("bool statisticsConservation::applyWindow(int _halfWindow) ");

    // Calculate the MDK array if it has not been calculated previously
    if (MDK == nullptr)
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
            delete[] MDK_Window;

        MDK_Window = nullptr;
        return true;
    }

    // Initialize the values used in the calculation
    int i, j, window;

    // Initialize the MDK window array if it's null
    if (MDK_Window == nullptr)
        MDK_Window = new float[residues + 1];


    halfWindowApplied = _halfWindow;
    window = 2 * halfWindowApplied + 1;

    // Do the average window calculations 
    for (i = 0; i < residues; i++) {
        MDK_Window[i] = 0.F;
        for (j = i - halfWindowApplied; j <= i + halfWindowApplied; j++) {
            if (j < 0) MDK_Window[i] += MDK[-j];
            else if (j >= residues) MDK_Window[i] += MDK[((2 * residues - j) - 2)];
            else MDK_Window[i] += MDK[j];
        }

        // Calculate the average value, by dividing the values
        MDK_Window[i] = MDK_Window[i] / (float) window;
    }

    return true;
}

bool statisticsConservation::isDefinedWindow() {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("bool statisticsConservation::isDefinedWindow(void) ");

    return (halfWindowRequested > 0);
}

float *statisticsConservation::getMdkWindowedVector() {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("float *statisticsConservation::getMdkWindowedVector(void) ");

    // If a window is defined
    if (isDefinedWindow()) {
        // Check if the window has been applied
        if (halfWindowRequested != halfWindowApplied)
            applyWindow(halfWindowRequested);
        // Return the windowed value
        return MDK_Window;
    }
    // Return the original values
    else return MDK;
}

bool statisticsConservation::setSimilarityMatrix(similarityMatrix *sm) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("bool statisticsConservation::setSimilarityMatrix(similarityMatrix *sm) ");

    // Checks if a similarity matrix is already being used.
    if (sm == nullptr)
        return false;

    if (simMatrix == sm)
        return true;

    delete simMatrix;

    simMatrix = sm;
    return true;
}

bool statisticsConservation::isSimMatrixDef() {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("bool statisticsConservation::isSimMatrixDef(void) ");

    return simMatrix != nullptr;
}

double statisticsConservation::calcCutPoint(float baseLine, float conservationPct) {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("double statisticsConservation::calcCutPoint(float baseLine, float conservationPct) ");
    // It computes the cutting point based on alignment's conservation values -
    // the so-called 'similarity'. It also takes into account the minimum percentage
    // from the input alignment to be kept. Depending on those two values, the
    // method will select a different cutting-point.

    double cuttingPoint_MinimumConserv, cuttingPoint_SimilThreshold;
    int i, highestPos;
    float *vectAux;

    vectAux = new float[residues];

    // Sort a copy of the vector containing the similarity values after applying
    // any windows methods. Take the columns value that it lower than the minimum
    // similarity threshold set by the user
    utils::copyVect(getMdkWindowedVector(), vectAux, residues);
    utils::quicksort(vectAux, 0, residues - 1);

    for (i = residues - 1; i >= 0; i--)
        if (vectAux[i] < conservationPct)
            break;
    cuttingPoint_SimilThreshold = vectAux[i];

    // It is possible that due to number casting, we get a number out of the
    // vector containing the similarity values - it is not reporting an overflow
    // situation but giving back a 0 when it should be a number equal (or closer)
    // to 1.
    highestPos = (int) ((double) (residues - 1) * (100.0 - baseLine) / 100.0);
    highestPos = highestPos < (residues - 1) ? highestPos : residues - 1;
    cuttingPoint_MinimumConserv = vectAux[highestPos];

    delete[] vectAux;

    // Return the minimum cutting point between the one set by the threshold and
    // the one set by the minimum percentage of the input alignment to be kept
    return (cuttingPoint_MinimumConserv < cuttingPoint_SimilThreshold ?
            cuttingPoint_MinimumConserv : cuttingPoint_SimilThreshold);
}

void statisticsConservation::printConservationColumns() {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("void statisticsConservation::printConservationColumns(void) ");

    int i, size = 20;

    std::string fname = _alignment->filename.substr(6, _alignment->filename.size() - 7);


    cout << std::fixed
         << std::setw(fname.length() + 7)
         << std::setfill(' ')
         << std::left << "" << endl;

    cout << "#\33[0;31m File :\33[0;1m" << fname << "\33[0m";

    fname = std::to_string(size);

    cout
            << std::setw(fname.length() + 7)
            << std::setfill(' ')
            << std::left << "" << endl;

    cout << "#\33[0;36m BlockSize : \33[0;1m" << fname << "\33[0m" << endl;

    fname = " Similarity per Column";

    cout << "#\33[0;32m Statistic :\33[0;1m" << fname << "\33[0m" << endl;

    cout << std::setw(_alignment->filename.substr(6, _alignment->filename.size() - 7).length() + 7)
         << std::setfill('-')
         << std::left << ""
         << std::setfill(' ')
         << endl;

    cout << "\33[0;33;1m"
         << std::setw(size) << std::left << " Residue" << std::left << " Similarity" << endl
         << std::setw(size) << std::left << " Number" << std::left << " Value" << endl
         << std::setfill('-')
         << "\33[0;m"
         << std::setw(size) << std::right << "  "
         << std::setw(size) << std::right << "  " << endl
         << std::setfill(' ');

    cout.precision(10);

    float *values;

    // If MDK_Window vector is defined, we use it to print the conservation's values. 
    if (MDK_Window != nullptr)
        values = MDK_Window;
        // In others cases, we uses the MDK vector to print the conservation's vlaues.
    else
        values = MDK;

    for (i = 0; i < residues; i++)
        cout << setw(size) << std::left << i << values[i] << endl;
}

void statisticsConservation::printConservationAcl() {
    // Create a timer that will report times upon its destruction
    //	which means the end of the current scope.
    StartTiming("void statisticsConservation::printConservationAcl(void) ");

    float refer, *vectAux;
    int i, num, acm;
    int size = 20;
    // Allocate memory 
    vectAux = new float[residues];

    // Select the conservation's value source and copy that vector in a auxiliar vector 
    if (MDK_Window != nullptr) utils::copyVect(MDK_Window, vectAux, residues);
    else utils::copyVect(MDK, vectAux, residues);

    // Sort the auxiliar vector. 
    utils::quicksort(vectAux, 0, residues - 1);

    // Print filename
    std::string fname = _alignment->filename.substr(6, _alignment->filename.size() - 7);

    cout << std::fixed
         << std::setw(fname.length() + 7)
         << std::setfill(' ')
         << std::left << "" << endl;

    cout << "#\33[0;31m File :\33[0;1m" << fname << "\33[0m";

    fname = std::to_string(size);

    cout
            << std::setw(fname.length() + 7)
            << std::setfill(' ')
            << std::left << "" << endl;

    cout << "#\33[0;36m BlockSize : \33[0;1m" << fname << "\33[0m" << endl;

    fname = " Similarity Total";

    cout << "#\33[0;32m Statistic :\33[0;1m" << fname << "\33[0m" << endl;

    cout << std::setw(_alignment->filename.substr(6, _alignment->filename.size() - 7).length() + 7)
         << std::setfill('-')
         << std::left << ""
         << std::setfill(' ')
         << endl;


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

    cout << "\33[0;33;1m"
         << firstLine.rdbuf() << endl
         << secondLine.rdbuf() << endl
         << thirdLine.rdbuf() << endl
         << "\33[0;m"
         << std::setfill('-');

    for (i = 0; i < 5; i++)
        cout << setw(size) << std::right << "   ";

    cout << endl << setfill(' ');
    cout.precision(10);


    // Initializate some values 
    refer = vectAux[residues - 1];
    acm = 0;
    num = 1;

    // Count the columns with the same conservation's value and compute this information to shows the accumulative
    // statistics in the alignment. 
    for (i = residues - 2; i >= 0; i--) {
        acm++;

        if (refer != vectAux[i]) {

            cout
                    << setw(size) << std::left << num

                    << setw(size) << std::left
                    << setw(size - 6) << std::right << ((float) num / residues * 100.0F)
                    << setw(6) << std::right << " "

                    << setw(size) << std::left << acm

                    << setw(size) << std::left
                    << setw(size - 6) << std::right << ((float) acm / residues * 100.0F)
                    << setw(6) << std::right << " "

                    << setw(size) << std::left << refer

                    << endl;
            refer = vectAux[i];
            num = 1;
        } else num++;
    }
    acm++;

    cout
            << setw(size) << std::left << num

            << setw(size) << std::left
            << setw(size - 6) << std::right << ((float) num / residues * 100.0F)
            << setw(6) << std::right << " "

            << setw(size) << std::left << acm

            << setw(size) << std::left
            << setw(size - 6) << std::right << ((float) acm / residues * 100.0F)
            << setw(6) << std::right << " "

            << setw(size) << std::left << refer

            << endl;

    // Deallocate the reserved memory. 
    delete[] vectAux;
}
