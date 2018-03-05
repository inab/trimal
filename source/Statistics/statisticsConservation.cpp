#include <TimerFactory.h>
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

    // Initializate the similarity matrix to nullptr.
    simMatrix = nullptr;

}

statisticsConservation::~statisticsConservation(void) {
	 // Create a timer that will report times upon its destruction
	 //	which means the end of the current scope.
	StartTiming("statisticsConservation::~statisticsConservation() ");

    // Deallocate memory, if it have been allocated previously.
    if (Q != nullptr) {
        delete[] Q;
        delete[] MDK;

        if (halfWindow > 0)
            delete[] MDK_Window;

        if (matrixIdentity != nullptr)
        for (int i = 0; i < _alignment->sequenNumber; i++)
            delete[] matrixIdentity[i];
        delete[] matrixIdentity;
    }
}

void statisticsConservation::calculateMatrixIdentity() {
	 // Create a timer that will report times upon its destruction
	 //	which means the end of the current scope.
	StartTiming("void statisticsConservation::calculateMatrixIdentity() ");

    char indet;
    int i, ii, j, jj, k, sum, length;

    matrixIdentity = new float *[sequences];
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

            // Calculate the value of matrixidn for columns j and i 
            matrixIdentity[jj][ii] = (100.0F - ((float) sum / length) * 100.0F);
            matrixIdentity[ii][jj] = matrixIdentity[jj][ii];

        }
    }

}

bool statisticsConservation::calculateVectors(int *gaps) {
	 // Create a timer that will report times upon its destruction
	 //	which means the end of the current scope.
	StartTiming("bool statisticsConservation::calculateVectors(int *gaps) ");

    // Calculation methods call
    calculateMatrixIdentity();

    char indet;
    int i, ii, j, jj, k, kk;
    float num, den;

    // Depending on alignment type, indetermination symbol will be one or other 
    indet = (_alignment->getAlignmentType() & SequenceTypes::AA) ? 'X' : 'N';
    // A conservation matrix must be defined. If not, return false 
    if (simMatrix == nullptr)
        return false;

    // For each column calculate the Q value and the MD value using an equation 
    for (i = 0, ii = -1; i < _alignment->originalResidNumber; i++) {
        if (_alignment->saveResidues[i] == -1) continue;
        ii++;
        // For each AAs/Nucleotides' pair in the column we compute its distance 
        if (gaps != nullptr)
            if (((float) gaps[ii] / sequences) >= 0.8F) {
                MDK[ii] = 0.F;
                continue;
            }

        for (j = 0, jj = -1, num = 0, den = 0; j < _alignment->originalSequenNumber; j++) {
            if (_alignment->saveSequences[j] == -1) continue;
            jj++;
            // We don't compute the distant if the first element is a indeterminate (X) or a gap (-) element. 
            if ((_alignment->sequences[j][i] != '-') && (_alignment->sequences[j][i] != indet))
                for (k = j + 1, kk = jj; k < _alignment->originalSequenNumber; k++) {
                    if (_alignment->saveSequences[k] == -1) continue;
                    kk++;
                    // We don't compute the distant between the pair if the second element is a indeterminate or a gap element 
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

        // If we are procesing a column with only one AA/nucleotide, the denominator is 0 and we don't execute the division
        // and we set the Q[i] value to 0. 
        Q[ii] = (den == 0) ? 0 : num / den;
        MDK[ii] = exp(-Q[ii]);

        // If the column has 80% or more gaps then we set its conservation value to 0 

        // If the MDK value is more than 1, we normalized this value to 1. 
        if (MDK[ii] > 1) MDK[ii] = 1;

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
    if (_halfWindow == 0)
    {
        if (halfWindow != 0)
            delete[] MDK_Window;

        MDK_Window = nullptr;
        return true;
    }
    if ((halfWindow == _halfWindow)) return true;

    if (_halfWindow > residues / 4)
    {
        debug.report(ErrorCode::SimilarityWindowTooBig);
        return false;
    }

    int i, j, window;
    delete[] MDK_Window;
    MDK_Window = new float[residues];

    // If one of this conditions is true, we return FALSE:                         
    //    .- If already exists a previously calculated vector for this window size 
    //    .- If mediumWinSize value is greater than 1/4 of alignment length

    halfWindow = _halfWindow;
    window = 2 * halfWindow + 1;

    // Do the average window calculations 
    for (i = 0; i < residues; i++) {
        MDK_Window[i] = 0.F;
        for (j = i - halfWindow; j <= i + halfWindow; j++) {
            if (j < 0) MDK_Window[i] += MDK[-j];
            else if (j >= residues) MDK_Window[i] += MDK[((2 * residues - j) - 2)];
            else MDK_Window[i] += MDK[j];
        }

        // Calculate the similiraty value for the i column 
        MDK_Window[i] = MDK_Window[i] / (float) window;
    }

    return true;
}

bool statisticsConservation::isDefinedWindow(void) {
	 // Create a timer that will report times upon its destruction
	 //	which means the end of the current scope.
	StartTiming("bool statisticsConservation::isDefinedWindow(void) ");

    return (halfWindow != -1);
}

float *statisticsConservation::getMdkwVector(void) {
	 // Create a timer that will report times upon its destruction
	 //	which means the end of the current scope.
	StartTiming("float *statisticsConservation::getMdkwVector(void) ");
    if (isDefinedWindow())
        return MDK_Window;
    else return MDK;
}

bool statisticsConservation::setSimilarityMatrix(similarityMatrix *sm) {
	 // Create a timer that will report times upon its destruction
	 //	which means the end of the current scope.
	StartTiming("bool statisticsConservation::setSimilarityMatrix(similarityMatrix *sm) ");

    // Checks if a similarity matrix is being used. 
    if (sm == nullptr)
        return false;

    // if a similarity matrix isn't being used, we associate a pointer gives as input parameter to object simMatrix's
    // pointer and return true. 
    simMatrix = sm;
    return true;
}

bool statisticsConservation::isSimMatrixDef(void) {
	 // Create a timer that will report times upon its destruction
	 //	which means the end of the current scope.
	StartTiming("bool statisticsConservation::isSimMatrixDef(void) ");

    return simMatrix != nullptr;
}

double statisticsConservation::calcCutPoint(float baseLine, float conservationPct) {
	 // Create a timer that will report times upon its destruction
	 //	which means the end of the current scope.
	StartTiming("double statisticsConservation::calcCutPoint(float baseLine, float conservationPct) ");
    /* It computes the cutting point based on alignment's conservation values -
      * the so-called 'similarity'. It also takes into account the minimum percentage
      * from the input alignment to be kept. Depending on those two values, the
      * method will select a different cutting-point. */

    double cuttingPoint_MinimumConserv, cuttingPoint_SimilThreshold;
    int i, highestPos;
    float *vectAux;

    /* Allocate memory */
    vectAux = new float[residues];

    /* Sort a copy of the vector containing the similarity values after applying
     * any windows methods. Take the columns value that it lower than the minimum
     * similarity threshold set by the user */
    utils::copyVect(getMdkwVector(), vectAux, residues);
    utils::quicksort(vectAux, 0, residues-1);

    for(i = residues - 1; i >= 0; i--)
        if(vectAux[i] < conservationPct)
            break;
    cuttingPoint_SimilThreshold = vectAux[i];

    /* It is possible that due to number casting, we get a number out of the
     * vector containing the similarity values - it is not reporting an overflow
     * situation but giving back a 0 when it should be a number equal (or closer)
     * to 1. */
    highestPos = (int) ((double)(residues - 1) * (100.0 - baseLine)/100.0);
    highestPos = highestPos < (residues - 1) ? highestPos : residues - 1;
    cuttingPoint_MinimumConserv = vectAux[highestPos];

    /* Deallocate memory */
    delete[] vectAux;

    /* Return the minimum cutting point between the one set by the threshold and
     * the one set by the minimum percentage of the input alignment to be kept */
    return (cuttingPoint_MinimumConserv < cuttingPoint_SimilThreshold ?
            cuttingPoint_MinimumConserv : cuttingPoint_SimilThreshold);
}

void statisticsConservation::printConservationColumns(void) {
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

//     cout << setiosflags(std::ios_base::width(14));

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

void statisticsConservation::printConservationAcl(void) {
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
