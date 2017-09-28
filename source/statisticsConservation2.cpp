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

#include "../include/statisticsConservation2.h"
#include "../include/defines.h"
#include "../include/newAlignment.h"
#include "../include/reportsystem.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  statisticsConservation2::statisticsConservation2(char **, int, int)                                                   |
|                                                                                                                      |
|       Class constructor. This method uses the inputs parameters to put the information in the new object that        |
|       has been created.                                                                                              |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

// statisticsConservation2::statisticsConservation2(string *alignmentMatrix, int species, int aminos, int dataType_) {
//
//   /* Initializate values to its corresponds values */
//   columns = aminos;
//   sequences = species;
//   dataType = dataType_;
//   halfWindow = -1;
//
//   /* Allocate memory to the structures and initializates it */
//   Q = new float[columns];
//   utils::initlVect(Q, columns, 0);
//
//   MDK = new float[columns];
//   utils::initlVect(MDK, columns, 0);
//
//   MDK_Window = new float[columns];
//   utils::initlVect(MDK_Window, columns, 0);
//
//   matrixIdentity = new float*[sequences];
//   for(int i = 0; i < sequences; i++){
//     matrixIdentity[i] = new float[sequences];
//     utils::initlVect(matrixIdentity[i], sequences, 0);
//   }
//
//   /* Initializate the similarity matrix to NULL. */
//   simMatrix = NULL;
//
//   /* Calculation methods call */
//   calculateMatrixIdentity(alignmentMatrix);
// }

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  statisticsConservation2::statisticsConservation2(void)                                                                |
|                                                                                                                      |
|       Class constructor.                                                                                             |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

statisticsConservation2::statisticsConservation2(newAlignment * parentAlignment) {

    _alignment = parentAlignment;

    residues = _alignment->residNumber;

    Q = new float[residues];
    utils::initlVect(Q, residues, 0);

    MDK = new float[residues];
    utils::initlVect(MDK, residues, 0);

    MDK_Window = new float[residues];
    utils::initlVect(MDK_Window, residues, 0);

    sequences = _alignment->sequenNumber;

    matrixIdentity = new float*[sequences];
    for(int i = 0; i < sequences; i++) {
        matrixIdentity[i] = new float[sequences];
        utils::initlVect(matrixIdentity[i], sequences, 0);
    }

    /* Initializate the similarity matrix to NULL. */
    simMatrix = NULL;

    /* Calculation methods call */
    calculateMatrixIdentity();

//     calculateVectors(NULL);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  statisticsConservation2::~statisticsConservation2(void)                                                               |
|                                                                                                                      |
|       Class destroyer.                                                                                               |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

statisticsConservation2::~statisticsConservation2(void) {

//   /* Deallocate memory, if it have been allocated previously. */
    if(Q != NULL) {
        delete[] Q;
        delete[] MDK;
        delete[] MDK_Window;

        for(int i = 0; i < _alignment->sequenNumber; i++)
            delete[] matrixIdentity[i];
        delete[] matrixIdentity;
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  void statisticsConservation2::calculateMatrixIdentity(char **, int, int)                                             |
|                                                                                                                      |
|       This method computes the matrix identity between all the sequences in the alignment.                           |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void statisticsConservation2::calculateMatrixIdentity() {

    char indet;
    int i, ii, j, jj, k, sum, length;

    /* Depending on alignment type, indetermination symbol will be one or other */
    indet = (_alignment->getAlignmentType() & SequenceTypes::AA) ? 'X' : 'N';

    /* For each sequences' pair */
    for(i = 0, ii = -1; i < _alignment->originalSequenNumber; i++) {
        if (_alignment->saveSequences[i] == -1) continue;
        ii++;
        for(j = i + 1, jj = ii; j < _alignment->originalSequenNumber; j++) {
            if (_alignment->saveSequences[j] == -1) continue;
            jj++;
            /* For each position in the alignment of that pair than we are processing */
            for(k = 0, sum = 0, length = 0; k < _alignment->originalResidNumber; k++) {
                if (_alignment->saveResidues[k] == -1) continue;

                /* If we find a element that is not a gap or an X aminoacid in the first sequence of the pair */
                if((_alignment->sequences[i][k] != '-') && (_alignment->sequences[i][k] != indet)) {

                    /* If we also find a valid element in the second sequence  */
                    if((_alignment->sequences[j][k] != '-') && (_alignment->sequences[j][k] != indet))

                        /* If the two valid elements are the same increase the sum */
                        if(_alignment->sequences[j][k] ==  _alignment->sequences[i][k])
                            sum++;

                    /* Increase the length of the sequence free of gaps and X elements */
                    length++;
                }

                /* If the first processed element is invalid and in the second we find a valid element increase the length of
                   the sequence free of gaps and X elements */
                else if((_alignment->sequences[j][k] != '-') && (_alignment->sequences[j][k] != indet))
                    length++;
            }

            /* Calculate the value of matrixidn for columns j and i */
            matrixIdentity[jj][ii] = (100.0 - ((float) sum / length) * 100.0);
            matrixIdentity[ii][jj] = matrixIdentity[jj][ii];

        }
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  bool statisticsConservation2::calculateVectors(char **, int *)                                                       |
|                                                                                                                      |
|       This method computes the distance between pairs for each column in the alignment.                              |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

bool statisticsConservation2::calculateVectors(int *gaps) {

    char indet;
    int i, ii, j, jj, k, kk;
    float num, den;

    /* Depending on alignment type, indetermination symbol will be one or other */
    indet = (_alignment->getAlignmentType() & SequenceTypes::AA) ? 'X' : 'N';
    /* A conservation matrix must be defined. If not, return false */
    if(simMatrix == NULL)
        return false;

    /* For each column calculate the Q value and the MD value using an equation */
    for(i = 0, ii = -1; i < _alignment->originalResidNumber; i++) {
        if (_alignment->saveResidues[i] == -1) continue;
        ii ++;
        /* For each AAs/Nucleotides' pair in the column we compute its distance */
        if(_alignment->sgaps->getGapsWindow() != NULL)
            if(((float) _alignment->sgaps->gapsWindow[ii] / sequences) >= 0.8) 
            {
                MDK[ii] = 0;
                continue;
            }
            
        for(j = 0, jj = -1, num = 0, den = 0; j < _alignment->originalSequenNumber; j++) {
            if (_alignment->saveSequences[j] == -1) continue;
            jj++;
            /* We don't compute the distant if the first element is a indeterminate (X) or a gap (-) element. */
            if((_alignment->sequences[j][i] != '-') && (_alignment->sequences[j][i] != indet))
                for(k = j + 1, kk = jj; k < _alignment->originalSequenNumber; k++)
                {
                    if (_alignment->saveSequences[k] == -1) continue;
                    kk++;
                    /* We don't compute the distant between the pair if the second element is a indeterminate or a gap element */
                    if((_alignment->sequences[k][i] != '-') && (_alignment->sequences[k][i] != indet)) {
                        /* We use the identity value for the two pairs and its distance based on similarity matrix's value. */
                        num += matrixIdentity[jj][kk] * simMatrix -> getDistance(_alignment->sequences[j][i], _alignment->sequences[k][i]);;
                        den += matrixIdentity[jj][kk];
                    }
                }
        }
        
        /* If we are procesing a column with only one AA/nucleotide, the denominator is 0 and we don't execute the division
           and we set the Q[i] value to 0. */
        Q[ii] = (den == 0) ? 0 : num / den;
        MDK[ii] = (float) exp(-Q[ii]);

        /* If the column has 80% or more gaps then we set its conservation value to 0 */

        /* If the MDK value is more than 1, we normalized this value to 1. */
        if(MDK[ii] > 1) MDK[ii] = 1;
//         Debug  << setw(20) << left << i << " "  << setw(20) << left << Q[i] << " "  << setw(20) << left << num << " "  << setw(20) << left << den << endl;

    }

    return true;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  bool statisticsConservation2::applyWindow(int)                                                                       |
|                                                                                                                      |
|       This method computes for each column's alignment its conservationwindows' value. For this purpose, the method  |
|       uses the values that previously has been calculated and the window's size value.                               |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

bool statisticsConservation2::applyWindow(int _halfWindow) {

    int i, j, window;

    /* If one of this conditions is true, we return FALSE:                         */
    /*    .- If already exists a previously calculated vector for this window size */
    /*    .- If mediumWinSize value is greater than 1/4 of alignment length        */
    if((halfWindow == _halfWindow) || (_halfWindow > residues/4))
        return false;

    halfWindow = _halfWindow;
    window = 2 * halfWindow + 1;

    /* Do the average window calculations */
    for(i = 0; i < residues; i++) {
        for(j = i - halfWindow; j <= i + halfWindow; j++) {
            if(j < 0) MDK_Window[i] += MDK[-j];
            else if(j >= residues) MDK_Window[i] += MDK[((2 * residues - j) - 2)];
            else MDK_Window[i] += MDK[j];
        }

        /* Calculate the similiraty value for the i column */
        MDK_Window[i] = MDK_Window[i] / (float) window;
    }
    return true;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  bool statisticsConservation2::isDefinedWindow(void)                                                                  |
|                                                                                                                      |
|       This method returns true if a similarity matrix has been defined and false in others cases.                    |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

bool statisticsConservation2::isDefinedWindow(void) {

    return (halfWindow != -1);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  bool statisticsConservation2::getMdkwVector(void)                                                                    |
|                                                                                                                      |
|       This method returns a pointer to conservation values' vector.                                                  |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

float *statisticsConservation2::getMdkwVector(void) {

    return MDK_Window;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  bool statisticsConservation2::setSimilarityMatrix(similarityMatrix *)                                                |
|                                                                                                                      |
|       This method associated a pointer to similarity matrix gives as input parameter. If a conservation matrix is    |
|       being used the methods return false and doesn't do anything.                                                   |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

bool statisticsConservation2::setSimilarityMatrix(similarityMatrix *sm) {

    /* Checks if a similarity matrix is being used. */
    if(sm == NULL)
        return false;

    /* if a similarity matrix isn't being used, we associate a pointer gives as input parameter to object simMatrix's
       pointer and return true. */
    simMatrix = sm;
    return true;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  bool statisticsConservation2::isSimMatrixDef(void)                                                                   |
|                                                                                                                      |
|       This method returns true if a similarity matrix is being used and false in others cases.                       |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

bool statisticsConservation2::isSimMatrixDef(void) {

    return (simMatrix != NULL);
}

double statisticsConservation2::calcCutPoint(float baseLine, float conservationPct) {
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

    utils::copyVect(MDK_Window, vectAux, residues);

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
    return std::min(cuttingPoint_MinimumConserv , cuttingPoint_SimilThreshold);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  void statisticsConservation2::printConservationColumns(void)                                                         |
|                                                                                                                      |
|       This method prints the conservation's value for each column in the alignment.                                  |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void statisticsConservation2::printConservationColumns(void) {

    int i;

    /* We set the output precision and print the header. */
    cout << "| Residue\t Similarity  |" << endl;
    cout << "| Number \t    Value    |" << endl;
    cout << "+----------------------------+" << endl;
    cout.precision(10);

    /* If MDK_Window vector is defined, we use it to print the conservation's values. */
    if(MDK_Window != NULL)
        for(i = 0; i < residues; i++)
            cout << "  " << setw(5) << i << "\t\t" << setw(7) << MDK_Window[i] << endl;

    /* In others cases, we uses the MDK vector to print the conservation's vlaues. */
    else
        for(i = 0; i < residues; i++)
            cout << "  " << setw(5) << i << "\t\t" << setw(7) << MDK[i] << endl;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
|  void statisticsConservation2::printConservationAcl(void)                                                             |
|                                                                                                                      |
|       This method prints the accumulative statistics related to conservation in the alignment.                       |
|                                                                                                                      |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void statisticsConservation2::printConservationAcl(void) {

    float refer, *vectAux;
    int i, num, acm;

    /* Allocate memory */
    vectAux = new float[residues];

    /* Select the conservation's value source and copy that vector in a auxiliar vector */
    if(MDK_Window != NULL) utils::copyVect(MDK_Window, vectAux, residues);
    else utils::copyVect(MDK, vectAux, residues);

    /* Sort the auxiliar vector. */
    utils::quicksort(vectAux, 0, residues-1);

    /* We set the output precision and print the header. */
    cout << "| Number of\t        \t|\t Cumulative \t% Cumulative\t|   Similarity   |" << endl;
    cout << "| Residues \t% Length\t|\tNumberResid.\t   Length   \t|     Value      |" << endl;
    cout << "+-------------------------------+---------------------------------------+----------------+" << endl;
    cout.precision(10);


    /* Initializate some values */
    refer = vectAux[residues-1];
    acm = 0;
    num = 1;

    /* Count the columns with the same conservation's value and compute this information to shows the accunulative
       statistics in the alignment. */
    for(i = residues-2; i >= 0; i--) {
        acm++;

        if(refer != vectAux[i]) {
            cout << "  " << num << "\t\t" << setw(10) << ((float) num/residues * 100.0)
                 << "\t\t" << acm << "\t\t" << setw(10) << ((float) acm/residues * 100.0) << "\t"
                 << setw(15) << refer << endl;
            refer = vectAux[i];
            num = 1;
        }
        else num++;
    }
    acm++;
    cout << "  " << num << "\t\t" << setw(10) << ((float) num/residues * 100.0)
         << "\t\t" << acm << "\t\t" << setw(10) << ((float) acm/residues * 100.0) << "\t"
         << setw(15) << refer << endl;

    /* Deallocate the reserved memory. */
    delete [] vectAux;
}
