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

#include "../include/sequencesMatrix.h"
#include "../include/TimerFactory.h"
#include "../include/newAlignment.h"
#include "../include/utils.h"

sequencesMatrix::sequencesMatrix(void) {
	 // Create a timer that will report times upon its destruction
	 //	which means the end of the current scope.
	StartTiming("sequencesMatrix::sequencesMatrix(void) ");

    resNumber = 0;
    seqsNumber = 0;

    seqsName = nullptr;
    matrix = nullptr;

}

sequencesMatrix::sequencesMatrix(newAlignment *parent) {
	 // Create a timer that will report times upon its destruction
	 //	which means the end of the current scope.
	StartTiming("sequencesMatrix::sequencesMatrix(newAlignment *parent) ");
    _alignment = parent;
    int i, j, k;

    seqsNumber = _alignment->originalSequenNumber;
    resNumber = _alignment->originalResidNumber;

    seqsName = _alignment->seqsName;

    matrix = new int *[seqsNumber];
    for (i = 0; i < seqsNumber; i++) {
        matrix[i] = new int[resNumber];
        utils::initlVect(matrix[i], resNumber, 0);
    }

    // Determinate the sequence for each alignment specie
    for (i = 0, k = 1; i < seqsNumber; i++, k = 1) {
        for (j = 0; j < resNumber; j++) {
            if (_alignment->sequences[i][j] != '-') {
                matrix[i][j] = k;
                k++;
            }
        }
    }

}


sequencesMatrix::sequencesMatrix(string *alignmentMatrix, string *alignmentSeqsName, int sequences, int residues) {
	 // Create a timer that will report times upon its destruction
	 //	which means the end of the current scope.
	StartTiming("sequencesMatrix::sequencesMatrix(string *alignmentMatrix, string *alignmentSeqsName, int sequences, int residues) ");
    int i, j, k;

    seqsNumber = sequences;
    resNumber = residues;

    seqsName = new string[seqsNumber];
    for (i = 0; i < seqsNumber; i++)
        seqsName[i] = alignmentSeqsName[i];


    matrix = new int *[seqsNumber];
    for (i = 0; i < seqsNumber; i++) {
        matrix[i] = new int[resNumber];
        utils::initlVect(matrix[i], resNumber, 0);
    }

    // Determinate the sequence for each alignment specie
    for (i = 0, k = 1; i < seqsNumber; i++, k = 1) {
        for (j = 0; j < resNumber; j++) {
            if (alignmentMatrix[i][j] != '-') {
                matrix[i][j] = k;
                k++;
            }
        }
    }
}

sequencesMatrix &sequencesMatrix::operator=(const sequencesMatrix &old) {
    int i, j;

    if (this != &old) {

        seqsNumber = old.seqsNumber;
        resNumber = old.resNumber;

        seqsName = new string[seqsNumber];
        for (i = 0; i < seqsNumber; i++)
            seqsName[i] = old.seqsName[i];

        matrix = new int *[seqsNumber];
        for (i = 0; i < seqsNumber; i++) {
            matrix[i] = new int[resNumber];
            for (j = 0; j < resNumber; j++)
                matrix[i][j] = matrix[i][j];
        }

    }
    return *this;
}

sequencesMatrix::~sequencesMatrix() {
	 // Create a timer that will report times upon its destruction
	 //	which means the end of the current scope.
	StartTiming("sequencesMatrix::~sequencesMatrix(void) ");
    int i;

    if (matrix != nullptr) {
        for (i = 0; i < seqsNumber; i++)
            delete matrix[i];
        delete[] matrix;
    }

    seqsNumber = 0;
    resNumber = 0;

    matrix = nullptr;
    seqsName = nullptr;
}

void sequencesMatrix::printMatrix(void) {
	 // Create a timer that will report times upon its destruction
	 //	which means the end of the current scope.
	StartTiming("void sequencesMatrix::printMatrix(void) ");
    int i, j, k;

    for (i = 0; i < resNumber; i += 20) {
        for (j = 0; j < seqsNumber; j++) {
            for (k = i; k < (20 + i) && k < resNumber; k++) {
                cout << setw(4) << matrix[j][k] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
}

void sequencesMatrix::getColumn(int column, int *columnSeqMatrix) {
	 // Create a timer that will report times upon its destruction
	 //	which means the end of the current scope.
	StartTiming("void sequencesMatrix::getColumn(int column, int *columnSeqMatrix) ");
    int i;

    if (column < resNumber)
        for (i = 0; i < seqsNumber; i++)
            columnSeqMatrix[i] = matrix[i][column];

    else
        for (i = 0; i < seqsNumber; i++)
            columnSeqMatrix[i] = 0;

}

void sequencesMatrix::getColumn(int value, int row, int *columnSeqMatrix) {
	 // Create a timer that will report times upon its destruction
	 //	which means the end of the current scope.
	StartTiming("void sequencesMatrix::getColumn(int value, int row, int *columnSeqMatrix) ");
    int i, j;

    for (i = 0; i < resNumber; i++)
        if (matrix[row][i] == value) break;

    if (i < resNumber)
        for (j = 0; j < seqsNumber; j++)
            columnSeqMatrix[j] = matrix[j][i];

    else
        for (j = 0; j < seqsNumber; j++)
            columnSeqMatrix[j] = -1;
}

void sequencesMatrix::setOrder(int *order) {
	 // Create a timer that will report times upon its destruction
	 //	which means the end of the current scope.
	StartTiming("void sequencesMatrix::setOrder(int *order) ");
    int i, j, **resg;

    resg = new int *[seqsNumber];
    for (i = 0; i < seqsNumber; i++)
        resg[i] = new int[resNumber];

    for (i = 0; i < seqsNumber; i++)
        for (j = 0; j < resNumber; j++)
            resg[i][j] = matrix[order[i]][j];

    for (i = 0; i < seqsNumber; i++) {
        for (j = 0; j < resNumber; j++)
            matrix[i][j] = resg[i][j];
        delete[] resg[i];
    }
    delete[] resg;
}

bool sequencesMatrix::getSequence(string seqName, int *sequence) {
	 // Create a timer that will report times upon its destruction
	 //	which means the end of the current scope.
	StartTiming("bool sequencesMatrix::getSequence(string seqName, int *sequence) ");
    int i, pos;

    for (pos = 0; pos < seqsNumber; pos++)
        if (seqsName[pos].compare(seqName) == 0)
            break;

    if (pos == seqsNumber)
        return false;

    for (i = 0; i < resNumber; i++)
        sequence[i] = matrix[pos][i];

    return true;
}

int sequencesMatrix::getSeqNumber(void) {
	 // Create a timer that will report times upon its destruction
	 //	which means the end of the current scope.
	StartTiming("int sequencesMatrix::getSeqNumber(void) ");
    return seqsNumber;
}

int sequencesMatrix::getResidNumber(void) {
	 // Create a timer that will report times upon its destruction
	 //	which means the end of the current scope.
	StartTiming("int sequencesMatrix::getResidNumber(void) ");
    return resNumber;
}


