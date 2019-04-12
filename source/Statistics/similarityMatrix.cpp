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

#ifndef SIMMatrix
#define SIMMatrix

#define NUMAMINOS 20
#define TAMABC 28
#define LINE_LENGTH 256
#define REFER 65

#include "defines.h"

#endif

#include "InternalBenchmarker.h"
#include "Statistics/similarityMatrix.h"
#include "reportsystem.h"
#include "values.h"
#include "utils.h"

namespace statistics {

    similarityMatrix::similarityMatrix() {
        // Create a timer that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("similarityMatrix::similarityMatrix() ");
        numPositions = 0;
        vhash = nullptr;
        simMat = nullptr;
        distMat = nullptr;
    }

    void similarityMatrix::memoryAllocation(int nPos) {
        // Create a timer that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("void similarityMatrix::memoryAllocation(int nPos) ");
        int i, j;

        // Initialize square table dimension to store the distances
        // and to store the similarity matrix.
        if (numPositions != 0) memoryDeletion();
        numPositions = nPos;

        // Reserve memory for all structures
        vhash = new int[TAMABC];

        simMat = new float *[nPos];
        distMat = new float *[nPos];

        for (i = 0; i < nPos; i++) {
            simMat[i] = new float[nPos];
            distMat[i] = new float[nPos];

            for (j = 0; j < nPos; j++) {
                distMat[i][j] = 0.0;
                simMat[i][j] = 0.0;
            }
        }
    }

    similarityMatrix::~similarityMatrix() {
        // Create a timer that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("similarityMatrix::~similarityMatrix() ");

        if (numPositions != 0) memoryDeletion();

    }

    void similarityMatrix::memoryDeletion() {
        // Create a timer that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("void similarityMatrix::memoryDeletion() ");
        int i;

        for (i = 0; i < numPositions; i++) {
            delete[] simMat[i];
            delete[] distMat[i];
        }

        delete[] distMat;
        delete[] simMat;
        delete[] vhash;

        numPositions = 0;
        distMat = nullptr;
        simMat = nullptr;
        vhash = nullptr;
    }

    bool similarityMatrix::loadSimMatrix(char *filename) {
        // Create a timer that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("bool similarityMatrix::loadSimMatrix(char *filename) ");
        char aux[LINE_LENGTH + 1],
                first[LINE_LENGTH],
                listSym[LINE_LENGTH + 1];
        int i, j, k;
        float sum;
        bool firstColumn = true;
        std::ifstream file;

        // We try to open the file, if we can't open the file
        // we return false.
        file.open(filename);
        if (file.fail())
            return false;

        // Read the first line of the file and, depending on the
        // line length (free of spaces an tabulators), we allocate
        // memory for the object structures
        file.getline(aux, LINE_LENGTH);
        utils::removeSpaces(aux, listSym);
        memoryAllocation(strlen(listSym));

        utils::initlVect(vhash, TAMABC, -1);
        // We create the hashing vector
        for (i = 0; i < numPositions; i++) {
            listSym[i] = (char) toupper((int) listSym[i]);

            if ((listSym[i] >= 'A') && (listSym[i] <= 'Z')) {
                if ((vhash[listSym[i] - 'A']) != -1) {
                    memoryDeletion();
                    return false;
                }
                vhash[listSym[i] - 'A'] = i;

            } else {
                memoryDeletion();
                return false;
            }
        }

        for (i = 0; i < numPositions; i++) {
            // Read the first symbol of the line
            j = 0;
            file >> first;

            // If the format includes the first aminoacid in the line
            if (firstColumn) {
                first[0] = (char) toupper((int) first[0]);

                // Format checking. The first token must not be a valid number
                if (((first[0] >= '0' && first[0] <= '9') ||
                     (first[0] == '-' && (first[1] >= '0' && first[1] <= '9')))
                    && i > 0) {
                    memoryDeletion();
                    return false;
                }

                // If in the token is a character, there is "first column"
                // in the format of the alignment
                if ((first[0] >= 'A') && (first[0] <= 'Z')) {
                    firstColumn = true;

                    if ((vhash[first[0] - 'A']) == -1) {
                        memoryDeletion();
                        return false;
                    }
                }

                    // If we have read a number there is no "first column"
                    // in the format of the alignment
                else if ((first[0] >= '0' && first[0] <= '9') || (first[0] == '-' && (first[1] >= '0' && first[1] <= '9'))) {
                    firstColumn = false;
                    j = 1;

                    simMat[i][0] = atof(first);
                    first[0] = listSym[i];
                }

            } else {
                j = 1;

                // Do some checkings
                if ((first[0] >= 'A') && (first[0] <= 'Z') && (i > 0)) {
                    memoryDeletion();
                    return false;
                }

                simMat[i][0] = atof(first);
                first[0] = listSym[i];
            }

            // Read the corresponding number row
            for (; j < numPositions; j++)
                file >> simMat[vhash[first[0] - 'A']][j];
        }

        // Calculate the average between two simmetric positions
        // respect to the diagonal of the matrix (between [i][j] and [j][i]
        // If the input is a non-symmetric matrix, the output will be a
        // symmetric matrix

        for (i = 0; i < numPositions; i++) {
            for (j = i + 1; j < numPositions; j++) {
                if (simMat[i][j] != simMat[j][i]) {
                    float value = (simMat[i][j] + simMat[j][i]) / 2.0;
                    simMat[i][j] = value;
                    simMat[j][i] = value;
                }
            }
        }

        // Calculate the distances between aminoacids
        // based on Euclidean distance

        for (j = 0; j < numPositions; j++) {
            for (i = 0; i < numPositions; i++) {
                if ((i != j) && (distMat[i][j] == 0.0)) {
                    for (k = 0, sum = 0; k < numPositions; k++)
                        sum += ((simMat[k][j] - simMat[k][i]) * (simMat[k][j] - simMat[k][i]));
                    sum = (float) sqrt(sum);
                    distMat[i][j] = sum;
                    distMat[j][i] = sum;
                }
            }
        }

        file.close();
        return true;
    }

    void similarityMatrix::defaultAASimMatrix(void) {
        // Create a timer that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("void similarityMatrix::defaultAASimMatrix(void) ");

        int i, j, k;
        float sum;

        memoryAllocation(20);
        for (i = 0; i < TAMABC; i++)
            vhash[i] = -1;

        // We create the hashing vector
        for (i = 0; i < numPositions; i++)
            vhash[listAASym[i] - 'A'] = i;

        for (i = 0; i < numPositions; i++)
            for (j = 0; j < numPositions; j++)
                simMat[i][j] = defaultAAMatrix[i][j];

        // Calculate the distances between aminoacids
        // based on Euclidean distance
        for (j = 0; j < numPositions; j++) {
            for (i = 0; i < numPositions; i++) {
                if ((i != j) && (distMat[i][j] == 0.0)) {
                    for (k = 0, sum = 0; k < numPositions; k++)
                        sum += ((simMat[k][j] - simMat[k][i]) * (simMat[k][j] - simMat[k][i]));
                    sum = (float) sqrt(sum);
                    distMat[i][j] = sum;
                    distMat[j][i] = sum;
                }
            }
        }
    }

    void similarityMatrix::defaultNTSimMatrix(void) {
        // Create a timer that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("void similarityMatrix::defaultNTSimMatrix(void) ");
        int i, j, k;
        float sum;

        memoryAllocation(5);
        for (i = 0; i < TAMABC; i++)
            vhash[i] = -1;

        // We create the hashing vector
        for (i = 0; i < numPositions; i++)
            vhash[listNTSym[i] - 'A'] = i;

        for (i = 0; i < numPositions; i++)
            for (j = 0; j < numPositions; j++)
                simMat[i][j] = defaultNTMatrix[i][j];

        // Calculate the distances between aminoacids
        // based on Euclidean distance
        for (j = 0; j < numPositions; j++) {
            for (i = 0; i < numPositions; i++) {
                if ((i != j) && (distMat[i][j] == 0.0)) {
                    for (k = 0, sum = 0; k < numPositions; k++)
                        sum += ((simMat[k][j] - simMat[k][i]) * (simMat[k][j] - simMat[k][i]));
                    sum = (float) sqrt(sum);
                    distMat[i][j] = sum;
                    distMat[j][i] = sum;
                }
            }
        }
    }

    void similarityMatrix::defaultNTDegeneratedSimMatrix(void) {
        // Create a timer that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("void similarityMatrix::defaultNTDegeneratedSimMatrix(void) ");
        int i, j, k;
        float sum;

        memoryAllocation(15);
        for (i = 0; i < TAMABC; i++)
            vhash[i] = -1;

        // We create the hashing vector
        for (i = 0; i < numPositions; i++)
            vhash[listNTDegenerateSym[i] - 'A'] = i;

        for (i = 0; i < numPositions; i++)
            for (j = 0; j < numPositions; j++)
                simMat[i][j] = defaultNTDegeneratedMatrix[i][j];

        // Calculate the distances between nucleotides based on Euclidean distance
        for (j = 0; j < numPositions; j++) {
            for (i = 0; i < numPositions; i++) {
                if ((i != j) && (distMat[i][j] == 0.0)) {
                    for (k = 0, sum = 0; k < numPositions; k++)
                        sum += ((simMat[k][j] - simMat[k][i]) * (simMat[k][j] - simMat[k][i]));
                    sum = (float) sqrt(sum);
                    distMat[i][j] = sum;
                    distMat[j][i] = sum;
                }
            }
        }
    }

    void similarityMatrix::alternativeSimilarityMatrices(int matrix_code, \
        int datatype) {
        int i, j, k;
        float sum;

        // Allocate memory depending on the input datatype
        switch (datatype) {
            case SequenceTypes::AA:
                memoryAllocation(20);
                break;
            case SequenceTypes::DNA:
            case SequenceTypes::RNA:
                memoryAllocation(5);
                break;
            case SequenceTypes::DNA | SequenceTypes::DEG:
            case SequenceTypes::RNA | SequenceTypes::DEG:
                memoryAllocation(15);
                break;
        }

        for (i = 0; i < TAMABC; i++)
            vhash[i] = -1;

        // We create the hashing vector taking into account the input datatype
        for (i = 0; i < numPositions; i++) {
            switch (datatype) {
                case SequenceTypes::AA:
                    vhash[listAASym[i] - 'A'] = i;
                    break;
                case SequenceTypes::DNA:
                case SequenceTypes::RNA:
                    vhash[listNTSym[i] - 'A'] = i;
                    break;
                case SequenceTypes::DNA | SequenceTypes::DEG:
                case SequenceTypes::RNA | SequenceTypes::DEG:
                    vhash[listNTDegenerateSym[i] - 'A'] = i;
                    break;
            }
        }

        // Working similarity matrix is set depending on the pre loaded matrices
        for (i = 0; i < numPositions; i++) {
            for (j = 0; j < numPositions; j++) {
                switch (matrix_code) {
                    case 1:
                        simMat[i][j] = alternative_1_NTDegeneratedMatrix[i][j];
                        break;
                }
            }
        }

        // Calculate the distances between residues based on Euclidean distance
        for (j = 0; j < numPositions; j++) {
            for (i = 0; i < numPositions; i++) {
                if ((i != j) && (distMat[i][j] == 0.0)) {
                    for (k = 0, sum = 0; k < numPositions; k++)
                        sum += ((simMat[k][j] - simMat[k][i]) * (simMat[k][j] - simMat[k][i]));
                    sum = (float) sqrt(sum);
                    distMat[i][j] = sum;
                    distMat[j][i] = sum;
                }
            }
        }
    }

    void similarityMatrix::printMatrix() {
        // Create a timer that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("void similarityMatrix::printMatrix() ");

        for (int i = 0; i < numPositions; i++) {
            for (int j = 0; j < numPositions; j++)
                std::cerr << std::setw(8) << std::setprecision(4) << std::right << simMat[i][j];
            std::cerr << std::endl;
        }
    }

    float similarityMatrix::getDistance(char a, char b) {
        // Create a timer that will report times upon its destruction
        //	which means the end of the current scope.
        StartTiming("float similarityMatrix::getDistance(char a, char b) ");
        int numa = -1, numb = -1;

        // Search the first character position
        if ((a >= 'A') && (a <= 'Z')) numa = vhash[a - 'A'];
        else {
            debug.report(ErrorCode::IncorrectSymbol, new std::string[1]{std::string(1, a)});
            return -1;
        }

        // Search the second character position
        if ((b >= 'A') && (b <= 'Z')) numb = vhash[b - 'A'];
        else {
            debug.report(ErrorCode::IncorrectSymbol, new std::string[1]{std::string(1, b)});
            return -1;
        }

        // We check if the two character positions are valid positions
        if (numa == -1) {
            debug.report(ErrorCode::UndefinedSymbol, new std::string[1]{std::string(1, a)});
            return -1;
        }

        if (numb == -1) {
            debug.report(ErrorCode::UndefinedSymbol, new std::string[1]{std::string(1, b)});
            return -1;
        }

        // Return the distance value between a and b
        return distMat[numa][numb];
    }
}