/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
   ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

    trimAl v1.4: a tool for automated alignment trimming in large-scale
                 phylogenetics analyses.

    readAl v1.4: a tool for automated alignment conversion among different
                 formats.

    2009-2011 Capella-Gutierrez S. and Gabaldon, T.
              [scapella, tgabaldon]@crg.es

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

***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

#ifndef NEW_ALIGNMENT_H
#define NEW_ALIGNMENT_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <memory>
#include <cmath>
#include <ctime>

// Forward declarations
class sequencesMatrix;
class StatisticsManager;
class Cleaner;

class newAlignment {

    int dataType = 0;

public:
// BEGIN of Submodules

    Cleaner *Cleaning = nullptr;

    StatisticsManager *Statistics = nullptr;

    sequencesMatrix *SequencesMatrix = nullptr;

// END of Submodules

    int *SeqRef = nullptr;

    int originalSequenNumber = 0;

    int sequenNumber = 0;

    int originalResidNumber = 0;

    int residNumber = 0;

    bool isAligned = false;

    std::string *sequences = nullptr;

    std::string *seqsName = nullptr;

    std::string *seqsInfo = nullptr;

    std::string filename;

    std::string aligInfo;

    /* Sequence Identities */
    float **identities = nullptr;

    /* Sequence Overlaps */
    float **overlaps = nullptr;

    /* New Info */
    int *saveResidues = nullptr;
    int *saveSequences = nullptr;

    bool fillMatrices(bool aligned);

public:

    newAlignment();

    newAlignment(newAlignment &);

    ~newAlignment();

    int getNumSpecies();

    void getSequences(std::string *names);

    void getSequences(std::string *names, int *lenghts);

    void getSequences(std::string *names, std::string *sequences, int *lenghts);

    bool getSequenceNameOrder(std::string *names, int *orderVector);

    int getNumAminos();

    void setWindowsSize(int ghWindow, int shWindow);

    void setBlockSize(int blockSize);

    void calculateSeqOverlap();

    void printSeqIdentity();

    void printSeqOverlap();

    int getAlignmentType();

    bool isFileAligned();

    newAlignment *getTranslationCDS(newAlignment *proteinAlignment);

    bool checkCorrespondence(std::string *names, int *lenghts, int totalInputSequences, int multiple);

    void calculateColIdentity(float *columnIdentity);

    void setKeepSequencesFlag(bool newFlagValue);

    void printAlignmentInfo(std::ostream &output);

    bool prepareCodingSequence(bool splitByStopCodon, bool ignoreStopCodon, newAlignment *proteinAlignment);

    bool alignmentSummaryHTML(newAlignment &_trimmedAlignment, char *destFile);

    bool alignmentSummarySVG(newAlignment &_trimmedAlignment, char *destFile, int blocks);

    void updateSequencesAndResiduesNums(bool countSequences = true, bool countResidues = true);
};

#endif
