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

#ifndef NEW_ALIGNMENT_H
#define NEW_ALIGNMENT_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <memory>
#include <cmath>
#include <ctime>

#include "defines.h"

// Forward declarations
class Cleaner;
namespace statistics {
    class Manager;
}

class Alignment {

    int dataType = SequenceTypes::NotDefined;

public:
// BEGIN of Submodules

    class sequencesMatrix;

    Cleaner *Cleaning = nullptr;

    statistics::Manager *Statistics = nullptr;

    sequencesMatrix *SequencesMatrix = nullptr;

// END of Submodules

    int *SeqRef = nullptr;

    int originalNumberOfSequences = 0;

    int numberOfSequences = 0;

    int originalNumberOfResidues = 0;

    int numberOfResidues = 0;

    bool isAligned = false;

    std::string *sequences = nullptr;

    std::string *seqsName = nullptr;

    std::string *seqsInfo = nullptr;

    std::string filename;

    std::string alignmentInfo;

    /* Sequence Identities */
    float **identities = nullptr;

    /* Sequence Overlaps */
    float **overlaps = nullptr;

    /* New Info */
    int *saveResidues = nullptr;
    int *saveSequences = nullptr;

    bool fillMatrices(bool aligned, bool checkInvalidChars = true);

public:

    Alignment();

    Alignment(Alignment &);

    ~Alignment();

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

    int getAlignmentType() const;

    bool isFileAligned();

    Alignment *getTranslationCDS(Alignment *proteinAlignment);

    bool checkCorrespondence(std::string *names, int *lenghts, int totalInputSequences, int multiple);

    void calculateColIdentity(float *columnIdentity);

    void setKeepSequencesFlag(bool newFlagValue);

    void printAlignmentInfo(std::ostream &output);

    bool prepareCodingSequence(bool splitByStopCodon, bool ignoreStopCodon, Alignment *proteinAlignment);

    bool alignmentSummaryHTML(const Alignment &trimmedAlig, const char * const destFile);

    bool statSVG(const char *const destFile);

    bool alignmentSummarySVG(Alignment &trimmedAlig, const char *destFile, int blocks);

    void updateSequencesAndResiduesNums(bool countSequences = true, bool countResidues = true);
};

#endif
