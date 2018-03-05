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

#include <fstream>
#include <iostream>
#include <memory>

#include <ctime>

#include <cstring>
#include <cstdlib>
#include <cmath>

#include "sequencesMatrix.h"
#include "similarityMatrix.h"
#include "utils.h"

#include "Cleaner.h"

#include "../include/Statistics/StatisticsManager.h"
using namespace std;


class newAlignment {

    int dataType = 0;
    
public:
// BEGIN of Submodules

    Cleaner * Cleaning                  = nullptr;

    StatisticsManager * Statistics      = nullptr;

    sequencesMatrix * SequencesMatrix   = nullptr;

// END of Submodules
    
    int * SeqRef                        = nullptr;
    
    int originalSequenNumber;
    
    int sequenNumber;
    
    int originalResidNumber;

    int residNumber;

    bool isAligned;

    string *sequences           = nullptr;

    string *seqsName            = nullptr;

    string *seqsInfo            = nullptr;

    string filename;

    string aligInfo;

    /* Sequence Identities */
    float **identities          = nullptr;
    
    /* Sequence Overlaps */
    float **overlaps            = nullptr;

    /* New Info */
    int *saveResidues           = nullptr;
    int *saveSequences          = nullptr;

    bool fillMatrices(bool aligned);
  
 public:

    newAlignment(void);

    newAlignment(newAlignment&);

    ~newAlignment(void);

    int getNumSpecies(void);

    void getSequences(string * names);

    void getSequences(string * names, int * lengths);

    void getSequences(string * names, string * sequences, int * lenghts);

    bool getSequenceNameOrder(string * names, int * orderVector);

    int getNumAminos(void);

    void setWindowsSize(int ghWindow, int shWindow);
    
    void setBlockSize(int blockSize);

    int getBlockSize(void);
    
    void calculateSeqIdentity(void);
      
    void calculateRelaxedSeqIdentity(void);
    
    void calculateSeqOverlap(void);

    void printSeqIdentity(void);
    
    void printSeqOverlap(void);

    int getAlignmentType(void);

    int *getCorrespResidues(void);

    int *getCorrespSequences(void);

    bool isFileAligned(void);

    newAlignment * getTranslationCDS(newAlignment * proteinAlignment);

    bool checkCorrespondence(string * names, int * lenghts, int totalInputSequences, int multiple);

    void fillNewDataStructure(string * newMatrix, string * newNames);

    void calculateColIdentity(float * columnIdentity);

    void printColumnsIdentity_DescriptiveStats(void);

    void setKeepSequencesFlag(bool newFlagValue);

    void setKeepSeqsHeaderFlag(bool newFlagValue);

    void printAlignmentInfo(ostream & output);

    bool prepareCodingSequence(bool splitByStopCodon, bool ignoreStopCodon, newAlignment * proteinAlignment);

    bool alignmentSummaryHTML(newAlignment &_trimmedAlignment, char *destFile, statisticsConsistency *consistency);
    
    bool alignmentSummarySVG(newAlignment &_trimmedAlignment, char *destFile, int blocks);

    bool alignmentColourHTML(ostream &file);
    
    void updateSequencesAndResiduesNums(bool countSequences = true, bool countResidues = true);

};

#endif
