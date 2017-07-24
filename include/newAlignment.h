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

#include <time.h>

#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "statisticsGaps.h"
#include "sequencesMatrix.h"
#include "statisticsConservation.h"
#include "similarityMatrix.h"
#include "utils.h"

#ifndef ALIGNMENT_H

#define SINGLE  1
#define MULTI   2

#define DNAType 1
#define RNAType 2
#define AAType  3


struct newValues {
  int residues;
  int sequences;
  string *matrix;
  string *seqsName;
};

#endif

using namespace std;

#include "Cleaner.h"
#include "StatisticsManager.h"
#include "ReadWriteManager.h"

/** \brief Class containing an alignment
 *
 * This class stores the alignment. It provides methods
 * to \b clean the alignment and generate the clean alignment.
 * It also provides methods for \b statistics \b calculation and
 * \b statistics \b printing.
 */
class newAlignment {

public:
    Cleaner* Cleaning;
    StatisticsManager* Statistics;
    ReadWriteManager* ReadWrite;
    sequencesMatrix *SequencesMatrix;

  int sequenNumber;
  int residNumber;

  bool isAligned;
  bool reverse;

  int dataType;

  string *sequences;
  string *seqsName;
  string *seqsInfo;

    string filename;
    string aligInfo;

  /* Statistics */
  statisticsGaps *sgaps;
  
  statisticsConservation *scons;

  /* Sequence Identities */
  float **identities;

  /* New Info */
  bool oldAlignment;
  int *residuesNumber;
  int *saveResidues;
  int *saveSequences;

  /* ***** Fill the matrices from the input alignment ***** */
  bool fillMatrices(bool aligned);
  
 public:

  /* Constructor */
  newAlignment(void);

  /* Copy Constructor */
  newAlignment(newAlignment&);

  /* Destructor */
  ~newAlignment(void);

  int getNumSpecies(void);

  void getSequences(string * names);

  void getSequences(string * names, int * lengths);

  void getSequences(string * names, string * sequences, int * lenghts);

  bool getSequenceNameOrder(string * names, int * orderVector);

  int getNumAminos(void);

  void setWindowsSize(int ghWindow, int shWindow);

  void setBlockSize(int blockSize);

//   void setOutputFormat(int, bool);

  void setReverseFlag(bool newFlagValue);

//   int getShortNames(void);

  int getReverseFlag(void);

  int getBlockSize(void);

  void printSeqIdentity(void);

  int getAlignmentType(void);

  int *getCorrespResidues(void);

  int *getCorrespSequences(void);

  bool isFileAligned(void);

  newAlignment * getTranslationCDS(int newResidues, int newSequences, int * columnsToKeep, string * oldSequencesNames, sequencesMatrix * sequenceMatrix, newAlignment * proteinAlignment);

  bool checkCorrespondence(string * names, int * lenghts, int totalInputSequences, int multiple);

  void fillNewDataStructure(string * newMatrix, string * newNames);
  
  void fillNewDataStructure(newValues * data);

  void calculateColIdentity(float * columnIdentity);
  
  void printColumnsIdentity_DescriptiveStats(void);

  void setKeepSequencesFlag(bool newFlagValue);

  void setKeepSeqsHeaderFlag(bool newFlagValue);

  void printAlignmentInfo(ostream & output);

  bool prepareCodingSequence(bool splitByStopCodon, bool ignoreStopCodon, newAlignment * proteinAlignment);

};

#endif
