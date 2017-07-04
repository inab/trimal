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
//
//class Cleaner;
//class StatisticsManager;
//class ReadWriteManager;

#include <Cleaner.h>
#include <StatisticsManager.h>
#include <ReadWriteManager.h>

/** \brief Class containing an alignment
 *
 * This class stores the alignment. It provides methods
 * to \b clean the alignment and generate the clean alignment.
 * It also provides methods for \b statistics \b calculation and
 * \b statistics \b printing.
 */
class newAlignment {
    friend class Cleaner;
    friend class StatisticsManager;
    friend class ReadWriteManager;
    friend class sequencesMatrix;

public:
    Cleaner* Cleaning;
    StatisticsManager* Statistics;
    ReadWriteManager* ReadWrite;
    sequencesMatrix *SequencesMatrix;

private:

    bool oldnewAlignment;

  int sequenNumber;
  int residNumber;

  bool isAligned;
  bool reverse;

  bool terminalGapOnly;

  int iformat;
  int oformat;
  bool shortNames;

  bool forceCaps;
  bool upperCase;
  bool lowerCase;

  bool keepSequences;
  bool keepHeader;

  string gapSymbol;

  int dataType;

  int ghWindow;
  int shWindow;

  int blockSize;

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

 private:
  /* ***** Fill the matrices from the input alignment ***** */
  bool fillMatrices(bool aligned);
 public:

  /* Constructors */
  newAlignment(void);

    newAlignment(string, string, string *, string *, string *, int, int, int, int,
    bool, int, int, bool, bool, bool, bool, int, int, int *, int *, int *, int,
     int, int, float **);

  /* Overlap the operator = to use it as a constructor */
  newAlignment &operator=(const newAlignment &);

  /* Destructor */
  ~newAlignment(void);

   /* Alignment's Info */

  /** \brief Gets alignment's sequenNumber number.
   * \return the alignment's sequenNumber number.
   *
   * This method returns the alignment's sequenNumber number.
   */
  int getNumSpecies(void);

  /** \brief Gets alignment's sequenNumber names.
   * \param characters' matrix used to storage sequenNumber names.
   *
   * This method returns the alignment's sequenNumber names.
   */
  void getSequences(string *);

  void getSequences(string *, int *);

  void getSequences(string *, string *, int *);

  bool getSeqNameOrder(string *, int *);

  /** \brief Gets alignment's amino acids number.
   * \return the alignment's amino acids number.
   *
   * This method returns the alignment's amino acids number.
   */
  int getNumAminos(void);


   /* Alignments' Compare */

//  /** \brief Building alignment's sequence matrix method.
//   *
//   * This method builds an alignment's sequence matrix.
//   */
//  void sequenMatrix(void);
//
//  void destroySequenMatrix(void);
//
//  /** \brief Printing alignment's sequence matrix method.
//   *
//   * This method prints an alignment's sequence matrix.
//   */
//  void printSequenMatrix(void);
//
//  /** \brief Returns a column from alignment's sequence matrix.
//   * \param colum, sequence matrix index
//   * \param columnSeqMatrix, vector used to storage a column from alignment sequence matrix.
//   *
//   * This method returns a column from alignment sequence matrix.
//   */
//  void getColumnSeqMatrix(int, int *);
//
//  /** \brief Returns a column from alignment's sequence matrix.
//   * \param value to look in a sequence matrix row.
//   * \param sequence matrix row where look for a value.
//   * \param columnSeqMatrix, vector used to storage a column from alignment sequence matrix.
//   *
//   * Method that returns a column from the aligment's sequence matrix with the same value that
//   * "value" at matrix's position (row, i)
//   */
//  void getColumnSeqMatrix(int, int, int *);
//
//  void setSeqMatrixOrder(int *);
//
//  sequencesMatrix *getSeqMatrix(void);




  void setWindowsSize(int, int);

  void setBlockSize(int);

  void setOutputFormat(int, bool);

  void setReverse(void);



  int getShortNames(void);

  int getReverse(void);

  int getBlockSize(void);

  void calculateSeqIdentity(void);


    // New
  void calculateRelaxedSeqIdentity(void);

  int selectMethod(void);

  void printSeqIdentity(void);


//  void checkTypeAlignment(void);

  int getTypeAlignment(void);

  int *getCorrespResidues(void);

  int *getCorrespSequences(void);

  bool isFileAligned(void);

    newAlignment * getTranslationCDS(int, int, int *, string *, sequencesMatrix *, newAlignment *);

  bool checkCorrespondence(string *, int *, int, int);

  int *calculateRepresentativeSeq(float maximumIdent);

  /* New code: version 1.4 */

  void computeComplementaryAlig(bool, bool);



//  void removeCols_SeqsAllGaps(newValues *);

  void fillNewDataStructure(string *, string *);
  void fillNewDataStructure(newValues *);

  // New Code: February/2012
  void calculateColIdentity(float *);
  void printColumnsIdentity_DescriptiveStats(void);

  // New Code: May/2012
  void setKeepSequencesFlag(bool);

  // New Code: Mar/2013
  void setKeepSeqsHeaderFlag(bool);

  void printAlignmentInfo(ostream &);

  // Updated: June/2013
  bool prepareCodingSequence(bool, bool, newAlignment *);

};

#endif
