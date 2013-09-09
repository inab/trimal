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

#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <fstream>
#include <iostream>

#include <time.h>

#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "statisticsGaps.h"
#include "sequencesMatrix.h"
#include "statisticsConservation.h"
#include "similarityMatrix.h"
#include "utils.h"

#define SINGLE  1
#define MULTI   2

#define DNAType 1
#define RNAType 2
#define AAType  3

using namespace std;

struct newValues {
  int residues;
  int sequences;
  string *matrix;
  string *seqsName;
};

/** \brief Class containing an alignment
 *
 * This class stores the alignment. It provides methods
 * to \b clean the alignment and generate the clean alignment.
 * It also provides methods for \b statistics \b calculation and
 * \b statistics \b printing.
 */

class alignment {

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

  /* Sequences */
  sequencesMatrix *seqMatrix;

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
  /* ***** ***** ***** ***** ***** ***** ***** ***** ***** */

  /* Alignment cleaning */

  alignment *cleanByCutValue(double, float, const int *, bool);

  alignment *cleanByCutValue(float, float, const float *, bool);

  alignment *cleanByCutValue(double, const int *, float, float, const float *, bool);

  alignment *cleanStrict(int, const int *, float, const float *, bool, bool);

  alignment *cleanOverlapSeq(float, float *, bool);
  /* ***** ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** ***** */

 public:

  /* Constructors */
  alignment(void);

  alignment(string, string, string *, string *, string *, int, int, int, int,
    bool, int, int, bool, bool, bool, bool, int, int, int *, int *, int *, int,
     int, int, float **);

  /* Overlap the operator = to use it as a constructor */
  alignment &operator=(const alignment &);

  /* Destructor */
  ~alignment(void);

  /* Basic operations  */

  /** \brief Alignment load method.
   * \param alignmentFile Alignment file name.
   * \return \e true if the load is ok, \e false if the load was wrong (i.e. the file doesn't exists).
   *
   * Method that loads an alignment from a file.
   */
  bool loadAlignment(char *alignmentFile);

  /** \brief Alignment storing method.
   * \param destFile Destination file name of the alignment.
   * \return \e true if all is ok, \e false if there were errors (i.e. the file doesn't exists).
   *
   * Method that saves an alignment to a file.
   */
  bool saveAlignment(char *destFile);

  /** \brief Alignment printing method.
   * \return \e true if all is ok, \e false if there were errors.
   *
   * Method that prints an alignment to the standard output.
   */
  bool printAlignment(void);


  /* Alignment trimming. */

  /** \brief Alignment trimming using gap method.
   * \param baseLine base line, security, percentage of colums in the alignment.
   * \param threshold maximum percentage of gaps permitted per column.
   * \return the clean alignment if all is ok, \e NULL if there were errors
   *
   * Method that cleans the alignment using gap method, and generates a new,
   * clean alignment.
   */

  /* ****************************************************************************************************************** */
  /* ********************************************************************************************************************/
  alignment *cleanGaps(float, float, bool);
  /* ****************************************************************************************************************** */
  /* ****************************************************************************************************************** */

  /** \brief Alignment trimming using conservation method.
   * \param baseLine base line, security, percentage of colums in the alignment.
   * \return the clean alignment if all is ok, \e NULL if there were errors
   *
   * Method that cleans the alignment using conservation method, and generates a new,
   * clean alignment.
   */

  /* *************************************************************************************************************** */
  /* *************************************************************************************************************** */
  alignment *cleanConservation(float, float, bool);
  /* *************************************************************************************************************** */
  /* *************************************************************************************************************** */


  /* *************************************************************************************************************** */
  /* *************************************************************************************************************** */
  alignment *cleanCombMethods(bool, bool);
  /* *************************************************************************************************************** */
  /* *************************************************************************************************************** */

  /* *************************************************************************************************************** */
  /* *************************************************************************************************************** */
  alignment *cleanNoAllGaps(bool);
  /* *************************************************************************************************************** */
  /* *************************************************************************************************************** */

  /** \brief Alignment trimming using gap method and conservation method.
   * \param baseLine base line, security, percentage of colums in the alignment.
   * \param threshold maximum percentage of gaps permitted per column.
   * \param conservationPct minimum conservatoin percentage permitted in the clean alignment.
   * \return the clean alignment if all is ok, \e NULL if there were errors.
   *
   * Method that cleans the alignment using both gaps method and conservation method, and
   * generates a new, clean alignment.
   */

  /* *************************************************************************************************************** */
  /* *************************************************************************************************************** */
  alignment *clean(float, float, float, bool);
  /* *************************************************************************************************************** */
  /* *************************************************************************************************************** */


  /** \brief Alignment trimming using a vector of comparison values.
   * \param threshold minimal comparison value permitted per column.
   * \param vectValue comparison values vector from the alignment to clean.
   * \return the clean alignment if all is ok, \e NULL if there were errors
   *
   * Method that cleans the alignment using vector of comparison values, and generates
   *  a new clean alignment.
   */

  /* *************************************************************************************************************** */
  /* *************************************************************************************************************** */
  alignment *cleanCompareFile(float, float, float *, bool);
  /* *************************************************************************************************************** */
  /* *************************************************************************************************************** */

  /** \brief Alignment automatic trimming relaxed method.
   * \return the clean alignment if all is ok, \e NULL if there were errors.
   *
   * Method that cleans, automaticly, the alignment using a second slope method between alignment's lenght
   * and gaps percentage ratio and, genenerates a new clean alignment.
   */


  /* *************************************************************************************************************** */
  /* *************************************************************************************************************** */
  alignment *clean2ndSlope(bool);
  /* *************************************************************************************************************** */
  /* *************************************************************************************************************** */

  alignment *cleanSpuriousSeq(float, float, bool);

  alignment *removeColumns(int *, int, int, bool);

  alignment *removeSequences(int *, int, int, bool);

  alignment *getClustering(float);

  float getCutPointClusters(int);

   /* Statistics calculation */

  /** \brief Basic conservation statistics calculation.
   * \return \e true if all is ok, \e false if there were errors (i.e. there is no similarity matrix defined
   * in conservation statistics).
   *
   * This method calculates conservation statistics with the previously defined similarity matrix.
   */
  bool calculateConservationStats(void);

  /** \brief Conservation statistics calculation.
   * \param sm similarity matrix used for the statistics calculation.
   * \return \e true if all is ok, \e false if there were errors.
   *
   * This method calculates conservation statistics using the \b sm similarity matrix.
   */
  bool setSimilarityMatrix(similarityMatrix *sm);

  /** \brief Gap statistics calculation.
   * \return \e true if all is ok, \e false if there were errors.
   *
   * This method calculates gap statistics without window calculation (half window value = 0).
   */
  bool calculateGapStats(void);

  /* Output Statistics */

  /** \brief Printing normal gap statistics method.
   *
   * This method prints gap statistics for each column.
   */
  void printStatisticsGapsColumns(void);

  /** \brief Printing accumulated gap statistics method.
   *
   * This method prints accumulated gap statistics.
   */
  void printStatisticsGapsTotal(void);

  /** \brief Printing conservation values method.
   *
   * This method prints conservation value for each column of the alignment.
   */
  void printStatisticsConservationColumns(void);

  /** \brief Printing accumulated conservation statistics method.
   * This method prints the accumulated number of columns for each conservation
   * value from the the alignment.
   */
  void printStatisticsConservationTotal(void);

  void printCorrespondence(void);

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

  /** \brief Building alignment's sequence matrix method.
   *
   * This method builds an alignment's sequence matrix.
   */
  void sequenMatrix(void);

  void destroySequenMatrix(void);

  /** \brief Printing alignment's sequence matrix method.
   *
   * This method prints an alignment's sequence matrix.
   */
  void printSequenMatrix(void);

  /** \brief Returns a column from alignment's sequence matrix.
   * \param colum, sequence matrix index
   * \param columnSeqMatrix, vector used to storage a column from alignment sequence matrix.
   *
   * This method returns a column from alignment sequence matrix.
   */
  void getColumnSeqMatrix(int, int *);

  /** \brief Returns a column from alignment's sequence matrix.
   * \param value to look in a sequence matrix row.
   * \param sequence matrix row where look for a value.
   * \param columnSeqMatrix, vector used to storage a column from alignment sequence matrix.
   *
   * Method that returns a column from the aligment's sequence matrix with the same value that
   * "value" at matrix's position (row, i)
   */
  void getColumnSeqMatrix(int, int, int *);

  void setSeqMatrixOrder(int *);

  sequencesMatrix *getSeqMatrix(void);

  /* ********** NEW CODE ********** */
  /* ********** ******** ********** */
  int formatInputAlignment(char *);

  int formatInputFile(void);

  int typeInputFile(void);

  bool loadPhylipAlignment(char *);

  bool loadFastaAlignment(char *);

  bool loadClustalAlignment(char *);

  bool loadNexusAlignment(char *);

  bool loadMegaInterleavedAlignment(char *);

  bool loadMegaNonInterleavedAlignment(char *);

  bool loadNBRF_PirAlignment(char *);

  bool loadPhylip3_2Alignment(char *);
  /* ********** ******** ********** */
  /* ********** ******** ********** */

  /* Alignment to a stream */
  void alignmentClustalToFile(ostream &);

  void alignmentNBRF_PirToFile(ostream &);

  void alignmentFastaToFile(ostream &);

  void alignmentPhylip3_2ToFile(ostream &);

  void alignmentPhylipToFile(ostream &);

  void alignmentPhylip_PamlToFile(ostream &);

  void alignmentNexusToFile(ostream &);

  void alignmentMegaToFile(ostream &);

  bool alignmentSummaryHTML(char *, int, int, int *, int *, float *);

  bool alignmentColourHTML(ostream &);

  void getSequences(ostream &);
  /* ********** ******** ********** */

  void saveStatistics(similarityMatrix *sm);

  void saveStatistics(similarityMatrix *, int, int);

  void trimTerminalGaps(bool);

  void setWindowsSize(int, int);

  void setBlockSize(int);

  void setOutputFormat(int, bool);

  void setReverse(void);

  int getInputFormat(void);

  int getOutputFormat(void);

  int getShortNames(void);

  int getReverse(void);

  int getBlockSize(void);

  void calculateSeqIdentity(void);

  // New
  void calculateRelaxedSeqIdentity(void);

  int selectMethod(void);

  void printSeqIdentity(void);

  bool calculateSpuriousVector(float, float *);

  void checkTypeAlignment(void);

  int getTypeAlignment(void);

  int *getCorrespResidues(void);

  int *getCorrespSequences(void);

  bool isFileAligned(void);

  alignment * getTranslationCDS(int, int, int *, string *, sequencesMatrix *, alignment *);

  bool checkCorrespondence(string *, int *, int, int);

  int *calculateRepresentativeSeq(float maximumIdent);

  /* New code: version 1.4 */

  void computeComplementaryAlig(bool, bool);

  void removeSmallerBlocks(int);

  bool removeOnlyTerminal(void);

  newValues removeCols_SeqsAllGaps(void);

  void removeCols_SeqsAllGaps(newValues *);

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
  bool prepareCodingSequence(bool, bool, alignment *);
};

#endif
