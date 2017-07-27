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
// #include "ReadWriteManager.h"

// //  * This class stores the alignment. It provides methods
// //  * to \b clean the alignment and generate the clean alignment.
// //  * It also provides methods for \b calculation and
// //  * \b statistics \b printing.

/** \brief Class containing an alignment\n
 * This class stores the alignment sequences with it's names, residues and extra information.\n
 * It contains multiple methods regarding the sequences.\n
 * It also contains submodules that provide methods for<b> Calculating statistics</b>,
 * <b> Cleaning the alignment</b> and<b> Printing sequences information</b>.
 */
class newAlignment {

public:
/* Submodules */
    /**
    * \brief Cleaning submodule.\n
    * It contains methods and variables related to trimming.
    */
    Cleaner* Cleaning;
    /**
     * \brief Statistics submodule.\n
     * It contains methods and variables related to Statistics calculation and reporting.
     */
    StatisticsManager * Statistics;
    /**
     * \brief SequencesMatrix submodule\n
     * \todo Give a good explanation of the module.
     */
    sequencesMatrix * SequencesMatrix;
/* Statistics */
    /**
     * \brief Gaps Statistics submodule.
     */
  statisticsGaps * sgaps;
  /**
   * \brief Conservation Statistics submodule.
   */
  statisticsConservation * scons;

    /**
     * \brief Number of sequences present on the alignment.
     */
  int sequenNumber;
    /**
    * \brief Number of residues present on the alignment if it is aligned.
    */
  int residNumber;
    /**
     * \brief Flag that indicates if all sequences on the alignment have the same length.
     */
  bool isAligned;
    /**
     * \brief Flag indicating that this alignment should be reversed.
     */
//   bool reverse;
    /**
     * \brief Number that represents the type of the alignment.
     * \todo We should change this from a define to a enum.
     */
  int dataType;

    /**
     * \brief Sequences vector containing the sequences residues.
     */
  string *sequences;
    /**
     * \brief Sequences names vector containing the sequences names.
     */
  string *seqsName;
    /**
     * \brief Sequences information vector containing extra information about each sequence.
     */
  string *seqsInfo;

  string filename;
    /**
     * \brief String containing information of the alignment as a whole.
     */
  string aligInfo;


  /* Sequence Identities */
  float **identities;

  /* New Info */
  int *residuesNumber;
  int *saveResidues;
  int *saveSequences;

  /**
   * \brief Method to intialize data that has been modified on the alignment.\n
   * The method checks if the sequences have been correctly loaded and are free of errors.\n
   * It checks if sequences contain unknown characters, sets the isAligned flag and initializes 'saveResidues' and 'saveSequences', depending on the sizes of the sequences (wether or not they have the same length).
   * \param aligned Flag to make the method check if the alignment is aligned or not.
   * \note Even with the aligned flag set to false, if the alignment is aligned it will initialize the variables 'saveResidues' and 'saveSequences', allowing the alignment to be trimmed.
   * \return \b True if the alignment information was ok. \n \b False if there was a problem.\n
   * This could happen if the sequences contain any unknown character or if aligned flag is set up to corrent, if the sequences have variable sizes.
   */
  bool fillMatrices(bool aligned);
  
 public:

  /** 
   * \brief
   * Constructor
   */
  newAlignment(void);

  /** \brief
   * Copy Constructor 
   */
  newAlignment(newAlignment&);

  /** \brief
   * Destructor 
   */
  ~newAlignment(void);
    /**
     * \brief Getter for the number of species.
     * \return Number of sequences in the alignment.
     */
  int getNumSpecies(void);
    /**
     * \brief Getter for the sequences names.
     * \param[out] names Vector of sequences names to fill.
     */
  void getSequences(string * names);
    /**
     * \brief Getter for the sequences names and its lenghts.
     * \param[out] names Vector of sequences names to fill.
     * \param[out] lengths Vector of lenghts to fill.
     */
  void getSequences(string * names, int * lengths);
    /**
     * \brief Getter for the sequences, its names and lenghts.
     * \param[out] names Vector of sequences names to fill.
     * \param[out] sequences Vector of sequences to fill.
     * \param[out] lenghts Vector of lenghts to fill.
     */
  void getSequences(string * names, string * sequences, int * lenghts);
    /**
     * \brief Method to map two sets of sequences, own and external.\n
     * The method accepts a vector of names, and test them against the alignment owned sequences.\n
     * If the sequence is present, its order on the alignment sequences vector will be stored in orderVector, in the same position it's name was in the names vector.\n
     * \param names Vector of sequences names to map.
     * \param[out] orderVector Vector of orders to fill.
     * \return \b True if both sets have the same sequences.\n \b False otherwise.
     */
  bool getSequenceNameOrder(string * names, int * orderVector);
    /**
     * \brief Residues getter.
     * \return Number of residues in the alignment.
     */
  int getNumAminos(void);
    /**
     * \brief Windows setter
     * \param ghWindow Half the Gap Window.
     * \param shWindow Half the similarity Window.
     */
  void setWindowsSize(int ghWindow, int shWindow);
    /**
     * \brief blockSize Setter.
     * \param blockSize New value.
     */
  void setBlockSize(int blockSize);
  
    /**
     * \brief blockSize getter.
     * \return blockSize.
     */
  int getBlockSize(void);
    /**
     * \brief Method to print different identity values computed from the alignment.\n
     * In this method we assess the identity values matrix, as well as diferent average values. \n
     * Moreover, the method computes which one is the most similar sequence, in term of identity values, for each one on this alignment
     */
  void printSeqIdentity(void);
    /**
     * \brief Alignment type getter.
     * \return Int representing the alignment type.\n
     *  - DNAType 1
     *  - RNAType 2
     *  - AAType  3
     */
  int getAlignmentType(void);
    /**
     * \brief Method that returns the correspondence between old and new alignment residues.
     * \return Vector that represents which residues are kept (!= -1) and which are rejected (== -1) 
     */
  int *getCorrespResidues(void);
    /**
     * \brief Method that returns the correspondence between old and new alignment sequences.
     * \return Vector that represents which sequences are kept (!= -1) and which are rejected (== -1) 
     */
  int *getCorrespSequences(void);
    /**
     * \brief isAligned getter.
     * \return Wheter the alignment is aligned or not.
     */
  bool isFileAligned(void);
    /**
     * \brief Method to obtain a DNA alignment from a Protein alignment.
     * \todo Check the need for so many arguments.
     * \param newResidues Number of residues of the new alignment, it should be 3 x old.residNumber
     * \param newSequences Number of residues of the new alignment, it should be the same as old.sequenNumber
     * \param columnsToKeep Columns to keep ? IDK really.
     * \param oldSequencesNames Sequences names.
     * \param sequenceMatrix sequence Matrix.
     * \param proteinAlignment original alignment.
     * \return Pointer to new alignment.
     */
  newAlignment * getTranslationCDS(int newResidues, int newSequences, int * columnsToKeep, string * oldSequencesNames, sequencesMatrix * sequenceMatrix, newAlignment * proteinAlignment);
    /**
     * \brief Function to check CDS file.\n
     * \todo multiple argument is needed?
     * It checks if sequences of input alignment are all present on the CDS file.\n
     * If nucleotide sequence is larger than protein sequence length * 3, warns about it and cuts the nucleotide sequence\n
     * If sequence has indetermination symbols, it warns about it.\n
     * If nucleotide sequence is smaller than protein sequence * 3, it adds some 'N' at the end of the nucleotide sequence.
     * \param names Vector containing the names to check.
     * \param lenghts Vector containing the length of each sequence.
     * \param totalInputSequences Number of sequences present.
     * \param multiple Multiplies the length of each sequence by this number.
     * \return \b True if all went right.\n \b False if there is a sequence in the alignment but not in the names vector.
     */
  bool checkCorrespondence(string * names, int * lenghts, int totalInputSequences, int multiple);
    /**
     * \brief Method that copies information of the alignment to the pointers given.
     * \param[out] newMatrix Sequences vector to fill.
     * \param[out] newNames Names vector to fill.
     */
  void fillNewDataStructure(string * newMatrix, string * newNames);
  /**
     * \brief Method that copies information of the alignment to the pointers given.
     * \param[out] data newValues struct to fill with sequences and sequences names.
     */
  void fillNewDataStructure(newValues * data);
    /**
     * \brief Method that calculates the columns identity value.\n
     * This is, the frequency of the most present element in the column, being residue, indetermination and gap allowed.
     * \param[out] columnIdentity Vector to fill with identities for each column.
     */
  void calculateColIdentity(float * columnIdentity);
  /**
   * \brief Method to print Indentity stats: Min, Max, Avg and Std.
   */
  void printColumnsIdentity_DescriptiveStats(void);
    /**
     * \brief Keep Sequences setter.
     * \param newFlagValue New value
     */
  void setKeepSequencesFlag(bool newFlagValue);
    /**
     * \brief Keep Header setter.
     * \param newFlagValue New value
     */
  void setKeepSeqsHeaderFlag(bool newFlagValue);
    /**
     * \brief Print information about sequences number, average sequence length, maximum and minimum sequences length, etc 
     * \param output Output stream.
     */
  void printAlignmentInfo(ostream & output);
    /**
     * \brief Method to check if the CDS file is correct.\n
     * Based on nature of residues: DNA/RNA (Most of the residues)\n
     * There is no gaps on the whole dataset.\n
     * Each sequence is a multiple of 3.\n
     * \n
     * It will also remove or split stop codons depending on the flags passed.
     * \param splitByStopCodon Flag that informs the method to split sequences if it finds any stop codon.
     * \param ignoreStopCodon Flag that informs the method to stop reading sequence if it finds any stop codon.
     * \param proteinAlignment Alignment containing protein sequences which names contains all names in the current alignment.
     */
  bool prepareCodingSequence(bool splitByStopCodon, bool ignoreStopCodon, newAlignment * proteinAlignment);
    /**
     * \brief Method to report the trimming results in HTML.\n
     * It contains information about the original and trimmed alignments, showing which residues and sequences have been kept and which ones have been removed on the output.
     * \param destFile Filename where to save the results.
     * \param residues Number of residues present on the trimmed alignment.
     * \param seqs Number of sequences present on the trimmed alignment.
     * \param selectedRes Vector containing the indexes of the original alignment residues that are kept on the trimmed alignment.
     * \param selectedSeq Vector containing the indexes of the original alignment sequences that are kept on the trimmed alignment.
     * \param consValues Vector containing the consistency values of columns.\n The consistency may be a null pointer, meaning we don't want to report consistency values.
     * \return \b True if everything went ok. \n \b False if file couldn't be open or alignment is not aligned.
     */
  bool alignmentSummaryHTML(char * destFile, int residues, int seqs, int * selectedRes, int * selectedSeq, float * consValues);
    /**
     * \brief Method that saves a clustal-color based visualization of the alignment.\n
     * It doesn't report any trimming result, meaning this only shows information of an alignment.
     * \param file Where to save the visualization.
     * \return \b True if everything went ok. \n \b False if alignment is not aligned.
     */
  bool alignmentColourHTML(ostream &file);
    
//     /**
//     * \brief Seems to be not implemented.
//     */
// //   void getSequences(ostream &);
};

#endif
