/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
   ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

    trimAl v1.4: a tool for automated alignment trimming in large-scale
                 phylogenetics analyses.

    readAl v1.4: a tool for automated alignment conversion among different
                 formats.

    statAl v1.4: a tool for getting stats about multiple sequence alignments.


    2009-2012 Capella-Gutierrez S. and Gabaldon, T.
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
using namespace std;

#include <float.h>
#include "alignment.h"
#include "rwAlignment.cpp"
#include "autAlignment.cpp"

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Class constructor */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

alignment::alignment(void) {

  /* Alignment parameter */
  sequenNumber = 0;
  residNumber =  0;

  /* Are the input sequences aligned? */
  isAligned = false;

  /* Should the output file be reversed? */
  reverse   = false;

  /* Should be trimmed only terminal gaps? */
  terminalGapOnly = false;

  /* Input and output formats */
  iformat = 0;
  oformat = 0;
  shortNames = false;

  forceCaps = false;
  upperCase = false;
  lowerCase = false;

  /* Indicate whether sequences composed only by gaps should be kept or not */
  keepSequences = false;

  /* Indicate whether original header, they may include non-alphanumerical
   * characters, should be dumped into output stream without any preprocessing
   * step */
  keepHeader = false;

  gapSymbol = "-";

  /* Sequence datatype: DNA, RNA or Protein */
  dataType = 0;

  /* Window sizes to trim the input alignment */
  ghWindow = 0;
  shWindow = 0;

  /* Minimum block size in the new alignment */
  blockSize = 0;

  /* Is this alignmnet new? */
  oldAlignment  = false;

  /* Sequence residues number */
  residuesNumber = NULL;

  /* Columns and sequences that have been previously selected */
  saveResidues  = NULL;
  saveSequences = NULL;

  /* Input sequences as well other information such as sequences name, etc */
  sequences = NULL;
  seqsName  = NULL;
  seqsInfo  = NULL;

  /* Information about input alignment */
  filename = "";
  aligInfo = "";

  /* Information computed from alignment */
  sgaps =     NULL;
  scons =     NULL;
  seqMatrix = NULL;
  identities = NULL;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Class constructor */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

alignment::alignment(string o_filename, string o_aligInfo, string *o_sequences, string *o_seqsName,
                     string *o_seqsInfo, int o_sequenNumber, int o_residNumber, int o_iformat, int o_oformat,
                     bool o_shortNames, int o_dataType, int o_isAligned, bool o_reverse, bool o_terminalGapOnly,
                     bool o_keepSeqs, bool o_keepHeader, int OldSequences, int OldResidues, int *o_residuesNumber,
                     int *o_saveResidues, int *o_saveSequences, int o_ghWindow, int o_shWindow, int o_blockSize,
                     float **o_identities) {

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  int i, j, k, ll;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  oldAlignment = true;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Assign the parameter values to the variables */
  sequenNumber = o_sequenNumber;
  residNumber  = o_residNumber;

  iformat = o_iformat;
  oformat = o_oformat;
  shortNames = o_shortNames;

  dataType = o_dataType;

  ghWindow = o_ghWindow;
  shWindow = o_shWindow;

  blockSize = o_blockSize;

  isAligned = o_isAligned;
  reverse   = o_reverse;

  terminalGapOnly = o_terminalGapOnly;

  filename = o_filename;
  aligInfo = o_aligInfo;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  keepSequences = o_keepSeqs;
  keepHeader = o_keepHeader;

  forceCaps = false;
  upperCase = false;
  lowerCase = false;

  gapSymbol = "-";


  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Basic information for the new alignment */
  sequences = new string[sequenNumber];
  for(i = 0; i < sequenNumber; i++)
    sequences[i] = o_sequences[i];

  seqsName = new string[sequenNumber];
  for(i = 0; i < sequenNumber; i++)
    seqsName[i] = o_seqsName[i];
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  residuesNumber = new int[sequenNumber];
  if((isAligned) || (o_residuesNumber != NULL)) {
    for(i = 0; i < sequenNumber; i++)
      residuesNumber[i] = residNumber;
  } else {
    for(i = 0; i < sequenNumber; i++)
      residuesNumber[i] = o_residuesNumber[i];
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  if(o_seqsInfo != NULL) {
    seqsInfo = new string[sequenNumber];
    for(i = 0; i < sequenNumber; i++)
      seqsInfo[i] = o_seqsInfo[i];
  } else seqsInfo = NULL;

  saveResidues  = NULL;
  if(o_saveResidues != NULL) {
    saveResidues = new int[residNumber];
    for(i = 0, j = 0; i < OldResidues; i++)
      if(o_saveResidues[i] != -1) {
        saveResidues[j] = o_saveResidues[i];
        j++;
      }
  }

  saveSequences = NULL;
  if(o_saveSequences != NULL) {
    saveSequences = new int[sequenNumber];
    for(i = 0, j = 0; i < OldSequences; i++)
      if(o_saveSequences[i] != -1) {
        saveSequences[j] = o_saveSequences[i];
        j++;
      }
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  identities = NULL;
  if(o_identities != NULL) {
    identities = new float*[sequenNumber];
    for(i = 0, j = 0; i < OldSequences; i++) {
      if(o_saveSequences[i] != -1) {
        identities[j] = new float[sequenNumber];
        for(k = 0, ll = 0; k < OldSequences; k++) {
          if(o_saveSequences[k] != -1) {
            identities[j][ll] = o_identities[i][k];
            ll++;
          }
        }
        j++;
      }
    }
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Any structure associated to the new alignment is
   * initialize to NULL. In this way, these structure,
   * if it will be necessary, has to be computed */
  sgaps  =     NULL;
  scons  =     NULL;
  seqMatrix =  NULL;
  identities = NULL;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Overlapping operator = to use it as a kind of class constructor */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

alignment &alignment::operator=(const alignment &old) {
  int i, j;

  if(this != &old) {

    oldAlignment = true;

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Assign the parameter values to the variables */
    sequenNumber = old.sequenNumber;
    residNumber =  old.residNumber;

    isAligned =  old.isAligned;
    reverse   =  old.reverse;

    terminalGapOnly = old.terminalGapOnly;

    iformat = old.iformat;
    oformat = old.oformat;
    shortNames = old.shortNames;

    dataType = old.dataType;

    ghWindow = old.ghWindow;
    shWindow = old.shWindow;

    blockSize = old.blockSize;

    filename = old.filename;
    aligInfo = old.aligInfo;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    keepSequences = old.keepSequences;
    keepHeader = old.keepHeader;

    forceCaps = old.forceCaps;
    upperCase = old.upperCase;
    lowerCase = old.lowerCase;

    gapSymbol = old.gapSymbol;

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    sequences = new string[sequenNumber];
    for(i = 0; i < sequenNumber; i++)
      sequences[i] = old.sequences[i];

    seqsName = new string[sequenNumber];
    for(i = 0; i < sequenNumber; i++)
      seqsName[i] = old.seqsName[i];
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    delete [] residuesNumber;
    if(old.residuesNumber) {
      residuesNumber = new int[sequenNumber];
      for(i = 0; i < sequenNumber; i++)
        residuesNumber[i] = old.residuesNumber[i];
    } residuesNumber = NULL;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    delete [] seqsInfo;
    if(old.seqsInfo) {
      seqsInfo = new string[sequenNumber];
      for(i = 0; i < sequenNumber; i++)
        seqsInfo[i] = old.seqsInfo[i];
    } else seqsInfo = NULL;

    delete [] saveResidues;
    if(old.saveResidues) {
      saveResidues = new int[residNumber];
      for(i = 0; i < residNumber; i++)
        saveResidues[i] = old.saveResidues[i];
    } else saveResidues = NULL;

    delete [] saveSequences;
    if(old.saveSequences) {
      saveSequences = new int[sequenNumber];
      for(i = 0; i < sequenNumber; i++)
        saveSequences[i] = old.saveSequences[i];
    } else saveSequences = NULL;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    delete [] identities;
    if(old.identities) {
      identities = new float*[sequenNumber];
      for(i = 0; i < sequenNumber; i++) {
        identities[i] = new float[sequenNumber];
        for(j = 0; j < sequenNumber; j++)
          identities[i][j] = old.identities[i][j];
      }
    } else identities = NULL;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    delete sgaps;
    sgaps = NULL;

    delete scons;
    scons = NULL;

    delete seqMatrix;
    seqMatrix = old.seqMatrix;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  return *this;
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Class destructor */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

alignment::~alignment(void) {
  int i;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  if(sequences != NULL)
    delete [] sequences;
  sequences = NULL;

  if(seqsName != NULL)
    delete [] seqsName;
  seqsName = NULL;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  if(residuesNumber != NULL)
    delete [] residuesNumber;
  residuesNumber = NULL;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  if(seqsInfo != NULL)
    delete [] seqsInfo;
  seqsInfo = NULL;

  if(saveResidues != NULL)
    delete[] saveResidues;
  saveResidues = NULL;

  if(saveSequences != NULL)
    delete[] saveSequences;
  saveSequences = NULL;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  if(identities != NULL) {
    for(i = 0; i < sequenNumber; i++)
      delete [] identities[i];
    delete [] identities;
  }
  identities = NULL;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  if(sgaps != NULL)
    delete sgaps;
  sgaps = NULL;

  if(scons != NULL)
    delete scons;
  scons = NULL;

  if(seqMatrix != NULL)
    delete seqMatrix;
  seqMatrix = NULL;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  oldAlignment = false;

  sequenNumber = 0;
  residNumber  = 0;

  isAligned = false;
  reverse   = false;

  iformat =  0;
  oformat =  0;
  shortNames = false;

  dataType = 0;

  ghWindow = 0;
  shWindow = 0;

  blockSize = 0;

  filename.clear();
  aligInfo.clear();
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This function is useful to detect the format from a given alignment */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
bool alignment::loadAlignment(char *alignmentFile) {

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Analyze the input alignment to detect its format */

  iformat = formatInputAlignment(alignmentFile);
  oformat = iformat;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Depending on the input format, the program use an
   * appropiate function to read that alignment */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  switch(iformat) {
    case 1:
      return loadClustalAlignment(alignmentFile);
    case 3:
      return loadNBRF_PirAlignment(alignmentFile);
    case 8:
      return loadFastaAlignment(alignmentFile);
    case 11:
      return loadPhylip3_2Alignment(alignmentFile);
    case 12:
      return loadPhylipAlignment(alignmentFile);
    case 17:
      return loadNexusAlignment(alignmentFile);
    case 21:
      return loadMegaInterleavedAlignment(alignmentFile);
    case 22:
      return loadMegaNonInterleavedAlignment(alignmentFile);
    default:
      return false;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
  }
  return false;
}
/* ***** ***** ***** ***** ***** ***** ***** ***** */

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This function returns the input alignment format */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
int alignment::formatInputFile(void) {

  return iformat;
}
/* ***** ***** ***** ***** ***** ***** ***** ***** */

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This function returns the input alignment type between a single or a
 *  multialignment. */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
int alignment::typeInputFile(void) {

  return SINGLE;
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This function prints the alignmnet to the standard output using an
 * appropiate function depending on the output format */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
bool alignment::printAlignment(void){

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  if(sequences == NULL)
    return false;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  if((residNumber == 0)  || (sequenNumber == 0)) {
    cerr << endl << "WARNING: Output alignment has not been generated. "
      << "It is empty" << endl << endl;
    return true;
  }

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  switch(oformat) {
    case 1:
      alignmentClustalToFile(cout);
      break;
    case 3:
      alignmentNBRF_PirToFile(cout);
      break;
    case 8:
      alignmentFastaToFile(cout);
      break;
    case 11:
      alignmentPhylip3_2ToFile(cout);
      break;
    case 12:
      alignmentPhylipToFile(cout);
      break;
    case 13:
      alignmentPhylip_PamlToFile(cout);
      break;
    case 17:
      alignmentNexusToFile(cout);
      break;
    case 21: case 22:
      alignmentMegaToFile(cout);
      break;
    case 99:
      getSequences(cout);
      break;
    case 100:
      alignmentColourHTML(cout);
      break;
    default:
      return false;
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  return true;
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This function puts the alignment in a given file */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
bool alignment::saveAlignment(char *destFile) {

  ofstream file;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  if(sequences == NULL)
    return false;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  if((residNumber == 0)  || (sequenNumber == 0)) {
    cerr << endl << "WARNING: Output alignment has not been generated. "
      << "It is empty." << endl << endl;
    return true;
  }

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* File open and correct open check */
  file.open(destFile);
  if(!file) return false;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Depending on the output format, we call to the
   * appropiate function */
  switch(oformat) {
    case 1:
      alignmentClustalToFile(file);
      break;
    case 3:
      alignmentNBRF_PirToFile(file);
      break;
    case 8:
      alignmentFastaToFile(file);
      break;
    case 11:
      alignmentPhylip3_2ToFile(file);
      break;
    case 12:
      alignmentPhylipToFile(file);
      break;
    case 13:
      alignmentPhylip_PamlToFile(file);
      break;
    case 17:
      alignmentNexusToFile(file);
      break;
    case 21: case 22:
      alignmentMegaToFile(file);
      break;
    case 99:
      getSequences(file);
      break;
    case 100:
      alignmentColourHTML(file);
      break;
    default:
      return false;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  }

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Close the output file */
  file.close();

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* All is OK, return true */
  return true;
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This function trimms a given alignment based on the gap distribution
 * values. To trim the alignment, this function uses two parameters:
 *
 *      baseLine:   Minimum percentage of columns that should
 *                  be conserved in the new alignment.
 *      gapsPct:    Maximum percentage of gaps per column that
 *                  should be in the new alignment.
 *
 * The function selects that combination of parameters that maximize the
 * final number of columns in the new alignment. There is an extra parameter
 * to ask for the complementary alignment. A complementary alignment consits
 * of those columns that would delete applying the standard approach. */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
alignment *alignment::cleanGaps(float baseLine, float gapsPct, bool complementary) {

  alignment *ret;
  double cut;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* If gaps statistics are not calculated, we
   * calculate them */
  if(calculateGapStats() != true)
    return NULL;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Obtain the cut point using the given parameters */
  cut = sgaps -> calcCutPoint(baseLine, gapsPct);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Once we have the cut value proposed, we call the
   * appropiate method to clean the alignment and, then,
   * generate the new alignment. */
  ret = cleanByCutValue(cut, baseLine, sgaps -> getGapsWindow(), complementary);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Return a reference of the new alignment */
  return ret;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This function trimms a given alignment based on the similarity distribution
 * values. To trim the alignment, this function uses two parameters:
 *
 *      baseLine:   Minimum percentage of columns that should
 *                  be conserved in the new alignment.
 *      conservat:  Minimum value of similarity per column that
 *                  should be in the new alignment.
 *
 * The function selects that combination of parameters that maximize the
 * final number of columns in the new alignment. There is an extra parameter
 * to ask for the complementary alignment. A complementary alignment consits
 * of those columns that would delete applying the standard approach. */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
alignment *alignment::cleanConservation(float baseLine, float conservationPct, bool complementary) {

  alignment *ret;
  float cut;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* If conservation's statistics are not calculated,
   * we calculate them */
  if(calculateConservationStats() != true)
    return NULL;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Calculate the cut point using the given parameters */
  cut = scons -> calcCutPoint(baseLine, conservationPct);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Once we have the cut value, we call the appropiate
   * method to clean the alignment and, then, generate
     the new alignment */
  ret = cleanByCutValue(cut, baseLine, scons -> getMdkwVector(), complementary);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Return a reference of the new alignment */
  return ret;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This function trimms a given alignment based on the similarity and gaps
 * distribution values. To trim the alignment, this function uses three
 * parameters:
 *
 *      baseLine:   Minimum percentage of columns that should
 *                  be conserved in the new alignment.
 *      gapsPct:    Maximum percentage of gaps per column that
 *                  should be in the new alignment.
 *      conservat:  Minimum value of similarity per column that
 *                  should be in the new alignment.
 *
 * The function selects that combination of parameters that maximize the
 * final number of columns in the new alignment. If the baseLine parameter
 * is the most strict, the other two ones would be relaxed to keep columns
 * that would not be selected. There is an extra parameter to ask for the
 * complementary alignment. A complementary alignment consits of those
 * columns that would delete applying the standard approach. */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
alignment *alignment::clean(float baseLine, float GapsPct, float conservationPct, bool complementary) {

  alignment *ret;
  float cutCons;
  double cutGaps;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* If gaps statistics are not calculated, we calculate
   *  them */
  if(calculateGapStats() != true)
    return NULL;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* If conservation's statistics are not calculated,
   * we calculate them */
  if(calculateConservationStats() != true)
    return NULL;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Calculate the two cut points using the parameters */
  cutGaps = sgaps->calcCutPoint(baseLine, GapsPct);
  cutCons = scons->calcCutPoint(baseLine, conservationPct);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Clean the alingment using the two cut values, the
   * gapsWindow and MDK_Windows vectors and the baseline
   * value */
  ret = cleanByCutValue(cutGaps, sgaps -> getGapsWindow(), baseLine, cutCons, scons -> getMdkwVector(), complementary);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Return a reference of the clean alignment object */
  return ret;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This function cleans a given alignment based on the consistency values
 * asses from a dataset of alignments. The method takes three parameters:
 *
 *      cutpoint:   Lower limit (0-1) of comparefile value admits in
 *                  the new alignment.
 *      baseline:   Minimum percentage of columns to have in the new alignment
 *      vectValues: A vector with alignment's consistency values.
 *
 * The function computes the optimal parameter combination value to trim
 * an alignment based on the consistency value from the comparison among
 * a dataset of alignments with the same sequences. */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
alignment *alignment::cleanCompareFile(float cutpoint, float baseLine, float *vectValues, bool complementary) {

  alignment *ret;
  float cut, *vectAux;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Allocate memory */
  vectAux = new float[residNumber];
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Sort a copy of the vectValues vector, and take the
   *  value at 100% - baseline position. */
  utils::copyVect((float *) vectValues, vectAux, residNumber);
  utils::quicksort(vectAux, 0, residNumber-1);
  cut = vectAux[(int) ((float)(residNumber - 1) * (100.0 - baseLine)/100.0)];
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* We have to decide which is the smallest value
   * between the cutpoint value and the value from
   * the minimum percentage threshold */
  cut = cutpoint < cut ? cutpoint : cut;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Deallocate memory */
  delete [] vectAux;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Clean the selected alignment using the input parameters. */
  ret = cleanByCutValue(cut, baseLine, vectValues, complementary);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Return a refernce of the new alignment */
  return ret;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This function removes those sequences that are misaligned with the rest
 * of sequences in the alignment at given sequence and residue overlaps.
 *
 *      overlapColumn:  Set the minimum similarity fraction value that has
 *                      to have a position from a given sequence to be
 *                      considered as an hit.
 *      minimumOverlap: Set the minimum proportion of hits that has to be
 *                      a sequence to keep it in the new alignment.
 *
 * At the same time, it's possible to get an alignment with those columns
 * that this method should remove setting for that purpose the complementary
 * parameter to True */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
alignment *alignment::cleanSpuriousSeq(float overlapColumn, float minimumOverlap, bool complementary) {

  float *overlapVector;
  alignment *newAlig;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Allocate local memory */
  overlapVector = new float[sequenNumber];
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Compute the overlap's vector using the overlap
   * column's value */
  if(!calculateSpuriousVector(overlapColumn, overlapVector))
    return NULL;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Select and remove the sequences with a overlap less
   * than threshold's overlap and create a new alignemnt */
  newAlig = cleanOverlapSeq(minimumOverlap, overlapVector, complementary);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Deallocate local memory */
  delete [] overlapVector;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Return a reference of the clean alignment object */
  return newAlig;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This function implements the gappyout approach. To trim the alignment, this
 * method computes the slope from the gaps distribution from the input alig.
 * Using this information, the method asses the most abrupt change between
 * three consecutive points from that distribution and selects the first of
 * this point as aan cut-off point.
 *
 * The method also can return the complementary alignment that consists of
 * those columns that would be delete applying this algorithm. */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
alignment *alignment::clean2ndSlope(bool complementarity) {

  alignment *ret;
  int cut;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* If gaps statistics are not calculated, we calculate
   *  them */
  if(calculateGapStats() != true)
    return NULL;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* We get the cut point using a automatic method for
   * this purpose. */
  cut = sgaps -> calcCutPoint2ndSlope();
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Using the cut point calculates in last steps, we
   * clean the alignment and generate a new alignment */
  ret = cleanByCutValue(cut, 0, sgaps->getGapsWindow(), complementarity);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Returns the new alignment. */
  return ret;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This method implements the strict and strictplus approaches, to apply
 * these algorithms, the program computes the gaps and similarity distribution
 * values to decide which ones are the optimal cutoff points. To compute these
 * cutoff points, the method applies those steps described in the manual. */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
alignment *alignment::cleanCombMethods(bool complementarity, bool variable) {

  float simCut, first20Point, last80Point, *simil, *vectAux;
  int i, j, acm, gapCut, *positions, *gaps;
  double inic, fin, vlr;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* If gaps statistics are not calculated, we calculate
   *  them */
  if(calculateGapStats() != true)
    return NULL;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Computes the gap cut point using a automatic method
   * and at the same time, we get the gaps values from
   * the alignment. */
  gapCut = sgaps -> calcCutPoint2ndSlope();
  gaps = sgaps -> getGapsWindow();
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* If conservation's statistics are not calculated,
   * we calculate them */
  if(calculateConservationStats() != true)
    return NULL;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Computes the conservations value for each column
   * in the alignment. At the same time, the method get
   * the vector with those values. */
  scons -> calculateVectors(sequences, sgaps -> getGapsWindow());
  simil = scons -> getMdkwVector();
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Allocate local memory and initializate it to -1 */
  positions = new int[residNumber];
  utils::initlVect(positions, residNumber, -1);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* The method only selects columns with gaps number
   * less or equal than the gap's cut point. Counts the
   * number of columns that have been selected */
  for(i = 0, acm = 0; i < residNumber; i++) {
    if(gaps[i] <= gapCut) {
      positions[i] = i;
      acm++;
    }
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Allocate local memory and save the similaritys
   * values for the columns that have been selected */
  vectAux = new float[acm];
  for(i = 0, j = 0; i < residNumber; i++)
    if(positions[i] != -1)
      vectAux[j++] = simil[i];
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Sort the conservation's value vector. */
  utils::quicksort(vectAux, 0, acm-1);

  /* ...and search for the vector points at the 20 and
   * 80% of length. */
  first20Point = 0;
  last80Point  = 0;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  for(i = acm - 1, j = 1; i >= 0; i--, j++) {
    if((((float) j/acm) * 100.0) <= 20.0)
      first20Point = vectAux[i];
    if((((float) j/acm) * 100.0) <= 80.0)
      last80Point = vectAux[i];
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Computes the logaritmic's values for those points.
   * Finally the method computes the similarity cut
   * point using these values. */
  inic = log10(first20Point);
  fin  = log10(last80Point);
  vlr  = ((inic - fin) / 10) + fin;
  simCut = (float) pow(10, vlr);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Clean the alignment and generate a new alignment
   * object using the gaps cut and the similaritys cut
   *  values */
  alignment *ret = cleanStrict(gapCut, sgaps -> getGapsWindow(),
        simCut, scons -> getMdkwVector(), complementarity, variable);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Deallocate local memory */
  delete [] vectAux;
  delete [] positions;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Return a reference of the new alignment */
  return ret;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This function allows to delete those columns consistents of only gaps.
 * This function is specially useful when we remove misaligned sequences from
 * a given alignmnet. As the rest of functions, this function can return the
 * complementary alignment but this alignment would have only columns with
 * all gaps. */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
alignment *alignment::cleanNoAllGaps(bool complementarity) {

  alignment *ret;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* If gaps statistics are not calculated, we calculate
   *  them */
  if(calculateGapStats() != true)
    return NULL;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* We want to conserve the columns with gaps' number
   * less or equal than sequences' number - 1  */
  ret = cleanByCutValue((sequenNumber - 1), 0, sgaps->getGapsWindow(), complementarity);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Returns the new alignment. */
  return ret;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This method computes the overlap for each sequence given a overlap value.
 * This overlap set the minimum fraction that has to be a position for the
 * selected sequence to count as an hit. This proportion measures how much
 * is similar, (in terms of residues, indetermination and gaps), a given
 * position respect to the element on its column. */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
bool alignment::calculateSpuriousVector(float overlap, float *spuriousVector) {

  int i, j, k, seqValue, ovrlap, hit;
  float floatOverlap;
  char indet;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Compute the overlap */
  floatOverlap = overlap * float(sequenNumber-1);
  ovrlap = int(overlap * (sequenNumber-1));

  if(floatOverlap > float(ovrlap))
    ovrlap++;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* If the spurious vectos is NULL, returns false. */
  if(spuriousVector == NULL)
    return false;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Depending on the kind of alignment, we have
   * different indetermination symbol */
  if(getTypeAlignment() == AAType)
    indet = 'X';
  else
    indet = 'N';
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* For each alignment's sequence, computes its overlap */
  for(i = 0, seqValue = 0; i < sequenNumber; i++, seqValue = 0) {

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* For each alignment's column, computes the overlap
     * between the selected sequence and the other ones */
    for(j = 0; j < residNumber; j++) {

      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      /* For sequences are before the sequence selected */
      for(k = 0, hit = 0; k < i; k++) {

        /* ***** ***** ***** ***** ***** ***** ***** ***** */
        /* If the element of sequence selected is the same
         * that the element of sequence considered, computes
         * a hit */
        if(sequences[i][j] == sequences[k][j])
          hit++;

        /* If the element of sequence selected isn't a 'X' nor
         * 'N' (indetermination) or a '-' (gap) and the element
         * of sequence considered isn't a  a 'X' nor 'N'
         * (indetermination) or a '-' (gap), computes a hit */
        else if((sequences[i][j] != indet) && (sequences[i][j] != '-')
          && (sequences[k][j] != indet) && (sequences[k][j] != '-'))
          hit++;
      }
      /* ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      /* For sequences are after the sequence selected */
      for(k = (i + 1); k < sequenNumber; k++) {

        /* ***** ***** ***** ***** ***** ***** ***** ***** */
        /* If the element of sequence selected is the same
         * that the element of sequence considered, computes
         * a hit */
        if(sequences[i][j] == sequences[k][j])
          hit++;

        /* If the element of sequence selected isn't a 'X' nor
         * 'N' (indetermination) or a '-' (gap) and the element
         * of sequence considered isn't a  a 'X' nor 'N'
         * (indetermination) or a '-' (gap), computes a hit */
        else if((sequences[i][j] != indet) && (sequences[i][j] != '-')
          && (sequences[k][j] != indet) && (sequences[k][j] != '-'))
          hit++;
      }
      /* ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      /* Finally, if the hit's number divided by number of
       * sequences minus one is greater or equal than
       * overlap's value, computes a column's hit. */
      if(hit >= ovrlap)
        seqValue++;
      /* ***** ***** ***** ***** ***** ***** ***** ***** */
    }

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* For each alignment's sequence, computes its spurious's
     * or overlap's value as the column's hits -for that
       sequence- divided by column's number. */
    spuriousVector[i] = ((float) seqValue / residNumber);
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
  }

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* If there is not problem in the method, return true */
  return true;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This method compute the cluster belongs to a given alignment at certain
 * identity threshold. To know the clusters from the alignment, those
 * clusters are calculated as follow: 1) Select the longest sequences. 2)
 * Using the previous computed identity values, compare if the second
 * longest sequence belongs, because the identity value with the cluster
 * representative is equal or greater than a given threshold, to this cluster
 * or not. 3) If the sequence belongs to this cluster, we add it, if not
 * belongs, we create a new cluster and fix this sequence as its representative
 * 4) Continue with the rest of sequences. In the case that a given sequence
 * can belong to more than one cluster, we choose the cluster which one
 * maximize the identity value respect to its representative sequence */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
int *alignment::calculateRepresentativeSeq(float maximumIdent) {

  int i, j, pos, clusterNum, **seqs;
  int *cluster;
  static int *repres;
  float max;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Ask for the sequence identities assesment */
  if(identities == NULL)
    calculateSeqIdentity();
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  seqs = new int*[sequenNumber];
  for(i = 0; i < sequenNumber; i++) {
    seqs[i] = new int[2];
    seqs[i][0] = utils::removeCharacter('-', sequences[i]).size();
    seqs[i][1] = i;
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  utils::quicksort(seqs, 0, sequenNumber-1);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  cluster = new int[sequenNumber];
  cluster[0] = seqs[sequenNumber - 1][1];
  clusterNum = 1;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  for(i = sequenNumber - 2; i >= 0; i--) {

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    for(j = 0, max = 0, pos = -1; j < clusterNum; j++) {
      if(identities[seqs[i][1]][cluster[j]] > maximumIdent) {
        if(identities[seqs[i][1]][cluster[j]] > max) {
          max = identities[seqs[i][1]][cluster[j]];
          pos = j;
        }
      }
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    if(pos == -1) {
    cluster[j] = seqs[i][1];
      clusterNum++;
    }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  }

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  repres = new int[clusterNum + 1];
  repres[0] = clusterNum;
  for(i = 0; i < clusterNum; i++)
    repres[i+1] = cluster[i];
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Deallocate dinamic memory */
  for(i = 0; i < sequenNumber; i++)
  delete [] seqs[i];

  delete [] cluster;
  delete [] seqs;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  return repres;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This method looks for the optimal cut point for a given clusters number.
 * The idea is to find a identity cut-off point that can be used to get a
 * number of representative sequences similar to the input parameter */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
float alignment::getCutPointClusters(int clusterNumber) {

  float max, min, avg, gMax, gMin, startingPoint, prevValue = 0, iter = 0;
  int i, j, clusterNum, *cluster, **seqs;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* If the user wants only one cluster means that all
   * of sequences have to be in the same clauster.
   * Otherwise, if the users wants the maximum number of
   * clusters means that each sequence have to be in their
   * own cluster */
  if(clusterNumber == sequenNumber) return 1;
  else if(clusterNumber == 1) return 0;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Ask for the sequence identities assesment */
  if(identities == NULL)
    calculateSeqIdentity();
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Compute the maximum, the minimum and the average
   * identity values from the sequences */
  for(i = 0,gMax = 0, gMin = 1, startingPoint = 0; i < sequenNumber; i++) {
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    for(j = 0, max = 0, avg = 0, min = 1; j < i; j++) {
      if(max < identities[i][j])
        max  = identities[i][j];
      if(min > identities[i][j])
        min  = identities[i][j];
      avg += identities[i][j];
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    for(j = i + 1; j < sequenNumber; j++) {
      if(max < identities[i][j])
        max  = identities[i][j];
      if(min > identities[i][j])
        min  = identities[i][j];
      avg += identities[i][j];
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    startingPoint += avg / (sequenNumber - 1);
    if(max > gMax) gMax = max;
    if(min < gMin) gMin = min;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Take the starting point as the average value */
  startingPoint /= sequenNumber;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Compute and sort the sequence length */
  seqs = new int*[sequenNumber];
  for(i = 0; i < sequenNumber; i++) {
    seqs[i] = new int[2];
    seqs[i][0] = utils::removeCharacter('-', sequences[i]).size();
    seqs[i][1] = i;
  }
  utils::quicksort(seqs, 0, sequenNumber-1);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Create the data structure to store the different
   * clusters for a given thresholds */
  cluster    = new int[sequenNumber];
  cluster[0] = seqs[sequenNumber - 1][1];
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /**** ***** ***** ***** ***** ***** ***** ***** */
  /* Look for the optimal identity value to get the
   * number of cluster set by the user. To do that, the
   * method starts in the average identity value and moves
   * to the right or lefe side of this value depending on
   * if this value is so strict or not. We set an flag to
   * avoid enter in an infinite loop, if we get the same
   * value for more than 10 consecutive point without get
   * the given cluster number means that it's impossible
   * to get this cut-off and we need to stop the search */
  do {
    clusterNum = 1;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Start the search */
    for(i = sequenNumber - 2; i >= 0; i--) {
      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      for(j = 0; j < clusterNum; j++)
        if(identities[seqs[i][1]][cluster[j]] > startingPoint)
          break;
      /* ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      if(j == clusterNum) {
        cluster[j] = seqs[i][1];
        clusterNum++;
      }
      /* ***** ***** ***** ***** ***** ***** ***** ***** */
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Given the cutoff point, if we get the same cluster
     * number or we have been getting for more than 10
     * consecutive times the same cutoff point, we stop */
    if((clusterNum == clusterNumber) || (iter > 10))
      break;
    else {
      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      /* We have to move to the left side from the average
       * to decrease the cutoff value and get a smaller number
       * of clusters */
      if(clusterNum > clusterNumber) {
        gMax = startingPoint;
        startingPoint = (gMax + gMin)/2;
      } else {
      /* In the opposite side, we have to move to the right
       * side from the cutting point to be more strict and get
       * a bigger number of clusters */
        gMin = startingPoint;
        startingPoint = (gMax + gMin)/2;
      }
      /* ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      /* If the clusters number from the previous iteration
       * is different from the current number, we put the
       * iteration number to 0 and store this new value */
      if(prevValue != clusterNum) {
        iter = 0;
        prevValue = clusterNum;
      /* Otherwise, we increase the iteration number to
       * avoid enter in an infinitve loop */
      } else iter++;
      /* ***** ***** ***** ***** ***** ***** ***** ***** */
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
  } while(true);

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Deallocate dinamic memory */
  for(i = 0; i < sequenNumber; i++) delete [] seqs[i];
  delete [] seqs;
  delete [] cluster;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  return startingPoint;
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This method select one representative sequence (the longest one) per each
 * cluster from the input alignment and generate a new alignment */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
alignment *alignment::getClustering(float identityThreshold) {

  string *matrixAux, *newSeqsName;
  int i, j, *clustering;
  alignment *newAlig;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Get the representative member for each cluster
   * given a maximum identity threshold */
  clustering = calculateRepresentativeSeq(identityThreshold);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Put all sequences to be deleted and get back those
   * sequences that are representative for each cluster
   * */
  for(i = 0; i < sequenNumber; i ++)
    saveSequences[i] = -1;
  for(i = 1; i <= clustering[0]; i ++)
    saveSequences[clustering[i]] = clustering[i];
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* We allocate memory to save the sequences selected */
  matrixAux = new string[clustering[0]];
  newSeqsName = new string[clustering[0]];

  /* Copy to new structures the information that have
   * been selected previously. */
  for(i = 0, j = 0; i < sequenNumber; i++)
    if(saveSequences[i] != -1) {
       newSeqsName[j] = seqsName[i];
       matrixAux[j] = sequences[i];
       j++;
    }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* When we have all parameters, we create the new
   * alignment */
  newAlig = new alignment(filename, aligInfo, matrixAux, newSeqsName, seqsInfo,
    clustering[0], residNumber, iformat, oformat, shortNames, dataType, isAligned,
    reverse, terminalGapOnly, keepSequences, keepHeader, sequenNumber, residNumber, residuesNumber,
    saveResidues, saveSequences, ghWindow, shWindow, blockSize, identities);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Deallocated auxiliar memory */
  delete [] matrixAux;
  delete [] newSeqsName;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Return the new alignment reference */
  return newAlig;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This function returns the backtranslation for a given protein processed
 * alignment into its CDS alignment. To do this convertion, the function needs
 * the Coding sequences as well the original composition of each protein
 * sequence. Also, the function needs to know which columns/sequences will be
 * in the final alignment to carray out the conversion. */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
alignment *alignment::getTranslationCDS(int newResidues, int newSequences, int *ColumnsToKeep, string *oldSeqsName, sequencesMatrix *seqMatrix, alignment *ProtAlig) {

  string *matrixAux;
  alignment *newAlig;
  int i, j, k, l, oldResidues;
  int *mappedSeqs, *tmpSequence, *selectedRes;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Map the selected protein sequences to the input
   * coding sequences */
  mappedSeqs = new int[newSequences];
  for(i = 0; i < sequenNumber; i++)
    for(j = 0; j < newSequences; j++)
    if(!seqsName[i].compare(oldSeqsName[j])) {
      mappedSeqs[j] = i;
      break;
    }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Get the original alignment size as well the original
   * selected/non-selected columns */
  oldResidues = seqMatrix -> getResidNumber();
  selectedRes = new int[oldResidues];

  for(i = 0; i < oldResidues; i++)
    selectedRes[i] = i;

  for(j = 0; j < ColumnsToKeep[0]; j++)
    selectedRes[j] = -1;

  for(i = 0; i < newResidues - 1; i++)
    for(j = ColumnsToKeep[i] + 1; j < ColumnsToKeep[i+1]; j++)
      selectedRes[j] = -1;

  for(j = ColumnsToKeep[newResidues - 1] + 1; j < oldResidues; j++)
    selectedRes[j] = -1;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* We allocate memory to save the columns selected */
  matrixAux = new string[newSequences];
  tmpSequence = new int[oldResidues];

  /* Using the information about which residues for each
   * sequence was selected by others function, we process
   * these residues to recover the corresponding codons */
  for(i = 0; i < newSequences; i++)
  if(seqMatrix -> getSequence(oldSeqsName[i], tmpSequence)) {
    for(j = 0; j < oldResidues; j++) {
    if((selectedRes[j] != -1) && (tmpSequence[j] != 0)) {

      for(k = 3 * (tmpSequence[j] - 1), l = 0; l < 3; k++, l++) {
        /* Check whether the nucleotide sequences end has been reached or not.
         * If it has been reached, complete backtranslation using indetermination
         * symbols 'N' */
        if((int) sequences[mappedSeqs[i]].length() > k)
          matrixAux[i].resize(matrixAux[i].size() + 1, sequences[mappedSeqs[i]][k]);
        else
          matrixAux[i].resize(matrixAux[i].size() + 1, 'N');
    }
    } else if(selectedRes[j] != -1) {
          matrixAux[i].resize(matrixAux[i].size() + 3, '-');
    }
    }
  /* If there is any problems with a sequence then
   * the function returns an error */
  } else {
    return NULL;
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */


  /* Remove columns/sequences composed only by gaps before making the retrotranslation */


  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* When we have all parameters, we create the new
   * alignment */
  newAlig = new alignment(filename, "", matrixAux, oldSeqsName, NULL, newSequences,
    newResidues * 3, ProtAlig -> getInputFormat(), ProtAlig -> getOutputFormat(),
    ProtAlig -> getShortNames(), DNAType, true, ProtAlig -> getReverse(),
    terminalGapOnly, keepSequences, keepHeader, sequenNumber, oldResidues * 3, NULL, NULL,
    NULL, 0, 0, ProtAlig -> getBlockSize(), NULL);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Deallocated auxiliar memory */
  delete [] matrixAux;
  delete mappedSeqs;
  delete tmpSequence;
  delete selectedRes;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Return the new alignment reference */
  return newAlig;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This function computes the gaps statistics for the input alignment. */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
bool alignment::calculateGapStats(void) {

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* If alignment matrix is not created, return false */
  if(sequences == NULL)
    return false;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* If sgaps object is not created, we create them
     and calculate the statistics */
  if(sgaps == NULL) {
    sgaps = new statisticsGaps(sequences, sequenNumber, residNumber, dataType);
    sgaps -> applyWindow(ghWindow);
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  return true;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Print the gaps value for each column in the alignment */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void alignment::printStatisticsGapsColumns(void) {

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Check if there is computed the gaps statistics */
  if(calculateGapStats())
    /* then prints the information */
    sgaps -> printGapsColumns();
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Print the acumulative gaps distribution value from the input alignment */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void alignment::printStatisticsGapsTotal(void) {

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Check it there is computed the gaps statistics */
  if(calculateGapStats())
    /* then prints the information */
    sgaps -> printGapsAcl();
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Set the similarity matrix. This matrix is necessary for some methods in
 * the program */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
bool alignment::setSimilarityMatrix(similarityMatrix *sm) {

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* If scons object is not created, we create them */
  if(scons == NULL)
    scons = new statisticsConservation(sequences, sequenNumber, residNumber, dataType);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Associate the matrix to the similarity statistics
   * object */
  if(!scons -> setSimilarityMatrix(sm))
    return false;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* If it's OK, we return true */
  return true;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This method compute the similarity values from the input alignment if it
 * has not been computed before */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
bool alignment::calculateConservationStats(void) {

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* It the gaps statistics object has not been created
   * we create it */
  if(calculateGapStats() != true)
    return false;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* It the similarity statistics object has not been
   * created we create it */
  if(scons == NULL)
    return false;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Ask for the similarity matrix */
  if(scons -> isSimMatrixDef() != true)
    return false;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Compute the similarity statistics from the input
   * alignment */
  if(!scons -> calculateVectors(sequences, sgaps->getGapsWindow()))
    return false;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Ask to know if it is necessary to apply any window
   * method. If it's necessary, we apply it */
  if(scons->isDefinedWindow())
    return true;
  else
    return scons->applyWindow(shWindow);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Print the similarity value for each column from the alignment */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void alignment::printStatisticsConservationColumns(void) {

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Check if the similarity statistics object has been
   * created */
  if(calculateConservationStats())
    /* then prints the information */
    scons -> printConservationColumns();
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Print the accumulative similarity distribution values from the alignment */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void alignment::printStatisticsConservationTotal(void) {

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Check if the similarity statistics object has been
   * created */
  if(calculateConservationStats())
    /* then prints the information */
    scons -> printConservationAcl();
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This method prints the correspondece between the columns in the original
 * and in the trimmed alignment */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void alignment::printCorrespondence(void) {
  int i;

  cout << "#ColumnsMap\t";
  /* Print the saveResidues relathionship */
  for(i = 0; i < residNumber - 1; i++)
    cout << saveResidues[i] << ", ";
  cout << saveResidues[i] << endl;

}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Return the number of the sequences in the alignment */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
int alignment::getNumSpecies(void) {
  return sequenNumber;
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Return the number of residues in the alignment */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
int alignment::getNumAminos(void) {
  return residNumber;
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Sets the condition to remove only terminal gaps after applying any
 * trimming method or not.   */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void alignment::trimTerminalGaps(bool terminalOnly_) {
  terminalGapOnly = terminalOnly_;
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This method stores the diferent windows values in the alignment object */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void alignment::setWindowsSize(int ghWindow_, int shWindow_) {
  ghWindow = ghWindow_;
  shWindow = shWindow_;
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This method lets to change the output format */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void alignment::setOutputFormat(int format_, bool shortNames_) {
  oformat    = format_;
  shortNames = shortNames_;
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This method sets a new block size value */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void alignment::setBlockSize(int blockSize_) {
  blockSize = blockSize_;
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Return the input format aligment */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
int alignment::getInputFormat(void) {
  return iformat;
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Return the output format alignment */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
int alignment::getOutputFormat(void) {
  return oformat;
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Return true if the shortNames flag has been setted */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
int alignment::getShortNames(void) {
  return shortNames;
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Return true if the reverse flag has been setted. */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
int alignment::getReverse(void) {
  return reverse;
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This method lets to change the output alignment orientation */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void alignment::setReverse(void) {
  reverse = true;
}


/* Set appropiate flag to decide whether sequences composed only by gaps should
 * be kept or not */
void alignment::setKeepSequencesFlag(bool flag) {
  keepSequences = flag;
}

/* Set appropiate flag to decide whether original sequences header should be
 * dumped into the output stream with or without any preprocessing step */
void alignment::setKeepSeqsHeaderFlag(bool flag) {
  keepHeader = flag;
}
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Return the block size value */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
int alignment::getBlockSize(void) {
  return blockSize;
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This method returns the sequences name */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void alignment::getSequences(string *Names) {
  for(int i = 0; i < sequenNumber; i++)
    Names[i] = seqsName[i];
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void alignment::getSequences(string *Names, int *lengths) {
  for(int i = 0; i < sequenNumber; i++) {
    lengths[i] = (int) utils::removeCharacter('-', sequences[i]).length();
    Names[i] = seqsName[i];
  }
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void alignment::getSequences(string *Names, string *Sequences, int *Lengths) {
  for(int i = 0; i < sequenNumber; i++) {
    Names[i] = seqsName[i];
    Sequences[i] = utils::removeCharacter('-', sequences[i]);
    Lengths[i] = (int) Sequences[i].length();
  }
}


/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* For a given set of sequences name, check for its correspondence with
 * the current sequences name. If there is not a complete correspondence
 * between both sets, the method return false, otherwise, return true */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
bool alignment::getSeqNameOrder(string *names, int *orderVector) {
  int i, j, numNames;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* For each name in the current alignment, we look
   * for its correspondence in the input set */
  for(i = 0, numNames = 0; i < sequenNumber; i++) {
    for(j = 0; j < sequenNumber; j++) {
      if(seqsName[i] == names[j]) {
        orderVector[i] = j;
        numNames++;
        break;
      }
    }
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Depending on if we get the complete correspondence
   * between both sets of names, we return true or not */
  if(numNames == sequenNumber) return true;
  else return false;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This method lets to build a sequence matrix. A sequence matrix contains
 * the residue position for each sequence without taking into account the
 * gaps in the sequence. This means that at position 10 we have the residue
 * 1 for that sequence */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void alignment::sequenMatrix(void) {
  if(seqMatrix == NULL)
    seqMatrix = new sequencesMatrix(sequences, seqsName, sequenNumber, residNumber);
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Use this method to destroy a given sequence matrix */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void alignment::destroySequenMatrix(void) {
  if(seqMatrix != NULL)
    delete seqMatrix;
  seqMatrix = NULL;
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Print the alignment's sequence matrix */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void alignment::printSequenMatrix(void) {
  if(seqMatrix == NULL)
    seqMatrix = new sequencesMatrix(sequences, seqsName, sequenNumber, residNumber);
  seqMatrix -> printMatrix();
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Given an index, this method returns the sequence matrix column for that
 * index */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void alignment::getColumnSeqMatrix(int column, int *columnSeqMatrix) {
  if(seqMatrix == NULL)
    seqMatrix = new sequencesMatrix(sequences, seqsName, sequenNumber, residNumber);
  seqMatrix -> getColumn(column, columnSeqMatrix);
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This function looks if a given value belongs a given row. In the affirmative
 * case, returns the columns value for that row and value combination */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void alignment::getColumnSeqMatrix(int value, int row, int *columnSeqMatrix) {
  if(seqMatrix == NULL)
    seqMatrix = new sequencesMatrix(sequences, seqsName, sequenNumber, residNumber);
  seqMatrix -> getColumn(value, row, columnSeqMatrix);
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This function change the rows order in the sequence matrix */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
void alignment::setSeqMatrixOrder(int *order) {
  if(seqMatrix == NULL)
    seqMatrix = new sequencesMatrix(sequences, seqsName, sequenNumber, residNumber);
  seqMatrix -> setOrder(order);
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* This function change the rows order in the sequence matrix */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
sequencesMatrix *alignment::getSeqMatrix(void) {
  if(seqMatrix == NULL)
    seqMatrix = new sequencesMatrix(sequences, seqsName, sequenNumber, residNumber);
  return seqMatrix;
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Computes, if it's necessary, and return the alignment's type */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
int alignment::getTypeAlignment(void) {
  if(dataType == 0)
    dataType = utils::checkTypeAlignment(sequenNumber, residNumber, sequences);
  return dataType;
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Returns the correspondence between the columns in the original and in the
 * trimmed alignment */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
int *alignment::getCorrespResidues(void) {
  return saveResidues;
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Returns the correspondence between the sequences in the original and in
 * the trimmed alignment */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
int *alignment::getCorrespSequences(void) {
  return saveSequences;
}

/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
/* Returns if the alignment is aligned or not */
/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
bool alignment::isFileAligned(void) {
  return isAligned;
}

/* *****************************************************************************
 * *****************************************************************************
 * *****************************************************************************
 * ************************************************************************** */

/* This method removes those columns that exceed a given threshold. If the
 * number of columns in the trimmed alignment is lower than a given percentage
 * of the original alignment. The program relaxes the threshold until to add
 * enough columns to achieve this percentage, to add those columns the program
 * start in the middle of the original alignment and make, alternatively,
 * movements to the left and right side from that point. This method also can
 * return the complementary alignment. */
alignment *alignment::cleanByCutValue(double cut, float baseLine,
  const int *gInCol, bool complementary) {

  int i, j, k, jn, oth, pos, block, *vectAux;
  string *matrixAux, *newSeqsName;
  alignment *newAlig;
  newValues counter;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Select the columns with a gaps value less or equal
   * than the cut point. */
  for(i = 0, counter.residues = 0; i < residNumber; i++)
    if(gInCol[i] <= cut) counter.residues++;
    else saveResidues[i] = -1;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Compute, if it's necessary, the number of columns
   * necessary to achieve the minimum number of columns
   * fixed by coverage parameter. */
  oth = utils::roundInt((((baseLine/100.0) - (float) counter.residues/residNumber)) * residNumber);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* if it's necessary to recover some columns, we
   * applied this instructions to recover it */
  if(oth > 0) {
    counter.residues += oth;

    /* Allocate memory */
    vectAux = new int[residNumber];

    /* Sort a copy of the gInCol vector, and take the value of the column that marks the % baseline */
    utils::copyVect((int *) gInCol, vectAux, residNumber);
    utils::quicksort(vectAux, 0, residNumber-1);
    cut = vectAux[(int) ((float)(residNumber - 1) * (baseLine)/100.0)];

    /* Deallocate memory */
    delete [] vectAux;
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Fixed the initial size of blocks as 0.5% of
   * alignment's length */
  for(k = utils::roundInt(0.005 * residNumber); (k >= 0) && (oth > 0); k--) {

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* We start in the alignment middle then we move on
     * right and left side at the same time. */
    for(i = (residNumber/2), j = (i + 1); (((i > 0) || (j < (residNumber - 1))) && (oth > 0)); i--, j++) {

      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      /* Left side. Here, we compute the block's size. */
      for(jn = i; ((saveResidues[jn] != -1) && (jn >= 0) && (oth > 0)); jn--) ;

      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      /* if block's size is greater or equal than the fixed
       * size then we save all columns that have not been
       * saved previously. */
      if((i - jn) >= k) {
        for( ; ((saveResidues[jn] == -1) && (jn >= 0) && (oth > 0)); jn--) {
          if(gInCol[jn] <= cut) {
            saveResidues[jn] = jn;
            oth--;
          } else
            break;
        }
      }
      i = jn;
      /* ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      /* Right side. Here, we compute the block's size. */
      for(jn = j; ((saveResidues[jn] != -1) && (jn < residNumber) && (oth > 0)); jn++) ;

      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      /* if block's size is greater or equal than the fixed
       * size then we save all columns that have not been
       * saved previously. */
      if((jn - j) >= k) {
        for( ; ((saveResidues[jn] == -1) && (jn < residNumber) && (oth > 0)); jn++) {
          if(gInCol[jn] <= cut) {
            saveResidues[jn] = jn;
            oth--;
          } else
            break;
        }
      }
      j = jn;
      /* ***** ***** ***** ***** ***** ***** ***** ***** */
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
  }

  /* Keep only columns blocks bigger than an input columns block size */
  if(blockSize != 0) {

    /* Traverse all alignment looking for columns blocks greater than LONGBLOCK,
     * everytime than a column is not selected by the trimming method, check
     * whether the current block size until that point is big enough to be kept
     * or it should be removed from the final alignment */
    for(i = 0, pos = 0, block = 0; i < residNumber; i++) {
      if(saveResidues[i] != -1)
        block++;
      else {
        /* Remove columns from blocks smaller than input blocks size */
        if(block < blockSize)
          for(j = pos; j <= i; j++)
            saveResidues[j] = -1;
        pos = i + 1;
        block = 0;
      }
    }
    /* Check final block separately since it could happen than last block is not
     * big enough but because the loop end could provoke to ignore it */
    if(block < blockSize)
      for(j = pos; j < i; j++)
        saveResidues[j] = -1;
  }

  /* If the flag -terminalony is set, apply a method to look for internal
   * boundaries and get back columns inbetween them, if they exist */
  if(terminalGapOnly == true)
    if(!removeOnlyTerminal())
      return NULL;

  /* Once the columns/sequences selection is done, turn it around
   * if complementary flag is active */
  if(complementary == true)
    computeComplementaryAlig(true, false);

  /* Check for any additional column/sequence to be removed */
  /* Compute new sequences and columns numbers */
  counter = removeCols_SeqsAllGaps();

  /* Allocate memory  for selected sequences/columns */
  matrixAux = new string[counter.sequences];
  newSeqsName = new string[counter.sequences];

  /* Fill local allocated memory with previously selected data */
  fillNewDataStructure(matrixAux, newSeqsName);

  /* When we have all parameters, we create the new alignment */
  newAlig = new alignment(filename, aligInfo, matrixAux, newSeqsName, seqsInfo,
    counter.sequences, counter.residues, iformat, oformat, shortNames, dataType,
    isAligned, reverse, terminalGapOnly, keepSequences, keepHeader, sequenNumber, residNumber,
    residuesNumber, saveResidues, saveSequences, ghWindow, shWindow, blockSize,
    identities);

  /* Deallocate local memory */
  delete[] matrixAux;
  delete[] newSeqsName;

  /* Return the new alignment reference */
  return newAlig;
}

/* This method removes those columns that not achieve a given threshold. If the
 * number of columns in the trimmed alignment is lower than a given percentage
 * of the original alignment. The program relaxes the threshold until to add
 * enough columns to achieve this percentage, to add those columns the program
 * start in the middle of the original alignment and make, alternatively,
 * movements to the left and right side from that point. This method also can
 * return the complementary alignment. */
alignment *alignment::cleanByCutValue(float cut, float baseLine,
  const float *ValueVect, bool complementary) {

  int i, j, k, jn, oth, pos, block;
  string *matrixAux, *newSeqsName;
  alignment *newAlig;
  newValues counter;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Select the columns with a conservation's value
   * greater than the cut point. */
  for(i = 0, counter.residues = 0; i < residNumber; i++)
    if(ValueVect[i] > cut) counter.residues++;
    else saveResidues[i] = -1;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Compute, if it's necessary, the number of columns
   * necessary to achieve the minimum number of columns
   * fixed by coverage value. */
  oth = utils::roundInt((((baseLine/100.0) - (float) counter.residues/residNumber)) * residNumber);
  if(oth > 0) counter.residues += oth;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Fixed the initial size of blocks as 0.5% of
   * alignment's length */
  for(k = utils::roundInt(0.005 * residNumber); (k >= 0) && (oth > 0); k--) {

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* We start in the alignment middle then we move on
     * right and left side at the same time. */
    for(i = (residNumber/2), j = (i + 1); (((i > 0) || (j < (residNumber - 1))) && (oth > 0)); i--, j++) {

      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      /* Left side. Here, we compute the block's size. */
      for(jn = i; ((saveResidues[jn] != -1) && (jn >= 0) && (oth > 0)); jn--) ;

       /* ***** ***** ***** ***** ***** ***** ***** ***** */
      /* if block's size is greater or equal than the fixed
       * size then we save all columns that have not been
       * saved previously. */
      /* Here, we only accept column with a conservation's
       * value equal to conservation cut point. */
      if((i - jn) >= k) {
        for( ; ((saveResidues[jn] == -1) && (jn >= 0) && (oth > 0)); jn--) {
          if(ValueVect[jn] == cut) {
            saveResidues[jn] = jn;
            oth--;
          } else
            break;
        }
      }
      i = jn;
      /* ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      /* Right side. Here, we compute the block's size. */
      for(jn = j; ((saveResidues[jn] != -1) && (jn < residNumber) && (oth > 0)); jn++) ;

      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      /* if block's size is greater or equal than the fixed
       * size then we select the column and save the block's
       * size for the next iteraction. it's obvius that we
       * decrease the column's number needed to finish. */
     /* Here, we only accept column with a conservation's
      * value equal to conservation cut point. */
      if((jn - j) >= k) {
        for( ; ((saveResidues[jn] == -1) && (jn < residNumber) && (oth > 0)); jn++) {
          if(ValueVect[jn] == cut) {
            saveResidues[jn] = jn;
            oth--;
          } else
            break;
        }
      }
      j = jn;
      /* ***** ***** ***** ***** ***** ***** ***** ***** */
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
  }

  /* Keep only columns blocks bigger than an input columns block size */
  if(blockSize != 0) {

    /* Traverse all alignment looking for columns blocks greater than LONGBLOCK,
     * everytime than a column is not selected by the trimming method, check
     * whether the current block size until that point is big enough to be kept
     * or it should be removed from the final alignment */
    for(i = 0, pos = 0, block = 0; i < residNumber; i++) {
      if(saveResidues[i] != -1)
        block++;
      else {
        /* Remove columns from blocks smaller than input blocks size */
        if(block < blockSize)
          for(j = pos; j <= i; j++)
            saveResidues[j] = -1;
        pos = i + 1;
        block = 0;
      }
    }
    /* Check final block separately since it could happen than last block is not
     * big enough but because the loop end could provoke to ignore it */
    if(block < blockSize)
      for(j = pos; j < i; j++)
        saveResidues[j] = -1;
  }

  /* If the flag -terminalony is set, apply a method to look for internal
   * boundaries and get back columns inbetween them, if they exist */
  if(terminalGapOnly == true)
    if(!removeOnlyTerminal())
      return NULL;

  /* Once the columns/sequences selection is done, turn it around
   * if complementary flag is active */
  if(complementary == true)
    computeComplementaryAlig(true, false);

  /* Check for any additional column/sequence to be removed */
  /* Compute new sequences and columns numbers */
  counter = removeCols_SeqsAllGaps();

  /* Allocate memory  for selected sequences/columns */
  matrixAux = new string[counter.sequences];
  newSeqsName = new string[counter.sequences];

  /* Fill local allocated memory with previously selected data */
  fillNewDataStructure(matrixAux, newSeqsName);

  /* When we have all parameters, we create the new alignment */
  newAlig = new alignment(filename, aligInfo, matrixAux, newSeqsName, seqsInfo,
    counter.sequences, counter.residues, iformat, oformat, shortNames, dataType,
    isAligned, reverse, terminalGapOnly, keepSequences, keepHeader, sequenNumber, residNumber,
    residuesNumber, saveResidues, saveSequences, ghWindow, shWindow, blockSize,
    identities);

  /* Deallocate local memory */
  delete[] matrixAux;
  delete[] newSeqsName;

  /* Return the new alignment reference */
  return newAlig;
}

/* This method removes those columns that not achieve the similarity threshond,
 * by one hand, and exceed the gaps threshold. If the number of columns that
 * have been selected is lower that a given coverage, the program relaxes both
 * thresholds in order to get back some columns and achieve this minimum
 * number of columns in the new alignment. For that purpose, the program
 * starts at the middle of the original alignment and makes, alternatively,
 * movememts to the right and left sides looking for those columns necessary
 * to achieve the minimum number of columns set by the coverage parameter.
 * This method also can return the complementary alignmnet consists of those
 * columns that will be deleted from the original alignment applying the
 * standard method. */
alignment *alignment::cleanByCutValue(double cutGaps, const int *gInCol,
  float baseLine, float cutCons, const float *MDK_Win, bool complementary) {

  int i, j, k, oth, pos, block, jn, blGaps, *vectAuxGaps;
  string *matrixAux, *newSeqsName;
  float blCons, *vectAuxCons;
  alignment *newAlig;
  newValues counter;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Select the columns with a conservation's value
   * greater than the conservation cut point AND
   * less or equal than the gap cut point. */
  for(i = 0, counter.residues = 0; i < residNumber; i++)
    if((MDK_Win[i] > cutCons) && (gInCol[i] <= cutGaps)) counter.residues++;
    else saveResidues[i] = -1;
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Compute, if it's necessary, the number of columns
   * necessary to achieve the minimum number of it fixed
   * by the coverage parameter. */
  oth = utils::roundInt((((baseLine/100.0) - (float) counter.residues/residNumber)) * residNumber);
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* If it's needed to add new columns, we compute the
   * news thresholds */
  if(oth > 0) {
    counter.residues += oth;

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Allocate memory */
    vectAuxCons = new float[residNumber];
    vectAuxGaps = new int[residNumber];
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Sort a copy of the MDK_Win vector and of the gInCol
     * vector, and take the value of the column that marks
     * the % baseline */
    utils::copyVect((float *) MDK_Win, vectAuxCons, residNumber);
    utils::copyVect((int *) gInCol, vectAuxGaps, residNumber);

    utils::quicksort(vectAuxCons, 0, residNumber-1);
    utils::quicksort(vectAuxGaps, 0, residNumber-1);

    blCons = vectAuxCons[(int) ((float)(residNumber - 1) * (100.0 - baseLine)/100.0)];
    blGaps = vectAuxGaps[(int) ((float)(residNumber - 1) * (baseLine)/100.0)];
    /* ***** ***** ***** ***** ***** ***** ***** ***** */

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* Deallocate memory */
    delete [] vectAuxCons;
    delete [] vectAuxGaps;
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
  }

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  /* Fixed the initial size of blocks as 0.5% of
   * alignment's length */
  for(k = utils::roundInt(0.005 * residNumber); (k >= 0) && (oth > 0); k--) {

    /* ***** ***** ***** ***** ***** ***** ***** ***** */
    /* We start in the alignment middle then we move on
     * right and left side at the same time. */
    for(i = (residNumber/2), j = (i + 1); (((i > 0) || (j < (residNumber - 1))) && (oth > 0)); i--, j++) {

      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      /* Left side. Here, we compute the block's size. */
      for(jn = i; ((saveResidues[jn] != -1) && (jn >= 0) && (oth > 0)); jn--) ;

      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      /* if block's size is greater or equal than the fixed
       * size then we select the column and save the block's
       * size for the next iteraction. it's obvius that we
       * decrease the column's number needed to finish. */
      /* Here, we accept column with a conservation's value
       * greater or equal than the conservation cut point OR
       * less or equal than the gap cut point. */
      if((i - jn) >= k) {
        for( ; ((saveResidues[jn] == -1) && (jn >= 0) && (oth > 0)); jn--) {
          if((MDK_Win[jn] >= blCons) || (gInCol[jn] <= blGaps)) {
            saveResidues[jn] = jn;
            oth--;
          } else
            break;
        }
      }
      i = jn;
      /* ***** ***** ***** ***** ***** ***** ***** ***** */

      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      /* Right side. Here, we compute the block's size. */
      for(jn = j; ((saveResidues[jn] != -1) && (jn < residNumber) && (oth > 0)); jn++) ;

      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      /* if block's size is greater or equal than the fixed
       * size then we select the column and save the block's
       * size for the next iteraction. it's obvius that we
       * decrease the column's number needed to finish. */
      /* Here, we accept column with a conservation's value
       * greater or equal than the conservation cut point OR
       * less or equal than the gap cut point. */
      if((jn - j) >= k) {
        for( ; ((saveResidues[jn] == -1) && (jn < residNumber) && (oth > 0)); jn++) {
          if((MDK_Win[jn] >= blCons) || (gInCol[jn] <= blGaps)) {
            saveResidues[jn] = jn;
            oth--;
          } else
            break;
        }
      }
      j = jn;
      /* ***** ***** ***** ***** ***** ***** ***** ***** */
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** */
  }

  /* Keep only columns blocks bigger than an input columns block size */
  if(blockSize != 0) {

    /* Traverse all alignment looking for columns blocks greater than LONGBLOCK,
     * everytime than a column is not selected by the trimming method, check
     * whether the current block size until that point is big enough to be kept
     * or it should be removed from the final alignment */
    for(i = 0, pos = 0, block = 0; i < residNumber; i++) {
      if(saveResidues[i] != -1)
        block++;
      else {
        /* Remove columns from blocks smaller than input blocks size */
        if(block < blockSize)
          for(j = pos; j <= i; j++)
            saveResidues[j] = -1;
        pos = i + 1;
        block = 0;
      }
    }
    /* Check final block separately since it could happen than last block is not
     * big enough but because the loop end could provoke to ignore it */
    if(block < blockSize)
      for(j = pos; j < i; j++)
        saveResidues[j] = -1;
  }

  /* If the flag -terminalony is set, apply a method to look for internal
   * boundaries and get back columns inbetween them, if they exist */
  if(terminalGapOnly == true)
    if(!removeOnlyTerminal())
      return NULL;

  /* Once the columns/sequences selection is done, turn it around
   * if complementary flag is active */
  if(complementary == true)
    computeComplementaryAlig(true, false);

  /* Check for any additional column/sequence to be removed */
  /* Compute new sequences and columns numbers */
  counter = removeCols_SeqsAllGaps();

  /* Allocate memory  for selected sequences/columns */
  matrixAux = new string[counter.sequences];
  newSeqsName = new string[counter.sequences];

  /* Fill local allocated memory with previously selected data */
  fillNewDataStructure(matrixAux, newSeqsName);

  /* When we have all parameters, we create the new alignment */
  newAlig = new alignment(filename, aligInfo, matrixAux, newSeqsName, seqsInfo,
    counter.sequences, counter.residues, iformat, oformat, shortNames, dataType,
    isAligned, reverse, terminalGapOnly, keepSequences, keepHeader, sequenNumber, residNumber,
    residuesNumber, saveResidues, saveSequences, ghWindow, shWindow, blockSize,
    identities);

  /* Deallocate local memory */
  delete[] matrixAux;
  delete[] newSeqsName;

  /* Return the new alignment reference */
  return newAlig;
}

/* This method carries out the strict and strict plus method. To trim the
 * alignment in an automated way, the method uses as an input a given gaps
 * thresholds over a gaps vector values, a given similarity threshold over a
 * similarity vector values. With a flag, we can decide which method the
 * program shall apply. With another one, the user can ask for the
 * complementary alignment. In this method, those columns that has been marked
 * to be deleted but has, at least, 3 of 4 surronding neighbour selected are
 * get back to the new alignmnet. At other part of the method, those column
 * blocks that do not have a minimum size fix by the method will be deleted
 * from the trimmed alignment. */
alignment *alignment::cleanStrict(int gapCut, const int *gInCol, float simCut,
  const float *MDK_W, bool complementary, bool variable) {

  int i, num, lenBlock;
  alignment *newAlig;

  /* Reject columns with gaps number greater than the gap threshold. */
  for(i = 0; i < residNumber; i++)
    if(gInCol[i] > gapCut)
      saveResidues[i] = -1;

  /* Reject columns with similarity score under the threshold. */
  for(i = 0; i < residNumber; i++)
    if(MDK_W[i] < simCut)
      saveResidues[i] = -1;

  /* Search for those columns that have been removed but have,
   * at least, 3 adjacent columns selected. */

  /* Special cases:
   * Second column */
  if((saveResidues[0] != -1) && (saveResidues[2] != -1)
    && (saveResidues[3] != -1))
    saveResidues[1] = 1;
  else
    saveResidues[1] = -1;

  /* Second last column  */
  if((saveResidues[residNumber-1] != -1) && (saveResidues[residNumber-3] != -1)
    && (saveResidues[residNumber-4] != -1))
    saveResidues[(residNumber - 2)] = (residNumber - 2);
  else
    saveResidues[(residNumber - 2)] = -1;

  /* Normal cases */
  for(i = 2, num = 0; i < (residNumber - 2); i++, num = 0)
    if(saveResidues[i] == -1) {
      if(saveResidues[(i - 2)] != -1)
        num++;
      if(saveResidues[(i - 1)] != -1)
        num++;
      if(saveResidues[(i + 1)] != -1)
        num++;
      if(saveResidues[(i + 2)] != -1)
        num++;
      saveResidues[i] = (num >= 3) ? i : -1;
    }

  /* Select blocks size based on user input. It can be set either to 5 or to a
   * variable number between 3 and 12 depending on alignment length (1% alig) */
  if(!variable)
    lenBlock = 5;
  else {
    lenBlock = utils::roundInt(residNumber * 0.01);
    lenBlock = lenBlock >  3 ? (lenBlock < 12 ? lenBlock : 12) : 3;
  }

  /* Allow to change minimal block size */
  blockSize = blockSize > 0 ? blockSize : lenBlock;

  /* Keep only columns blocks bigger than either a computed dinamically or
   * set by the user block size */
  removeSmallerBlocks(blockSize);

  /* If the flag -terminalony is set, apply a method to look for internal
   * boundaries and get back columns inbetween them, if they exist */
  if(terminalGapOnly == true)
    if(!removeOnlyTerminal())
      return NULL;

  /* Once the columns/sequences selection is done, turn it around
   * if complementary flag is active */
  if(complementary == true)
    computeComplementaryAlig(true, false);

  /* Check for any additional column/sequence to be removed */
  /* Compute new sequences and columns numbers */

  newValues counter;
  counter = removeCols_SeqsAllGaps();

  /* Allocate memory  for selected sequences/columns */
  //~ matrixAux = new string[counter.sequences];
  //~ newSeqsName = new string[counter.sequences];

  /* Fill local allocated memory with previously selected data */
  //~ fillNewDataStructure(matrixAux, newSeqsName);

  fillNewDataStructure(&counter);

  //~ cerr << counter -> seqsName[counter -> sequences - 1] << endl;

  /* When we have all parameters, we create the new alignment */
  //~ newAlig = new alignment(filename, aligInfo, matrixAux, newSeqsName, seqsInfo, counter.sequences, counter.residues,
  newAlig = new alignment(filename, aligInfo, counter.matrix, counter.seqsName,
    seqsInfo, counter.sequences, counter.residues, iformat, oformat, shortNames,
    dataType, isAligned, reverse, terminalGapOnly, keepSequences, keepHeader, sequenNumber,
    residNumber, residuesNumber, saveResidues, saveSequences, ghWindow, shWindow,
    blockSize, identities);

  /* Deallocate local memory */
  //~ delete[] matrixAux;
  //~ delete[] newSeqsName;

  /* Return the new alignment reference */
  return newAlig;
}

/* Remove those sequences with an overlap less than a given threshold. It can
 * return the complementary alignment if appropiate flag is set. */
alignment *alignment::cleanOverlapSeq(float minimumOverlap, float *overlapSeq,
  bool complementary) {

  string *matrixAux, *newSeqsName;
  alignment *newAlig;
  newValues counter;
  int i;

  /* Keep only those sequences with an overlap value equal or greater than
   * the minimum overlap value set by the user.  */
  for(i = 0; i < sequenNumber; i++)
    if(overlapSeq[i] < minimumOverlap)
      saveSequences[i] = -1;

  /* Once the columns/sequences selection is done, turn it around
   * if complementary flag is active */
  if(complementary == true)
    computeComplementaryAlig(false, true);

  /* Check for any additional column/sequence to be removed */
  /* Compute new sequences and columns numbers */
  counter = removeCols_SeqsAllGaps();

  /* Allocate memory  for selected sequences/columns */
  matrixAux = new string[counter.sequences];
  newSeqsName = new string[counter.sequences];

  /* Fill local allocated memory with previously selected data */
  fillNewDataStructure(matrixAux, newSeqsName);

  /* When we have all parameters, we create the new alignment */
  newAlig = new alignment(filename, aligInfo, matrixAux, newSeqsName, seqsInfo,
    counter.sequences, counter.residues, iformat, oformat, shortNames, dataType,
    isAligned, reverse, terminalGapOnly, keepSequences, keepHeader, sequenNumber, residNumber,
    residuesNumber, saveResidues, saveSequences, ghWindow, shWindow, blockSize,
    identities);

  /* Deallocate local memory */
  delete [] matrixAux;
  delete [] newSeqsName;

  /* Return the new alignment reference */
  return newAlig;
}

/* Remove those columns, expressed as range, set by the user. It can return
 * the complementary alignmnet if appropiate flags is set. */
alignment *alignment::removeColumns(int *columns, int init, int size,
  bool complementary) {

  string *matrixAux, *newSeqsName;
  alignment *newAlig;
  newValues counter;
  int i, j;

  /* Delete those range columns defines in the columns vector */
  for(i = init; i < size + init; i += 2)
    for(j = columns[i]; j <= columns[i+1]; j++)
      saveResidues[j] = -1;

  /* Once the columns/sequences selection is done, turn it around
   * if complementary flag is active */
  if(complementary == true)
    computeComplementaryAlig(true, false);

  /* Check for any additional column/sequence to be removed */
  /* Compute new sequences and columns numbers */
  counter = removeCols_SeqsAllGaps();

  /* Allocate memory  for selected sequences/columns */
  matrixAux = new string[counter.sequences];
  newSeqsName = new string[counter.sequences];

  /* Fill local allocated memory with previously selected data */
  fillNewDataStructure(matrixAux, newSeqsName);

  /* When we have all parameters, we create the new alignment */
  newAlig = new alignment(filename, aligInfo, matrixAux, newSeqsName, seqsInfo,
    counter.sequences, counter.residues, iformat, oformat, shortNames, dataType,
    isAligned, reverse, terminalGapOnly, keepSequences, keepHeader, sequenNumber, residNumber,
    residuesNumber, saveResidues, saveSequences, ghWindow, shWindow, blockSize,
    identities);

  /* Deallocate local memory */
  delete[] matrixAux;
  delete[] newSeqsName;

  /* Return the new alignment reference */
  return newAlig;
}

/* This method removes those sequences, expressed as range of sequences, set by
 * the user. The method also can return the complementary alignment, it means,
 * those sequences that originally should be removed from the input alignment */
alignment *alignment::removeSequences(int *seqs, int init, int size,
  bool complementary) {

  string *matrixAux, *newSeqsName;
  alignment *newAlig;
  newValues counter;
  int i, j;

  /* Delete those range of sequences defines by the seqs vector */
  for(i = init; i < size + init; i += 2)
    for(j = seqs[i]; j <= seqs[i+1]; j++)
      saveSequences[j] = -1;

  /* Once the columns/sequences selection is done, turn it around
   * if complementary flag is active */
  if(complementary == true)
    computeComplementaryAlig(false, true);

  /* Check for any additional column/sequence to be removed */
  /* Compute new sequences and columns numbers */
  counter = removeCols_SeqsAllGaps();

  /* Allocate memory  for selected sequences/columns */
  matrixAux = new string[counter.sequences];
  newSeqsName = new string[counter.sequences];

  /* Fill local allocated memory with previously selected data */
  fillNewDataStructure(matrixAux, newSeqsName);

  /* When we have all parameters, we create the new alignment */
  newAlig = new alignment(filename, aligInfo, matrixAux, newSeqsName, seqsInfo,
    counter.sequences, counter.residues, iformat, oformat, shortNames, dataType,
    isAligned, reverse, terminalGapOnly, keepSequences, keepHeader, sequenNumber, residNumber,
    residuesNumber, saveResidues, saveSequences, ghWindow, shWindow, blockSize,
    identities);

  /* Free local memory */
  delete [] matrixAux;
  delete [] newSeqsName;

  /* Return the new alignment reference */
  return newAlig;
}

/* Function for computing the complementary alignment. It just turn around the
 * current columns/sequences selection */
void alignment::computeComplementaryAlig(bool residues, bool sequences) {
  int i;

  for(i = 0; i < residNumber && residues; i++)
    saveResidues[i] = (saveResidues[i] == -1) ? i : -1;

  for(i = 0; i < sequenNumber && sequences; i++)
    saveSequences[i] = (saveSequences[i] == -1) ? i : -1;
}

/* Function designed for identifying right and left borders between central
 * and terminal regions in the alignment. The borders are those columns, first
 * and last, composed by only residues. Everything inbetween left and right
 * borders are keept independently of the applied methods */
bool alignment::removeOnlyTerminal(void) {

  int i, left_boundary, right_boundary;
  const int *gInCol;

  /* Get alignments gaps stats and copy it */
  if(calculateGapStats() != true) {
    cerr << endl << "WARNING: Impossible to apply 'terminal-only' method"
      << endl << endl;
    return false;
  }
  gInCol = sgaps -> getGapsWindow();

  /* Identify left and right borders. First and last columns with no gaps */
  for(i = 0; i < residNumber && gInCol[i] != 0; i++) ;
  left_boundary = i;

  for(i = residNumber - 1; i > -1 && gInCol[i] != 0; i--) ;
  right_boundary = i + 1;

  /* Once the interal boundaries have been established, if these limits exist
   * then retrieved all columns inbetween both boundaries. Columns out of these
   * limits will remain selected or not depending on the algorithm applied */
  for(i = left_boundary; i < right_boundary; i++)
    if(saveResidues[i] == -1)
      saveResidues[i] = i;

  return true;
}

/* Function designed to identify and remove those columns blocks smaller than
 * a given size */
void alignment::removeSmallerBlocks(int blockSize) {

  int i, j, pos, block;

  if(blockSize == 0)
    return ;

  /* Traverse the alignment looking for blocks greater than BLOCKSIZE, everytime
   * than a column hasn't been selected, check whether the current block is big
   * enough to be kept or it should be removed from the final alignment */
  for(i = 0, pos = 0, block = 0; i < residNumber; i++) {
    if(saveResidues[i] != -1)
      block++;
    else {
      /* Remove columns from blocks smaller than input blocks size */
      if(block < blockSize)
        for(j = pos; j <= i; j++)
          saveResidues[j] = -1;
      pos = i + 1;
      block = 0;
    }
  }

  /* Check final block separately since it could happen than last block is not
   * big enough but because the loop end could provoke to ignore it */
  if(block < blockSize)
    for(j = pos; j < i; j++)
      saveResidues[j] = -1;
  return ;
}

/* Function designed to identify columns and sequences composed only by gaps.
 * Once these columns/sequences have been identified, they are removed from
 * final alignment. */
newValues alignment::removeCols_SeqsAllGaps(void) {
  int i, j, valid, gaps;
  bool warnings = false;
  newValues counter;

  /* Check all valid columns looking for those composed by only gaps */
  for(i = 0, counter.residues = 0; i < residNumber; i++) {
    if(saveResidues[i] == -1)
      continue;

    for(j = 0, valid = 0, gaps = 0; j < sequenNumber; j++) {
      if (saveSequences[j] == -1)
        continue;
      if (sequences[j][i] == '-')
        gaps ++;
      valid ++;
    }
    /* Once a column has been identified, warm about it and remove it */
    if(gaps == valid) {
      if(!warnings)
        cerr << endl;
      warnings = true;
      cerr << "WARNING: Removing column '" << i << "' composed only by gaps"
        << endl;
      saveResidues[i] = -1;
    } else {
      counter.residues ++;
    }
  }

  /* Check for those selected sequences to see whether there is anyone with
   * only gaps */
  for(i = 0, counter.sequences = 0; i < sequenNumber; i++) {
    if(saveSequences[i] == -1)
      continue;

    for(j = 0, valid = 0, gaps = 0; j < residNumber; j++) {
      if(saveResidues[j] == -1)
        continue;
      if (sequences[i][j] == '-')
        gaps ++;
      valid ++;
    }
    /* Warm about it and remove each sequence composed only by gaps */
    if(gaps == valid) {
      if(!warnings)
        cerr << endl;
      warnings = true;

      if(keepSequences) {
        cerr << "WARNING: Keeping sequence '" << seqsName[i]
          << "' composed only by gaps" << endl;
        counter.sequences ++;
      } else {
        cerr << "WARNING: Removing sequence '" << seqsName[i]
          << "' composed only by gaps" << endl;
        saveSequences[i] = -1;
      }
    } else {
      counter.sequences ++;
    }
  }
  if(warnings)
    cerr << endl;

  counter.matrix = new string[counter.sequences];
  counter.seqsName = new string[counter.sequences];

  return counter;
}

/* Function for copying to previously allocated memory those data selected
 * for being in the final alignment */
void alignment::fillNewDataStructure(string *newMatrix, string *newNames) {
  int i, j, k;

  /* Copy only those sequences/columns selected */
  for(i = 0, j = 0; i < sequenNumber; i++) {
    if(saveSequences[i] == -1)
      continue;

    newNames[j] = seqsName[i];
    for(k = 0; k < residNumber; k++) {
      if(saveResidues[k] == -1)
        continue;
      newMatrix[j].resize(newMatrix[j].size() + 1, sequences[i][k]);
    }
    j++;
  }
}

/* Function for copying to previously allocated memory those data selected
 * for being in the final alignment */
void alignment::fillNewDataStructure(newValues *data) {
  int i, j, k;

  /* Copy only those sequences/columns selected */
  for(i = 0, j = 0; i < sequenNumber; i++) {
    if(saveSequences[i] == -1)
      continue;

    data -> seqsName[j] = seqsName[i];
    for(k = 0; k < residNumber; k++) {
      if(saveResidues[k] == -1)
        continue;
      data -> matrix[j].resize(data -> matrix[j].size() + 1, sequences[i][k]);
    }
    j++;
  }
  // cerr << data -> seqsName[j-1] << endl;
}

/* Check if CDS file is correct based on: Residues are DNA/RNA (at most). There
 * is not gaps in the whole dataset. Each sequence is multiple of 3. At the same
 * time, the function will remove stop codons if appropiate flags are used */
bool alignment::prepareCodingSequence(bool splitByStopCodon, bool ignStopCodon,\
  alignment *proteinAlig) {

  bool warning = false;
  size_t found;
  int i;

  /* New code: We now care about the presence of wildcards/indeterminations
   * characters such as 'X' or 'B' into protein sequences as well as about the
   * presence of Selenocysteines ('U') or Pyrrolysines ('O'). It only works, by
   * now, with the universal genetic code */
  int *protSeqsLengths, numbProtSeqs, current_prot;
  string *protSeqsNames, *protSequences;
  char aminoAcid;

  numbProtSeqs =  proteinAlig -> getNumSpecies();
  protSeqsNames = new string[numbProtSeqs];
  protSequences = new string[numbProtSeqs];
  protSeqsLengths = new int[numbProtSeqs];

  proteinAlig -> getSequences(protSeqsNames, protSequences, protSeqsLengths);

  /* Check read sequences are real DNA/RNA */
  if (getTypeAlignment() == AAType) {
    cerr << endl << "ERROR: Check input CDS file. It seems to content protein "
      << "residues." << endl << endl;
      return false;
  }

  for(i = 0; i < sequenNumber; i++) {

    /* Get protein sequence to compare against any potential stop codon in the
     * coding sequence. If there is not protein sequence for current coding
     * sequence, skip its analysis */
    for(current_prot = 0; current_prot < numbProtSeqs; current_prot++)
      if(protSeqsNames[current_prot] == seqsName[i])
        break;
    if(current_prot == numbProtSeqs)
      continue;

    if(sequences[i].find("-") != string::npos) {
      if (!warning)
        cerr << endl;
      cerr << "ERROR: Sequence \"" << seqsName[i] << "\" has, at least, one gap"
        << endl << endl;
      return false;
    }

    if((sequences[i].length() % 3) != 0) {
      if (!warning)
        cerr << endl;
      warning = true;
      cerr << "WARNING: Sequence length \"" << seqsName[i] << "\" is not "
        << "multiple of 3 (length: " << sequences[i].length() << ")" << endl;
    }

    /* Ignore stop codons from the CDS if set by the user */
    if (ignStopCodon)
      continue;

    /* Detect universal stop codons in the CDS. Then, compare those stop codons
     * against the protein sequence to see whether they are real stop codons or
     * are representing rare amino-acids such as Selenocysteines or Pyrrolysines
     * It also allows stop-codons when there are wildcards/indet characters in
     * the protein sequence. CDS sequences could be splitted using stop codons
     * from the sequence itself */

    /* Initialize first appearence of a given stop codon to -1.
     * That means that it has not been found yet */
    found = -1;
    do {
      found = sequences[i].find("TGA", found + 1);

      /* If a stop codon has been found and its position is multiple of 3.
       * Analize it */
      if((found != string::npos) && (((int) found % 3) == 0)) {

        aminoAcid = (char) toupper(protSequences[current_prot][(int) found/3]);
        /* It may be a Selenocysteine ('TGA') which should be represented as 'U'
         * or wildcard/indet characters such as 'X' or 'B' */
        //~ if ((aminoAcid == 'U') || (aminoAcid == 'X') || (aminoAcid == 'B'))
        /* If a rare amino-acids such as 'U'/'O' or a wildcard/indet character
         * such as 'B'/'X' is present, skip current stop codon */
        if((aminoAcid == 'U') || (aminoAcid == 'O') || (aminoAcid == 'X') || \
          (aminoAcid == 'B'))
          continue;

        /* If split_by_stop_codon flag is activated then cut input CDS sequence
         * up to first appearance of a stop codon */
        else if(splitByStopCodon) {
          if (!warning)
            cerr << endl;
          warning = true;
          cerr << "WARNING: Cutting sequence \"" << seqsName[i] << "\" at first"
            << " appearance of stop codon \"TGA\" (residue \"" << aminoAcid
            << "\") at position " << (int) found + 1 << " (length: "
            << sequences[i].length() << ")" << endl;
          sequences[i].resize((int) found);
        }
        /* Otherwise, warn about it and return an error */
        else {
          if (!warning)
            cerr << endl;
          cerr << "ERROR: Sequence \"" << seqsName[i] << "\" has stop codon \""
            << "TGA\" (residue \"" << aminoAcid << "\") at position "
            << (int) found + 1 << " (length: " << sequences[i].length() << ")"
            << endl << endl;
          return false;
        }
      }
    /* Iterate over the CDS until not stop codon is found */
    } while(found != string::npos);

    /* Initialize first appearence of a given stop codon to -1.
     * That means that it has not been found yet */
    found = -1;
    do {
      found = sequences[i].find("TAA", found + 1);

      /* If a stop codon has been found and its position is multiple of 3.
       * Analize it */
      if((found != string::npos) && (((int) found % 3) == 0)) {

        aminoAcid = (char) toupper(protSequences[current_prot][(int) found/3]);
        /* Check if there is any wildcard/indet characters such as 'X' or 'B' */
        //~ if ((aminoAcid == 'X') || (aminoAcid == 'B'))
        /* If a rare amino-acids such as 'U'/'O' or a wildcard/indet character
         * such as 'B'/'X' is present, skip current stop codon */
        if((aminoAcid == 'U') || (aminoAcid == 'O') || (aminoAcid == 'X') || \
          (aminoAcid == 'B'))
          continue;

        /* If split_by_stop_codon flag is activated then cut input CDS sequence
         * up to first appearance of a stop codon */
        else if(splitByStopCodon) {
          if (!warning)
            cerr << endl;
          warning = true;
          cerr << "WARNING: Cutting sequence \"" << seqsName[i] << "\" at first"
            << " appearance of stop codon \"TAA\" (residue \"" << aminoAcid
            << "\") at position " << (int) found + 1 << " (length: "
            << sequences[i].length() << ")" << endl;
          sequences[i].resize((int) found);
        }
        /* Otherwise, warn about it and return an error */
        else {
          if (!warning)
            cerr << endl;
          cerr << "ERROR: Sequence \"" << seqsName[i] << "\" has stop codon \""
            << "TAA\" (residue \"" << aminoAcid << "\") at position "
            << (int) found + 1 << " (length: " << sequences[i].length() << ")"
            << endl << endl;
          return false;
        }
      }
    /* Iterate over the CDS until not stop codon is found */
    } while(found != string::npos);

    /* Initialize first appearence of a given stop codon to -1.
     * That means that it has not been found yet */
    found = -1;
    do {
      found = sequences[i].find("TAG", found + 1);
      /* If a stop codon has been found and its position is multiple of 3.
       * Analize it */
      if((found != string::npos) && (((int) found % 3) == 0)) {

        aminoAcid = (char) toupper(protSequences[current_prot][(int) found/3]);
        /* It may be a Pyrrolysine ('TAG') which should be represented as 'O'
         * or wildcard/indet characters such as 'X' or 'B' */
        //~ if ((aminoAcid == 'O') || (aminoAcid == 'X') || (aminoAcid == 'B'))
        /* If a rare amino-acids such as 'U'/'O' or a wildcard/indet character
         * such as 'B'/'X' is present, skip current stop codon */
        if((aminoAcid == 'U') || (aminoAcid == 'O') || (aminoAcid == 'X') || \
          (aminoAcid == 'B'))
          continue;

        /* If split_by_stop_codon flag is activated then cut input CDS sequence
         * up to first appearance of a stop codon */
        else if(splitByStopCodon) {
          if (!warning)
            cerr << endl;
          warning = true;
          cerr << "WARNING: Cutting sequence \"" << seqsName[i] << "\" at first"
            << " appearance of stop codon \"TAG\" (residue \"" << aminoAcid
            << "\") at position " << (int) found + 1 << " (length: "
            << sequences[i].length() << ")" << endl;
          sequences[i].resize((int) found);
        }
        /* Otherwise, warn about it and return an error */
        else {
          if (!warning)
            cerr << endl;
          cerr << "ERROR: Sequence \"" << seqsName[i] << "\" has stop codon \""
            << "TAG\" (residue \"" << aminoAcid << "\") at position "
            << (int) found + 1 << " (length: " << sequences[i].length() << ")"
            << endl << endl;
          return false;
        }
      }
    /* Iterate over the CDS until not stop codon is found */
    } while(found != string::npos);
  }

  /* If everything was return an OK to informat about it. */
  return true;
}

/* Function designed to check whether input CDS file is correct or not based on
 * some features: Sequences are in both files (it could be more on CDS file),
 * they have (more or less) same ength. Otherwise, some nucleotides could be
 * excluded or some 'N's added to fit protein length. */
bool alignment::checkCorrespondence(string *names, int *lengths, int \
  totalInputSeqs, int multiple = 1) {

  int i, j, seqLength, indet;
  bool warnings = false;
  string tmp;

  /* For each sequence in the current protein alignment, look for its coding
   * DNA sequence checking that they have the same size. */
  for(i = 0; i < sequenNumber; i++) {

    /* Get protein sequence length removing any possible gap. Get as well last
     * residue from current sequence */

    tmp = utils::removeCharacter('-', sequences[i]);
    seqLength = tmp.length() * multiple;
    indet = ((int) tmp.length() - utils::min((int) tmp.find_last_not_of("X"), \
      (int) tmp.find_last_not_of("x"))) - 1;

    /* Go through all available CDS looking for the one with the same ID */
    for(j = 0; j < totalInputSeqs; j++) {

      /* Once both ID matchs, compare its lengths */
      if(seqsName[i] == names[j]) {

        /* If both sequences have the same length, stop the search */
        if(seqLength == lengths[j])
          break;

        /* If nucleotide sequence is larger than protein sequence, warn about
         * it and continue the verification process. It will used the 'Nth'
         * first nucleotides for the conversion */
        else if(seqLength < lengths[j]) {
          if (!warnings)
            cerr << endl;
          warnings = true;
          cerr << "WARNING: Sequence \"" << seqsName[i] << "\" will be cutted "
            << "at position " << seqLength << " (length: "<< lengths[j] << ")"
            << endl;
          break;
        }

        /* It has been detected some indeterminations at the end of the protein
         * sequence. That issue could be cause by some incomplete codons in the
         * nucleotide sequences. This issue is solved adding as much 'N' symbols
         * as it is needed to preserve the backtranslated alignment */
        else if((indet > 0) && (indet > (seqLength - lengths[j])/3)) {
          if (!warnings)
            cerr << endl;
          warnings = true;
          cerr << "WARNING: Sequence \"" << seqsName[i] << "\" has some inde"
            << "termination symbols 'X' at the end of sequence. They will be"
            << " included in the final alignment." << endl;
          break;
        }

        /* If nucleotide sequence is shorter than protein sequence, return an
         * error since it is not feasible to cut the input protein aligment to
         * fit it into CDNA sequences size */
        else {
          if (!warnings)
            cerr << endl;
          warnings = true;
          cerr << "WARNING: Sequence \"" << seqsName[i] << "\" has less nucleo"
            << "tides (" << lengths[j] << ") than expected (" << seqLength
            << "). It will be added N's to complete the sequence"  << endl;
          break;
        }
      }
    }

    /* Warn about a mismatch a sequences name level */
    if(j == totalInputSeqs) {
      cerr << endl << "ERROR: Sequence \"" << seqsName[i] << "\" is not in "
        << "CDS file." << endl << endl;
      return false;
    }
  }

  /* If everything is OK, return an appropiate flag */
  return true;
}
