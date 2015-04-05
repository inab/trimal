/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
   ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

    statAl v1.4: a tool for getting descriptive alignment features/scores.

    2009-2013 Capella-Gutierrez S. and Gabaldon, T.
              [scapella, tgabaldon]@crg.es

    This file is part of statAl.

    statAl is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, the last available version.

    statAl is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with trimAl. If not, see <http://www.gnu.org/licenses/>.

***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

#include <fstream>
#include <iostream>
#include <iomanip>

#include <stdlib.h>
#include <string.h>

#include "compareFiles.h"
#include "compareFiles.h"
#include "alignment.h"
#include "defines.h"
#include "utils.h"

void show_menu(void);
void show_examples(void);

int main(int argc, char *argv[]){

  /* Input values */
  char *inFile = NULL, *forceFile = NULL, *setAlignments = NULL, *matrix = NULL;
  int windowSize = -1, gapWindow = -1, simWindow = -1, conWindow = -1;
  bool stats_gaps_columns = 0, stats_gaps_dist = 0, stats_simil_columns = 0,
    stats_simil_dist = 0, stats_seqs_ident = 0, stats_col_ident_gen = 0,
    stats_file_columns = 0, stats_file_dist = 0;
  alignment *origAlig = NULL, **compAlig  = NULL;

  /* Internal variables */
  int i = 1, numFiles = 0, maxResidues = 0, referFile = 0, alignDataType = -1;
  similarityMatrix *similMatrix = NULL;
  char **filesToCompare = NULL;
  bool appearErrors = false;
  float *compareVect = NULL;
  ifstream algsPaths;
  string line;

  /* ***** ***** ***** ***** ***** Help functions ***** ***** ***** ***** *** */
  /* Show help and exit either help flag is set or not arguments are provided */
  if((argc == 1) || ((!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) &&
    (i+1 == argc))) {
    show_menu();
    show_examples();
    return 0;
  }

  /* Show version and revision if it is asked for */
  if((!strcmp(argv[i], "-v") || !strcmp(argv[i], "--version")) &&
    (i+1 == argc)) {
    cout << endl << "statAl v" << VERSION << ".rev" << REVISION
         << " build[" << BUILD << "]" << endl << endl;
    return 0;
  }

  /* Allocate local memory for input alignment */
  origAlig = new alignment;

  /* ***** ***** ***** ***** Input parameters loop ***** ***** ***** ***** ** */
  while(i < argc) {

    /* Get input alignment */
    if((!strcmp(argv[i], "-i") || !strcmp(argv[i], "-in")) && (i+1 != argc) &&
      (inFile == NULL)) {

      /* Try to load input alignment */
      inFile = new char[(int) strlen(argv[++i]) + 1];
      strcpy(inFile, argv[i]);

      if(!origAlig -> loadAlignment(inFile)) {
        cerr << endl << "ERROR: Alignment not loaded: \"" << inFile
          << "\" Check input aligmment content." << endl << endl;
        appearErrors = true;
      }
    }

    /* Force selection of a specific input alignment as a reference to perform
     * alignment comparisons */
    else if(!strcmp(argv[i], "-forceselect") && (i+1 != argc) &&
    (forceFile == NULL)) {

      forceFile = new char[(int) strlen(argv[++i])  + 1];
      strcpy(forceFile, argv[i]);

      if(!origAlig -> loadAlignment(forceFile)) {
        cerr << endl << "ERROR: Alignment not loaded: \"" << forceFile
          << "\" Check input aligmment content." << endl << endl;
        appearErrors = true;
      }
    }

    /* Load a set of alignments to perform alignments comparisos among them */
    else if(!strcmp(argv[i], "-compareset") && (i+1 != argc) &&
      (setAlignments == NULL)) {
      setAlignments = new char[(int) strlen(argv[++i])  + 1];
      strcpy(setAlignments, argv[i]);
    }

    /* Load a similarity matrix */
    else if(!strcmp(argv[i], "-matrix") && (i+1 != argc) && (matrix == NULL)) {
      matrix = new char[int(strlen(argv[++i])) + 1];
      strcpy(matrix, argv[i]);

      similMatrix = new similarityMatrix();
      if(!similMatrix -> loadSimMatrix(matrix)) {
        cerr << endl << "ERROR: Similarity Matrix not loaded: \"" << matrix
          << "\" Check input file content." << endl << endl;
        appearErrors = true;
      }
    }

    /* Set flags for computing different stats from input alignments */
    else if((!strcmp(argv[i], "-sgc")) && (!stats_gaps_columns))
      stats_gaps_columns = true;
    else if((!strcmp(argv[i], "-sgt")) && (!stats_gaps_dist))
      stats_gaps_dist = true;
    else if((!strcmp(argv[i], "-ssc")) && (!stats_simil_columns))
      stats_simil_columns = true;
    else if((!strcmp(argv[i], "-sst")) && (!stats_simil_dist))
      stats_simil_dist = true;
    else if((!strcmp(argv[i], "-sident")) && (!stats_seqs_ident))
      stats_seqs_ident = true;
    else if((!strcmp(argv[i], "-scolidentt")) && (!stats_col_ident_gen))
      stats_col_ident_gen = true;
    else if((!strcmp(argv[i], "-sfc")) && (!stats_file_columns))
      stats_file_columns = true;
    else if((!strcmp(argv[i], "-sft")) && (!stats_file_dist))
      stats_file_dist = true;

    /* Windows size parameters. Setting others values than default could affect
     * stats computed from selected alignment */
    else if(!strcmp(argv[i], "-w") && (i+1 != argc) && (windowSize == -1)) {
      if(utils::isNumber(argv[i+1])) {
        windowSize = atoi(argv[++i]);
      }
      else {
        cerr << endl << "ERROR: Window size value should be a number\n\n.";
        appearErrors = true;
      }
      if(windowSize < 1) {
        cerr << endl << "ERROR: Windows size value should be equal or greater "
          << "than 1. Check your command-line parameter" << endl << endl;
        appearErrors = true;
      }
    }
    /* Window size value only affecting stats related to gaps */
    else if(!strcmp(argv[i], "-gw") && (i+1 != argc) && (gapWindow == -1)){
      if(utils::isNumber(argv[i+1])) {
        gapWindow = atoi(argv[++i]);
      }
      else {
        cerr << endl << "ERROR: Gap Window size value should be a number\n\n.";
        appearErrors = true;
      }
      if(gapWindow < 1) {
        cerr << endl << "ERROR: Gap Windows size should be equal or greater "
          << "than 1. Check your command-line parameter" << endl << endl;
        appearErrors = true;
      }
    }
    /* Window size value only affecting stats related to similarity */
    else if(!strcmp(argv[i], "-sw") && (i+1 != argc) && (simWindow == -1)) {
      if(utils::isNumber(argv[i+1])) {
        simWindow = atoi(argv[++i]);
      }
      else {
        cerr << endl << "ERROR: Similarity Window size should be a number\n\n.";
        appearErrors = true;
      }
      if(simWindow < 1) {
        cerr << endl << "ERROR: Similarity Windows size should be equal or "
          << "gretaer than 1. Check your command-line parameter\n\n";
        appearErrors = true;
      }
    }
    /* Window size only affecting stats related to alignment comparisons */
    else if(!strcmp(argv[i], "-cw") && (i+1 != argc) && (conWindow == -1)) {
      if(utils::isNumber(argv[i+1])) {
        conWindow = atoi(argv[++i]);
      }
      else {
        cerr << endl << "ERROR: Consistency Window size should be a number\n\n";
        appearErrors = true;
      }
      if(conWindow < 1) {
        cerr << endl << "ERROR: Consistency Windows size should be equal or "
          << "gretaer than 1. Check your command-line parameter\n\n";
        appearErrors = true;
      }
    }

    else {
      cerr << endl << "ERROR: Parameter \"" << argv[i] << "\" not valid\n\n.";
      appearErrors = true;
    }
    i++;

    if(appearErrors)
      break;
  }

  /* ***** ***** ***** ** Control input parameters errors ***** ***** ***** * */
  if((!appearErrors) && (inFile != NULL)) {
    if((stats_file_columns) || (stats_file_dist)) {
      cerr << endl << "ERROR: No stats about alignments comparison can be "
        << "obtained from Single Input Alignment" << endl << endl;
      appearErrors = true;
    }
    if((forceFile != NULL) || (setAlignments != NULL)) {
      cerr << endl << "ERROR: Incompatibilities detected: Sigle Alignment "
        << "provided as well as Set of Alignments for being compared\n\n";
      appearErrors = true;
    }
  }

  if((!appearErrors) && (forceFile != NULL)) {
    if(setAlignments == NULL) {
      cerr << endl << "ERROR: An alignmnet has been set as reference but it has"
        << "not been defined a set to compare with" << endl << endl;
      appearErrors = true;
    }
  }

  if((!appearErrors) && (stats_file_columns || stats_file_dist)) {
    if(setAlignments == NULL) {
      cerr << endl << "ERROR: A set of alignments should be provided to compute"
        << "stats about its comparison" << endl << endl;
      appearErrors = true;
    }
  }

  if((!appearErrors) && (windowSize != -1)) {
    if((gapWindow != -1) || (simWindow != -1) || (conWindow != -1)) {
      cerr << endl << "ERROR: General and specific windows size specified. Just"
        << " select one of them" << endl << endl;
      appearErrors = true;
    }
  }

  if((!appearErrors) && (setAlignments != NULL)) {

    /* Open and check whether input file containing the path to alignments for
     * being compared is correct or not */
    algsPaths.open(setAlignments, ifstream::in);
    if(!algsPaths) {
      cerr << endl << "ERROR: Check input file \"" << setAlignments
        << "\" containing alignments paths" << endl << endl;
      appearErrors = true;
    }
    /* Count number of alignments */
    while(getline(algsPaths, line))
      numFiles++;
    algsPaths.close();
  }

  /* ***** ***** ***** ** Input set of alignments: Load it ***** ***** ***** */
  if((!appearErrors) && (numFiles != 0)) {

    /* allocate memory for input alignments */
    filesToCompare = new char*[numFiles];
    compAlig = new alignment*[numFiles];

    algsPaths.open(setAlignments, ifstream::in);
    for(i = 0; (i < numFiles) && (!appearErrors); i++) {

      /* Get alignment path */
      getline(algsPaths, line);
      if(line.size() == 0)
        continue;
      /* Store alignment path and try to load it */
      filesToCompare[i] = new char [line.size() + 1];
      strcpy(filesToCompare[i], line.c_str());

      /* Check currently load alignment */
      compAlig[i] = new alignment;
      if(!compAlig[i] -> loadAlignment(filesToCompare[i])) {
        cerr << endl << "ERROR: Alignment not loaded: \"" << filesToCompare[i]
          << "\". Check input file content" << endl << endl;
        appearErrors = true;
        break;
      }
      /* Check whether each individual alignment is aligned. It is impossible
       * to make any comparison with unaligned sequences */
      if(!compAlig[i] -> isFileAligned()) {
        cerr << endl << "ERROR: Input alignment \"" << filesToCompare[i] << "\""
          << "should have sequences aligned for computing any stats from them."
          << endl << endl;
        appearErrors = true;
        break;
      }
      /* Check all input alignments contain the same kind of biological data */
      if(alignDataType == -1)
        alignDataType = compAlig[i] -> getTypeAlignment();
      else if(compAlig[i] -> getTypeAlignment() != alignDataType) {
        cerr << endl << "ERROR: Check input alignment data-type. Alignment \""
          << filesToCompare[0] << "\" contains a different data-type different "
          << "than \"" << filesToCompare[i] << "\"" << endl << endl;
        appearErrors = true;
        break;
      }
      /* Construct residues positions matrix to make possible any comparison */
      compAlig[i] -> sequenMatrix();
      if(compAlig[i] -> getNumAminos() > maxResidues)
        maxResidues = compAlig[i] -> getNumAminos();
    }
    
    /* ***** ***** ***** ** Input set of alignments: Compare it ***** ***** * */
    if((!appearErrors) && (forceFile == NULL)) {
      /* If reference alignment has not been set, just compute the most
       * consistent one and choose it as the reference one */
      compareVect = new float[maxResidues];
      referFile = compareFiles::algorithm(compAlig, filesToCompare, compareVect,
        numFiles, true);
      origAlig -> loadAlignment(filesToCompare[referFile]);
    }
    else if((!appearErrors) && (forceFile != NULL)) {
      /* Compute consistency vector for the aligment set as the reference one */
      compareVect = new float[origAlig -> getNumAminos()];
      appearErrors = !(compareFiles::forceComparison(compAlig, numFiles,
        origAlig, compareVect));
    }
    /* Apply any possible windows for future stats */
    conWindow = (conWindow != -1) ? conWindow : windowSize;
    if((conWindow != -1) && (!appearErrors))
      compareFiles::applyWindow(origAlig -> getNumAminos(), conWindow,
        compareVect);
    /* Deallocate memory for the alignments used to make the comparisons */
    for(i = 0; i < numFiles; delete compAlig[i], i++)
      delete filesToCompare[i];
  }
  
  /* ***** ***** ***** ***** ***** Final Checks ***** ***** ***** ***** ***** */
  if((!appearErrors) && (!origAlig -> isFileAligned() &&  (stats_gaps_columns or
    stats_gaps_dist or stats_simil_columns or stats_simil_dist or stats_seqs_ident
    or stats_col_ident_gen or stats_file_columns or stats_file_dist))) {
    cerr << endl << "ERROR: Reference alignment should aligned to compute any "
      << "stats using it" << endl << endl;
    appearErrors = true;
  }

  /* It is mandatory to choose an option for processing input alignment */
  if((!appearErrors) && (!stats_gaps_columns and !stats_gaps_dist and
    !stats_simil_columns and !stats_simil_dist and !stats_seqs_ident and
    !stats_col_ident_gen and !stats_file_columns and !stats_file_dist)) {
    cerr << endl << "ERROR: An option has to be chosen." << endl << endl;
    appearErrors = true;
  }
  
  /* ***** ***** ***** ***** * Apply input parameters ***** ***** ***** ***** */
  if(!appearErrors) {
    /* Apply window size to either input single alignment or the most consistent
     * one */
    gapWindow = (windowSize != -1) ? windowSize : (gapWindow != -1) ?
      gapWindow : 0;
    simWindow = (windowSize != -1) ? windowSize : (simWindow != -1) ?
      simWindow : 0;
    origAlig -> setWindowsSize(gapWindow, simWindow);
  }
  /* ***** ***** ***** ***** * Load Similarity Matrix ***** ***** ***** ***** */
  if(!appearErrors) {
    if((stats_simil_columns) || (stats_simil_dist)) {
      /* Load predefined similarity matrix, if none has been set by the user */
      if(matrix == NULL){
        similMatrix = new similarityMatrix();

        if((origAlig -> getTypeAlignment()) == AAType)
          similMatrix -> defaultAASimMatrix();
        else
          similMatrix -> defaultNTSimMatrix();
      }
      /* Assign similarity matrix, either predefined one or supply by the user,
       * to the alignment that is used to compute stats */
      if(!origAlig -> setSimilarityMatrix(similMatrix)) {
        cerr << endl << "ERROR: Impossible to associate similarity matrix to "
          << "selected alignment. Check input matrix file content\n\n";
        appearErrors = true;
      }
    }
  }
  /* ***** ***** ***** ***** ** Compute/Show stats ***** ***** ***** ***** ** */
  if(!appearErrors) {
    if(stats_gaps_columns) {
      cout << endl << "## Gaps scores per column" << endl;
      origAlig -> printStatisticsGapsColumns();
    }
    if(stats_gaps_dist) {
      cout << endl << "## Gaps scores distribution" << endl;
      origAlig -> printStatisticsGapsTotal();
    }
    if(stats_simil_columns) {
      cout << endl << "## Similarity values per column" << endl;
      origAlig -> printStatisticsConservationColumns();
    }
    if(stats_simil_dist) {
      cout << endl << "## Similarity values distribution" << endl;
      origAlig -> printStatisticsConservationTotal();
    }
    if(stats_seqs_ident) {
      cout << endl << "## Sequence Identity Values" << endl;
      origAlig -> printSeqIdentity();
    }
    if(stats_col_ident_gen) {
      cout << endl << "## Columns Identity Descriptive Statistics" << endl;
      origAlig -> printColumnsIdentity_DescriptiveStats();
    }
    if(stats_file_columns) {
      cout << endl << "## Consistency values per column for selected align\n";
      compareFiles::printStatisticsFileColumns(origAlig -> getNumAminos(),
        compareVect);
    }
    if(stats_file_dist) {
      cout << endl << "## Consistency values distribution for selected align\n";
      compareFiles::printStatisticsFileAcl(origAlig -> getNumAminos(),
        compareVect);
    }
  }
  /* ***** ***** ***** ***** *** Deallocate memory ***** ***** ***** ***** ** */
  delete origAlig;
  delete[] compAlig;

  delete similMatrix;

  delete[] filesToCompare;
  delete[] compareVect;

  delete[] inFile;
  delete[] matrix;
  delete[] forceFile;

  return (appearErrors == true ? -1 : 0);
}

void show_menu(void) {

  cout << endl
    << "statAl v" << VERSION << ".rev" << REVISION << " build[" << BUILD
    << "]. " << AUTHORS << endl << endl;

  cout << "statAl/trimAl webpage: http://trimal.cgenomics.org" << endl << endl;

  cout << "This program is free software: you can redistribute it and/or modify"
    << "\nit under the terms of the GNU General Public License as published by "
    << "the\nFree Software Foundation, the last available version.\n\n";

  cout << "Please cite:" << endl
    << "\ttrimAl: a tool for automated alignment trimming in large-scale "
    << "phylogenetic analyses.\n\tSalvador Capella-Gutierrez; Jose M. Silla-"
    << "Martinez; Toni Gabaldon. Bioinformatics 2009, 25:1972-1973.\n\n";

  cout << "Basic usage:" << endl
       << "\tstatal -in <inputfile> (options)." << endl << endl;

  cout << "Available options:" << endl << endl;
  cout << "    -h --help                "
    << "Print this information and show some examples." << endl;
  cout << "    -v --version             "
    << "Print the trimAl version." << endl << endl;

  cout << "    -i -in <inputfile>       "
    << "Input file in several formats (clustal, fasta, nexus, phylip, etc)."
    << endl << endl;

  cout << "    -compareset <inputfile>  "
    << "Input list of paths for the alignments to compare." << endl;
  cout << "    -forceselect <inputfile> "
    << "Force selection of a given file as reference for being compare with "
    << "others." << endl << endl;

  cout << "    -matrix <inpufile>       "
    << "Input file for user-defined similarity matrix (default: Blosum62)."
    << endl << endl;

  cout << "    -sgc                     "
    << "Print gap score per column from input alignment." << endl;
  cout << "    -sgt                     "
    << "Print accumulated gap scores distribution from input alignment."
    << endl << endl;

  cout << "    -ssc                     "
    << "Print similarity score per column from input alignment." << endl;
  cout << "    -sst                     "
    << "Print accumulated similarity scores distribution for input alignment."
    << endl << endl;

  cout << "    -sfc                     "
    << "Print sum-of-pairs score per column for the selected alignment" << endl;
  cout << "    -sft                     "
    << "Print accumulated sum-of-pairs scores distribution for the selected "
    << "alignment" << endl << endl;

  cout << "    -sident                  "
    << "Print identity scores for sequences in the alignemnt." << endl;
  cout << "    -scolidentt              "
    << "Print general descriptive statistics for column identity scores from "
    << "input alignemnt." << endl << endl;

  cout << "    -w <n>                   "
    << "(half) Window size, score of position i is the average of the window "
    << "(i - n) to (i + n)." << endl;
  cout << "    -gw <n>                  "
    << "(half) Window size only applies to statistics based on Gaps." << endl;
  cout << "    -sw <n>                  "
    << "(half) Window size only applies to statistics based on Similarity.\n";
  cout << "    -cw <n>                  "
    << "(half) Window size only applies to statistics based on Consistency."
    << endl << endl;
}

void show_examples(void) {

  cout << "Some Examples:" << endl << endl;

  cout << "1) Get information about gaps distribution for input alignment"
       << endl << "   statal -in <inputfile> -sgt" << endl << endl;

  cout << "2) Get information about consistency score per column for the most "
       << "consistent input alignment, if more than one is provided"
       << endl << "   statal -in <inputfile> -sfc" << endl << endl;

  cout << "3) Get general descriptive statistics for columns identity"
       << endl << "   statal -in <inputfile> -scolidentt" << endl << endl;

  cout << "4) Change the windows size for computing similarity score per column"
       << endl << "   statal -in <inputfile> -sw 3 -ssc" << endl << endl;
}
